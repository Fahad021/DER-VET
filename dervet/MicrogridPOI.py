"""
Copyright (c) 2021, Electric Power Research Institute

 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification,
 are permitted provided that the following conditions are met:

     * Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.
     * Neither the name of DER-VET nor the names of its contributors
       may be used to endorse or promote products derived from this software
       without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
"""
MicrogridPOI.py

"""

import pandas as pd
from storagevet.POI import POI
import cvxpy as cvx
from storagevet.ErrorHandling import *
import numpy as np


class MicrogridPOI(POI):
    """
        This class holds the load data for the case described by the user defined model parameter.
        It will also impose any constraints that should be opposed at the microgrid's POI.
    """

    def __init__(self, params, technology_inputs_map, technology_class_map):
        super().__init__(params, technology_inputs_map, technology_class_map)
        self.is_sizing_optimization = self.check_if_sizing_ders()
        if self.is_sizing_optimization:
            self.error_checks_on_sizing()

        # add thermal site load time series
        for der in self.active_ders:
            try:
                self.site_steam_load = der.site_steam_load
            except AttributeError:
                self.site_steam_load = None
            try:
                self.site_hotwater_load = der.site_hotwater_load
            except AttributeError:
                self.site_hotwater_load = None
            try:
                self.site_cooling_load = der.site_cooling_load
            except AttributeError:
                self.site_cooling_load = None

    def check_if_sizing_ders(self):
        """ This method will iterate through the initialized DER instances and return a logical OR
        of all of their 'being_sized' methods.

        Returns: True if ANY DER is getting sized

        """
        for der_instance in self.der_list:
            try:
                solve_for_size = der_instance.being_sized()
            except AttributeError:
                solve_for_size = False
            if solve_for_size:
                return True
        return False

    def grab_active_ders(self, indx):
        """ drops DER that are not considered active in the optimization window's horizon

        """
        year = indx.year[0]
        active_ders = [der_inst for der_inst in self.der_list if der_inst.operational(year)]
        self.active_ders = active_ders

    def error_checks_on_sizing(self):
        # perform error checks on DERs that are being sized
        # collect errors and raise if any were found
        errors_found = [1 if der.sizing_error() else 0 for der in self.der_list]
        if sum(errors_found):
            raise ParameterError(f'Sizing of DERs has an error. Please check error log.')

    def is_any_sizable_der_missing_power_max(self):
        return bool(sum([1 if not der_inst.max_power_defined else 0
                         for der_inst in self.der_list]))

    def set_size(self, value_streams, start_year):
        """ part of Deferral's sizing module: TODO USE THIS INSTEAD OF set_size IN MICROGRID SERVICE AGGREGATOR
        iterates over a list of DER+DERExtension objects and sets their minimum size
        based on the P and E requirements set by MIN_YEAR objective.

        Args:
            value_streams:
            start_year:

        Returns: der_list with size minimums

        """
        deferral = value_streams.get('Deferral')
        min_year = deferral.min_years
        last_year_to_defer = start_year.year + min_year - 1
        p_e_req = deferral.deferral_df.loc[last_year_to_defer, :]
        min_power = p_e_req.loc['Power Capacity Requirement (kW)']
        min_energy = p_e_req.loc['Energy Capacity Requirement (kWh)']
        ess_inst = self.der_list[0]
        if len(value_streams.keys()) > 1:
            ess_inst.size_constraints += [cvx.NonPos(min_energy - ess_inst.ene_max_rated)]
            ess_inst.size_constraints += [cvx.NonPos(min_power - ess_inst.ch_max_rated)]
            ess_inst.size_constraints += [cvx.NonPos(min_power - ess_inst.dis_max_rated)]
        else:
            ess_inst.ch_max_rated = min_power
            ess_inst.dis_max_rated = min_power
            ess_inst.ene_max_rated = min_energy

    def is_dcp_error(self, is_binary_formulation):
        """ If trying to sizing power of batteries (or other DERs) AND using the binary formulation
         (of ESS) our linear model will not be linear anymore

        Args:
            is_binary_formulation (bool):

        Returns: a boolean

        """
        solve_for_size = False
        for der_instance in self.der_list:
            if der_instance.tag == 'Battery':
                solve_for_size = solve_for_size or \
                                 (der_instance.is_power_sizing() and is_binary_formulation)
        return solve_for_size

    def get_state_of_system(self, mask):
        """ POI method to measure the state of POI depending on available types of DERs. used in
        SET_UP_OPTIMIZATION. Extends StorageVET's method to take into account types of
        technologies added by DERVET, and thermal recovery

        Args:
            mask (DataFrame): DataFrame of booleans used, the same length as time_series. The value
             is true if the corresponding column in time_series is included in the data to be
             optimized.

        Returns:
            aggregation of loads
            aggregation of generation from variable resources
            aggregation of generation from other sources
            total net power from ESSs
            total state of energy stored in the system
            aggregation of all the power flows into the POI
            aggregation of all the power flows out if the POI

            aggregation of steam thermal heating power (heat recovered)
            aggregation of hotwater thermal heating power (heat recovered)
            aggregation of thermal cooling power (cold recovered)
        """
        # get values from storagevet/poi method
        load_sum, var_gen_sum, gen_sum, tot_net_ess, total_soe, agg_power_flows_in, \
            agg_power_flows_out, agg_steam_heating_power, agg_hotwater_heating_power, \
            agg_thermal_cooling_power = super().get_state_of_system(mask)

        # dervet-specific
        for der_inst in self.active_ders:
            # add to aggregate values for dervet-specific technology-types
            if der_inst.technology_type == 'Electric Vehicle':
                load_sum += der_inst.get_charge(mask)
                # total_soe += der_instance.get_state_of_energy(mask)

            # thermal power recovered: hot (steam/hotwater) and cold
            if der_inst.is_hot:
                if self.site_steam_load is None and self.site_hotwater_load is None:
                    TellUser.warning('A heat source technology is active ' +
                                     f'({der_inst.unique_tech_id()}), but you have set the ' +
                                     f'scenario parameter incl_thermal_load to False. Any ' +
                                     f'thermal load will be ignored.')
                else:
                    if self.site_steam_load is not None:
                        TellUser.debug(f'adding heat (steam) recovered from this DER: ' +
                                       f'{der_inst.unique_tech_id()}')
                        agg_steam_heating_power += der_inst.get_steam_recovered(mask)
                    if self.site_hotwater_load is not None:
                        TellUser.debug(f'adding heat (hotwater) recovered from this DER: ' +
                                       f'{der_inst.unique_tech_id()}')
                        agg_hotwater_heating_power += der_inst.get_hotwater_recovered(mask)
            if der_inst.is_cold:
                if self.site_cooling_load is None:
                    TellUser.warning(f'A cold source technology is active ' +
                                     f'({der_inst.unique_tech_id()}), but you have set the ' +
                                     f'scenario parameter incl_thermal_load to False. Any ' +
                                     f'thermal load will be ignored.')
                else:
                    TellUser.debug(f'adding cold recovered from this DER: ' +
                                   f'{der_inst.unique_tech_id()}')
                    agg_thermal_cooling_power += der_inst.get_cold_recovered(mask)

        return load_sum, var_gen_sum, gen_sum, tot_net_ess, total_soe, agg_power_flows_in, \
            agg_power_flows_out, agg_steam_heating_power, agg_hotwater_heating_power, \
            agg_thermal_cooling_power

    def optimization_problem(self, mask, power_in, power_out, steam_in, hotwater_in, cold_in,
                             annuity_scalar=1):
        """ Builds the master POI constraint list for the subset of time series data being
        optimized. Due to VS power reservations, control constraints, import/export constraints,
        and energy throughput requirements. Builds onto storagevet method to add thermal balance
        constraints for dervet

        Args:
            mask (DataFrame): A boolean array that is true for indices corresponding to
                time_series data included in the subs data set
            power_in (cvx.Expression):
            power_out (cvx.Expression):
            steam_in (cvx.Expression):
            hotwater_in (cvx.Expression):
            cold_in (cvx.Expression):
            annuity_scalar (float): a scalar value to be multiplied by any yearly cost or benefit
                that helps capture the cost/benefit over the entire project lifetime (only to be
                set iff sizing)

        Returns:
            A dictionary with the portion of the objective function that it affects, labeled by
                the expression's key.
            A list of constraints being set by the POI: power reservations, control constraints
                requirements, max import, max export, etc.
        """
        obj_expression, constraint_list = super().optimization_problem(mask, power_in, power_out,
                                                                       steam_in, hotwater_in,
                                                                       cold_in, annuity_scalar)

        # thermal power balance constraints
        if self.site_steam_load is not None:
            if steam_in.variables():
                TellUser.debug('adding steam thermal power balance constraint')
                constraint_list += [cvx.NonPos(-steam_in + self.site_steam_load)]
        if self.site_hotwater_load is not None:
            if hotwater_in.variables():
                TellUser.debug('adding hot water thermal power balance constraint')
                constraint_list += [cvx.NonPos(-hotwater_in + self.site_hotwater_load)]
        if self.site_cooling_load is not None:
            if cold_in.variables():
                TellUser.debug('adding thermal cooling power balance constraint')
                constraint_list += [cvx.NonPos(-cold_in + self.site_cooling_load)]

        return obj_expression, constraint_list

    def sizing_summary(self):
        rows = list(map(lambda der: der.sizing_summary(), self.der_list))
        sizing_df = pd.DataFrame(rows)
        sizing_df.set_index('DER')
        return sizing_df

    def merge_reports(self, is_dispatch_opt, index):
        """ Collects and merges the optimization results for all DERs into
        Builds extends StorageVET's method to take into account types of technologies added by
        DERVET

        Returns: A timeseries dataframe with user-friendly column headers that summarize the
            results pertaining to this instance

        """
        results = pd.DataFrame(index=index)
        monthly_data = pd.DataFrame()

        # initialize all the data columns that will ALWAYS be present in our results
        results.loc[:, 'Total Original Load (kW)'] = 0
        results.loc[:, 'Total Load (kW)'] = 0
        results.loc[:, 'Total Generation (kW)'] = 0
        results.loc[:, 'Total Storage Power (kW)'] = 0
        results.loc[:, 'Aggregated State of Energy (kWh)'] = 0

        for der in self.der_list:
            report_df = der.timeseries_report()
            results = pd.concat([report_df, results], axis=1)
            if is_dispatch_opt:
                if der.technology_type in ['Generator', 'Intermittent Resource']:
                    results.loc[:, 'Total Generation (kW)'] += \
                        results[f'{der.unique_tech_id()} Electric Generation (kW)']
                if der.technology_type == 'Energy Storage System':
                    results.loc[:, 'Total Storage Power (kW)'] += \
                        results[f'{der.unique_tech_id()} Power (kW)']
                    results.loc[:, 'Aggregated State of Energy (kWh)'] += \
                        results[f'{der.unique_tech_id()} State of Energy (kWh)']
                if der.technology_type == 'Load':
                    results.loc[:, 'Total Original Load (kW)'] += \
                        results[f'{der.unique_tech_id()} Original Load (kW)']
                    if der.tag == "ControllableLoad":
                        results.loc[:, 'Total Load (kW)'] += \
                            results[f'{der.unique_tech_id()} Load (kW)']
                    else:
                        results.loc[:, 'Total Load (kW)'] += \
                            results[f'{der.unique_tech_id()} Original Load (kW)']
                if der.technology_type == 'Electric Vehicle':
                    results.loc[:, 'Total Load (kW)'] += \
                        results[f'{der.unique_tech_id()} Charge (kW)']
                    if der.tag == 'ElectricVehicle1':
                        results.loc[:, 'Aggregated State of Energy (kWh)'] += \
                            results[f'{der.unique_tech_id()} State of Energy (kWh)']
            report = der.monthly_report()
            monthly_data = pd.concat([monthly_data, report], axis=1, sort=False)
        # assumes the orginal net load only does not contain the Storage system
        # check if Total Original Load and Total Load are the same.
        if np.all(results['Total Load (kW)'] == results['Total Original Load (kW)']):
            # Drop Total Original Load
            results.drop('Total Original Load (kW)', axis=1, inplace=True)
        # net load is the load see at the POI
        results.loc[:, 'Net Load (kW)'] = \
            results.loc[:, 'Total Load (kW)'] - results.loc[:, 'Total Generation (kW)'] - \
            results.loc[:, 'Total Storage Power (kW)']
        return results, monthly_data
