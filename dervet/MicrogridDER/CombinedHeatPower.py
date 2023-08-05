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
CHP Sizing class

This Python class contains methods and attributes specific for technology analysis within StorageVet.
"""
import cvxpy as cvx
from dervet.MicrogridDER.CombustionTurbine import CT
import storagevet.Library as Lib
from storagevet.ErrorHandling import *


class CHP(CT):
    """ Combined Heat and Power (CHP) Technology, with sizing optimization

    """

    def __init__(self, params):
        """ Initialize all technology with the following attributes.

        Args:
            params (dict): Dict of parameters for initialization
        """
        TellUser.debug(f"Initializing {__name__}")
        # base class is CT
        super().__init__(params)

        # overrides
        self.tag = 'CHP'
        self.is_hot = True

        self.electric_ramp_rate = params['electric_ramp_rate']      # MW/min # TODO this is not being used? --AE
        self.electric_heat_ratio = params['electric_heat_ratio']    # elec/heat (generation)
        self.max_steam_ratio = params['max_steam_ratio']           # steam/hotwater relative ratio
        # time series inputs
        self.site_steam_load = params.get('site_steam_load')    # BTU/hr
        self.site_hotwater_load = params.get('site_hotwater_load')    # BTU/hr

    def grow_drop_data(self, years, frequency, load_growth):
        if self.site_steam_load is not None:
            self.site_steam_load = Lib.fill_extra_data(self.site_steam_load, years, 0, frequency)
            # TODO use a non-zero growth rate of steam load? --AE
            self.site_steam_load = Lib.drop_extra_data(self.site_steam_load, years)
        if self.site_hotwater_load is not None:
            self.site_hotwater_load = Lib.fill_extra_data(self.site_hotwater_load, years, 0, frequency)
            # TODO use a non-zero growth rate of hotwater load? --AE
            self.site_hotwater_load = Lib.drop_extra_data(self.site_hotwater_load, years)

    def initialize_variables(self, size):
        # rotating generation
        super().initialize_variables(size)
        # plus heat (steam and hotwater)
        self.variables_dict.update({
            'steam': cvx.Variable(shape=size, name=f'{self.name}-steamP', nonneg=True),
            'hotwater': cvx.Variable(shape=size, name=f'{self.name}-hotwaterP', nonneg=True),
        })

    def constraints(self, mask):
        constraint_list = super().constraints(mask)
        elec = self.variables_dict['elec']
        steam = self.variables_dict['steam']
        hotwater = self.variables_dict['hotwater']

        # to ensure that CHP never produces more steam than it can
        constraint_list += [cvx.NonPos(steam - self.max_steam_ratio * hotwater)]

        constraint_list += [cvx.Zero((steam + hotwater) * self.electric_heat_ratio - elec)]

        return constraint_list

    def get_steam_recovered(self, mask):
        # thermal power is recovered in a CHP plant whenever electric power is being generated
        # it is proportional to the electric power generated at a given time
        return self.variables_dict['steam']

    def get_hotwater_recovered(self, mask):
        # thermal power is recovered in a CHP plant whenever electric power is being generated
        # it is proportional to the electric power generated at a given time
        return self.variables_dict['hotwater']

    def timeseries_report(self):

        tech_id = self.unique_tech_id()
        results = super().timeseries_report()

        results[f'{tech_id} Steam Generation (kW)'] = self.variables_df['steam']
        results[f'{tech_id} Hot Water Generation (kW)'] = self.variables_df['hotwater']
        if self.site_steam_load is not None:
            results[f'{tech_id} Site Steam Thermal Load (BTU/hr)'] = self.site_steam_load
        if self.site_hotwater_load is not None:
            results[
                f'{tech_id} Site Hot Water Thermal Load (BTU/hr)'
            ] = self.site_hotwater_load

        return results

    def objective_function(self, mask, annuity_scalar=1):

        return super().objective_function(mask, annuity_scalar)
