#!/usr/bin/env python

import numpy as np
from six.moves import range

from landlab import Component, CLOSED_BOUNDARY
from ...utils.decorators import use_file_name_or_kwds


class SoilInfiltrationGreenAmpt(Component):

    """Infiltrate surface water into a soil following the Green-Ampt method.

    This code is based on an overland flow model by Francis Rengers and
    colleagues, after Julien et al., 1995. The infiltration scheme follows the
    Green and Ampt equation.

    It was implemented in Landlab by DEJH, March '16. Please cite
    Rengers et al., 2016, Model Predictions of Water Runoff in Steep
    Catchments after Wildfire, WRR.

    Construction::

        SoilInfiltrationGreenAmpt(grid, hydraulic_conductivity=0.005,
                                  soil_bulk_density=1590, rock_density=2650,
                                  initial_soil_moisture_content=0.15,
                                  soil_type='sandy loam',
                                  volume_fraction_coarse_fragments=0.2,
                                  coarse_sed_flag=False,
                                  surface_water_minimum_depth=1.e-8,
                                  soil_pore_size_distribution_index=None,
                                  soil_bubbling_pressure=None,
                                  wetting_front_capillary_pressure_head=None)

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    hydraulic_conductivity : float, array, or field name (m/s)
        The soil effective hydraulic conductivity.
    soil_bulk_density : float (kg/m**3)
        The dry bulk density of the soil.
    rock_density : float (kg/m**3)
        The density of the soil constituent material (i.e., lacking porosity).
    initial_soil_moisture_content : float (m**3/m**3, 0. to 1.)
        The fraction of the initial pore space filled with water.
    soil_type : {'sand', loamy sand', 'sandy loam', 'loam', 'silt loam',
                 'sandy clay loam', 'clay loam', 'silty clay loam',
                 'sandy clay', 'silty clay', 'clay'}, or None
        A soil type to automatically set soil_pore_size_distribution_index
        and soil_bubbling_pressure, using mean values from Rawls et al.,
        1992.
    volume_fraction_coarse_fragments : float (m**3/m**3, 0. to 1.)
        The fraction of the soil made up of rocky fragments with very
        little porosity, with diameter > 2 mm.
    coarse_sed_flag : boolean, optional
        If this flag is set to true, the fraction of coarse material in the
        soil column with be used as a correction for phi, the porosity factor.
    surface_water_minimum_depth : float (m), optional
        A minimum water depth to stabilize the solutions for surface flood
        modelling. Leave as the default in most normal use cases.
    soil_pore_size_distribution_index : float, optional
        An index describing the distribution of pore sizes in the soil,
        and controlling effective hydraulic conductivity at varying water
        contents, following Brooks and Corey (1964). Can be set by
        soil_type. Typically denoted "lambda".
    soil_bubbling_pressure : float (m), optional
        The bubbling capillary pressure of the soil, controlling effective
        hydraulic conductivity at varying water contents, following Brooks
        and Corey (1964). Can be set by soil_type. Typically denoted "h_b".
    wetting_front_capillary_pressure_head : float (m), optional
        The effective head at the wetting front in the soil driven by
        capillary pressure in the soil pores. If not set, will be
        calculated by the component from the pore size distribution and
        bubbling pressure, following Brooks and Corey.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import SoilInfiltrationGreenAmpt
    >>> mg = RasterModelGrid((4,3), spacing=10.)
    >>> hydraulic_conductivity = mg.ones('node')*1.e-6
    >>> hydraulic_conductivity.reshape((4,3))[0:2,:] *= 10000.
    >>> h = mg.add_ones('node', 'surface_water__depth')
    >>> h *= 0.01
    >>> d = mg.add_ones('node', 'soil_water_infiltration__depth', dtype=float)
    >>> d *= 0.2
    >>> SI = SoilInfiltrationGreenAmpt(
    ...     mg,hydraulic_conductivity=hydraulic_conductivity)
    >>> for i in range(10):  # 100s total
    ...     SI.run_one_step(10.)
    >>> mg.at_node['surface_water__depth']
    array([  1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             1.00000000e-08,   1.00000000e-08,   1.00000000e-08,
             9.88530416e-03,   9.88530416e-03,   9.88530416e-03,
             9.88530416e-03,   9.88530416e-03,   9.88530416e-03])
    >>> mg.at_node['soil_water_infiltration__depth']
    array([ 0.20999999,  0.20999999,  0.20999999,  0.20999999,  0.20999999,
           0.20999999,  0.2001147 ,  0.2001147 ,  0.2001147 ,  0.2001147 ,
           0.2001147 ,  0.2001147 ])
    """

    _name = 'SoilInfiltrationGreenAmpt'

    _input_var_names = (
        'surface_water__depth',
        'soil_water_infiltration__depth'
        )

    _output_var_names = (
        'surface_water__depth',
        'soil_water_infiltration__depth',
        )

    _var_units = {
        'surface_water__depth': 'm',
        'soil_water_infiltration__depth': 'm'
        }

    _var_mapping = {
        'surface_water__depth': 'node',
        'soil_water_infiltration__depth': 'node'
        }

    _var_doc = {
        'surface_water__depth': 'Depth of water above the surface',
        'soil_water_infiltration__depth': (
            'Water column height above the surface previously absorbed into ' +
            'the soil. Note that this is NOT the actual depth of the wetted ' +
            'front, which also depends on the porosity.')
        }

    @use_file_name_or_kwds
    def __init__(self, grid, hydraulic_conductivity=0.005,
                 soil_bulk_density=1590., rock_density=2650.,
                 initial_soil_moisture_content=0.15, soil_type='sandy loam',
                 volume_fraction_coarse_fragments=0.2,
                 coarse_sed_flag=False,
                 surface_water_minimum_depth=1.e-8,
                 soil_pore_size_distribution_index=None,
                 soil_bubbling_pressure=None,
                 wetting_front_capillary_pressure_head=None, **kwds):
        """
        Initialize the kinematic wave approximation overland flow component.
        """

        self._grid = grid

        # This follows mean values from Rawls et al., 1992; lambda then h_b
        soil_props = {'sand': (0.694, 0.0726),
                      'loamy sand': (0.553, 0.0869),
                      'sandy loam': (0.378, 0.1466),
                      'loam': (0.252, 0.1115),
                      'silt loam': (0.234, 0.2076),
                      'sandy clay loam': (0.319, 0.2808),
                      'clay loam': (0.242, 0.2589),
                      'silty clay loam': (0.177, 0.3256),
                      'sandy clay': (0.223, 0.2917),
                      'silty clay': (0.150, 0.3419),
                      'clay': (0.165, 0.3730),
                      }

        # set up the soil properties. if a string is passed as a keyword
        # when the component is initialized, the component solves this 
        # automatically. if no soil_type is passed, lambda and h_b must be
        # set manually by the user using the soil_pore_size_distribution_index
        # (lambda) and soil_bubbling_pressure (h_b) keywords.
        assert (soil_type in soil_props.keys()) or (soil_type is None), \
            "Soil type must be a recognised soil type string, or None"
        if soil_pore_size_distribution_index is None:
            soil_pore_size_distribution_index = soil_props[soil_type][0]
        if soil_bubbling_pressure is None:
            soil_bubbling_pressure = soil_props[soil_type][1]
        if soil_type is None:
            assert (soil_pore_size_distribution_index is not None and
                    soil_bubbling_pressure is not None), \
                ("If you do not supply a soil_type, you must provide lambda " +
                 "and hb manually!")

        # sets hydraulic conductivity as a field or as a float.
        if type(hydraulic_conductivity) is str:
            self._Ks = self.grid.at_node[hydraulic_conductivity]
        else:
            self._Ks = hydraulic_conductivity
          
        # maintaining some layer of water on the surface at all times...
        self._lilwater = surface_water_minimum_depth
        
        # set the soil porosity. 
        phi = 1. - float(soil_bulk_density)/float(rock_density)
        assert 0. < phi < 1.
        
        # if the user wants to correct for coarse sediment, this flag
        # recalculate porosity
        if coarse_sed_flag == True:
            assert 0. <= volume_fraction_coarse_fragments < 1.
            coarse_frag_correction = 1. - volume_fraction_coarse_fragments
            phi = phi * coarse_frag_correction
        
        # setting the soil moisture deficit - making sure it is positive. 
        # code will crash if moisture deficit is <= 0.
        moisture_deficit = phi - initial_soil_moisture_content
        self._Md = moisture_deficit
        assert self._Md > 0., \
                    "Moisture deficit must be positive! Initial " \
                    "soil moisture content must be less than porosity."
        
        # Pressure head is set using lambda and h_b from above, using an
        # equation after Brooks-Corey, following Rawls et al., 1992
        if wetting_front_capillary_pressure_head is None:
            self._psi_f = ((2.+3.*soil_pore_size_distribution_index) /
                           (1.+3.*soil_pore_size_distribution_index) *
                           soil_bubbling_pressure/2.)
        else:
            self._psi_f = wetting_front_capillary_pressure_head

        # Setting water depth field, assuring that water depth is POSITIVE.
        self._water_depth = self.grid.at_node['surface_water__depth']
        assert np.sum(self._water_depth[
                          self.grid.core_nodes] < self._lilwater) == 0, \
            "Water depths must all be positive!"
        
        # Setting up the array of infiltration depths. 
        self._infiltration_depth = self.grid.at_node[
            'soil_water_infiltration__depth']
        # ^This guy is in "surface water equivalent", i.e., the actual depth
        # of penetration is greater because it only fills the pores

    def run_one_step(self, dt):
        """Update fields with current hydrologic conditions.

        Parameters
        ----------
        dt : float (s)
            The imposed timestep for the model.
        """
        # Actual infiltration capacity for a given time step.
        self.i_act = np.zeros(self.grid.number_of_nodes)
        
        # wetting front is a function of total infiltration depth and 
        # moisture deficit.
        self.wettingfront_depth = self._infiltration_depth/self._Md
        
        # Calculate infilitration capacity (m/s) if Ks is a float...
        try:
            self.infilt_cap = self._Ks * ((self.wettingfront_depth + 
                                          self._psi_f +
                                       self._water_depth) /
                                      self.wettingfront_depth)
        # Calculate infilitration capacity (m/s) if Ks is an array...    
        except ValueError:  
            self.infilt_cap = self._Ks*((self.wettingfront_depth+
                                            self._psi_f +
                                             self._water_depth) /
                                            self.wettingfront_depth)
            
        # Potential infiltration depth (m)
        self.potential_infilt = self.infilt_cap * dt  
        
        self.potential_infilt[np.where(self.potential_infilt < 0.0)] = 0.0

        # Where is water depth (per time, m/s) less than infil rate (m/s)?
        full_infil_indx = np.where(self._water_depth <= self.potential_infilt)
        
        # Where is water depth (per time, m/s) more than infil rate?
        notfull = np.where(self._water_depth > self.potential_infilt)
        
        # Set actual infiltration rate (m/s) to the calculated value (m/s)
        self.i_act = self.potential_infilt
        self.i_act[full_infil_indx] = (self._water_depth) - self._lilwater
        
        # Where water is completely infiltrated set to the minimum water value
        self._water_depth[full_infil_indx]= self._lilwater

        
        # Where water is not completely infiltrated, only subtract that 
        self._water_depth[notfull] = (self._water_depth[notfull] - 
                                                        (self.i_act[notfull]))
        
        # Assure water depths are all positive. Code will break if not!
        assert np.all(self._water_depth[
                          self.grid.core_nodes]) > 0, \
            "Water depths must all be positive!"
            
        # Update infiltration depth field. 
        self._infiltration_depth += self.i_act 
