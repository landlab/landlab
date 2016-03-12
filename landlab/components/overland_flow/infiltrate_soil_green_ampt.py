#!/usr/bin/env python

import numpy as np

from landlab import Component, CLOSED_BOUNDARY
from ...utils.decorators import use_file_name_or_kwds


class SoilInfiltration(Component):

    """
    This code is based on an overland flow model by Francis Rengers and
    colleagues, after Julien et al., 1995. The infiltration scheme follows the
    Green and Ampt equation.
    It was implemented in Landlab by DEJH, March '16. Please cite
    Rengers et al., in review, Model Predictions of Water Runoff in Steep
    Catchments after Wildfire, WRR.

    Construction::

        SoilInfiltration(grid, hydraulic_conductivity=0.005,
                         soil_bulk_density=1590, rock_density=2650,
                         initial_soil_moisture_content=0.15,
                         soil_type='sandy loam',
                         volume_fraction_coarse_fragments=0.2,
                         surface_water_minimum_depth=1.e-8,
                         soil_pore_size_distribution_index=None,
                         soil_bubbling_pressure=None,
                         wetting_front_capillary_pressure_head=None)

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.

    Examples
    --------
    >>> from landlab import RasterModelGrid

    """

    _name = 'SoilInfiltration'

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
                 soil_bulk_density=1590, rock_density=2650,
                 initial_soil_moisture_content=0.15, soil_type='sandy loam',
                 volume_fraction_coarse_fragments=0.2,
                 surface_water_minimum_depth=1.e-8,
                 soil_pore_size_distribution_index=None,
                 soil_bubbling_pressure=None,
                 wetting_front_capillary_pressure_head=None, **kwds):
        """Initialize the kinematic wave approximation overland flow component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        hydraulic_conductivity : float, array, or field_name (m/s)
            Hydraulic conductivity of the soil.
        soil_bulk_density : float (kg/m**3)
            Density of the dry soil, including porosity.
        rock_density : float (kg/m**3)
            Density of the soil component particles.
        initial_soil_moisture_content : XXXX
            XXXX
        soil_type : {'sand', loamy sand', 'sandy loam', 'loam', 'silt loam',
                     'sandy clay loam', 'clay loam', 'silty clay loam',
                     'sandy clay', 'silty clay', 'clay'}, or None
            A soil type to automatically set soil_pore_size_distribution_index
            and soil_bubbling_pressure, using mean values from Rawls et al.,
            1992.
        volume_fraction_coarse_fragments : float, 0. to 1.
            The fraction of the soil made up of rocky fragments with very
            little porosity, with diameter > 2 mm.
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

        if type(hydraulic_conductivity) is str:
            self._Ks = mg.at_node[hydraulic_conductivity]
        else:
            self._Ks = hydraulic_conductivity
        self._lilwater = surface_water_minimum_depth
        # build the key model params from the inputs
        phi = 1. - soil_bulk_density/rock_density
        assert 0. < phi < 1.
        assert 0. <= volume_fraction_coarse_fragments < 1.
        coarse_frag_correction = 1. - volume_fraction_coarse_fragments
        self._phi_c = phi * coarse_frag_correction
        moisture_deficit = phi - initial_soil_moisture_content
        self._Md = moisture_deficit
        if wetting_front_capillary_pressure_head is None:
            self._psi_f = ((2.+3.*soil_pore_size_distribution_index) /
                           (1.+3.*soil_pore_size_distribution_index) *
                           soil_bubbling_pressure/2.)
            # equ. after Brooks-Corey, following Rawls et al., 1992
        else:
            self._psi_f = wetting_front_capillary_pressure_head

        self._water_depth = mg.at_node['surface_water__depth']
        assert np.sum(self._water_depth[
                          self.grid.core_nodes] < self._lilwater) == 0, \
            "Water depths must all be positive!"
        self._infiltration_depth = mg.at_node['soil_water_infiltration__depth']
        # ^This guy is in "surface water equivalent", i.e., the actual depth
        # of penetration is greater because it only fills the pores

    def update_one_timestep(self, dt):
        """Update fields with current hydrologic conditions.

        Parameters
        ----------

        """
        active = self.grid.status_at_node != CLOSED_BOUNDARY
        active_IDs = np.where(active)[0]
        wettingfront_depth = self._infiltration_depth[active]/self._Md
        try:
            infilt_cap = self._Ks * ((wettingfront_depth + self._psi_f +
                                      self._water_depth[active]) /
                                     wettingfront_depth)
        except ValueError:  # broadcasting error...
            infilt_cap = self._Ks[active]*((wettingfront_depth+self._psi_f +
                                            self._water_depth[active]) /
                                           wettingfront_depth)
        # where is depth > Ic?
        potential_infilt = infilt_cap*dt
        all_water_drawn = self._water_depth[active] <= potential_infilt
        not_all_water_drawn = np.logical_not(all_water_drawn)
        actual_infilt = potential_infilt
        actual_infilt[all_water_drawn] = self._water_depth[active][
                                            all_water_drawn]
        self._water_depth[active_IDs[all_water_drawn]] = self._lilwater
        self._water_depth[active_IDs[not_all_water_drawn]] -= actual_infilt[
            not_all_water_drawn]
        self._infiltration_depth[active] += actual_infilt
