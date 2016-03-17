"""
This metacomponent combines a pit filling operation with a kinematic wave
shallow water flow routing algorithm and soil infiltration to produce a
simple bare-landscape hydrological model, after Rengers et al., in review,
and Julien et al., 1995.
"""

import numpy as np

from landlab.components import SinkFiller
from landlab.components import KinematicWave, SoilInfiltrationGreenAmpt
from ...utils.decorators import use_file_name_or_kwds


class FillInfiltrateKinematicWave(Component):
    """

    """

    _name = 'FillInfiltrateKinematicWave'

    _input_var_names = (
        'topographic__elevation',
        'surface_water__depth',
        'soil_water_infiltration__depth'
    )

    _output_var_names = (
        'topographic__elevation',
        'sediment_fill__depth',
        'surface_water__depth',
        'soil_water_infiltration__depth'
        'surface_water__discharge',
        'surface_water__velocity'
    )

    _var_units = {
        'topographic__elevation': 'm',
        'sediment_fill__depth': 'm',
        'surface_water__depth': 'm',
        'soil_water_infiltration__depth': 'm',
        'surface_water__discharge': 'm**3/s',
        'surface_water__velocity': 'm/s'
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'sediment_fill__depth': 'node',
        'surface_water__depth': 'node',
        'soil_water_infiltration__depth': 'node'
        'surface_water__discharge': 'node',
        'surface_water__velocity': 'node'
    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'sediment_fill__depth': 'Depth of sediment added at each' +
                                'node',
        'surface_water__depth': 'Depth of water above the surface',
        'soil_water_infiltration__depth': (
            'Water column height above the surface previously absorbed into ' +
            'the soil. Note that this is NOT the actual depth of the wetted ' +
            'front, which also depends on the porosity.'),
        'surface_water__discharge': ('Magnitude of discharge of water above ' +
                                     'the surface'),
        'surface_water__velocity': 'Speed of water flow above the surface'
    }

    def __init__(self, grid, hydraulic_conductivity=0.005,
                 soil_bulk_density=1590., rock_density=2650.,
                 initial_soil_moisture_content=0.15, soil_type='sandy loam',
                 volume_fraction_coarse_fragments=0.2,
                 surface_water_minimum_depth=1.e-8,
                 soil_pore_size_distribution_index=None,
                 soil_bubbling_pressure=None,
                 wetting_front_capillary_pressure_head=None,
                 mannings_n=0.03, critical_flow_depth=0.003,
                 mannings_epsilon=0.33333333, dt_max=0.3, max_courant=0.2,
                 **kwds):
        """
        """
        self._grid = grid
        self._sinkfill = SinkFiller(grid, routing='D4', apply_slope=True,
                                    fill_slope=1.e-6)
        spsdi = soil_pore_size_distribution_index
        wfcph = wetting_front_capillary_pressure_head
        self._infilt = SoilInfiltrationGreenAmpt(
            grid, hydraulic_conductivity=hydraulic_conductivity,
            soil_bulk_density=soil_bulk_density, rock_density=rock_density,
            initial_soil_moisture_content=initial_soil_moisture_content,
            soil_type=soil_type,
            volume_fraction_coarse_fragments=volume_fraction_coarse_fragments,
            surface_water_minimum_depth=surface_water_minimum_depth,
            soil_pore_size_distribution_index=spsdi,
            soil_bubbling_pressure=soil_bubbling_pressure,
            wetting_front_capillary_pressure_head=wfcph, **kwds)
        self._wave = KinematicWave(
            grid, mannings_n=mannings_n,
            critical_flow_depth=critical_flow_depth,
            mannings_epsilon=mannings_epsilon,
            dt_max=dt_max, max_courant=max_courant,
            min_surface_water_depth=surface_water_minimum_depth, **kwds)
