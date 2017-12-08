# -*- coding: utf-8 -*-
"""
Unit tests for landlab.components.overland_flow.KinwaveOverlandFlowModel

last updated: 3/14/16
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.overland_flow import KinwaveOverlandFlowModel

(_SHAPE, _SPACING, _ORIGIN) = ((10, 10), (25, 25), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((10, 10), spacing=25)
    grid.add_zeros('node', 'topographic__elevation', dtype=float)
    grid.add_zeros('node', 'topographic__gradient')

    globals().update({
        'KinWaveOF': KinwaveOverlandFlowModel(grid)})


@with_setup(setup_grid)
def test_KinWaveOF_name():
    assert_equal(KinWaveOF.name, 'KinwaveOverlandFlowModel')


@with_setup(setup_grid)
def test_KinWaveOF_input_var_names():
    assert_equal(KinWaveOF.input_var_names,  ('topographic__elevation',
                                          'topographic__gradient',))


@with_setup(setup_grid)
def test_KinWaveOF_output_var_names():
    assert_equal(KinWaveOF.output_var_names, ('surface_water__depth',
                                          'water__velocity',
                                          'water__specific_discharge'))

@with_setup(setup_grid)
def test_KinWaveOF_var_units():
    assert_equal(set(KinWaveOF.input_var_names) |
                 set(KinWaveOF.output_var_names),
                 set(dict(KinWaveOF.units).keys()))

    assert_equal(KinWaveOF.var_units('topographic__elevation'), 'm')
    assert_equal(KinWaveOF.var_units('topographic__gradient'), 'm/m')
    assert_equal(KinWaveOF.var_units('surface_water__depth'), 'm')
    assert_equal(KinWaveOF.var_units('water__velocity'), 'm/s')
    assert_equal(KinWaveOF.var_units('water__specific_discharge'), 'm2/s')

@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(KinWaveOF.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(KinWaveOF.grid.number_of_node_columns, _SHAPE[1])


# @with_setup(setup_grid)
# def test_calc_soil_pressure():
#     np.testing.assert_almost_equal(KinWaveOF.calc_soil_pressure('KinWaveOFlt loam'),
#                             0.1647870740305523, decimal=6)
# @with_setup(setup_grid)
# def test_calc_soil_head():
#     soil_props = SoilInfiltrationGreenAmpt.SOIL_PROPS['loam']
#     np.testing.assert_almost_equal(KinWaveOF.calc_pressure_head(soil_props[0],
#             soil_props[1]), 0.087498292, decimal=6)
#
# @with_setup(setup_grid)
# def test_calc_moisture_deficit():
#     np.testing.assert_almost_equal(KinWaveOF.calc_moisture_deficit(
#         soil_bulk_denKinWaveOFty=1700., rock_denKinWaveOFty=2650.,
#         volume_fraction_coarse_fragments=0.,
#         soil_moisture_content=0.2), 0.15849056603, decimal=6)
#
# @with_setup(setup_grid)
# def test_run_one_step():
#     pass
#
# def test_negative_infiltration():
#     pass
