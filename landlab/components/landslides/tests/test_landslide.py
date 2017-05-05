"""
Unit tests for landlab.components.landslides.landslide_probability
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
from numpy.testing import assert_array_almost_equal
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components import LandslideProbability


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    grid = RasterModelGrid((20, 20), spacing=10e0)
    grid.at_node['topographic__slope'] = (
        np.zeros(grid.number_of_nodes, dtype=float))
    LS_prob = LandslideProbability(grid)
    globals().update({
        'LS_prob': LandslideProbability(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert_equal(LS_prob.name, 'Landslide Probability')


@with_setup(setup_grid)
def test_input_var_names():
    assert_equal(sorted(LS_prob.input_var_names),
                 ['soil__density',
                  'soil__internal_friction_angle',
                  'soil__maximum_total_cohesion',
                  'soil__minimum_total_cohesion',
                  'soil__mode_total_cohesion',
                  'soil__thickness',
                  'soil__transmissivity',
                  'topographic__slope',
                  'topographic__specific_contributing_area'])


@with_setup(setup_grid)
def test_output_var_names():
    assert_equal(sorted(LS_prob.output_var_names),
                 ['landslide__probability_of_failure',
                  'soil__mean_relative_wetness'])


@with_setup(setup_grid)
def test_var_units():
    assert_equal(set(LS_prob.input_var_names) |
                 set(LS_prob.output_var_names),
                 set(dict(LS_prob.units).keys()))

    assert_equal(LS_prob.var_units(
        'topographic__specific_contributing_area'), 'm')
    assert_equal(LS_prob.var_units('topographic__slope'), 'tan theta')
    assert_equal(LS_prob.var_units('soil__transmissivity'), 'm2/day')
    assert_equal(LS_prob.var_units(
        'soil__mode_total_cohesion'), 'Pa or kg/m-s2')
    assert_equal(LS_prob.var_units(
        'soil__minimum_total_cohesion'), 'Pa or kg/m-s2')
    assert_equal(LS_prob.var_units(
        'soil__maximum_total_cohesion'), 'Pa or kg/m-s2')
    assert_equal(LS_prob.var_units(
        'soil__internal_friction_angle'), 'degrees')
    assert_equal(LS_prob.var_units('soil__density'), 'kg/m3')
    assert_equal(LS_prob.var_units('soil__thickness'), 'm')
    assert_equal(LS_prob.var_units('soil__mean_relative_wetness'), 'None')
    assert_equal(LS_prob.var_units('landslide__probability_of_failure'), 'None')


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(LS_prob.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(LS_prob.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_grid_x_extent():
    assert_equal(LS_prob.grid.extent[1], (_SHAPE[1] - 1) * _SPACING[1])


@with_setup(setup_grid)
def test_grid_y_extent():
    assert_equal(LS_prob.grid.extent[0], (_SHAPE[0] - 1) * _SPACING[0])


@with_setup(setup_grid)
def test_field_getters():
    for name in LS_prob.grid['node']:
        field = LS_prob.grid['node'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (LS_prob.grid.number_of_node_rows *
                      LS_prob.grid.number_of_node_columns, ))

    assert_raises(KeyError, lambda: LS_prob.grid['not_a_var_name'])


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in LS_prob.grid['node']:
        field = LS_prob.grid['node'][name]
        assert_array_almost_equal(field, np.zeros(
           LS_prob.grid.number_of_nodes))


def test_calculate_landslide_probability():
    grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    gridnum = grid.number_of_nodes
    np.random.seed(seed=0)
    grid.at_node['topographic__slope'] = np.random.rand(gridnum)
    scatter_dat = np.random.randint(1, 10, gridnum)
    grid.at_node['topographic__specific_contributing_area']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid.at_node['soil__transmissivity']= (
             np.sort(np.random.randint(5, 20, gridnum),-1))
    grid.at_node['soil__mode_total_cohesion']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid.at_node['soil__minimum_total_cohesion']= (
             grid.at_node['soil__mode_total_cohesion'] - scatter_dat)
    grid.at_node['soil__maximum_total_cohesion']= (
             grid.at_node['soil__mode_total_cohesion'] + scatter_dat)
    grid.at_node['soil__internal_friction_angle']= (
             np.sort(np.random.randint(26, 40, gridnum)))
    grid.at_node['soil__thickness']= (
             np.sort(np.random.randint(1, 10, gridnum)))
    grid.at_node['soil__density']= (
             2000. * np.ones(grid.number_of_nodes))

    LS_prob = LandslideProbability(grid, number_of_iterations=10,
        groundwater__recharge_distribution='uniform',
        groundwater__recharge_min_value=20.,
        groundwater__recharge_max_value=120.)
    LS_prob.calculate_landslide_probability()
    np.testing.assert_almost_equal(
        grid.at_node['landslide__probability_of_failure'][5], 1.)
    np.testing.assert_almost_equal(
        grid.at_node['landslide__probability_of_failure'][9], 0.6)
    