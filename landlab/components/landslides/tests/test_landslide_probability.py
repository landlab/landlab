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
    """Setting up test raster grid.
    """
    grid = RasterModelGrid((20, 20), spacing=10e0)
    grid.at_node['topographic__slope'] = (
        np.zeros(grid.number_of_nodes, dtype=float))
    ls_prob = LandslideProbability(grid)
    globals().update({
        'ls_prob': LandslideProbability(grid)
    })


@with_setup(setup_grid)
def test_name():
    """Testing if the name is right.
    """
    assert_equal(ls_prob.name, 'Landslide Probability')


@with_setup(setup_grid)
def test_input_var_names():
    """Testing if the input_var_names outputs the right list.
    """
    assert_equal(sorted(ls_prob.input_var_names),
                 ['soil__density',
                  'soil__internal_friction_angle',
                  'soil__maximum_total_cohesion',
                  'soil__minimum_total_cohesion',
                  'soil__mode_total_cohesion',
                  'soil__saturated_hydraulic_conductivity',
                  'soil__thickness',
                  'soil__transmissivity',
                  'topographic__slope',
                  'topographic__specific_contributing_area'])


@with_setup(setup_grid)
def test_output_var_names():
    """Testing if output_var_names outputs the right list.
    """
    assert_equal(sorted(ls_prob.output_var_names),
                 ['landslide__probability_of_failure',
                  'soil__mean_relative_wetness',
                  'soil__probability_of_saturation'])


@with_setup(setup_grid)
def test_var_units():
    """Testing if units are right.
    """
    assert_equal(set(ls_prob.input_var_names) |
                 set(ls_prob.output_var_names),
                 set(dict(ls_prob.units).keys()))

    assert_equal(ls_prob.var_units(
        'topographic__specific_contributing_area'), 'm')
    assert_equal(ls_prob.var_units('topographic__slope'), 'tan theta')
    assert_equal(ls_prob.var_units('soil__transmissivity'), 'm2/day')
    assert_equal(ls_prob.var_units('soil__saturated_hydraulic_conductivity'),
                 'm/day')
    assert_equal(ls_prob.var_units(
        'soil__mode_total_cohesion'), 'Pa or kg/m-s2')
    assert_equal(ls_prob.var_units(
        'soil__minimum_total_cohesion'), 'Pa or kg/m-s2')
    assert_equal(ls_prob.var_units(
        'soil__maximum_total_cohesion'), 'Pa or kg/m-s2')
    assert_equal(ls_prob.var_units(
        'soil__internal_friction_angle'), 'degrees')
    assert_equal(ls_prob.var_units('soil__density'), 'kg/m3')
    assert_equal(ls_prob.var_units('soil__thickness'), 'm')
    assert_equal(ls_prob.var_units('soil__mean_relative_wetness'), 'None')
    assert_equal(ls_prob.var_units('landslide__probability_of_failure'),
                                   'None')
    assert_equal(ls_prob.var_units('soil__probability_of_saturation'),
                 'None')


@with_setup(setup_grid)
def test_grid_shape():
    """Testing if the grid shape matches the inputs.
    """
    assert_equal(ls_prob.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(ls_prob.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_grid_x_extent():
    """Testing if x extent is right.
    """
    assert_equal(ls_prob.grid.extent[1], (_SHAPE[1] - 1) * _SPACING[1])


@with_setup(setup_grid)
def test_grid_y_extent():
    """Testing if y extent is right.
    """
    assert_equal(ls_prob.grid.extent[0], (_SHAPE[0] - 1) * _SPACING[0])


@with_setup(setup_grid)
def test_field_getters():
    """Testing if the right field is called.
    """
    for name in ls_prob.grid['node']:
        field = ls_prob.grid['node'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (ls_prob.grid.number_of_node_rows *
                      ls_prob.grid.number_of_node_columns, ))

    assert_raises(KeyError, lambda: ls_prob.grid['not_a_var_name'])


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    """Testing if the fields are initialized with zeros.
    """
    for name in ls_prob.grid['node']:
        field = ls_prob.grid['node'][name]
        assert_array_almost_equal(field, np.zeros(
           ls_prob.grid.number_of_nodes))


def test_calculate_landslide_probability_uniform_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'uniform' method. 
    """
    grid_1 = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    gridnum = grid_1.number_of_nodes
    np.random.seed(seed=5)
    grid_1.at_node['topographic__slope'] = np.random.rand(gridnum)
    scatter_dat = np.random.randint(1, 10, gridnum)
    grid_1.at_node['topographic__specific_contributing_area']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid_1.at_node['soil__transmissivity']= (
             np.sort(np.random.randint(5, 20, gridnum), -1))
    grid_1.at_node['soil__mode_total_cohesion']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid_1.at_node['soil__minimum_total_cohesion']= (
             grid_1.at_node['soil__mode_total_cohesion'] - scatter_dat)
    grid_1.at_node['soil__maximum_total_cohesion']= (
             grid_1.at_node['soil__mode_total_cohesion'] + scatter_dat)
    grid_1.at_node['soil__internal_friction_angle']= (
             np.sort(np.random.randint(26, 37, gridnum)))
    grid_1.at_node['soil__thickness']= (
             np.sort(np.random.randint(1, 10, gridnum)))
    grid_1.at_node['soil__density']= (2000. * np.ones(gridnum))

    ls_prob_uniform = LandslideProbability(grid_1, number_of_iterations=10,
        groundwater__recharge_distribution='uniform',
        groundwater__recharge_min_value=20.,
        groundwater__recharge_max_value=120.,
        seed=5)
    ls_prob_uniform.calculate_landslide_probability()
    np.testing.assert_almost_equal(
        grid_1.at_node['landslide__probability_of_failure'][5], 1.)
    np.testing.assert_almost_equal(
        grid_1.at_node['landslide__probability_of_failure'][9], 0.)


def test_calculate_landslide_probability_lognormal_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal' method. 
    """
    grid_2 = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    gridnum = grid_2.number_of_nodes
    np.random.seed(seed=6)
    grid_2.at_node['topographic__slope'] = np.random.rand(gridnum)
    scatter_dat = np.random.randint(1, 10, gridnum)
    grid_2.at_node['topographic__specific_contributing_area']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid_2.at_node['soil__transmissivity']= (
             np.sort(np.random.randint(5, 20, gridnum), -1))
    grid_2.at_node['soil__mode_total_cohesion']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid_2.at_node['soil__minimum_total_cohesion']= (
             grid_2.at_node['soil__mode_total_cohesion'] - scatter_dat)
    grid_2.at_node['soil__maximum_total_cohesion']= (
             grid_2.at_node['soil__mode_total_cohesion'] + scatter_dat)
    grid_2.at_node['soil__internal_friction_angle']= (
             np.sort(np.random.randint(26, 37, gridnum)))
    grid_2.at_node['soil__thickness']= (
             np.sort(np.random.randint(1, 10, gridnum)))
    grid_2.at_node['soil__density']= (2000. * np.ones(gridnum))

    ls_prob_lognormal = LandslideProbability(grid_2, number_of_iterations=10,
        groundwater__recharge_distribution='lognormal',
        groundwater__recharge_mean=5.,
        groundwater__recharge_standard_deviation=0.25,
        seed=6)
    ls_prob_lognormal.calculate_landslide_probability()
    np.testing.assert_almost_equal(
        grid_2.at_node['landslide__probability_of_failure'][5], 0.8)
    np.testing.assert_almost_equal(
        grid_2.at_node['landslide__probability_of_failure'][9], 0.4)


def test_calculate_landslide_probability_lognormal_spatial_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal_spatial' method. 
    """
    grid_3 = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    gridnum = grid_3.number_of_nodes
    np.random.seed(seed=7)
    grid_3.at_node['topographic__slope'] = np.random.rand(gridnum)
    scatter_dat = np.random.randint(1, 10, gridnum)
    grid_3.at_node['topographic__specific_contributing_area']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid_3.at_node['soil__transmissivity']= (
             np.sort(np.random.randint(5, 20, gridnum), -1))
    grid_3.at_node['soil__mode_total_cohesion']= (
             np.sort(np.random.randint(30, 900, gridnum)))
    grid_3.at_node['soil__minimum_total_cohesion']= (
             grid_3.at_node['soil__mode_total_cohesion'] - scatter_dat)
    grid_3.at_node['soil__maximum_total_cohesion']= (
             grid_3.at_node['soil__mode_total_cohesion'] + scatter_dat)
    grid_3.at_node['soil__internal_friction_angle']= (
             np.sort(np.random.randint(26, 37, gridnum)))
    grid_3.at_node['soil__thickness']= (
             np.sort(np.random.randint(1, 10, gridnum)))
    grid_3.at_node['soil__density']= (2000. * np.ones(gridnum))

    ls_prob_lognormal_spatial = LandslideProbability(grid_3,
        number_of_iterations=10,
        groundwater__recharge_distribution='lognormal_spatial',
        groundwater__recharge_mean=np.random.randint(2,7, gridnum),
        groundwater__recharge_standard_deviation=np.random.rand(gridnum),
        seed=7)
    ls_prob_lognormal_spatial.calculate_landslide_probability()
    np.testing.assert_almost_equal(
        grid_3.at_node['landslide__probability_of_failure'][5], 0.4)
    np.testing.assert_almost_equal(
        grid_3.at_node['landslide__probability_of_failure'][9], 0.29999999)
