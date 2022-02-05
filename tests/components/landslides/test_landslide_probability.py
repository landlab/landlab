"""
Unit tests for landlab.components.landslides.landslide_probability
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import LandslideProbability

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(ls_prob):
    """Testing if the name is right."""
    assert ls_prob.name == "Landslide Probability"


def test_input_var_names(ls_prob):
    """Testing if the input_var_names outputs the right list."""
    assert sorted(ls_prob.input_var_names) == [
        "soil__density",
        "soil__internal_friction_angle",
        "soil__maximum_total_cohesion",
        "soil__minimum_total_cohesion",
        "soil__mode_total_cohesion",
        "soil__thickness",
        "topographic__slope",
        "topographic__specific_contributing_area",
    ]


def test_output_var_names(ls_prob):
    """Testing if output_var_names outputs the right list."""
    assert sorted(ls_prob.output_var_names) == [
        "landslide__probability_of_failure",
        "soil__mean_relative_wetness",
        "soil__probability_of_saturation",
    ]


def test_var_units(ls_prob):
    """Testing if units are right."""
    assert set(ls_prob.input_var_names) | set(ls_prob.output_var_names), set(
        dict(ls_prob.units).keys()
    )

    assert ls_prob.var_units("topographic__specific_contributing_area") == "m"
    assert ls_prob.var_units("topographic__slope") == "tan theta"
    assert ls_prob.var_units("soil__transmissivity") == "m2/day"
    assert ls_prob.var_units("soil__saturated_hydraulic_conductivity") == "m/day"
    assert ls_prob.var_units("soil__mode_total_cohesion") == "Pa or kg/m-s2"
    assert ls_prob.var_units("soil__minimum_total_cohesion") == "Pa or kg/m-s2"
    assert ls_prob.var_units("soil__maximum_total_cohesion") == "Pa or kg/m-s2"
    assert ls_prob.var_units("soil__internal_friction_angle") == "degrees"
    assert ls_prob.var_units("soil__density") == "kg/m3"
    assert ls_prob.var_units("soil__thickness") == "m"
    assert ls_prob.var_units("soil__mean_relative_wetness") == "None"
    assert ls_prob.var_units("landslide__probability_of_failure") == "None"
    assert ls_prob.var_units("soil__probability_of_saturation") == "None"


def test_grid_shape(ls_prob):
    """Testing if the grid shape matches the inputs."""
    assert ls_prob.grid.number_of_node_rows == _SHAPE[0]
    assert ls_prob.grid.number_of_node_columns == _SHAPE[1]


def test_grid_x_extent(ls_prob):
    """Testing if x extent is right."""
    assert ls_prob.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


def test_grid_y_extent(ls_prob):
    """Testing if y extent is right."""
    assert ls_prob.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


def test_field_getters(ls_prob):
    """Testing if the right field is called."""
    for name in ls_prob.grid["node"]:
        field = ls_prob.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            ls_prob.grid.number_of_node_rows * ls_prob.grid.number_of_node_columns,
        )

    with pytest.raises(KeyError):
        ls_prob.grid["not_a_var_name"]


def test_field_initialized_to_zero(ls_prob):
    """Testing if the fields are initialized with zeros."""
    for name in ls_prob.grid["node"]:
        field = ls_prob.grid["node"][name]
        assert_array_almost_equal(field, np.zeros(ls_prob.grid.number_of_nodes))


class Test_init(object):

    def test_both_groundwater_and_saturated_thickness(self, example_raster_model_grid):
        """User provides saturated thickness but does not set groundwater__recharge_distribution 
        to None, check that LandslideProbability runs using groundwater__recharge_distribution """
        
        mg = example_raster_model_grid       
        ls_prob1 = LandslideProbability(mg,
                                        groundwater__recharge_distribution="uniform")
        ls_prob1.calculate_landslide_probability()
        LP_e = mg.at_node["landslide__probability_of_failure"].copy()
        ls_prob2 = LandslideProbability(mg,saturated__thickness_distribution="event")
        ls_prob2.calculate_landslide_probability()
        LP = mg.at_node["landslide__probability_of_failure"].copy()        
        np.testing.assert_almost_equal(LP,LP_e)

    def test_wrong_recharge(self, example_raster_model_grid):
        """wrong recharge method, should raise value error message"""        
        mg = example_raster_model_grid
        with pytest.raises(ValueError) as exc_info:
            ls_prob1 = LandslideProbability(mg,
                                            groundwater__recharge_distribution="lots_of_water")
        msg_e = "not a recharge distribution option"
        assert exc_info.match(msg_e)

    def test_wrong_saturated_thickness(self, example_raster_model_grid):
        """wrong saturated thickness method, should raise value error message"""        
        mg = example_raster_model_grid     
        with pytest.raises(ValueError) as exc_info:
            ls_prob1 = LandslideProbability(mg,
                                            groundwater__recharge_distribution=None,
                                            saturated__thickness_distribution="dry_as_a_bone")
        msg_e = "not a saturated thickness distribution option"
        assert exc_info.match(msg_e)
        
    def test_no_soil_transmissivity_or_conductivity(self, example_raster_model_grid):
        """testing exception is raised if raster model grid does not have a 
        soid__transmissivity or soil__saturated_hydraulic_conductivity field"""
        mg = example_raster_model_grid
        mg.delete_field(name = "soil__saturated_hydraulic_conductivity", loc = "node")
        mg.delete_field(name = "soil__transmissivity", loc = "node")
        with pytest.raises(ValueError) as exc_info: 
            LandslideProbability(mg)    
        msg_e = "no 'soil__transmissivity' or 'soil__saturated_hydraulic_conductivity' field" 
        assert exc_info.match(msg_e)



def test_calculate_landslide_probability_uniform_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'uniform' method.
    """
    grid_1 = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    gridnum = grid_1.number_of_nodes
    np.random.seed(seed=5)
    grid_1.at_node["topographic__slope"] = np.random.rand(gridnum)
    scatter_dat = np.random.randint(1, 10, gridnum)
    # grid_1.add_zeros("soil__saturated_hydraulic_conductivity", at="node")
    grid_1.at_node["topographic__specific_contributing_area"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid_1.at_node["soil__transmissivity"] = np.sort(
        np.random.randint(5, 20, gridnum).astype(float), -1
    )
    grid_1.at_node["soil__mode_total_cohesion"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid_1.at_node["soil__minimum_total_cohesion"] = (
        grid_1.at_node["soil__mode_total_cohesion"] - scatter_dat
    )
    grid_1.at_node["soil__maximum_total_cohesion"] = (
        grid_1.at_node["soil__mode_total_cohesion"] + scatter_dat
    )
    grid_1.at_node["soil__internal_friction_angle"] = np.sort(
        np.random.randint(26, 37, gridnum).astype(float)
    )
    grid_1.at_node["soil__thickness"] = np.sort(
        np.random.randint(1, 10, gridnum).astype(float)
    )
    grid_1.at_node["soil__density"] = 2000.0 * np.ones(gridnum)

    ls_prob_uniform = LandslideProbability(
        grid_1,
        number_of_iterations=10,
        groundwater__recharge_distribution="uniform",
        groundwater__recharge_min_value=20.0,
        groundwater__recharge_max_value=120.0,
        seed=5,
    )
    ls_prob_uniform.calculate_landslide_probability()
    np.testing.assert_almost_equal(
        grid_1.at_node["landslide__probability_of_failure"][5], 1.0
    )
    np.testing.assert_almost_equal(
        grid_1.at_node["landslide__probability_of_failure"][9], 0.0
    )


def test_calculate_landslide_probability_lognormal_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal' method.
    """
    grid_2 = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    gridnum = grid_2.number_of_nodes
    np.random.seed(seed=6)
    grid_2.at_node["topographic__slope"] = np.random.rand(gridnum)
    # grid_2.add_zeros("soil__saturated_hydraulic_conductivity", at="node")
    scatter_dat = np.random.randint(1, 10, gridnum).astype(float)
    grid_2.at_node["topographic__specific_contributing_area"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid_2.at_node["soil__transmissivity"] = np.sort(
        np.random.randint(5, 20, gridnum).astype(float), -1
    )
    grid_2.at_node["soil__mode_total_cohesion"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid_2.at_node["soil__minimum_total_cohesion"] = (
        grid_2.at_node["soil__mode_total_cohesion"] - scatter_dat
    )
    grid_2.at_node["soil__maximum_total_cohesion"] = (
        grid_2.at_node["soil__mode_total_cohesion"] + scatter_dat
    )
    grid_2.at_node["soil__internal_friction_angle"] = np.sort(
        np.random.randint(26, 37, gridnum).astype(float)
    )
    grid_2.at_node["soil__thickness"] = np.sort(
        np.random.randint(1, 10, gridnum).astype(float)
    )
    grid_2.at_node["soil__density"] = 2000.0 * np.ones(gridnum)

    ls_prob_lognormal = LandslideProbability(
        grid_2,
        number_of_iterations=10,
        groundwater__recharge_distribution="lognormal",
        groundwater__recharge_mean=5.0,
        groundwater__recharge_standard_deviation=0.25,
        seed=6,
    )
    ls_prob_lognormal.calculate_landslide_probability()
    np.testing.assert_almost_equal(
        grid_2.at_node["landslide__probability_of_failure"][5], 0.8
    )
    np.testing.assert_almost_equal(
        grid_2.at_node["landslide__probability_of_failure"][9], 0.4
    )


def test_calculate_landslide_probability_lognormal_spatial_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal_spatial' method.
    """
    grid_3 = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    gridnum = grid_3.number_of_nodes
    np.random.seed(seed=7)
    # grid_3.add_zeros("soil__saturated_hydraulic_conductivity", at="node")
    grid_3.at_node["topographic__slope"] = np.random.rand(gridnum)
    scatter_dat = np.random.randint(1, 10, gridnum).astype(float)
    grid_3.at_node["topographic__specific_contributing_area"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid_3.at_node["soil__transmissivity"] = np.sort(
        np.random.randint(5, 20, gridnum).astype(float), -1
    )
    grid_3.at_node["soil__mode_total_cohesion"] = np.sort(
        np.random.randint(30, 900, gridnum).astype(float)
    )
    grid_3.at_node["soil__minimum_total_cohesion"] = (
        grid_3.at_node["soil__mode_total_cohesion"] - scatter_dat
    )
    grid_3.at_node["soil__maximum_total_cohesion"] = (
        grid_3.at_node["soil__mode_total_cohesion"] + scatter_dat
    )
    grid_3.at_node["soil__internal_friction_angle"] = np.sort(
        np.random.randint(26, 37, gridnum).astype(float)
    )
    grid_3.at_node["soil__thickness"] = np.sort(
        np.random.randint(1, 10, gridnum).astype(float)
    )
    grid_3.at_node["soil__density"] = 2000.0 * np.ones(gridnum)

    ls_prob_lognormal_spatial = LandslideProbability(
        grid_3,
        number_of_iterations=10,
        groundwater__recharge_distribution="lognormal_spatial",
        groundwater__recharge_mean=np.random.randint(2, 7, gridnum).astype(float),
        groundwater__recharge_standard_deviation=np.random.rand(gridnum),
        seed=7,
    )
    ls_prob_lognormal_spatial.calculate_landslide_probability()
    np.testing.assert_almost_equal(
        grid_3.at_node["landslide__probability_of_failure"][5], 0.4
    )
    np.testing.assert_almost_equal(
        grid_3.at_node["landslide__probability_of_failure"][9], 0.29999999
    )



class Test_saturated_thickness_event(object):
    def test_calculate_landslide_probability_normal_1(self, example_raster_model_grid):
        """test probability is correct when using modeled saturated thickness
        for a specific event - normal value test"""
        
        gridnum = example_raster_model_grid.number_of_nodes
    
        ls_prob = LandslideProbability(
        example_raster_model_grid,
        number_of_iterations=10,
        groundwater__recharge_distribution = None,
        saturated__thickness_distribution = "event",
        seed=7,
        )
        
        ls_prob.calculate_landslide_probability()
        
        expected_probability = np.array([ 1,  1,  1,  1,  0,  0])
        np.testing.assert_almost_equal(
            example_raster_model_grid.at_node["landslide__probability_of_failure"]\
                [example_raster_model_grid.core_nodes], expected_probability
        )
         
    def test_calculate_landslide_probability_boundary_1(self, example_raster_model_grid):
        """test probability is correct when using modeled saturated thickness
        for a specific dry event, saturated__thickness = 0, - boundary value test"""
        
        gridnum = example_raster_model_grid.number_of_nodes
        example_raster_model_grid.at_node["saturated__thickness"] = np.zeros(gridnum)
        ls_prob = LandslideProbability(
        example_raster_model_grid,
        number_of_iterations=10,
        groundwater__recharge_distribution = None,
        saturated__thickness_distribution = "event",
        seed=7,
        )
        
        ls_prob.calculate_landslide_probability()    
        expected_probability = np.array([0.2, 0.1, 0, 0.4, 0, 0])
        np.testing.assert_almost_equal(
            example_raster_model_grid.at_node["landslide__probability_of_failure"]\
                [example_raster_model_grid.core_nodes], expected_probability
        )
            
    def test_instantiate_LP_bad_saturated_thickness_1(self, example_raster_model_grid):
        """test a ValueError exception is raised if saturated__thickness is negative
        - bad value test"""
        
        with pytest.raises(ValueError) as exc_info:
            
            example_raster_model_grid.at_node["saturated__thickness"] = \
            np.array([ 0.03,  0.03,  0.03, -0.11,  0.14,  0.06,  0.  ,  0.03, -0.34,
           -0.4 , -0.51, -0.43, -0.49, -0.56, -0.64, -0.73, -0.76, -0.82,
           -0.84, -0.93])
            gridnum = example_raster_model_grid.number_of_nodes        
            ls_prob = LandslideProbability(
            example_raster_model_grid,
            number_of_iterations=10,
            groundwater__recharge_distribution = None,
            saturated__thickness_distribution = "event",
            seed=7,
            )
            
        assert exc_info.match("saturated__thickness cannot be negative")
    
    
    def test_instantiate_LP_bad_saturated_thickness_2(self, example_raster_model_grid):
        """test a ValueError exception is raised when dtw field is has the wrong
        number of values, e.g. each node is a list of values rather than a single 
        float value- bad value test"""
        
        with pytest.raises(ValueError) as exc_info:
            
            gridnum = example_raster_model_grid.number_of_nodes
            example_raster_model_grid.at_node["saturated__thickness"] = np.ones((gridnum,2))       
            ls_prob = LandslideProbability(
            example_raster_model_grid,
            number_of_iterations=10,
            groundwater__recharge_distribution = None,
            saturated__thickness_distribution = "event",
            seed=7,
            )
            
        assert exc_info.match("saturated__thickness should be a 1-d array")

class Test_saturated_thickness_lognormal_spatial(object):
    def test_calculate_landslide_probability_normal_1(self, example_raster_model_grid):
        """test probability matches expected when relative wetness is randomly 
        determined from a lognormal distribution of saturated thickness - normal value test"""
        
        gridnum = example_raster_model_grid.number_of_nodes
    
        ls_prob = LandslideProbability(
        example_raster_model_grid,
        number_of_iterations=10,
        groundwater__recharge_distribution = None,
        saturated__thickness_distribution = "lognormal_spatial",
        saturated__thickness_mean = example_raster_model_grid.at_node["mean__saturated_thickness"],
        saturated__thickness_standard_deviation = example_raster_model_grid.at_node["stdev__saturated_thickness"],
        seed=7,
        )
        
        ls_prob.calculate_landslide_probability()
        
        expected_probability = np.array([1, 0.9, 0.8, 1, 0, 0])
        
        np.testing.assert_almost_equal(
            example_raster_model_grid.at_node["landslide__probability_of_failure"]\
                [example_raster_model_grid.core_nodes], expected_probability
        )
     
        
    def test_calculate_landslide_probability_boundary_1(self, example_raster_model_grid):
        """test probability matches expected when mean and standard deviation of the
        saturated thickness equal zero - boundary value test"""
        
        gridnum = example_raster_model_grid.number_of_nodes
        
        example_raster_model_grid.at_node["mean__saturated_thickness"] = np.zeros(gridnum)
        example_raster_model_grid.at_node["stdev__saturated_thickness"] = np.zeros(gridnum)
        
        ls_prob = LandslideProbability(
        example_raster_model_grid,
        number_of_iterations=10,
        groundwater__recharge_distribution = None,
        saturated__thickness_distribution = "lognormal_spatial",
        saturated__thickness_mean = example_raster_model_grid.at_node["mean__saturated_thickness"],
        saturated__thickness_standard_deviation = example_raster_model_grid.at_node["stdev__saturated_thickness"],
        seed=7,
        )
        
        ls_prob.calculate_landslide_probability()
        
        expected_probability = np.array([0.2, 0.1, 0, 0.4, 0, 0])
        
        np.testing.assert_almost_equal(
            example_raster_model_grid.at_node["landslide__probability_of_failure"]\
                [example_raster_model_grid.core_nodes], expected_probability
        )
            
    def test_calculate_landslide_probability_boundary_2(self, example_raster_model_grid):
        """test probability matches expected when mean saturated thickness is equal 
        to soil thickness and standard deviation is zero - boundary value test"""
        
        gridnum = example_raster_model_grid.number_of_nodes
        
        example_raster_model_grid.at_node["mean__saturated_thickness"] = \
            example_raster_model_grid.at_node["soil__thickness"]
        
        example_raster_model_grid.at_node["stdev__saturated_thickness"] = np.zeros(gridnum)
        
        ls_prob = LandslideProbability(
        example_raster_model_grid,
        number_of_iterations=10,
        groundwater__recharge_distribution = None,
        saturated__thickness_distribution = "lognormal_spatial",
        saturated__thickness_mean = example_raster_model_grid.at_node["mean__saturated_thickness"],
        saturated__thickness_standard_deviation = example_raster_model_grid.at_node["stdev__saturated_thickness"],
        seed=7,
        )
        
        ls_prob.calculate_landslide_probability()
        
        expected_probability = np.array([1, 1, 1, 1, 0, 0])
        
        np.testing.assert_almost_equal(
            example_raster_model_grid.at_node["landslide__probability_of_failure"]\
                [example_raster_model_grid.core_nodes], expected_probability
        )
            
    def test_instantiate_LP_bad_saturated_thickness_mean_1(self, example_raster_model_grid):
        """test a ValueError exception is raised if mean__saturated_thickness is negative
        - bad value test"""
        
        with pytest.raises(ValueError) as exc_info:
            
            example_raster_model_grid.at_node["mean__saturated_thickness"] = \
            np.array([ 0.03,  0.03,  0.03, -0.11,  0.14,  0.06,  0.  ,  0.03, -0.34,
           -0.4 , -0.51, -0.43, -0.49, -0.56, -0.64, -0.73, -0.76, -0.82,
           -0.84, -0.93])
            gridnum = example_raster_model_grid.number_of_nodes        
            ls_prob = LandslideProbability(
            example_raster_model_grid,
            number_of_iterations=10,
            groundwater__recharge_distribution = None,
            saturated__thickness_distribution = "lognormal_spatial",
            saturated__thickness_mean = example_raster_model_grid.at_node["mean__saturated_thickness"],
            saturated__thickness_standard_deviation = example_raster_model_grid.at_node["stdev__saturated_thickness"],
            seed=7,
            )
            
        assert exc_info.match("negative mean__saturated_thickness and/or stdev__saturated_thicknes")
                                  
    def test_instantiate_LP_bad_saturated_thickness_stdev_1(self, example_raster_model_grid):
        """test a ValueError exception is raised if stdev__saturated_thickness is negative
        - bad value test"""
        
        with pytest.raises(ValueError) as exc_info:
            
            example_raster_model_grid.at_node["stdev__saturated_thickness"] = \
            np.array([ 0.03,  0.03,  0.03, -0.11,  0.14,  0.06,  0.  ,  0.03, -0.34,
           -0.4 , -0.51, -0.43, -0.49, -0.56, -0.64, -0.73, -0.76, -0.82,
           -0.84, -0.93])
            gridnum = example_raster_model_grid.number_of_nodes        
            ls_prob = LandslideProbability(
            example_raster_model_grid,
            number_of_iterations=10,
            groundwater__recharge_distribution = None,
            saturated__thickness_distribution = "lognormal_spatial",
            saturated__thickness_mean = example_raster_model_grid.at_node["mean__saturated_thickness"],
            saturated__thickness_standard_deviation = example_raster_model_grid.at_node["stdev__saturated_thickness"],
            seed=7,
            )
            
        assert exc_info.match("negative mean__saturated_thickness and/or stdev__saturated_thicknes")
    
    def test_instantiate_LP_bad_saturated_thickness_mean_2(self, example_raster_model_grid):
        """test a ValueError exception is raised when mean__saturated_thickness field 
        has the wrong number of values, i.e. each node is a list of values rather 
        than a single value- bad value test"""
        
        with pytest.raises(ValueError) as exc_info:
    
            gridnum = example_raster_model_grid.number_of_nodes        
            example_raster_model_grid.at_node["mean__saturated_thickness"] = np.ones((gridnum,2))       
            ls_prob = LandslideProbability(
            example_raster_model_grid,
            number_of_iterations=10,
            groundwater__recharge_distribution = None,
            saturated__thickness_distribution = "lognormal_spatial",
            saturated__thickness_mean = example_raster_model_grid.at_node["mean__saturated_thickness"],
            saturated__thickness_standard_deviation = example_raster_model_grid.at_node["stdev__saturated_thickness"],
            seed=7,
            )
            
        assert exc_info.match("mean__saturated_thickness and/or stdev__saturated_thickness not a 1-d array")
                        
    def test_instantiate_LP_bad_saturated_thickness_stdev_2(self, example_raster_model_grid):
        """test a ValueError exception is raised when stdev__saturated_thickness field 
        has the wrong number of values, i.e. each node is a list of values rather 
        than a single value- bad value test"""
        
        with pytest.raises(ValueError) as exc_info:
            
            gridnum = example_raster_model_grid.number_of_nodes
            example_raster_model_grid.at_node["stdev__saturated_thickness"] = np.ones((gridnum,2))       
            ls_prob = LandslideProbability(
            example_raster_model_grid,
            number_of_iterations=10,
            groundwater__recharge_distribution = None,
            saturated__thickness_distribution = "lognormal_spatial",
            saturated__thickness_mean = example_raster_model_grid.at_node["mean__saturated_thickness"],
            saturated__thickness_standard_deviation = example_raster_model_grid.at_node["stdev__saturated_thickness"],
            seed=7,
            )
            
        assert exc_info.match("mean__saturated_thickness and/or stdev__saturated_thickness not a 1-d array")