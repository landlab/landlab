"""
Unit tests for landlab.components.landslides.landslide_probability
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components.landslides import LandslideProbability

relative_wetness = 0.75


def unitgrid(dimensions):
    shape = (5, 4)
    spacing = (10e0, 10e0)
    coordinates = (0.0, 0.0)
    grid = RasterModelGrid(shape, spacing, coordinates)
    gridnum = grid.number_of_nodes

    assert gridnum == dimensions[0][0] * dimensions[0][1]

    return grid, gridnum


def calc_unit_default_value():

    unit_default_value = {}
    scatter = np.random.randint(
        1, 10
    )  # smallest min is the largest scatter, minus the mode cohesion
    unit_default_value["topographic__slope"] = (0.1, 0.8)
    unit_default_value["topographic__specific_contributing_area"] = (10, 100)
    unit_default_value["soil__transmissivity"] = (0.1, 10)
    unit_default_value["soil__mode_total_cohesion"] = (30, 900)
    mode = unit_default_value["soil__mode_total_cohesion"]
    # code used to run the model uses a random scatter <cohesion_scatter = np.random.randint(1,10)>
    scatter = np.random.randint(
        1, 10
    )  # smallest min is the largest scatter, minus the mode cohesion
    unit_default_value["soil__minimum_total_cohesion"] = (
        (mode[0] - scatter),
        mode[1] - scatter,
    )
    # large max is the largest scatter, plus the mode cohesion
    unit_default_value["soil__maximum_total_cohesion"] = (
        (mode[0] + scatter),
        mode[1] + scatter,
    )
    unit_default_value["soil__internal_friction_angle"] = (26, 37)
    unit_default_value["soil__thickness"] = (
        0.1,
        1,
    )  # functionally setting to uniform uniform 1 meter default; coded for the option of a range of values
    unit_default_value["soil__density"] = (
        1999,
        2001,
    )  # functionally setting to uniform 2000
    # Note on Ksat default range calculation
    # reverse order: (see above) lowest Transmissivity value are in upper right of synthethic grid,
    # not reversed: (see above) highest soil thickness in upper right of synthethic grid,
    # caution creating Ksat since parameterization of T and hs will have inverse relationship on Ksat
    # Kmin  = Tmin/hs_max   -slow- think of a deep thick clay layer
    kmin = (
        unit_default_value["soil__transmissivity"][0]
        / unit_default_value["soil__thickness"][1]
    )
    # Kmax = Tmax/hs_min    -fast- think of a sieve or a parking lot
    kmax = (
        unit_default_value["soil__transmissivity"][1]
        / unit_default_value["soil__thickness"][0]
    )
    unit_default_value["soil__saturated_hydraulic_conductivity"] = (kmin, kmax)

    # Output value range used in tests will change as a function of input default values
    unit_default_value["landslide__probability_of_failure"] = (0, 1)
    unit_default_value["soil__mean_relative_wetness"] = (0, 1)
    unit_default_value["soil__mean_watertable_depth"] = (0, 1000)
    unit_default_value["soil__mean_recharge"] = (0, 1000)
    unit_default_value["soil__probability_of_saturation"] = (0, 1)
    return unit_default_value


def get_default_nodevalues():
    """Print out the input range for
    all landslide component values (min/max).
    ----------
    Input: Landlab landslide model instance  """
    unit_default_value = calc_unit_default_value()

    return unit_default_value


def print_all_nodevalues(ls_prob, unit_default_value):
    """Print out the input range and resulting node values for 
    all landslide component values on the node with default parameters (min/max).
    ----------
    Input: Landlab landslide model instance  """

    unit_default_value = calc_unit_default_value()

    for name in ls_prob.grid["node"]:

        print_inputrange = print(
            "Input: (Min, Max) range of {value} default model inputs: {pars}".format(
                value=name, pars=unit_default_value[name]
            )
        )
        print_outputrange = print(
            "Ouput: {value} default value initialized for each node:".format(value=name)
        )
        field = ls_prob.grid["node"][name]
        print(field)
        field_min = ls_prob.grid["node"][name].min()
        field_max = ls_prob.grid["node"][name].max()
        unit_min = unit_default_value[name][0]
        unit_max = unit_default_value[name][1]

        assert field_min >= unit_min
        assert field_max <= unit_max


def print_list_nodevalues(ls_prob, value_list):
    """Print out the input range and resulting node values for 
    a list of landslide component values on the node.
    ----------
    Input: Landlab landslide model instance,   """

    unit_default_value = calc_unit_default_value()

    for name in value_list:

        print(
            "Input: (Min, Max) range of default {value} model inputs: {pars}".format(
                value=name, pars=unit_default_value[name]
            )
        )
        print(
            "Ouput: {value} default value initialized for each node:".format(value=name)
        )
        field = ls_prob.grid["node"][name]
        print(field)
        print("")


def print_failure_nodevalues(ls_prob):
    """Print out the input range and resulting node values for 
    a list of landslide component values on the node.
    ----------
    Input: Landlab landslide model instance,   """

    name = "landslide__probability_of_failure"

    print("Ouput Failure: {value} value calculated for each node:".format(value=name))
    field = ls_prob.grid["node"][name]
    print(field)
    print("")


def print_list_partest_nodevalues(ls_prob, value_list, unit_par_value):
    """Print out the input range and resulting node values for 
    a list of landslide component values on the node.  Call this to see results
    after running the function 
    ----------
    Input: Landlab landslide model instance,   """

    unit_default_value = calc_unit_default_value()

    for name in value_list:

        print(
            "Input: (Min, Max) range of default {value} model inputs: {pars}".format(
                value=name, pars=unit_default_value[name]
            )
        )
        print(
            "Input: (Min, Max) range of test parameter range for {value} model inputs: {pars}".format(
                value=name, pars=unit_par_value[name]
            )
        )
        print(
            "Ouput: {value} default value initialized for each node:".format(value=name)
        )
        field = ls_prob.grid["node"][name]
        print(field)
        print("")


def test_name(ls_prob):
    """Testing if the name is right.
    """
    assert ls_prob.name == "Landslide Probability"


def test_input_var_names(ls_prob):
    """Testing if the input_var_names outputs the right list.
    """
    assert sorted(ls_prob.input_var_names) == [
        "soil__density",
        "soil__internal_friction_angle",
        "soil__maximum_total_cohesion",
        "soil__minimum_total_cohesion",
        "soil__mode_total_cohesion",
        "soil__saturated_hydraulic_conductivity",
        "soil__thickness",
        "soil__transmissivity",
        "topographic__slope",
        "topographic__specific_contributing_area",
    ]


def test_output_var_names(ls_prob):
    """Testing if output_var_names outputs the right list.
    """
    assert sorted(ls_prob.output_var_names) == [
        "landslide__probability_of_failure",
        "soil__mean_recharge",
        "soil__mean_relative_wetness",
        "soil__mean_watertable_depth",
        "soil__probability_of_saturation",
    ]


def test_var_units(ls_prob):
    """Testing if units are right.
    """
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
    assert ls_prob.var_units("soil__mean_watertable_depth") == "m"
    assert ls_prob.var_units("soil__mean_recharge") == "mm/day"


def test_grid_shape(ls_prob):
    """Testing if the grid shape matches the inputs.
    """

    assert ls_prob.grid.number_of_node_rows == ls_prob.shape[0]
    assert ls_prob.grid.number_of_node_columns == ls_prob.shape[1]


def test_grid_x_extent(ls_prob):
    """Testing if x extent is right.
    """
    (shape, spacing, origin) = ((5, 4), (10e0, 10e0), (0.0, 0.0))

    assert ls_prob.grid.extent[1] == (shape[1] - 1) * spacing[1]


def test_grid_y_extent(ls_prob):
    """Testing if y extent is right.
    """
    (shape, spacing, origin) = ((5, 4), (10e0, 10e0), (0.0, 0.0))

    assert ls_prob.grid.extent[0] == (shape[0] - 1) * spacing[0]


def test_field_getters(ls_prob):
    """Testing if the right field is called.
    """
    for name in ls_prob.grid["node"]:
        field = ls_prob.grid["node"][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (
            ls_prob.grid.number_of_node_rows * ls_prob.grid.number_of_node_columns,
        )

    with pytest.raises(KeyError):
        ls_prob.grid["not_a_var_name"]


def test_field_initialized_to_zero(ls_prob):
    """Testing if the output fields are initialized with zeros.
    """
    for name in ls_prob.output_var_names:
        field = ls_prob.grid["node"][name]
        assert_array_almost_equal(field, np.zeros(ls_prob.grid.number_of_nodes))


def landslide_pars_ongrid(grid, gridnum, unit_default_value):

    parameter_range = unit_default_value

    scatter = np.random.randint(
        1, 10
    )  # smallest min is the largest scatter, minus the mode cohesion

    # Slope
    grid.at_node["topographic__slope"] = np.sort(
        np.random.uniform(
            parameter_range["topographic__slope"][0],
            parameter_range["topographic__slope"][1],
            gridnum,
        ).astype(float)
    )
    # reverse order: lowest slopes are in upper right of synthethic grid
    grid["node"]["topographic__slope"] = grid["node"]["topographic__slope"][::-1]

    # Area
    grid.at_node["topographic__specific_contributing_area"] = np.sort(
        np.random.randint(
            parameter_range["topographic__specific_contributing_area"][0],
            parameter_range["topographic__specific_contributing_area"][1],
            gridnum,
        ).astype(int)
        # no reverse order: highest contributing areas are in upper right of synthethic grid
    )

    # Transmissivity
    grid.at_node["soil__transmissivity"] = np.sort(
        np.random.uniform(
            parameter_range["soil__transmissivity"][0],
            parameter_range["soil__transmissivity"][1],
            gridnum,
        ).astype(float)
    )
    # reverse order: lowest Transmissivity value are in upper right of synthethic grid
    grid["node"]["soil__transmissivity"] = grid["node"]["soil__transmissivity"][::-1]

    # Cohesion
    grid.at_node["soil__mode_total_cohesion"] = np.random.randint(
        parameter_range["soil__mode_total_cohesion"][0],
        parameter_range["soil__mode_total_cohesion"][1],
        gridnum,
    ).astype(int)

    grid.at_node["soil__minimum_total_cohesion"] = (
        grid.at_node["soil__mode_total_cohesion"] - scatter
    ).astype(int)

    grid.at_node["soil__maximum_total_cohesion"] = (
        grid.at_node["soil__mode_total_cohesion"] + scatter
    ).astype(int)

    # Internal angle of friction
    grid.at_node["soil__internal_friction_angle"] = np.sort(
        np.random.randint(
            parameter_range["soil__internal_friction_angle"][0],
            parameter_range["soil__internal_friction_angle"][1],
            gridnum,
        ).astype(int)
    )

    # soil thickness
    grid.at_node["soil__thickness"] = np.sort(
        np.random.uniform(
            parameter_range["soil__thickness"][0],
            parameter_range["soil__thickness"][1],
            gridnum,
        ).astype(float)
        # no reverse order: greatest soil thickness is in upper right of synthethic grid
    )
    # soil density
    grid.at_node["soil__density"] = np.sort(
        np.random.randint(
            parameter_range["soil__density"][0],
            parameter_range["soil__density"][1],
            gridnum,
        ).astype(int)
        # no reverse order: greatest density is in upper right of synthethic grid
    )

    # K sat - #calculate T in driver with K * D; also usable as calibration multiplier value for T
    grid["node"]["soil__saturated_hydraulic_conductivity"] = (
        grid["node"]["soil__transmissivity"] / grid["node"]["soil__thickness"]
    )
    # reverse order: (see above) lowest Transmissivity value are in upper right of synthethic grid,
    # not reversed: (see above) highest soil thickness in upper right of synthethic grid,
    # caution creating Ksat since parameterization of T and hs will have inverse relationship on Ksat

    ls_grid = grid
    # ls_prob = LandslideProbability(grid)

    return ls_grid


def test_build_grid_unitarea():
    """ UNIT TEST. Testing the main method 'calculate_landslide_probability()' with default unit grid and default unit parameters 'uniform' method. No Input values for any input is required.  UNIT TEST.
    Input parameters:  None
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))

    grid, gridnum = unitgrid(dimensions)

    unit_default_value = calc_unit_default_value()

    ls_prob_defaultpar = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    return ls_prob_defaultpar


def test_build_unitgrid(shape, spacing):
    """ Customize UNIT TEST for bigger grids. Testing the main method 'calculate_landslide_probability()' with
    'uniform' method. Useful for testing run time and designing model experiment limits before GIS work.
    
    Input parameters:  dimensions of grid 
    """

    grid = RasterModelGrid(shape, spacing)
    gridnum = grid.number_of_nodes
    unit_default_value = calc_unit_default_value()
    ls_prob = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    return ls_prob


def test_field_initialized_to_range(ls_prob, unit_default_value):
    """Testing if the output fields are initialized within range of default.
    """
    for name in ls_prob.input_var_names:

        field_min = ls_prob.grid["node"][name].min()
        field_max = ls_prob.grid["node"][name].max()
        print(field_max)
        unit_min = unit_default_value[name][0]
        unit_max = unit_default_value[name][1]

        if name == "soil__maximum_total_cohesion":
            assert field_min >= unit_min
            assert (
                field_max <= unit_max + 10
            )  # set to max scatter applied to max cohesion
        else:
            assert field_min >= unit_min
            assert field_max <= unit_max


def test_testpars_unitlandslide(ls_prob):
    """Testing the main method 'calculate_landslide_probability()'. 
    To make values uniform, set min/max with a very small (0.001 difference)     
    Values are sorted to similate 'Upper Right' of the grid as a drainage outlet.
    
    Input parameters: LandslideProbability grid, dictionary of custom input parameter tuples
    """
    gridnum = ls_prob.grid.number_of_nodes
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        ls_prob.grid, gridnum, unit_default_value
    )  # see default values given in dictionary above
    ls_prob = LandslideProbability(ls_grid)

    # test_field_initialized_to_range() function to check output grid matches range given
    for name in ls_prob.input_var_names:

        field_min = ls_prob.grid["node"][name].min()
        field_max = ls_prob.grid["node"][name].max()

        unit_min = unit_default_value[name][0]
        unit_max = unit_default_value[name][1]

        if name == "soil__maximum_total_cohesion":
            assert field_min >= unit_min
            assert (
                field_max <= unit_max + 10
            )  # set to max scatter applied to max cohesion
        else:
            assert field_min >= unit_min
            assert field_max <= unit_max

    return ls_prob


# Recharge and depth unit values for a given relative wetness


def scenario_unit_explorer(rw, grid):
    theta_unit = grid.at_node["topographic__slope"]
    T_unit = grid.at_node["soil__transmissivity"]
    a_unit = grid["node"]["topographic__specific_contributing_area"]
    hs = grid["node"]["soil__thickness"]

    Recharge = ((rw * (T_unit * theta_unit)) / a_unit) * 1000  # mm/day

    Remin_value = Recharge.min()
    Remean = Recharge.mean()
    Restandard_deviation = Recharge.std()
    Remax_value = Recharge.max()
    rw_r = Recharge * a_unit / (T_unit * theta_unit)

    Default_R = [Remin_value, Remax_value, Remean, Restandard_deviation]

    De_sat_threshold = 0.001
    Depth = hs - rw * hs
    rw_d = (hs - Depth) / (hs - De_sat_threshold)
    hw = hs - Depth

    Demin_value = Depth.min()  # assumes unit test with soil thickness 1 m
    Demean = Depth.mean()
    Destandard_deviation = Depth.std()
    Demax_value = Depth.max()

    Default_D = [
        Demin_value.min(),
        Demax_value.max(),
        Demean.mean(),
        Destandard_deviation.mean(),
    ]

    return Default_R, Default_D


def print_testcore_nodevalues(ls_prob):
    name = "landslide__probability_of_failure"
    grid_test = ls_prob.grid
    print(
        "Unit test core node 6 = {value}".format(
            value=grid_test.at_node["landslide__probability_of_failure"][5]
        )
    )
    print(
        "Unit test core node 10 = {value}".format(
            value=grid_test.at_node["landslide__probability_of_failure"][9]
        )
    )


def test_calculate_landslide_probability_uniform_depth_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'uniform' method and depth to groundwater hydrologic forcings.
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__depth_distribution="uniform",
        groundwater__depth_min_value=Default_D[0],
        groundwater__depth_max_value=Default_D[1],
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 0.80
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.48
    )

    return ls_prob


def test_calculate_landslide_probability_uniform_recharge_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'uniform' method and recharge hydrologic forcing.
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__recharge_distribution="uniform",
        groundwater__recharge_min_value=Default_R[0],
        groundwater__recharge_max_value=Default_R[1],
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 0.96
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.84
    )

    return ls_prob


def test_calculate_landslide_probability_lognormal_depth_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal' method.
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__depth_distribution="lognormal",
        groundwater__depth_mean=Default_D[2],
        groundwater__depth_standard_deviation=Default_D[3],
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 0.88
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.32
    )

    return ls_prob


def test_calculate_landslide_probability_lognormal_recharge_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal' method.
    """

    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__recharge_distribution="lognormal",
        groundwater__recharge_mean=Default_R[2],
        groundwater__recharge_standard_deviation=Default_R[3],
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 1.0
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.80
    )

    return ls_prob


def test_calculate_landslide_probability_lognormal_spatial_depth_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal' method.
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    D_mu = np.random.uniform(Default_D[0], Default_D[1], gridnum)
    D_std = np.random.rand(gridnum)

    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__depth_distribution="lognormal_spatial",
        groundwater__depth_mean=D_mu,
        groundwater__depth_standard_deviation=D_std,
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 0.88
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.96
    )

    return ls_prob


def test_calculate_landslide_probability_lognormal_spatial_recharge_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal_spatial' method.
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    R_mu = np.random.uniform(Default_R[0], Default_R[1], gridnum)
    R_std = np.random.rand(gridnum)

    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__recharge_distribution="lognormal_spatial",
        groundwater__recharge_mean=R_mu,
        groundwater__recharge_standard_deviation=R_std,
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 1.0
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.24
    )

    return ls_prob


def test_calculate_landslide_probability_datadriven_spatial_depth_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal_spatial' method.
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    key = np.arange(gridnum)
    HSD_inputs_d = {}
    for i in range(gridnum):
        annualmaxDWT = np.random.uniform(Default_D[0], Default_D[1], 10)
        HSD_inputs_d[i] = annualmaxDWT

    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__depth_distribution="data_driven_spatial",
        groundwater__depth_HSD_inputs=HSD_inputs_d,
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 0.92
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.0
    )
    return ls_prob


def test_calculate_landslide_probability_datadriven_spatial_recharge_method():
    """Testing the main method 'calculate_landslide_probability()' with
    'lognormal_spatial' method.
    """
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    grid, gridnum = unitgrid(dimensions)
    unit_default_value = calc_unit_default_value()
    ls_grid = landslide_pars_ongrid(
        grid, gridnum, unit_default_value
    )  # see default values given in dictionary above

    Default_R, Default_D = scenario_unit_explorer(relative_wetness, ls_grid)
    ls_prob = {}
    HSD_dict = {}
    for vkey in range(2, 8):  # the first example vkey =2 is printed below
        HSD_dict[vkey] = np.random.uniform(
            Default_R[0], Default_R[1], 10
        )  # 20 - 120 mm/day
    HSD_id_dict = {}
    for ckey in grid.core_nodes:
        HSD_id_dict[ckey] = np.random.randint(2, 8, 2)
        fract_dict = {}
    for ckey in grid.core_nodes:
        fract_dict[ckey] = np.random.rand(2)

    HSD_inputs_r = [HSD_dict, HSD_id_dict, fract_dict]
    ls_prob = LandslideProbability(
        ls_grid,
        number_of_iterations=25,
        groundwater__recharge_distribution="data_driven_spatial",
        groundwater__recharge_HSD_inputs=HSD_inputs_r,
    )
    ls_prob.calculate_landslide_probability()

    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][5], 1.0
    )
    np.testing.assert_almost_equal(
        ls_prob.grid["node"]["landslide__probability_of_failure"][9], 0.84
    )
    return ls_prob
