"""

@author: YuvalShmilovitz


Doc tests and unit tests the SoilGrading component.
"""

import numpy as np
import pytest
from numpy import testing
from landlab import RasterModelGrid
from landlab.components.soil_grading import SoilGrading


def test_outputFields_soil():
    """
    SoilGrading should create the soil__depth field
    """
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("bedrock__elevation", at="node")
    SoilGrading(mg)

    assert "soil__depth" in mg.at_node


def test_inputFields_soil():
    """
    SoilGrading should raise an error if soil__depth field is already exists in the grid
    """
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")

    with pytest.raises(ValueError, match="Soil field already exists"):
        SoilGrading(mg)


def test_outputFields_bedrock():
    """
    SoilGrading should create the topographic__elevation field when not provided
    """

    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    SoilGrading(mg)

    assert "bedrock__elevation" in mg.at_node


def test_inputFields_topographic_elevation():
    """
    SoilGrading should create the bedrock__elevation field when not provided
    """

    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    SoilGrading(mg)

    assert "bedrock__elevation" in mg.at_node


def test_depth_and_grains_mass_match():
    """
    Test if the soil depth match the sum of grains mass
    """
    err_msg = "Soil dz do not match the sum of grains mass"
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    soil_grading = SoilGrading(mg)

    soil_dz_at_node = mg.at_node["soil__depth"][9]
    mass_converted_to_dz_at_node = np.sum(mg.at_node["grains__weight"][9]) / (
        mg.dx * mg.dy * soil_grading._soil_density * (1 - soil_grading._phi)
    )  # Take into account the porosity

    testing.assert_almost_equal(
        soil_dz_at_node, mass_converted_to_dz_at_node, decimal=5, err_msg=err_msg
    )


def test_bedrock_grains_proportions():
    """
    Test if summed proportions of bedrock grain sizes are summed to 1
    """
    err_msg = "Bedrock grains proportions not okay"
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    SoilGrading(mg)

    summed_proportions_at_node = np.sum(mg.at_node["bed_grains__proportions"][9])
    testing.assert_almost_equal(
        summed_proportions_at_node, 1, decimal=5, err_msg=err_msg
    )


def test_match_to_user_defined_median_size():
    """
    Test if the median grain size at node match the user-defined median size
    """

    err_msg = "Wrong median size at node"
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    meansizes = [0.001, 0.01, 0.1]
    initial_median_size = 0.01
    SoilGrading(
        mg, meansizes=meansizes, initial_median_size=initial_median_size
    )

    median_size_at_node = mg.at_node["median_size__weight"][0]

    testing.assert_array_equal(
        median_size_at_node, initial_median_size, err_msg=err_msg
    )


def test_match_to_user_defined_meansizes():
    """
    Test if the meansizes of grain size fractions match the user-defined median size
    """

    err_msg = "Meansizes of grains fractions do not match user input"
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    meansizes = [0.001, 0.01, 0.1]
    initial_median_size = 0.01

    SoilGrading(
        mg, meansizes=meansizes, initial_median_size=initial_median_size
    )

    testing.assert_array_equal(
        mg.at_node["grains_fractions__size"][9], meansizes, err_msg=err_msg
    )


@pytest.mark.filterwarnings(
    "ignore:Median size provided is larger than the distribution upper bound"
)
def test_beyond_range_user_defined_soil_median_size_larger():
    """
    Test if in the case that the user-defined median size in soil is larger
    than the size of the largest grain size fraction mean size,
    all the mass will be store in the largest grain fraction
    """

    err_msg = "Wrong proportion of grains fractions in soil"
    n_rows = 5
    n_columns = 5
    spacing = 1
    meansizes = [0.001, 0.01, 0.1]
    initial_median_size_soil = np.inf

    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    SoilGrading(
        mg, meansizes=meansizes, initial_median_size=initial_median_size_soil
    )

    proportion_in_larger_grainsize_fraction = mg.at_node["grains__weight"][
        9, -1
    ] / np.sum(mg.at_node["grains__weight"][9])

    testing.assert_approx_equal(
        proportion_in_larger_grainsize_fraction, 1, err_msg=err_msg
    )


@pytest.mark.filterwarnings(
    "ignore:Median size provided is larger than the distribution upper bound"
)
def test_beyond_range_user_defined_bedrock_median_size_larger():
    """
    Test if in the case that the user-defined median size in bedrock is larger
    than the size of the larger grain size fraction mean size,
    all the mass will be store in the largest grain fraction
    """

    err_msg = "Wrong proportion of grains fractions in bedrock"
    n_rows = 5
    n_columns = 5
    spacing = 1
    meansizes = [0.001, 0.01, 0.1]
    initial_median_size_soil = 0.01

    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    median_size_bedrock = np.inf

    soil_grading = SoilGrading(
        mg, meansizes=meansizes, initial_median_size=initial_median_size_soil
    )
    soil_grading.create_dist(
        is_bedrock_distribution_flag=True, median_size=median_size_bedrock
    )
    proportion_of_mass_in_the_larger_size_fraction = mg.at_node[
        "bed_grains__proportions"
    ][9, -1]

    testing.assert_approx_equal(
        proportion_of_mass_in_the_larger_size_fraction, 1, err_msg=err_msg
    )


@pytest.mark.filterwarnings(
    "ignore:Median size provided is smaller than the distribution lower bound"
)
def test_beyond_range_user_defined_soil_median_size_smaller():
    """
    Test if in the case that the user-defined median size in soil is smaller
    than the smallest size of grain size fraction mean size,
    all the mass will be store in the smallest grain fraction
    """

    err_msg = "Wrong proportion of grains fractions in soil"
    n_rows = 5
    n_columns = 5
    spacing = 1
    meansizes = [0.001, 0.01, 0.1]
    initial_median_size_soil = -np.inf

    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    SoilGrading(
        mg, meansizes=meansizes, initial_median_size=initial_median_size_soil
    )

    proportion_in_smallest_grainsize_fraction = mg.at_node["grains__weight"][
        9, 0
    ] / np.sum(mg.at_node["grains__weight"][9])
    testing.assert_approx_equal(
        proportion_in_smallest_grainsize_fraction, 1, err_msg=err_msg
    )


@pytest.mark.filterwarnings(
    "ignore:Median size provided is smaller than the distribution lower bound"
)
def test_beyond_range_user_defined_bedrock_median_size_smaller():
    """
    Test if in the case that the user-defined median size in bedrock is smaller
    than the size of the smallest grain size fraction mean size,
    all the mass will be store in the smallest grain fraction
    """

    err_msg = "Wrong proportion of grains fractions in bedrock"
    n_rows = 5
    n_columns = 5
    spacing = 1
    meansizes = [0.001, 0.01, 0.1]
    initial_median_size_soil = 0.01

    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    median_size_bedrock = -np.inf

    soil_grading = SoilGrading(
        mg, meansizes=meansizes, initial_median_size=initial_median_size_soil
    )
    soil_grading.create_dist(
        is_bedrock_distribution_flag=True, median_size=median_size_bedrock
    )
    proportion_of_mass_in_the_larger_size_fraction = mg.at_node[
        "bed_grains__proportions"
    ][9, 0]

    testing.assert_approx_equal(
        proportion_of_mass_in_the_larger_size_fraction, 1, err_msg=err_msg
    )
