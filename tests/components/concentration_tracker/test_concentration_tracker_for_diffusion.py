"""
Created on Wed Jul 12 12:25:27 2023

@author: LaurentRoberge
"""

import numpy as np
import pytest

from landlab import FieldError, RasterModelGrid
from landlab.components import ConcentrationTrackerForDiffusion

# %% Test input field errors


def test_input_soil_flux_from_diffuser():
    """
    ConcentrationTrackerForDiffusion should throw an error when output fields
    from a diffusion component do not exist (soil__flux)
    """
    # Make a raster model grid
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("soil_production__rate", at="node")
    mg.add_zeros("topographic__elevation", at="node")

    # Instantiate the component
    with pytest.raises(FieldError):
        ConcentrationTrackerForDiffusion(mg)


@pytest.mark.parametrize(
    "required_field", ["soil__depth", "soil_production__rate", "topographic__elevation"]
)
def test_input_fields_soil(required_field):
    """
    ConcentrationTrackerForDiffusion should throw an error when input fields
    are not provided (soil__depth, soil_production__rate, topographic__elevation)
    """
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")

    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("soil_production__rate", at="node")
    mg.add_zeros("topographic__elevation", at="node")

    mg.at_node.pop(required_field)
    with pytest.raises(FieldError):
        ConcentrationTrackerForDiffusion(mg)


# %% Test field instantiation


def test_field_instantiation():
    """
    ConcentrationTrackerForDiffusion should instantiate the following fields
    when they do not already exist ('bedrock_property__concentration' and
    'sediment_property__concentration')
    """
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("soil_production__rate", at="node")
    mg.add_zeros("topographic__elevation", at="node")

    ConcentrationTrackerForDiffusion(mg)

    missing_fields = {
        "bedrock_property__concentration",
        "sediment_property__concentration",
    } - set(mg.at_node)
    assert not missing_fields


# %% Test different user input options


# Test that default input produces correct fields with no pre-existing fields
def test_fields_for_default_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    ConcentrationTrackerForDiffusion(mg)

    node_fields = [
        mg.at_node["sediment_property__concentration"],
        mg.at_node["bedrock_property__concentration"],
    ]

    node_check = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    for node_field in node_fields:
        np.testing.assert_equal(node_field, node_check)


# Test that default input uses correct fields with pre-existing fields
def test_fields_for_default_input_with_preexisting_fields():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    mg.add_ones("sediment_property__concentration", at="node")
    mg.add_ones("bedrock_property__concentration", at="node")

    ConcentrationTrackerForDiffusion(mg)

    node_fields = [
        mg.at_node["sediment_property__concentration"],
        mg.at_node["bedrock_property__concentration"],
    ]

    node_check = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    for node_field in node_fields:
        np.testing.assert_equal(node_field, node_check)


# Test that user input of single values produces the correct fields
def test_fields_for_user_value_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    ConcentrationTrackerForDiffusion(
        mg,
        concentration_initial=1,
        concentration_in_bedrock=1,
    )

    node_fields = [
        mg.at_node["sediment_property__concentration"],
        mg.at_node["bedrock_property__concentration"],
    ]

    node_check = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    for node_field in node_fields:
        np.testing.assert_equal(node_field, node_check)


# Test that user input of arrays produces the correct fields
def test_fields_for_user_array_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    c_sed = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    c_br = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    ConcentrationTrackerForDiffusion(
        mg,
        concentration_initial=c_sed,
        concentration_in_bedrock=c_br,
    )

    node_fields = [
        mg.at_node["sediment_property__concentration"],
        mg.at_node["bedrock_property__concentration"],
    ]

    node_check = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    for node_field in node_fields:
        np.testing.assert_equal(node_field, node_check)


# Test that user input of grid fields produces the correct fields
def test_fields_for_user_field_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    c_sed = mg.add_ones("sediment_property__concentration", at="node")
    c_br = mg.add_ones("bedrock_property__concentration", at="node")

    ConcentrationTrackerForDiffusion(
        mg,
        concentration_initial=c_sed,
        concentration_in_bedrock=c_br,
    )

    node_fields = [
        mg.at_node["sediment_property__concentration"],
        mg.at_node["bedrock_property__concentration"],
    ]

    node_check = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    for node_field in node_fields:
        np.testing.assert_equal(node_field, node_check)


# %% Test concentration from weathering
# Test that default input produces C_w equal to C_br
def test_C_w_for_default_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    ct = ConcentrationTrackerForDiffusion(mg)

    C_w_check = ct.C_br

    np.testing.assert_equal(C_w_check, ct.C_w)


# Test that user input of a single value produces C_w array different to C_br
def test_C_w_for_user_value_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    ct = ConcentrationTrackerForDiffusion(
        mg,
        concentration_in_bedrock=0,
        concentration_from_weathering=1,
    )

    C_w_check = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    np.testing.assert_equal(C_w_check, ct.C_w)


# Test that user input of array produces C_w array different to C_br
def test_C_w_for_user_array_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    c_w = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    ct = ConcentrationTrackerForDiffusion(
        mg,
        concentration_in_bedrock=0,
        concentration_from_weathering=c_w,
    )

    C_w_check = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    np.testing.assert_equal(C_w_check, ct.C_w)


# Test that user input of grid fields produces C_w array different to C_br
def test_C_w_for_user_field_input():
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    c_w = mg.add_ones("weathering__concentration", at="node")

    ct = ConcentrationTrackerForDiffusion(
        mg,
        concentration_in_bedrock=0,
        concentration_from_weathering=c_w,
    )

    C_w_check = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    np.testing.assert_equal(C_w_check, ct.C_w)


# %% Test errors for physicaly impossible inputs
# Test that physically impossible inputs raise correct errors
def test_properties_concentrations():
    """
    ConcentrationTrackerForDiffusion should throw an error when input
    concentration values are negative.
    """
    mg = RasterModelGrid((3, 3))
    mg.add_zeros("soil__flux", at="link")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    # Instantiate the component
    with pytest.raises(ValueError):
        ConcentrationTrackerForDiffusion(mg, concentration_initial=-1)
    with pytest.raises(ValueError):
        ConcentrationTrackerForDiffusion(mg, concentration_in_bedrock=-1)
    with pytest.raises(ValueError):
        ConcentrationTrackerForDiffusion(mg, concentration_from_weathering=-1)


# %% Test against analytical solution


def test_concentration_from_soil_flux():
    """
    ConcentrationTrackerForDiffusion should correctly calculate concentration
    values based on known soil fluxes.
    """
    mg = RasterModelGrid((3, 5))
    mg.axis_units = ("m", "m")
    mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)
    mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE

    # Grid fields
    mg.add_ones("soil__depth", at="node")
    mg.add_zeros("bedrock__elevation", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.at_node["topographic__elevation"][:] += mg.at_node["bedrock__elevation"]
    mg.at_node["topographic__elevation"][:] += mg.at_node["soil__depth"]

    # Concentration field for soil
    mg.add_zeros("sediment_property__concentration", at="node")
    mg.at_node["sediment_property__concentration"][8] += 1

    # Forced soil flux and production rate fields
    mg.add_zeros("soil_production__rate", at="node")
    mg.add_zeros("soil__flux", at="link")
    # Soil flux for middle row of grid is negative 1
    middle_row_link_ids = [9, 10, 11, 12]
    mg.at_link["soil__flux"][middle_row_link_ids] -= 1

    # dx is 1 (by default) and soil depth is 1, so soil volume is 1 at each node.
    # Flux of -1 in the middle row should shift all C_sed values left by one node.

    ct = ConcentrationTrackerForDiffusion(mg)
    ct.run_one_step(1)

    C_sed_check = np.zeros(15)
    C_sed_check[7] = 1  # Node 7 should have all of the concentration from Node 8.

    np.testing.assert_equal(C_sed_check, mg.at_node["sediment_property__concentration"])


def test_concentration_from_weathering_without_C_w():
    """
    ConcentrationTrackerForDiffusion should correctly calculate concentration
    values based on known soil production rates from bedrock weathering when
    concentration_from_weathering input field is not defined.
    """
    mg = RasterModelGrid((3, 5))
    mg.axis_units = ("m", "m")
    mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)

    # Grid fields
    mg.add_ones("soil__depth", at="node")
    mg.add_zeros("bedrock__elevation", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.at_node["topographic__elevation"][:] += mg.at_node["bedrock__elevation"]
    mg.at_node["topographic__elevation"][:] += mg.at_node["soil__depth"]

    # Concentration fields
    mg.add_zeros("sediment_property__concentration", at="node")
    mg.at_node["sediment_property__concentration"][7] += 1
    mg.at_node["sediment_property__concentration"][8] += 1
    mg.add_zeros("bedrock_property__concentration", at="node")
    mg.at_node["bedrock_property__concentration"][7] += 1

    # Forced soil flux and production rate fields
    mg.add_zeros("soil__flux", at="link")
    # Soil production rate is 1
    mg.add_ones("soil_production__rate", at="node")

    ct = ConcentrationTrackerForDiffusion(mg)

    # Soil volume is 1 at each node. Soil production rate of 1 doubles volume.
    # This is normally done by the DepthDependentDiffuser. Here, it is forced.
    mg.at_node["soil__depth"] += 1

    # Node 7: C_sed remains 1 because parent bedrock had C_br of 1.
    # Node 8: C_sed is halved from 1 to 0.5 because parent bedrock had C_br = 0.

    ct.run_one_step(1)

    C_sed_check = np.zeros(15)
    C_sed_check[7] = 1  # Node 7 should have the same concentration as before.
    C_sed_check[8] = 0.5  # Node 8 should have half its previous concentration.

    np.testing.assert_equal(C_sed_check, mg.at_node["sediment_property__concentration"])


def test_concentration_from_weathering_with_C_w():
    """
    ConcentrationTrackerForDiffusion should correctly calculate concentration
    values based on known soil production rates from bedrock weathering when
    concentration_from_weathering input field is explicitly defined.
    """
    mg = RasterModelGrid((3, 5))
    mg.axis_units = ("m", "m")
    mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)

    # Grid fields
    mg.add_ones("soil__depth", at="node")
    mg.add_zeros("bedrock__elevation", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    mg.at_node["topographic__elevation"][:] += mg.at_node["bedrock__elevation"]
    mg.at_node["topographic__elevation"][:] += mg.at_node["soil__depth"]

    # Concentration fields
    mg.add_zeros("sediment_property__concentration", at="node")
    mg.at_node["sediment_property__concentration"][7] += 1
    mg.at_node["sediment_property__concentration"][8] += 1
    mg.add_zeros("bedrock_property__concentration", at="node")
    mg.at_node["bedrock_property__concentration"][7] += 1

    # Forced soil flux and production rate fields
    mg.add_zeros("soil__flux", at="link")
    # Soil production rate is 1
    mg.add_ones("soil_production__rate", at="node")

    ct = ConcentrationTrackerForDiffusion(mg, concentration_from_weathering=0)

    # Soil volume is 1 at each node. Soil production rate of 1 doubles volume.
    # This is normally done by the DepthDependentDiffuser. Here, it is forced.
    mg.at_node["soil__depth"] += 1

    # C_w overrides C_br values. In this case, no concentration is produced by
    # the weathering process, even at Node 7 where C_br = 1.

    # Node 7: C_sed is halved from 1 to 0.5 despite parent bedrock with C_br = 1.
    # Node 8: C_sed is halved from 1 to 0.5 because C_w = 0.

    ct.run_one_step(1)

    C_sed_check = np.zeros(15)
    C_sed_check[7] = 0.5  # Node 7 should have half its previous concentration.
    C_sed_check[8] = 0.5  # Node 8 should have half its previous concentration.

