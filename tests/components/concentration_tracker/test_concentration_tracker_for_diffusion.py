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

# Test concentration change for soil flux
mg = RasterModelGrid((3, 5), dx=1)
mg.axis_units = ('m', 'm')
mg.set_status_at_node_on_edges(right=4,
                               top=4,
                               left=4,
                               bottom=4)
mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE

mg.add_ones('soil__depth', at='node', units= ['m','m'])
mg.add_zeros('bedrock__elevation', at='node', units= ['m','m'])
mg.add_zeros('topographic__elevation', at='node', units= ['m','m'])
mg.at_node['topographic__elevation'][:] += mg.at_node['bedrock__elevation']
mg.at_node['topographic__elevation'][:] += mg.at_node['soil__depth']

# Colour band concentration fields
C_sed = mg.add_zeros('sediment_property__concentration', at='node', units= ['kg/m^3','kg/m^3'])
mg.at_node['sediment_property__concentration'][8] += 1
C_br = mg.add_zeros('bedrock_property__concentration', at='node', units= ['kg/m^3','kg/m^3'])
mg.at_node['bedrock_property__concentration'] += C_br_initial



# Test concentration change for bedrock weathering to soil
mg = RasterModelGrid((3, 5), dx=1)
mg.axis_units = ('m', 'm')
mg.set_status_at_node_on_edges(right=4,
                               top=4,
                               left=4,
                               bottom=4)
mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE

mg.add_ones('soil__depth', at='node', units= ['m','m'])
mg.add_zeros('bedrock__elevation', at='node', units= ['m','m'])
mg.add_zeros('topographic__elevation', at='node', units= ['m','m'])
mg.at_node['topographic__elevation'][:] += mg.at_node['bedrock__elevation']
mg.at_node['topographic__elevation'][:] += mg.at_node['soil__depth']

# Colour band concentration fields
C_sed = mg.add_zeros('sediment_property__concentration', at='node', units= ['kg/m^3','kg/m^3'])
mg.at_node['sediment_property__concentration'][8] += 1
C_br = mg.add_zeros('bedrock_property__concentration', at='node', units= ['kg/m^3','kg/m^3'])
mg.at_node['bedrock_property__concentration'] += C_br_initial

