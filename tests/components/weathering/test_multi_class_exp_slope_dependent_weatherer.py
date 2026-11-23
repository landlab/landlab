#!/usr/bin/env python3
"""

@author: YShmilovitz

Tests for the MultiClassExpSlopeDependentWeatherer component.
"""

import numpy as np
import pytest
from numpy import testing

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components.weathering import MultiClassExpSlopeDependentWeatherer

PHI = 0.4
RHO = 2650.0


def _make_grid(n_classes=2, mass_per_class=5.0, slope=0.0):
    mg = RasterModelGrid((3, 4), xy_spacing=10.0)
    mg.at_node["grains__mass"] = np.full(
        (mg.number_of_nodes, n_classes), mass_per_class
    )
    mg.at_node["bed_grains__proportions"] = np.full(
        (mg.number_of_nodes, n_classes), 1.0 / n_classes
    )
    soil = mg.add_field(
        "soil__depth",
        mg.at_node["grains__mass"].sum(axis=1) / (RHO * (1 - PHI)),
        at="node",
    )
    bedrock = mg.add_zeros("bedrock__elevation", at="node")
    mg.add_field("topographic__elevation", bedrock + soil, at="node")
    mg.add_zeros("topographic__steepest_slope", at="node")
    mg.at_node["topographic__steepest_slope"][mg.core_nodes] = slope
    return mg


def test_output_field_created():
    """
    The field soil_production__rate should be created by the component
    """
    mg = _make_grid()
    MultiClassExpSlopeDependentWeatherer(mg)

    assert "soil_production__rate" in mg.at_node


def test_missing_slope_field():
    """
    The component should raise a FieldError if the slope field does not exist
    """
    mg = _make_grid()
    with pytest.raises(FieldError):
        MultiClassExpSlopeDependentWeatherer(mg, slope_name="not_a_field")


def test_custom_slope_field_name():
    """
    The component should use a user-defined slope field
    """
    err_msg = "Custom slope field was not used"
    mg = _make_grid(slope=0.0)
    mg.add_zeros("slope", at="node")
    mg.at_node["slope"][mg.core_nodes] = 0.5
    p0, kappa = 1e-2, 1e-2
    h0 = mg.at_node["soil__depth"][mg.core_nodes].copy()
    MultiClassExpSlopeDependentWeatherer(
        mg,
        soil_production_maximum_rate=p0,
        slope_dependent_weathering_maximum_rate=kappa,
        slope_name="slope",
    ).run_one_step(dt=1.0)

    expected = (p0 + kappa * 0.5) * np.exp(-h0 / 1.0)
    testing.assert_allclose(
        mg.at_node["soil_production__rate"][mg.core_nodes], expected, err_msg=err_msg
    )


def test_slope_with_nan():
    """
    The weatherer should raise a ValueError if the slope has NaN at core nodes
    """
    mg = _make_grid()
    mg.at_node["topographic__steepest_slope"][mg.core_nodes[0]] = np.nan
    with pytest.raises(ValueError):
        MultiClassExpSlopeDependentWeatherer(mg)


@pytest.mark.parametrize("slope", [0.0, 0.1, 1.0])
def test_known_rate(slope):
    """
    The soil production rate should match (p0 + kappa * S) * exp(-h / h*)
    """
    err_msg = "Soil production rate does not match"
    p0, kappa, decay_depth = 1e-2, 1e-3, 0.5
    mg = _make_grid(slope=slope)
    h0 = mg.at_node["soil__depth"][mg.core_nodes].copy()

    MultiClassExpSlopeDependentWeatherer(
        mg,
        phi=PHI,
        rho=RHO,
        soil_production_maximum_rate=p0,
        slope_dependent_weathering_maximum_rate=kappa,
        soil_production_decay_depth=decay_depth,
    ).run_one_step(dt=1.0)

    expected = (p0 + kappa * slope) * np.exp(-h0 / decay_depth)
    testing.assert_allclose(
        mg.at_node["soil_production__rate"][mg.core_nodes], expected, err_msg=err_msg
    )


def test_depth_and_grains_mass_match():
    """
    Soil depth should still match the sum of grains mass
    """
    err_msg = "Soil depth does not match the sum of grains mass"
    mg = _make_grid(slope=0.1)
    MultiClassExpSlopeDependentWeatherer(mg, phi=PHI, rho=RHO).run_one_step(dt=10.0)

    mass = mg.at_node["grains__mass"][mg.core_nodes].sum(axis=1)
    testing.assert_allclose(
        mass / (RHO * (1 - PHI)),
        mg.at_node["soil__depth"][mg.core_nodes],
        err_msg=err_msg,
    )


def test_mass_split_by_bed_proportions():
    """
    New mass should be added to each grain class according to bed_grains__proportions
    """
    err_msg = "Weathered mass was not split according to bed proportions"
    mg = _make_grid(n_classes=3, mass_per_class=0.0)
    proportions = np.array([0.2, 0.3, 0.5])
    mg.at_node["bed_grains__proportions"][:] = proportions

    MultiClassExpSlopeDependentWeatherer(mg).run_one_step(dt=1.0)

    added = mg.at_node["grains__mass"][mg.core_nodes]
    testing.assert_allclose(
        added / added.sum(axis=1, keepdims=True),
        np.broadcast_to(proportions, added.shape),
        err_msg=err_msg,
    )
