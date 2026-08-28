# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 11:59:54 2026

@author: 12092
"""

import numpy as np
import pytest
from landlab import RasterModelGrid
from landlab.components import TaylorNonLinearDiffuser
from landlab.components import SoilDepthEvolver


class DummyDiffuser:
    """Simple diffuser used to test SoilDepthEvolver."""

    def __init__(self, grid, elevation_change=0.0):
        self.grid = grid
        self.elevation_change = elevation_change
        self.calls = 0

    def run_one_step(self, dt):
        self.grid.at_node["topographic__elevation"] += self.elevation_change
        self.calls += 1


@pytest.fixture
def grid():
    """Create a small Landlab grid for testing."""

    grid = RasterModelGrid((3, 3), xy_spacing=1.0)

    grid.add_zeros(
        "topographic__elevation",
        at="node",
    )

    soil_depth = grid.add_zeros(
        "soil__depth",
        at="node",
    )

    soil_depth[:] = 0.5

    return grid


def test_component_initializes(grid):
    """SoilDepthEvolver initializes with valid inputs."""

    diffuser = DummyDiffuser(grid)

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
    )

    assert component.current_time == 0.0


def test_negative_initial_soil_depth_raises(grid):
    """Negative initial soil depth should raise ValueError."""

    grid.at_node["soil__depth"][0] = -0.1

    diffuser = DummyDiffuser(grid)

    with pytest.raises(ValueError):
        SoilDepthEvolver(
            grid,
            diffuser=diffuser,
        )


def test_negative_soil_production_rate_raises(grid):
    """Negative soil-production rate should fail."""

    diffuser = DummyDiffuser(grid)

    with pytest.raises(ValueError):
        SoilDepthEvolver(
            grid,
            diffuser=diffuser,
            soil_production_rate=-0.0003,
        )


def test_zero_decay_depth_raises(grid):
    """Production decay depth must be positive."""

    diffuser = DummyDiffuser(grid)

    with pytest.raises(ValueError):
        SoilDepthEvolver(
            grid,
            diffuser=diffuser,
            soil_production_decay_depth=0.0,
        )


def test_zero_rock_density_raises(grid):
    """Rock density must be positive."""

    diffuser = DummyDiffuser(grid)

    with pytest.raises(ValueError):
        SoilDepthEvolver(
            grid,
            diffuser=diffuser,
            rock_density=0.0,
        )


def test_zero_soil_density_raises(grid):
    """Soil density must be positive."""

    diffuser = DummyDiffuser(grid)

    with pytest.raises(ValueError):
        SoilDepthEvolver(
            grid,
            diffuser=diffuser,
            soil_density=0.0,
        )


def test_diffuser_is_called(grid):
    """run_one_step should call the supplied diffuser once."""

    diffuser = DummyDiffuser(
        grid,
        elevation_change=0.01,
    )

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
    )

    component.run_one_step(10.0)

    assert diffuser.calls == 1


def test_current_time_advances(grid):
    """Current time should increase by dt."""

    diffuser = DummyDiffuser(grid)

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
    )

    component.run_one_step(25.0)

    assert component.current_time == 25.0

    component.run_one_step(25.0)

    assert component.current_time == 50.0


def test_soil_production_decreases_with_depth():
    """Deeper soil should have a lower production rate."""

    grid = RasterModelGrid((3, 3), xy_spacing=1.0)

    grid.add_zeros(
        "topographic__elevation",
        at="node",
    )

    soil_depth = grid.add_zeros(
        "soil__depth",
        at="node",
    )

    soil_depth[:] = np.array(
        [
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
        ]
    )

    diffuser = DummyDiffuser(grid)

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
    )

    result = component.run_one_step(1.0)

    production_rate = result["production_rate"]

    assert np.all(
        np.diff(production_rate) < 0.0
    )


def test_soil_production_equation(grid):
    """Production rate should match the exponential production law."""

    diffuser = DummyDiffuser(grid)

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
        soil_production_rate=0.0003,
        soil_production_decay_depth=0.5,
        rock_density=2000.0,
        soil_density=1600.0,
    )

    result = component.run_one_step(1.0)

    expected = (
        (2000.0 / 1600.0)
        * 0.0003
        * np.exp(-0.5 / 0.5)
    )

    np.testing.assert_allclose(
        result["production_rate"],
        expected,
    )


def test_positive_elevation_change_adds_soil(grid):
    """Deposition plus soil production should increase soil depth."""

    diffuser = DummyDiffuser(
        grid,
        elevation_change=0.1,
    )

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
        soil_production_rate=0.0,
    )

    result = component.run_one_step(1.0)

    np.testing.assert_allclose(
        result["elevation_change"],
        0.1,
    )

    np.testing.assert_allclose(
        result["soil_depth"],
        0.6,
    )

    np.testing.assert_allclose(
        result["soil_depth_change"],
        0.1,
    )


def test_negative_elevation_change_removes_soil(grid):
    """Surface lowering should reduce soil depth."""

    diffuser = DummyDiffuser(
        grid,
        elevation_change=-0.1,
    )

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
        soil_production_rate=0.0,
    )

    result = component.run_one_step(1.0)

    np.testing.assert_allclose(
        result["soil_depth"],
        0.4,
    )

    np.testing.assert_allclose(
        result["soil_depth_change"],
        -0.1,
    )


def test_erosion_cannot_make_soil_negative(grid):
    """Erosion exceeding available soil should not create negative depth."""

    diffuser = DummyDiffuser(
        grid,
        elevation_change=-1.0,
    )

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
        soil_production_rate=0.0,
    )

    result = component.run_one_step(1.0)

    np.testing.assert_allclose(
        result["soil_depth"],
        0.0,
    )


def test_exhausted_soil_retains_new_production(grid):
    """If old soil is exhausted, newly produced soil should remain."""

    diffuser = DummyDiffuser(
        grid,
        elevation_change=-1.0,
    )

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
        soil_production_rate=0.001,
        soil_production_decay_depth=0.5,
        rock_density=2000.0,
        soil_density=1600.0,
    )

    result = component.run_one_step(1.0)

    expected_production = (
        (2000.0 / 1600.0)
        * 0.001
        * np.exp(-0.5 / 0.5)
    )

    np.testing.assert_allclose(
        result["soil_depth"],
        expected_production,
    )


def test_run_one_step_returns_expected_keys(grid):
    """run_one_step should return all timestep diagnostics."""

    diffuser = DummyDiffuser(grid)

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
    )

    result = component.run_one_step(1.0)

    assert set(result) == {
        "soil_depth",
        "elevation_change",
        "soil_depth_change",
        "soil_produced",
        "production_rate",
    }


def test_two_timesteps_use_previous_soil_depth(grid):
    """Second timestep should use soil depth produced by the first."""

    diffuser = DummyDiffuser(
        grid,
        elevation_change=0.0,
    )

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
        soil_production_rate=0.001,
        soil_production_decay_depth=0.5,
    )

    first = component.run_one_step(1.0)

    depth_after_first = first[
        "soil_depth"
    ].copy()

    second = component.run_one_step(1.0)

    assert np.all(
        second["soil_depth"]
        > depth_after_first
    )

    assert np.all(
        second["production_rate"]
        < first["production_rate"]
    )

def test_with_taylor_nonlinear_diffuser():
    """SoilDepthEvolver should run with a real Landlab diffuser."""

    grid = RasterModelGrid((5, 5), xy_spacing=1.0)

    elevation = grid.add_zeros(
        "topographic__elevation",
        at="node",
    )

    # Give the grid a gentle slope so TNLD has something to act on.
    elevation[:] = grid.node_x * 0.1

    soil_depth = grid.add_zeros(
        "soil__depth",
        at="node",
    )

    soil_depth[:] = 0.5

    diffuser = TaylorNonLinearDiffuser(
        grid,
        linear_diffusivity=0.0042,
        slope_crit=1.25,
        nterms=2,
        dynamic_dt=True,
    )

    component = SoilDepthEvolver(
        grid,
        diffuser=diffuser,
        soil_production_rate=0.0003,
        soil_production_decay_depth=0.5,
        rock_density=2000.0,
        soil_density=1600.0,
    )

    result = component.run_one_step(1.0)

    assert component.current_time == 1.0

    assert result["soil_depth"].shape == (
        grid.number_of_nodes,
    )

    assert np.all(
        np.isfinite(result["soil_depth"])
    )

    assert np.all(
        result["soil_depth"] >= 0.0
    )