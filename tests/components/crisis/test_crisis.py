"""Test the CRISIS component.

"""

import io

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import CrisisComponent
from landlab.components import FlowAccumulator
from landlab.components.crisis.crisis_landlab_component import _node_to_2d
from landlab.components.crisis.crisis_landlab_component import check_params

GAMMA_W = 9.81


def _dem_header(nrows=3, ncols=100, cellsize=1.0, corner=False, x0=0.0, y0=0.0):
    """Return an in-memory ASCII DEM header (a new stream every call)."""
    kind = "corner" if corner else "center"
    return io.StringIO(
        f"ncols {ncols}\nnrows {nrows}\n"
        f"xll{kind} {x0}\nyll{kind} {y0}\ncellsize {cellsize}\n"
    )


def _sloped_grid(nrows=3, ncols=100, rise=0.9, dx=1.0):
    """A uniformly inclined slope draining to one outlet on the left edge."""
    grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
    grid.add_field("topographic__elevation", rise * grid.x_of_node, at="node")
    grid.add_field(
        "topographic__slope",
        np.degrees(grid.calc_slope_at_node("topographic__elevation")),
        at="node",
    )
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    grid.status_at_node[grid.nodes_at_left_edge[nrows // 2]] = (
        grid.BC_NODE_IS_FIXED_VALUE
    )
    FlowAccumulator(grid, runoff_rate=1).run_one_step()
    return grid


def _params(output_dir, **kwds):
    """constant-Ru parameter set."""
    params = {
        "dem_file": _dem_header(),
        "output_dir": str(output_dir),
        "Coordinate_Reference_System": "EPSG:32610",
        "Gamma": 22.0,
        "z_min": 0.125,
        "z_max": 1.875,
        "Depth_points": 8,
        "Min_LD_area": 12.0,
        "top_or_bottom_indicator": 1,
        "Triggering_Indicator": 1,
        "Hydrological_model_indicator": 1,
        "Ru": 0.8,
        "Strength_indicator": 1,
        "c": 20.0,
        "Phi": 37.0,
    }
    params.update(kwds)
    return params


@pytest.fixture
def grid():
    return _sloped_grid()



def test_node_to_2d():
    grid = RasterModelGrid((2, 3))
    grid.add_field("topographic__elevation", np.arange(6.0), at="node")

    assert_array_equal(
        _node_to_2d(grid, "topographic__elevation"),
        [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
    )


def test_check_params_complete(tmp_path):
    assert check_params(_params(tmp_path)) is None


@pytest.mark.parametrize("params", ([], "params", None, 1))
def test_check_params_not_a_dict(params):
    with pytest.raises(TypeError):
        check_params(params)


def test_check_params_ru_required(tmp_path):
    with pytest.raises(ValueError, match="Ru is missing"):
        check_params(_params(tmp_path, Ru=None))


def test_init_checks_params_first(tmp_path):
    """A bad params dict is reported."""
    grid = HexModelGrid((3, 3))

    with pytest.raises(ValueError, match="Missing or invalid entries"):
        CrisisComponent(grid, _params(tmp_path, Ru=None))


def test_init_requires_raster_grid(tmp_path):
    grid = HexModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")

    with pytest.raises(TypeError):
        CrisisComponent(grid, _params(tmp_path))


def test_init_depths_vector(grid, tmp_path):
    crisis = CrisisComponent(grid, _params(tmp_path))

    assert_array_almost_equal(crisis._Depths_vector, np.linspace(0.125, 1.875, 8))


def test_init_no_data_elevation(grid, tmp_path):
    grid.at_node["topographic__elevation"][0] = -9999.0
    crisis = CrisisComponent(grid, _params(tmp_path))

    assert np.isnan(crisis._Elev_Mat[0, 0])
    assert np.count_nonzero(np.isnan(crisis._Elev_Mat)) == 1


@pytest.mark.parametrize("ru", (0.0, 0.4, 0.8))
def test_constant_ru_pressure_head_field(grid, tmp_path, ru):
    crisis = CrisisComponent(grid, _params(tmp_path, Ru=ru))
    head = grid.at_node["soil_water__pressure_head"]

    assert head.shape == (grid.number_of_nodes, 8)
    expected = ru * 22.0 / GAMMA_W * crisis._Depths_vector
    assert_array_almost_equal(head, np.broadcast_to(expected, head.shape))


def test_constant_ru_calc_hydrology(grid, tmp_path):
    params = _params(tmp_path)
    crisis = CrisisComponent(grid, params)

    theta_3d, pressures_3d = crisis._calc_hydrology(0, params)

    assert theta_3d is None
    assert pressures_3d is None


@pytest.mark.parametrize("theta", (2, 3))
def test_constant_ru_ignores_theta(grid, tmp_path, theta):
    messages = []
    crisis = CrisisComponent(
        grid,
        _params(tmp_path, Strength_variation_with_theta=theta),
        log_callback=messages.append,
    )

    assert crisis._Strength_variation_with_theta == 1
    assert "soil_water__volume_fraction" not in grid.at_node
    assert len(messages) == 1
    assert "is ignored" in messages[0]