import numpy as np
import pytest
from numpy.testing import (
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_equal,
)

from landlab import HexModelGrid, RasterModelGrid
from landlab.components import AdvectionSolverTVD
from landlab.components.advection import (
    find_upwind_link_at_link,
    upwind_to_local_grad_ratio,
)
from landlab.grid.linkorientation import LinkOrientation


@pytest.mark.parametrize("u", [0.0, 1.0, 10.0, [1.0] * 17])
def test_upwind_link_at_link_raster_positive(u):
    grid = RasterModelGrid((3, 4))
    uwl = find_upwind_link_at_link(grid, u)

    assert_array_equal(
        uwl[grid.vertical_links].reshape((2, 4)),
        [
            [-1, -1, -1, -1],
            [3, 4, 5, 6],
        ],
    )
    assert_array_equal(
        uwl[grid.horizontal_links].reshape((3, 3)),
        [
            [-1, 0, 1],
            [-1, 7, 8],
            [-1, 14, 15],
        ],
    )


@pytest.mark.parametrize("u", [-1.0, -10.0, [-1.0] * 17])
def test_upwind_link_at_link_raster_negative(u):
    grid = RasterModelGrid((3, 4))
    uwl = find_upwind_link_at_link(grid, u)
    assert_array_equal(
        uwl[grid.vertical_links].reshape((2, 4)),
        [
            [10, 11, 12, 13],
            [-1, -1, -1, -1],
        ],
    )
    assert_array_equal(
        uwl[grid.horizontal_links].reshape((3, 3)),
        [
            [1, 2, -1],
            [8, 9, -1],
            [15, 16, -1],
        ],
    )


def test_upwind_link_at_link_raster_mixed():
    grid = RasterModelGrid((3, 4))
    u = np.zeros(grid.number_of_links)
    u[4:6] = -1
    u[7] = -1
    u[8:10] = 1
    u[11:13] = 1
    uwl = find_upwind_link_at_link(grid, u)
    assert_array_equal(
        uwl[grid.vertical_links].reshape((2, 4)),
        [
            [-1, 11, 12, -1],
            [3, 4, 5, 6],
        ],
    )
    assert_array_equal(
        uwl[grid.horizontal_links].reshape((3, 3)),
        [
            [-1, 0, 1],
            [8, 7, 8],
            [-1, 14, 15],
        ],
    )


def test_upwind_link_at_link_hex_horizontal():
    # Hex horizontal
    grid = HexModelGrid((3, 3), orientation="horizontal")
    uwl = find_upwind_link_at_link(grid, 1.0)

    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.E], [-1, 0, -1, 8, 9, -1, 17]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNE],
        [-1, -1, -1, -1, 3, 5],
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNW], [-1, -1, -1, 4, 6, -1]
    )

    uwl = find_upwind_link_at_link(grid, -1.0)
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.E],
        [1, -1, 9, 10, -1, 18, -1],
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNE],
        [13, 15, -1, -1, -1, -1],
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNW],
        [-1, 12, 14, -1, -1, -1],
    )

    u = np.zeros(grid.number_of_links)
    u[3:7] = -1
    u[8] = -1
    u[9:11] = 1
    u[12:16] = 1

    expected = np.choose(
        u >= 0, [find_upwind_link_at_link(grid, -1), find_upwind_link_at_link(grid, 1)]
    )
    assert_array_equal(find_upwind_link_at_link(grid, u), expected)


def test_upwind_link_at_link_hex_vertical():
    # Hex vertical
    grid = HexModelGrid((3, 3), orientation="vertical")
    uwl = find_upwind_link_at_link(grid, 1.0)

    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ENE], [-1, -1, 3, -1, 10, -1]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.N], [-1, -1, -1, 2, 5, 6, 9]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ESE], [-1, 7, -1, 14, -1, -1]
    )

    uwl = find_upwind_link_at_link(grid, -1.0)

    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ENE], [-1, 8, -1, 15, -1, -1]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.N], [9, 12, 13, 16, -1, -1, -1]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ESE], [-1, -1, 4, -1, 11, -1]
    )

    u = np.zeros(grid.number_of_links)
    u[2:4] = -1
    u[4] = 1
    u[7] = -1
    u[9] = 1
    u[10] = -1
    u[11] = 1
    u[14] = -1
    u[15] = 1

    expected = np.choose(
        u >= 0, [find_upwind_link_at_link(grid, -1), find_upwind_link_at_link(grid, 1)]
    )
    assert_array_equal(find_upwind_link_at_link(grid, u), expected)


@pytest.mark.parametrize("u", [-1.0, 0.0, 1.0])
def test_upwind_to_local_grad_ratio_zero_gradient(u):
    """Check that ratio is set to one when gradient is zero."""
    grid = HexModelGrid((4, 4))

    v = grid.ones(at="node")
    upwind_link_at_link = find_upwind_link_at_link(grid, u)
    r = upwind_to_local_grad_ratio(grid, v, upwind_link_at_link)

    assert_array_almost_equal(r, 1.0)


def test_upwind_to_local_grad_ratio_out_keyword():
    """Check that ratio is set to one when gradient is zero."""
    grid = HexModelGrid((4, 4))

    out = grid.empty(at="link")

    v = grid.ones(at="node")
    upwind_link_at_link = find_upwind_link_at_link(grid, -1.0)
    r = upwind_to_local_grad_ratio(grid, v, upwind_link_at_link, out=out)

    assert_array_almost_equal(out, 1.0)
    assert r is out


@pytest.mark.parametrize("cls", [RasterModelGrid, HexModelGrid])
@pytest.mark.parametrize("u", [-1.0, 0.0, 1.0])
def test_upwind_to_local_grad_ratio(cls, u):
    """
    Predicted upwind_to_local_grad_ratio for u>1
    """
    grid = cls((4, 4))
    v = np.arange(grid.number_of_nodes) ** 2

    upwind_link_at_link = find_upwind_link_at_link(grid, u)
    r = upwind_to_local_grad_ratio(grid, v, upwind_link_at_link)

    diff_at_link = grid.calc_grad_at_link(v)
    expected = diff_at_link[upwind_link_at_link] / diff_at_link

    assert_array_almost_equal(r[upwind_link_at_link == -1], 1.0)
    assert_array_almost_equal(
        r[upwind_link_at_link != -1], expected[upwind_link_at_link != -1]
    )


def test_step_function_with_no_quantity_named():
    grid = RasterModelGrid((3, 7))
    u = grid.add_zeros("advection__velocity", at="link")
    u[grid.horizontal_links] = 1.0

    advector = AdvectionSolverTVD(grid, advection_direction_is_steady=True)
    C = grid.at_node["advected__quantity"]
    C[grid.x_of_node < 2.5] = 1.0

    for _ in range(15):
        advector.run_one_step(0.2)

    assert_almost_equal(C[10], 1.0, decimal=2)

    # verify that "advection__flux" field has been created
    assert len(grid.at_link["advection__flux"]) == grid.number_of_links

    # test that we can handle existing field names
    advector = AdvectionSolverTVD(grid, advection_direction_is_steady=True)
    assert advector.grid.at_node["advected__quantity"][0] == 1.0


def test_step_function_with_single_named_quantity():
    grid = RasterModelGrid((3, 7))
    C = grid.add_zeros("one_advected__quantity", at="node")
    C[grid.x_of_node < 2.5] = 1.0
    u = grid.add_zeros("advection__velocity", at="link")
    u[grid.horizontal_links] = 1.0

    advector = AdvectionSolverTVD(
        grid,
        fields_to_advect="one_advected__quantity",
        advection_direction_is_steady=True,
    )

    for _ in range(15):
        advector.run_one_step(0.2)

    assert_almost_equal(C[10], 1.0, decimal=2)

    # verify that field "flux_of_one_advected__quantity" has been created
    assert len(grid.at_link["flux_of_one_advected__quantity"]) == grid.number_of_links


def test_with_two_advected_fields():
    grid = RasterModelGrid((3, 7))
    C1 = grid.add_zeros("first_advected__quantity", at="node")
    C1[grid.x_of_node < 2.5] = 10.0
    C2 = grid.add_zeros("second_advected__quantity", at="node")
    C2[grid.x_of_node < 2.5] = 1.0
    u = grid.add_zeros("advection__velocity", at="link")
    u[grid.horizontal_links] = 1.0

    advector = AdvectionSolverTVD(
        grid,
        fields_to_advect=["first_advected__quantity", "second_advected__quantity"],
        advection_direction_is_steady=True,
    )

    for _ in range(15):
        advector.run_one_step(0.2)

    assert_almost_equal(C1[10], 10.0, decimal=1)
    assert_almost_equal(C2[10], 1.0, decimal=2)

    # verify that flux fields have been created
    assert len(grid.at_link["flux_of_first_advected__quantity"]) == grid.number_of_links
    assert (
        len(grid.at_link["flux_of_second_advected__quantity"]) == grid.number_of_links
    )

    # verify the component can handle already-existing fields
    advector = AdvectionSolverTVD(
        grid,
        fields_to_advect=["first_advected__quantity", "second_advected__quantity"],
        advection_direction_is_steady=True,
    )
    assert (
        len(advector.grid.at_link["flux_of_first_advected__quantity"])
        == grid.number_of_links
    )
    assert (
        len(advector.grid.at_link["flux_of_second_advected__quantity"])
        == grid.number_of_links
    )


def test_with_arrays_instead_of_fields():
    grid = RasterModelGrid((3, 7))
    C1 = np.zeros(grid.number_of_nodes)
    C1[grid.x_of_node < 2.5] = 10.0
    C2 = np.zeros(grid.number_of_nodes)
    C2[grid.x_of_node < 2.5] = 1.0
    u = grid.add_zeros("advection__velocity", at="link")
    u[grid.horizontal_links] = 1.0

    advector = AdvectionSolverTVD(
        grid,
        fields_to_advect=[C1, C2],
        advection_direction_is_steady=True,
    )

    for _ in range(15):
        advector.run_one_step(0.2)

    assert_almost_equal(C1[10], 10.0, decimal=1)
    assert_almost_equal(C2[10], 1.0, decimal=2)

    # verify that flux fields have been created
    assert len(grid.at_link["advection__flux_0"]) == grid.number_of_links
    assert len(grid.at_link["advection__flux_1"]) == grid.number_of_links
