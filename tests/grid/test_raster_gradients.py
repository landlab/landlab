import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid.gradients import calc_diff_at_link as calc_diff_at_link_slow
from landlab.grid.gradients import calc_grad_at_link as calc_grad_at_link_slow
from landlab.grid.raster_gradients import calc_diff_at_link
from landlab.grid.raster_gradients import calc_grad_at_link


@pytest.mark.benchmark(group="calc_diff_at_link")
@pytest.mark.parametrize(
    "func",
    [calc_diff_at_link, calc_diff_at_link_slow],
    ids=["raster-specific", "general"],
)
def test_calc_diff_at_link_bench(benchmark, func):
    grid = RasterModelGrid((400, 5000), (1.0, 2.0))
    value_at_node = np.random.uniform(size=grid.number_of_links)
    out = grid.empty(at="link")

    benchmark(func, grid, value_at_node, out=out)


@pytest.mark.benchmark(group="calc_grad_at_link")
@pytest.mark.parametrize(
    "func",
    [calc_diff_at_link, calc_diff_at_link_slow],
    ids=["raster-specific", "general"],
)
def test_calc_grad_at_link_bench(benchmark, func):
    grid = RasterModelGrid((400, 5000), (1.0, 2.0))
    value_at_node = np.random.uniform(size=grid.number_of_links)
    out = grid.empty(at="link")

    benchmark(func, grid, value_at_node, out=out)


@pytest.mark.parametrize("shape", [(4, 5), (40, 50), (50, 40), (3, 3)])
@pytest.mark.parametrize("spacing", [(1.0, 3.0), (3.0, 1.0)])
def test_calc_diff_at_link_matches(shape, spacing):
    grid = RasterModelGrid(shape, xy_spacing=spacing)

    value_at_link = np.random.uniform(size=grid.number_of_links)
    actual = grid.empty(at="node")

    expected = calc_diff_at_link_slow(grid, value_at_link)
    actual = calc_diff_at_link(grid, value_at_link)

    assert_array_almost_equal(actual, expected)


@pytest.mark.parametrize("shape", [(4, 5), (40, 50), (50, 40), (3, 3)])
@pytest.mark.parametrize("spacing", [(1.0, 3.0), (3.0, 1.0)])
def test_calc_grad_at_link_matches(shape, spacing):
    grid = RasterModelGrid(shape, xy_spacing=spacing)

    value_at_link = np.random.uniform(size=grid.number_of_links)
    actual = grid.empty(at="node")

    expected = calc_grad_at_link_slow(grid, value_at_link)
    actual = calc_grad_at_link(grid, value_at_link)

    assert_array_almost_equal(actual, expected)


def test_calc_grad_at_active_d8():
    grid = RasterModelGrid((3, 4), xy_spacing=(4, 3))
    z = [[3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 0.0, 0.0], [3.0, 0.0, 0.0, 0.0]]
    grad_at_active_d8 = grid.calc_grad_at_d8(z)[grid.active_d8]

    assert_array_equal(
        grad_at_active_d8[: len(grid.active_links)],
        grid.calc_grad_at_link(z)[grid.active_links],
    )
    assert_array_equal(
        grad_at_active_d8[len(grid.active_links) :],
        grid.calc_grad_at_diagonal(z)[grid.active_diagonals],
    )


def test_calc_diff_at_d8():
    grid = RasterModelGrid((3, 4), xy_spacing=(4, 3))
    z = [[60.0, 60.0, 60.0, 60.0], [60.0, 60.0, 0.0, 0.0], [60.0, 0.0, 0.0, 0.0]]
    diff_at_d8 = grid.calc_diff_at_d8(z)

    assert_array_equal(diff_at_d8[: grid.number_of_links], grid.calc_diff_at_link(z))
    assert_array_equal(
        diff_at_d8[grid.number_of_links :], grid.calc_diff_at_diagonal(z)
    )


def test_calc_grad_at_d8():
    grid = RasterModelGrid((3, 4), xy_spacing=(4, 3))
    z = [[60.0, 60.0, 60.0, 60.0], [60.0, 60.0, 0.0, 0.0], [60.0, 0.0, 0.0, 0.0]]
    grad_at_d8 = grid.calc_grad_at_d8(z)

    assert_array_equal(grad_at_d8[: grid.number_of_links], grid.calc_grad_at_link(z))
    assert_array_equal(
        grad_at_d8[grid.number_of_links :], grid.calc_grad_at_diagonal(z)
    )


def test_calc_diff_at_link():
    grid = RasterModelGrid((3, 4), xy_spacing=3.0)
    z = [[0.0, 0.0, 0.0, 0.0], [0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    assert_array_equal(
        grid.calc_diff_at_link(z),
        [
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            1.0,
            0.0,
            -1.0,
            0.0,
            -1.0,
            -1.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )


def test_calc_grad_at_link():
    grid = RasterModelGrid((3, 4), xy_spacing=(2.0, 4.0))
    z = [[0.0, 0.0, 0.0, 0.0], [0.0, 4.0, 4.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    assert_array_equal(
        grid.calc_grad_at_link(z),
        [
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            2.0,
            0.0,
            -2.0,
            0.0,
            -1.0,
            -1.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )


def test_calc_diff_at_diagonal():
    grid = RasterModelGrid((3, 4), xy_spacing=3.0)
    z = [[0.0, 0.0, 0.0, 0.0], [0.0, 5.0, 5.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    assert_array_equal(
        grid.calc_diff_at_diagonal(z),
        np.array(
            [[5.0, 0.0, 5.0, 5.0, 0.0, 5.0], [0.0, -5.0, -5.0, -5.0, -5.0, 0.0]]
        ).flat,
    )


def test_calc_grad_at_diagonal():
    grid = RasterModelGrid((3, 4), xy_spacing=(3.0, 4.0))
    z = [[0.0, 0.0, 0.0, 0.0], [0.0, 5.0, 5.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    expected = np.array(
        [[1.0, 0.0, 1.0, 1.0, 0.0, 1.0], [0.0, -1.0, -1.0, -1.0, -1.0, 0.0]]
    ).flatten()
    assert_array_equal(grid.calc_grad_at_diagonal(z), expected)


def test_raster_calc_slope_at_node():
    grid = RasterModelGrid((4, 4))
    z = grid.x_of_node**2 + grid.y_of_node**2
    slopes, cmp = grid.calc_slope_at_node(z, method="Horn", return_components=True)
    assert_array_almost_equal(
        slopes[grid.core_nodes],
        [1.20591837, 1.3454815, 1.3454815, 1.39288142],
        decimal=1,
    )
    assert_array_almost_equal(
        cmp[0].reshape((4, 4))[:, 0], cmp[1].reshape((4, 4))[0, :]
    )


@pytest.mark.parametrize(
    "func",
    (
        "calc_diff_at_link",
        "calc_diff_at_diagonal",
        "calc_diff_at_d8",
        "calc_grad_at_link",
        "calc_grad_at_diagonal",
        "calc_grad_at_d8",
    ),
)
def test_out_keyword(func):
    grid = RasterModelGrid((3, 4), xy_spacing=(3.0, 4.0))
    z = [[0.0, 0.0, 0.0, 0.0], [0.0, 5.0, 5.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    expected = getattr(grid, func)(z)
    out = np.empty_like(expected)
    actual = getattr(grid, func)(z, out=out)

    assert actual is out
    assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    "func",
    (
        "calc_diff_at_link",
        "calc_diff_at_diagonal",
        "calc_diff_at_d8",
        "calc_grad_at_link",
        "calc_grad_at_diagonal",
        "calc_grad_at_d8",
    ),
)
def test_array_as_field(func):
    grid = RasterModelGrid((3, 4), xy_spacing=(3.0, 4.0))
    z = [[0.0, 0.0, 0.0, 0.0], [0.0, 5.0, 5.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    grid.add_field("elevation", z, at="node")

    assert_array_equal(getattr(grid, func)("elevation"), getattr(grid, func)(z))
