import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


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
