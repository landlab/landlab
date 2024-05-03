import numpy as np
import pytest
from pytest import approx

from landlab import HexModelGrid
from landlab import RadialModelGrid
from landlab import RasterModelGrid


def test_xy_of_reference_default_is_zero():
    grid = RasterModelGrid((9, 5))
    assert grid.xy_of_reference == approx((0.0, 0.0))


@pytest.mark.parametrize("to_iterable", [np.asarray, list, tuple])
def test_xy_of_reference_is_tuple(random_xy, to_iterable):
    grid = RasterModelGrid((9, 5), xy_of_reference=to_iterable(random_xy))
    assert isinstance(grid.xy_of_reference, tuple)
    assert grid.xy_of_reference == approx(random_xy)


def test_xy_of_reference_setter(random_xy):
    grid = RasterModelGrid((9, 5))
    grid.xy_of_reference = random_xy
    assert grid.xy_of_reference == approx(random_xy)


def test_move_reference_raster(random_xy):
    mg = RasterModelGrid((9, 5), xy_spacing=2.0, xy_of_lower_left=random_xy)
    assert mg.xy_of_lower_left == approx(random_xy)
    assert mg.x_of_node.min() == approx(random_xy[0])
    assert mg.y_of_node.min() == approx(random_xy[1])

    xy_of_new_lower_left = random_xy[0] * 10.0, random_xy[1] * 10.0

    mg.xy_of_lower_left = xy_of_new_lower_left
    assert mg.xy_of_lower_left == approx(xy_of_new_lower_left)
    assert mg.x_of_node.min() == approx(xy_of_new_lower_left[0])
    assert mg.y_of_node.min() == approx(xy_of_new_lower_left[1])


@pytest.mark.parametrize("to_iterable", [np.asarray, list, tuple])
def test_raster_lower_left_as_iterables(random_xy, to_iterable):
    expected = approx(tuple(random_xy))

    grid = RasterModelGrid((9, 5), xy_of_lower_left=to_iterable(random_xy))
    assert isinstance(grid.xy_of_lower_left, tuple)
    assert grid.xy_of_lower_left == expected


@pytest.mark.parametrize("to_iterable", [np.asarray, list, tuple])
def test_hex_lower_left_as_iterables(random_xy, to_iterable):
    expected = approx(tuple(random_xy))

    grid = HexModelGrid(
        (9, 5),
        xy_of_lower_left=to_iterable(random_xy),
        orientation="horizontal",
        node_layout="rect",
    )
    assert isinstance(grid.xy_of_lower_left, tuple)
    assert grid.xy_of_lower_left == expected


@pytest.mark.parametrize("to_iterable", [np.asarray, list, tuple])
def test_radial_center_as_iterables(random_xy, to_iterable):
    expected = approx(tuple(random_xy))

    grid = RadialModelGrid(
        n_rings=9, nodes_in_first_ring=8, xy_of_center=to_iterable(random_xy)
    )
    assert isinstance(grid.xy_of_center, tuple)
    assert grid.xy_of_center == expected


@pytest.mark.parametrize("orientation", ["horizontal", "vertical"])
@pytest.mark.parametrize("node_layout", ["rect", "hex"])
@pytest.mark.parametrize("n_cols", [12, 11, 10, 9])
@pytest.mark.parametrize("n_rows", [12, 11, 10, 9])
def test_move_reference_hex(random_xy, n_rows, n_cols, node_layout, orientation):
    grid = HexModelGrid(
        (n_rows, n_cols),
        spacing=2.0,
        xy_of_lower_left=random_xy,
        orientation=orientation,
        node_layout=node_layout,
    )

    assert grid.xy_of_lower_left == approx(random_xy)
    assert grid.x_of_node.min() == approx(random_xy[0])
    assert grid.y_of_node.min() == approx(random_xy[1])

    x, y = grid.x_of_node.copy(), grid.y_of_node.copy()

    grid.xy_of_lower_left = (30.0, 45.0)
    assert grid.xy_of_lower_left == (30.0, 45.0)

    assert grid.x_of_node - x == approx(30.0 - random_xy[0])
    assert grid.y_of_node - y == approx(45.0 - random_xy[1])


def test_move_reference_radial(random_xy):
    grid = RadialModelGrid(n_rings=9, nodes_in_first_ring=8, xy_of_center=random_xy)

    assert grid.xy_of_center == approx(random_xy)

    x, y = grid.x_of_node.copy(), grid.y_of_node.copy()

    grid.xy_of_center = (30.0, 45.0)
    assert grid.xy_of_center == (30.0, 45.0)

    assert grid.x_of_node - x == approx(30.0 - random_xy[0])
    assert grid.y_of_node - y == approx(45.0 - random_xy[1])


def test_bad_shape_xy_spacing():
    with pytest.raises(ValueError):
        RasterModelGrid((3, 3), xy_spacing=(4, 5, 5))


def test_bad_type_xy_spacing():
    with pytest.raises(ValueError):
        RasterModelGrid((3, 3), xy_spacing="spam and eggs")


def test_bad_origin():
    with pytest.raises(ValueError):
        RasterModelGrid((3, 3), xy_of_lower_left=(10, 13, 12))
