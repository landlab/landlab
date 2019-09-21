import numpy as np
import pytest
from numpy.testing import assert_array_equal
from pytest import approx

from landlab import HexModelGrid, RadialModelGrid, RasterModelGrid


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
        9,
        5,
        xy_of_lower_left=to_iterable(random_xy),
        orientation="horizontal",
        node_layout="rect",
    )
    assert isinstance(grid.xy_of_lower_left, tuple)
    assert grid.xy_of_lower_left == expected


@pytest.mark.parametrize("to_iterable", [np.asarray, list, tuple])
def test_radial_center_as_iterables(random_xy, to_iterable):
    expected = approx(tuple(random_xy))

    grid = RadialModelGrid(num_shells=9, dr=10.0, xy_of_center=to_iterable(random_xy))
    assert isinstance(grid.xy_of_center, tuple)
    assert grid.xy_of_center == expected


@pytest.mark.parametrize("orientation", ["horizontal", "vertical"])
@pytest.mark.parametrize("node_layout", ["rect"])
@pytest.mark.parametrize("n_cols", [12, 11, 10, 9])
@pytest.mark.parametrize("n_rows", [12, 11, 10, 9])
def test_move_reference_hex(random_xy, n_rows, n_cols, node_layout, orientation):
    mg = HexModelGrid(
        n_rows,
        n_cols,
        dx=2.0,
        xy_of_lower_left=random_xy,
        orientation=orientation,
        node_layout=node_layout,
    )

    assert mg.xy_of_lower_left == random_xy
    assert mg.x_of_node.min() == random_xy[0]
    assert mg.y_of_node.min() == random_xy[1]

    mg.xy_of_lower_left = (30.0, 45.0)
    assert mg._xy_of_lower_left == (30.0, 45.0)
    assert mg.x_of_node.min() == approx(30.0)
    assert mg.y_of_node.min() == approx(45.0)


def test_move_reference_radial(random_xy):
    mg = RadialModelGrid(num_shells=9, dr=10.0, xy_of_center=random_xy)

    assert mg.xy_of_center == approx(random_xy)

    pre_move_llc = (mg.x_of_node.min(), mg.y_of_node.min())

    new_xy_of_center = (30.0, 45.0)
    mg.xy_of_center = new_xy_of_center
    assert mg.xy_of_center == approx(new_xy_of_center)

    post_move_llc = (mg.x_of_node.min(), mg.y_of_node.min())

    actual_dydx = (
        post_move_llc[0] - pre_move_llc[0],
        post_move_llc[1] - pre_move_llc[1],
    )
    known_dydx = (
        new_xy_of_center[0] - random_xy[0],
        new_xy_of_center[1] - random_xy[1],
    )

    assert known_dydx == approx(actual_dydx)


def test_radial_deprecate_origin_x():
    with pytest.deprecated_call():
        mg = RadialModelGrid(num_shells=1, dr=1.0, origin_x=10)
    assert mg._xy_of_center == (10.0, 0.0)
    pts, npts = mg._create_radial_points(1, 1, xy_of_center=mg._xy_of_center)
    assert pts[0, 0] == approx(mg._xy_of_center[0])
    assert pts[0, 1] == approx(mg._xy_of_center[1])


def test_radial_deprecate_origin_y():
    with pytest.deprecated_call():
        mg = RadialModelGrid(num_shells=1, dr=1.0, origin_y=10)
    assert mg._xy_of_center == (0.0, 10.0)
    pts, npts = mg._create_radial_points(1, 1, xy_of_center=mg._xy_of_center)
    assert pts[0, 0] == approx(mg._xy_of_center[0])
    assert pts[0, 1] == approx(mg._xy_of_center[1])


def test_raster_with_args_and_shape():
    with pytest.deprecated_call():
        RasterModelGrid(3, 3, num_cols=3)


def test_raster_with_negative_shape():
    with pytest.raises(ValueError):
        with pytest.deprecated_call():
            RasterModelGrid(-2, 3)


def test_raise_deprecation_dx():
    with pytest.deprecated_call():
        mg = RasterModelGrid(3, 3, dx=3)
    X, Y = np.meshgrid([0.0, 3.0, 6.0], [0.0, 3.0, 6.0])
    assert_array_equal(mg.x_of_node, X.flatten())
    assert_array_equal(mg.y_of_node, Y.flatten())


def test_raise_deprecation_dxdx():
    with pytest.deprecated_call():
        mg = RasterModelGrid(3, 3, dx=(4, 5))
    X, Y = np.meshgrid([0.0, 5.0, 10.0], [0.0, 4.0, 8.0])
    assert_array_equal(mg.x_of_node, X.flatten())
    assert_array_equal(mg.y_of_node, Y.flatten())


def test_raise_deprecation_spacing():
    with pytest.deprecated_call():
        mg = RasterModelGrid(3, 3, spacing=(4, 5))
    X, Y = np.meshgrid([0.0, 5.0, 10.0], [0.0, 4.0, 8.0])
    assert_array_equal(mg.x_of_node, X.flatten())
    assert_array_equal(mg.y_of_node, Y.flatten())


def test_raise_deprecation_spacing_as_arg():
    with pytest.deprecated_call():
        mg = RasterModelGrid(3, 3, 3)
    X, Y = np.meshgrid([0.0, 3.0, 6.0], [0.0, 3.0, 6.0])
    assert_array_equal(mg.x_of_node, X.flatten())
    assert_array_equal(mg.y_of_node, Y.flatten())


def test_raise_deprecation_spacing2_as_arg():
    with pytest.deprecated_call():
        mg = RasterModelGrid(3, 3, (4, 5))
    X, Y = np.meshgrid([0.0, 5.0, 10.0], [0.0, 4.0, 8.0])
    assert_array_equal(mg.x_of_node, X.flatten())
    assert_array_equal(mg.y_of_node, Y.flatten())


def test_bad_shape_xy_spacing():
    with pytest.raises(ValueError):
        with pytest.deprecated_call():
            RasterModelGrid(3, 3, xy_spacing=(4, 5, 5))


def test_bad_type_xy_spacing():
    with pytest.raises(ValueError):
        with pytest.deprecated_call():
            RasterModelGrid(3, 3, xy_spacing="spam and eggs")


def test_deprecate_origin():
    with pytest.deprecated_call():
        mg = RasterModelGrid(3, 3, origin=(10, 13))
    X, Y = np.meshgrid([0.0, 1.0, 2.0], [0.0, 1.0, 2.0])
    assert_array_equal(mg.x_of_node, X.flatten() + 10.0)
    assert_array_equal(mg.y_of_node, Y.flatten() + 13.0)


def test_bad_origin():
    with pytest.raises(ValueError):
        with pytest.deprecated_call():
            RasterModelGrid(3, 3, xy_of_lower_left=(10, 13, 12))


def test_curent_vs_past_origin():
    with pytest.deprecated_call():
        mg1 = RasterModelGrid(3, 3, origin=(10, 13))
    with pytest.deprecated_call():
        mg2 = RasterModelGrid(3, 3, xy_of_lower_left=(10, 13))
    assert_array_equal(mg1.x_of_node, mg2.x_of_node)
    assert_array_equal(mg1.y_of_node, mg2.y_of_node)


def test_curent_vs_past_spacing():
    with pytest.deprecated_call():
        mg1 = RasterModelGrid(3, 3, spacing=(5, 4))
    with pytest.deprecated_call():
        mg2 = RasterModelGrid(3, 3, xy_spacing=(4, 5))
    assert_array_equal(mg1.x_of_node, mg2.x_of_node)
    assert_array_equal(mg1.y_of_node, mg2.y_of_node)
