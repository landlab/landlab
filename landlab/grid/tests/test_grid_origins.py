

from numpy.random import rand
from landlab import RasterModelGrid, HexModelGrid, RadialModelGrid, VoronoiDelaunayGrid


def test_move_origin_raster():
    mg = RasterModelGrid(9, 5, dx=2.0, origin=(10., 20.))
    assert mg._origin == (20., 10.)
    assert mg.x_of_node.min() == 20.
    assert mg.y_of_node.min() == 10.

    mg.origin == (0., 0.)
    assert mg._origin == (0., 0.)
    assert mg.x_of_node.min() == 0.
    assert mg.y_of_node.min() == 0.


def test_move_origin_hex():
    shapes = ["rect", "hex"]
    orientations = ["horizontal", "vertical"]

    for shape in shapes:
        for orientation in orientations:
            mg = HexModelGrid(
                9, 5, dx=2.0, origin=(10., 20.), orientation=orientation, shape=shape
            )
            assert mg._origin == (20., 10.)
            assert mg.x_of_node.min() == 20.
            assert mg.y_of_node.min() == 10.

            mg.origin == (0., 0.)
            assert mg._origin == (0., 0.)
            assert mg.x_of_node.min() == 0.
            assert mg.y_of_node.min() == 0.


def test_move_origin_radial():
    mg = RadialModelGrid(num_shells=9, dr=10., origin_x=20., origin_y=10.)
    assert mg._origin == (20., 10.)
    assert mg.x_of_node.min() == 20.
    assert mg.y_of_node.min() == 10.

    mg.origin = (0., 0.)
    assert mg._origin == (0., 0.)
    assert mg.x_of_node.min() == 0.
    assert mg.y_of_node.min() == 0.


def test_move_origin_voronoi():
    x, y = rand(25), rand(25)
    mg = VoronoiDelaunayGrid(x, y)

    assert mg._origin == (y.min(), x.min())
    mg.origin = (0., 0.)
    assert mg._origin == (0., 0.)
