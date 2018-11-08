

from numpy.random import rand
from landlab import RasterModelGrid, HexModelGrid, RadialModelGrid, VoronoiDelaunayGrid

_START_ORIGIN = (10., 20.)
_MOVE_ORIGIN = (30., 45.)
def test_move_origin_raster():
    mg = RasterModelGrid(9, 5, dx=2.0, origin=_START_ORIGIN)
    assert mg._origin == _START_ORIGIN
    assert mg.x_of_node.min() == _START_ORIGIN[1]
    assert mg.y_of_node.min() == _START_ORIGIN[0]

    mg.origin = _MOVE_ORIGIN
    assert mg._origin == _MOVE_ORIGIN
    assert mg.x_of_node.min() == _MOVE_ORIGIN[1]
    assert mg.y_of_node.min() == _MOVE_ORIGIN[0]


def test_move_origin_hex():
    shapes = ["rect", "hex"]
    orientations = ["horizontal", "vertical"]

    for shape in shapes:
        for orientation in orientations:
            mg = HexModelGrid(
                9, 5, dx=2.0, origin=_START_ORIGIN, orientation=orientation, shape=shape
            )
            assert mg._origin == _START_ORIGIN
            assert mg.x_of_node[0] == _START_ORIGIN[1]
            assert mg.y_of_node[0] == _START_ORIGIN[0]

            mg.origin = _MOVE_ORIGIN
            assert mg._origin == _MOVE_ORIGIN
            assert mg.x_of_node[0] == _MOVE_ORIGIN[1]
            assert mg.y_of_node[0] == _MOVE_ORIGIN[0]


def test_move_origin_radial():
    mg = RadialModelGrid(num_shells=9, dr=10., origin_x=_START_ORIGIN[1], origin_y=_START_ORIGIN[0])

    assert mg._origin == _START_ORIGIN

    pre_move_llc = (mg.y_of_node.min(),
                    mg.x_of_node.min())

    mg.origin = _MOVE_ORIGIN
    assert mg._origin == _MOVE_ORIGIN

    post_move_llc = (mg.y_of_node.min(),
                     mg.x_of_node.min())

    actual_dydx = (post_move_llc[0] - pre_move_llc[0],
                   post_move_llc[1] - pre_move_llc[1])
    known_dydx = (_MOVE_ORIGIN[0] - _START_ORIGIN[0],
                  _MOVE_ORIGIN[1] - _START_ORIGIN[1])

    assert  known_dydx == actual_dydx


def test_move_origin_voronoi():
    x, y = rand(25), rand(25)
    mg = VoronoiDelaunayGrid(x, y)

    assert mg._origin == (y.min(), x.min())
    mg.origin = _MOVE_ORIGIN
    assert mg._origin == _MOVE_ORIGIN
