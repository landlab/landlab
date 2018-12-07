

from numpy.random import rand
from landlab import RasterModelGrid, HexModelGrid, RadialModelGrid


_START_REFERENCE = (10., 20.)
_MOVE_REFERENCE = (30., 45.)


def test_move_reference_raster():
    mg = RasterModelGrid(9, 5,
                         xy_spacing=2.0,
                         xy_of_lower_left=_START_REFERENCE)
    assert mg._xy_of_lower_left == _START_REFERENCE
    assert mg.x_of_node.min() == _START_REFERENCE[0]
    assert mg.y_of_node.min() == _START_REFERENCE[1]

    mg.xy_of_lower_left = _MOVE_REFERENCE
    assert mg._xy_of_lower_left == _MOVE_REFERENCE
    assert mg.x_of_node.min() == _MOVE_REFERENCE[0]
    assert mg.y_of_node.min() == _MOVE_REFERENCE[1]


def test_move_reference_hex():
    shapes = ["rect", "hex"]
    orientations = ["horizontal", "vertical"]
    sizes = ((8, 7),
             (7, 8),
             (7, 11),
             (8, 12))
    for shape in shapes:
        for orientation in orientations:
            for size in sizes:
                mg = HexModelGrid(
                    size[0],
                    size[1],
                    dx=2.0,
                    xy_of_lower_left=_START_REFERENCE,
                    orientation=orientation,
                    shape=shape
                )
                assert mg._xy_of_lower_left == _START_REFERENCE
                assert mg.x_of_node[0] == _START_REFERENCE[0]
                assert mg.y_of_node[0] == _START_REFERENCE[1]

                mg.xy_of_lower_left = _MOVE_REFERENCE
                assert mg._xy_of_lower_left == _MOVE_REFERENCE
                assert mg.x_of_node[0] == _MOVE_REFERENCE[0]
                assert mg.y_of_node[0] == _MOVE_REFERENCE[1]


def test_move_reference_radial():
    mg = RadialModelGrid(num_shells=9,
                         dr=10.,
                         xy_of_center=_START_REFERENCE)

    assert mg._xy_of_center == _START_REFERENCE

    pre_move_llc = (mg.x_of_node.min(),
                    mg.y_of_node.min())

    mg.xy_of_center = _MOVE_REFERENCE
    assert mg._xy_of_center == _MOVE_REFERENCE

    post_move_llc = (mg.x_of_node.min(),
                     mg.y_of_node.min())

    actual_dydx = (post_move_llc[0] - pre_move_llc[0],
                   post_move_llc[1] - pre_move_llc[1])
    known_dydx = (_MOVE_REFERENCE[0] - _START_REFERENCE[0],
                  _MOVE_REFERENCE[1] - _START_REFERENCE[1])

    assert known_dydx == actual_dydx
