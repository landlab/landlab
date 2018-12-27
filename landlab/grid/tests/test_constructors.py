import numpy as np
import pytest


from landlab import CLOSED_BOUNDARY, HexModelGrid, RasterModelGrid, NetworkModelGrid, VoronoiDelaunayGrid, RadialModelGrid

from landlab.grid.raster import from_dict as raster_from_dict
from landlab.grid.hex import from_dict as hex_from_grid


def test_base_not_implemented():
    pass


def test_raster_old_from_dict_deprecated():
    pass


def test_hex_old_from_dict_deprecated():
    pass


def test_raster_from_file():
    pass


def test_raster_from_dict():
    params = {
        "shape": (10, 20),
        "xy_spacing": (25, 45),
        "bc": {
            "right": "closed",
            "top": "closed",
            "left": "closed",
            "bottom": "closed",
        },
        "xy_of_lower_left": (35, 55),
        "axis_name": ("spam", "eggs"),
        "axis_units": ("smoot", "parsec"),
    }

    mg = RasterModelGrid.from_dict(params)

    # assert things.
    assert mg.shape == mg.shape
    assert mg.dx == 25
    assert mg.dy == 45
    assert (mg.x_of_node.min(), mg.y_of_node.min()) == (35, 55)
    assert np.all(mg.status_at_node[mg.boundary_nodes] == CLOSED_BOUNDARY)
    assert mg.axis_units == ("smoot", "parsec")
    assert mg.axis_name == ("spam", "eggs")


def test_hex_from_dict():
    pass


def test_radial_from_dict():
    pass


def test_network_from_dict():
    pass


def test_voronoi_from_dict():
    pass
