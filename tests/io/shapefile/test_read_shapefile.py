import os

from io import BytesIO
import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx, raises
import shapefile
from shapefile import ShapefileException

from landlab.io.shapefile import read_shapefile

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_read_methow():
    file = os.path.join(_TEST_DATA_DIR, "Methow_Network.shp")
    grid = read_shapefile(file)
    assert grid.number_of_nodes == 721
    assert grid.number_of_links == 720

    assert grid.x_of_node[0] == approx(-1672349.0889982011)
    assert grid.y_of_node[0] == approx(1160800.240247)
    assert "x_of_polyline" in grid.at_link
    assert "y_of_polyline" in grid.at_link

    fields = ["Length_m", "ToLink", "usarea_km2", "uselev_m", "dselev_m", "Slope"]
    for field in fields:
        assert field in grid.at_link

        # todo add test for order.


def test_bad_file():
    file = os.path.join(_TEST_DATA_DIR, "bad_file.shp")
    with raises(ShapefileException):
        read_shapefile(file)


def test_points():
    file = os.path.join(_TEST_DATA_DIR, "points.shp")
    with raises(ValueError):
        read_shapefile(file)


def test_multipart():
    file = os.path.join(_TEST_DATA_DIR, "multipartpolyline.shp")
    with raises(ValueError):
        read_shapefile(file)


def test_simple_reorder():

    orders = [(0, 1, 2), (0, 2, 1),
              (1, 0, 2), (1, 2, 0),
              (2, 0, 1), (2, 1, 0)]

    lines = [[[[5,5],[10,10]]],
             [[[5,0],[5,5]]],
             [[[5,5],[0,10]]],
              ]

    records = [37, 100, 239]

    for order in orders:

        shp = BytesIO()
        shx = BytesIO()
        dbf = BytesIO()

        w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)

        w.shapeType = 3
        w.field("spam", "N")

        for o in order:
            w.line(lines[o])
            w.record(records[o])
        w.close()

        grid = read_shapefile(shp, dbf=dbf)

        assert_array_equal( grid.nodes, np.array([0, 1, 2, 3]))
        assert_array_equal( grid.x_of_node, np.array([  5.,   5.,   0.,  10.]))
        assert_array_equal( grid.y_of_node, np.array([  0.,   5.,  10.,  10.]))
        assert_array_equal( grid.nodes_at_link, np.array([[0, 1],
               [2, 1],
               [1, 3]]))
        assert "spam" in grid.at_link
        assert_array_equal(grid.at_link["spam"], np.array([100, 239, 37]))

        del grid, w, shp, shx, dbf
