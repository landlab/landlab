import os

from pytest import approx, raises
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
