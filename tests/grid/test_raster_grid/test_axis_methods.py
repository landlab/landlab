import pytest

from landlab import RasterModelGrid


def test_default_names():
    rmg = RasterModelGrid((4, 5))
    assert rmg.axis_name == ("x", "y")


def test_name_keyword():
    rmg = RasterModelGrid(
        (4, 5),
        xy_axis_name=["longitude", "latitude"],
        xy_axis_units=["degrees_east", "degrees_north"],
    )
    assert rmg.axis_name == ("longitude", "latitude")


def test_name_setter():
    rmg = RasterModelGrid(
        (4, 5),
        xy_axis_name=["longitude", "latitude"],
        xy_axis_units=["degrees_east", "degrees_north"],
    )
    rmg.axis_name = ("yyy", "xxx")
    assert rmg.axis_name == ("yyy", "xxx")


def test_name_setter_too_few_names():
    rmg = RasterModelGrid((4, 5))
    rmg.axis_name = "z"
    assert rmg.axis_name == ("z", "z")


def test_name_setter_too_many_names():
    rmg = RasterModelGrid((4, 5))
    with pytest.raises(ValueError):
        rmg.axis_name = ("z", "y", "x")


def test_default_units():
    rmg = RasterModelGrid((4, 5))
    assert rmg.axis_units == ("-", "-")


def test_axis_units_keyword():
    rmg = RasterModelGrid(
        (4, 5),
        xy_axis_name=["longitude", "latitude"],
        xy_axis_units=["degrees_east", "degrees_north"],
    )
    assert rmg.axis_units == ("degrees_east", "degrees_north")


def test_axis_units_setter():
    rmg = RasterModelGrid(
        (4, 5),
        xy_axis_name=["longitude", "latitude"],
        xy_axis_units=["degrees_east", "degrees_north"],
    )
    rmg.axis_units = ("mm", "cm")
    assert rmg.axis_units == ("mm", "cm")


def test_name_setter_too_few_units():
    grid = RasterModelGrid((4, 5))
    grid.axis_units = "m"
    assert grid.axis_units == ("m", "m")


def test_name_setter_too_many_units():
    rmg = RasterModelGrid((4, 5))
    with pytest.raises(ValueError):
        rmg.axis_units = ("m", "cm", "mm")
