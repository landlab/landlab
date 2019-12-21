#! /usr/bin/env python
import os

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.io import read_esri_ascii, write_esri_ascii

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_grid_with_no_fields(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    with tmpdir.as_cwd():
        with pytest.raises(ValueError):
            write_esri_ascii("test.asc", grid)


def test_grid_with_one_field(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    grid.add_field("node", "air__temperature", np.arange(20.0))
    with tmpdir.as_cwd():
        files = write_esri_ascii("test.asc", grid)
        assert files == ["test.asc"]
        for fname in files:
            assert os.path.isfile(fname)


def test_grid_with_two_fields(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    grid.add_field("node", "air__temperature", np.arange(20.0))
    grid.add_field("node", "land_surface__elevation", np.arange(20.0))
    with tmpdir.as_cwd():
        files = write_esri_ascii("test.asc", grid)
        files.sort()
        assert files == [
            "test_air__temperature.asc",
            "test_land_surface__elevation.asc",
        ]
        for fname in files:
            assert os.path.isfile(fname)


def test_names_keyword_as_str(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    grid.add_field("air__temperature", np.arange(20.0), at="node")
    grid.add_field("land_surface__elevation", np.arange(20.0), at="node")

    with tmpdir.as_cwd():
        files = write_esri_ascii("test.asc", grid, names="air__temperature")
        assert files == ["test.asc"]
        assert os.path.isfile("test.asc")


def test_names_keyword_as_list(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    grid.add_field("air__temperature", np.arange(20.0), at="node")
    grid.add_field("land_surface__elevation", np.arange(20.0), at="node")

    with tmpdir.as_cwd():
        files = write_esri_ascii("test.asc", grid, names=["air__temperature"])
        assert files == ["test.asc"]
        assert os.path.isfile("test.asc")


def test_names_keyword_multiple_names(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    grid.add_field("node", "air__temperature", np.arange(20.0))
    grid.add_field("node", "land_surface__elevation", np.arange(20.0))

    with tmpdir.as_cwd():
        files = write_esri_ascii(
            "test.asc", grid, names=["air__temperature", "land_surface__elevation"]
        )
        files.sort()
        assert files == [
            "test_air__temperature.asc",
            "test_land_surface__elevation.asc",
        ]
        for fname in files:
            assert os.path.isfile(fname)


def test_names_keyword_with_bad_name(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    grid.add_field("node", "air__temperature", np.arange(20.0))

    with tmpdir.as_cwd():
        with pytest.raises(ValueError):
            write_esri_ascii("test.asc", grid, names="not_a_name")


def test_clobber_keyword(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    grid.add_field("node", "air__temperature", np.arange(20.0))

    with tmpdir.as_cwd():
        write_esri_ascii("test.asc", grid)
        with pytest.raises(ValueError):
            write_esri_ascii("test.asc", grid)
        with pytest.raises(ValueError):
            write_esri_ascii("test.asc", grid, clobber=False)
        write_esri_ascii("test.asc", grid, clobber=True)


def test_write_then_read(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0), xy_of_lower_left=(15.0, 10.0))
    grid.add_field("node", "air__temperature", np.arange(20.0))

    with tmpdir.as_cwd():
        write_esri_ascii("test.asc", grid)
        new_grid, field = read_esri_ascii("test.asc")

    assert grid.number_of_node_columns == new_grid.number_of_node_columns
    assert grid.number_of_node_rows == new_grid.number_of_node_rows
    assert grid.dx == new_grid.dx
    assert (grid.x_of_node.min(), grid.y_of_node.min()) == (15.0, 10.0)
    assert_array_almost_equal(grid.node_x, new_grid.node_x)
    assert_array_almost_equal(grid.node_y, new_grid.node_y)
    assert_array_almost_equal(field, grid.at_node["air__temperature"])
