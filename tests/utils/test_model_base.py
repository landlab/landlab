from contextlib import chdir

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from requireit import ValidationError

from landlab import RasterModelGrid
from landlab.io.native_landlab import save_grid
from landlab.utils.model_base import read_arrays_from_files
from landlab.utils.model_base import setup_grid


def test_setup_grid_creates_grid():
    params = {
        "source": "create",
        "create_grid": {
            "RasterModelGrid": {
                "shape": (4, 5),
                "xy_spacing": (2.0, 4.0),
            }
        },
    }

    grid = setup_grid(params)

    assert isinstance(grid, RasterModelGrid)
    assert grid.shape == (4, 5)
    assert grid.dx == 2.0
    assert grid.dy == 4.0


def test_setup_grid_loads_grid(tmp_path):
    original = RasterModelGrid((3, 4), xy_spacing=(2.0, 4.0))
    original.add_zeros("topographic__elevation", at="node")
    path = tmp_path / "model.grid"
    save_grid(original, path)

    actual = setup_grid({"source": "file", "grid_file_name": path})

    assert isinstance(actual, RasterModelGrid)
    assert actual.shape == original.shape
    assert actual.spacing == original.spacing
    assert "topographic__elevation" in actual.at_node


def test_setup_grid_uses_existing_grid():
    expected = RasterModelGrid((3, 4))

    actual = setup_grid({"source": "grid_object", "grid_object": expected})

    assert actual is expected


def test_setup_grid_rejects_unknown_source():
    with pytest.raises(ValidationError, match="^source must be one of"):
        setup_grid({"source": "unknown"})


def test_setup_grid_rejects_non_grid_object():
    with pytest.raises(ValueError, match="^grid source must be"):
        setup_grid({"source": "grid_object", "grid_object": object()})


def test_read_arrays_from_files(tmp_path):
    expected_1d = np.arange(3) / 2
    expected_2d = np.arange(4).reshape((2, 2)) / 2
    expected_col = np.arange(10).reshape((-1, 1)) / 4

    np.save(tmp_path / "test1", expected_1d)
    np.save(tmp_path / "test2", expected_2d)
    np.save(tmp_path / "test3", expected_col)

    p = {
        "a": 123,
        "b": {"c": 456, "d": {"_filepath": "test1.npy"}},
        "e": {"_filepath": "test2.npy"},
        "f": {"_filepath": "test3.npy"},
    }
    with chdir(tmp_path):
        p = read_arrays_from_files(p)

    assert_array_equal(p["b"]["d"], expected_1d)
    assert_array_equal(p["e"], expected_2d)
    assert_array_equal(p["f"], expected_col)
