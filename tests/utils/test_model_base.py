from contextlib import chdir

import numpy as np
from numpy.testing import assert_array_equal

from landlab.utils.model_base import read_arrays_from_files


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
