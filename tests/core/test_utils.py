import numpy as np
import pytest
from requireit import ValidationError

from landlab.core.utils import require_id_array


@pytest.mark.parametrize("array", ([-1, 2, 3], (1, 2, 3), np.arange(3)))
def test_require_id_array_valid_returns_input(array):
    assert require_id_array(array) is array


@pytest.mark.parametrize("dtype", (int, np.intp, np.int32))
def test_require_id_array_with_integers(dtype):
    require_id_array(np.array([1, 2, 3], dtype=dtype))


@pytest.mark.parametrize("array", ([1.0, 2, 3], (True, False)))
def test_require_id_array_not_integer(array):
    with pytest.raises(ValidationError):
        require_id_array(array)


def test_require_id_array_with_bad_id():
    require_id_array([-2, 1, 2, 3], bad_id=-2)
    with pytest.raises(ValidationError, match="^array must"):
        require_id_array([-1, 1, 2, 3], bad_id=None)


def test_require_id_array_with_name():
    with pytest.raises(ValidationError, match="^foobar must"):
        require_id_array([1.0, 2, 3], name="foobar")


def test_require_id_array_with_shape():
    require_id_array([[1, 2, 3]], shape=("n", 3))
    with pytest.raises(ValidationError, match="^array must"):
        require_id_array([1, 2, 3], shape=("n", 3))


def test_require_id_array_empty():
    require_id_array(np.array([], dtype=int))
