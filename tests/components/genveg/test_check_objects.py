import numpy as np
import pytest

from landlab.components.genveg.check_objects import UnitTestChecks


def test_is_negative_present():
    utc = UnitTestChecks()
    # test a single negative value
    with pytest.raises(ValueError) as excinfo:
        utc.is_negative_present(-1.5, "shoot_sys_width")
    assert str(excinfo.value) == "shoot_sys_width is negative"

    # test an array with a single negative
    with pytest.raises(ValueError) as excinfo:
        utc.is_negative_present(np.array([0, 0.004, -0.678, 1.5, 3]), "shoot_sys_width")
    assert str(excinfo.value) == "A negative value was found in shoot_sys_width array"


def test_is_zero():
    utc = UnitTestChecks()
    with pytest.raises(ValueError) as excinfo:
        utc.is_zero(0, "max_n_stem")
    assert str(excinfo.value) == "max_n_stem is zero, cannot divide by zero"
