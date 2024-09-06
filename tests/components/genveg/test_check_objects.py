import pytest
import numpy as np

from landlab.components.genveg.check_objects import UnitTestChecks


def test_is_negative_present():
    utc = UnitTestChecks()
    with pytest.raises(ValueError):
        # test single negative value
        utc.is_negative_present(-1.5, 'shoot_sys_width')
        # test an array with a signle negative
        utc.is_negative_present(np.array([0, 0.004, -0.678, 1.5, 3]), 'shoot_sys_width')


def test_is_zero():
    utc = UnitTestChecks()
    with pytest.raises(ValueError):
        utc.is_zero(0, 'max_n_stem')
