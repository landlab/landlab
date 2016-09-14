"""Test the utils.find_shoreline function."""
import numpy as np
from nose.tools import raises, assert_almost_equal

from landlab.components.submarine_diffusion.utils import find_shoreline


# Test data.
x = np.arange(10.0)
z = 5.0 - x
sea_level = 0.25
expected_value = 5.0
expected_value_with_sea_level = 4.75
hi_value = 100.0
lo_value = -100.0


@raises(TypeError)
def test_find_shoreline_fails_with_no_args():
    """Test find_shoreline fails with no arguments"""
    find_shoreline()


@raises(TypeError)
def test_find_shoreline_fails_with_one_arg():
    """Test find_shoreline fails with one argument"""
    find_shoreline(x)


def test_find_shoreline_with_default_keywords():
    """Test find_shoreline with the keyword defaults"""
    find_shoreline(x, z)


@raises(NotImplementedError)
def test_find_shoreline_fails_with_unknown_kind():
    """Test find_shoreline fails with unknown interpolation"""
    find_shoreline(x, z, kind='foobarbaz')


def test_find_shoreline_return_value():
    """Test find_shoreline return value"""
    r = find_shoreline(x, z)
    assert_almost_equal(r, expected_value)


def test_find_shoreline_return_value_with_sea_level():
    """Test find_shoreline return value with sea level"""
    r = find_shoreline(x, z, sea_level=sea_level)
    assert_almost_equal(r, expected_value_with_sea_level)


def test_find_shoreline_return_value_with_hi_sea_level():
    """Test find_shoreline return value with high sea level"""
    r = find_shoreline(x, z, sea_level=hi_value)
    assert_almost_equal(r, x[0])


def test_find_shoreline_return_value_with_lo_sea_level():
    """Test find_shoreline return value with low sea level"""
    r = find_shoreline(x, z, sea_level=lo_value)
    assert_almost_equal(r, x[-1])
