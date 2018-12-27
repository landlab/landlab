

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from landlab import RasterModelGrid, CORE_NODE, CLOSED_NODE
from landlab.values import random, plane, constant, sine
from landlab.values.synthetic import _plane_function

_NORMAL = (1, 1, 1)
_POINT = (0, 0, 0)


def test_bad_grid_element_name():
    mg = RasterModelGrid((4, 4))
    with pytest.raises(KeyError):
        constant(mg, 'some_flux', 'not_a_grid_element', constant=1.0)


def test_bad_distribution_name():
    mg = RasterModelGrid((4, 4))
    with pytest.raises(ValueError):
        random(mg, 'values', 'node', distribution='not_a_distribution')


def test_vertical_plane():
    mg = RasterModelGrid((4, 4))
    with pytest.raises(ValueError):
        plane(mg, 'values', normal=(0., 1., 0.))


def test_no_xy_values():
    mg = RasterModelGrid((4, 4))
    with pytest.raises(ValueError):
        plane(mg, 'values', "patch")


def test_xy_node():
    mg = RasterModelGrid((4, 4))
    vals = plane(mg, 'values', "node", normal=_NORMAL)
    truth = _plane_function(mg.x_of_node, mg.y_of_node,
                            point=_POINT,
                            normal=_NORMAL)
    assert_array_equal(vals, truth)


def test_xy_cell():
    mg = RasterModelGrid((4, 4))
    vals = plane(mg, 'values', "cell", normal=_NORMAL)
    truth = _plane_function(mg.x_of_cell, mg.y_of_cell,
                            point=_POINT,
                            normal=_NORMAL)
    assert_array_equal(vals, truth)


def test_xy_link():
    mg = RasterModelGrid((4, 4))
    vals = plane(mg, 'values', "link", normal=_NORMAL)
    truth = _plane_function(mg.x_of_link, mg.y_of_link,
                            point=_POINT,
                            normal=_NORMAL)
    assert_array_equal(vals, truth)


def test_xy_patch():
    mg = RasterModelGrid((4, 4))
    vals = plane(mg, 'values', "face", normal=_NORMAL)
    truth = _plane_function(mg.x_of_face, mg.y_of_face,
                            point=_POINT,
                            normal=_NORMAL)
    assert_array_equal(vals, truth)


def test_where_status_with_patches():
    mg = RasterModelGrid((4, 4))
    with pytest.raises(ValueError):
        constant(mg, 'values', "patch", where=CORE_NODE)


def test_multiple_status_node():
    mg = RasterModelGrid((4, 4))
    mg.se

# test all grid elements


# where is size of elements
# where is list of types
# where is list but use patches (raises value error)


# x and y of multiple types of grid elements
# ask for a plane of something without X and y.
