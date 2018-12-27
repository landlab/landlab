import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import CLOSED_BOUNDARY, CORE_NODE
from landlab.values import constant, plane, random
from landlab.values.synthetic import _plane_function

_NORMAL = (1, 1, 1)
_POINT = (0, 0, 0)


def test_bad_grid_element_name(four_by_four_raster):
    with pytest.raises(KeyError):
        constant(four_by_four_raster, "some_flux", "not_a_grid_element", constant=1.0)


def test_bad_distribution_name(four_by_four_raster):
    with pytest.raises(ValueError):
        random(four_by_four_raster, "values", "node", distribution="not_a_distribution")


def test_vertical_plane(four_by_four_raster):
    with pytest.raises(ValueError):
        plane(four_by_four_raster, "values", normal=(0., 1., 0.))


def test_no_xy_values(four_by_four_raster):
    with pytest.raises(ValueError):
        plane(four_by_four_raster, "values", "patch")


def test_xy_node_raster(four_by_four_raster):
    vals = plane(four_by_four_raster, "values", "node", normal=_NORMAL)
    truth = _plane_function(
        four_by_four_raster.x_of_node,
        four_by_four_raster.y_of_node,
        point=_POINT,
        normal=_NORMAL,
    )
    assert_array_equal(vals, truth)


def test_xy_cell_raster(four_by_four_raster):
    vals = plane(four_by_four_raster, "values", "cell", normal=_NORMAL)
    truth = _plane_function(
        four_by_four_raster.x_of_cell,
        four_by_four_raster.y_of_cell,
        point=_POINT,
        normal=_NORMAL,
    )
    assert_array_equal(vals, truth)


def test_xy_link_raster(four_by_four_raster):
    vals = plane(four_by_four_raster, "values", "link", normal=_NORMAL)
    truth = _plane_function(
        four_by_four_raster.x_of_link,
        four_by_four_raster.y_of_link,
        point=_POINT,
        normal=_NORMAL,
    )
    assert_array_equal(vals, truth)


def test_xy_face_raster(four_by_four_raster):
    vals = plane(four_by_four_raster, "values", "face", normal=_NORMAL)
    truth = _plane_function(
        four_by_four_raster.x_of_face,
        four_by_four_raster.y_of_face,
        point=_POINT,
        normal=_NORMAL,
    )
    assert_array_equal(vals, truth)


def test_xy_node_network(simple_network):
    vals = plane(simple_network, "values", "node", normal=_NORMAL)
    truth = _plane_function(
        simple_network.x_of_node, simple_network.y_of_node, point=_POINT, normal=_NORMAL
    )
    assert_array_equal(vals, truth)


def test_xy_face_network(simple_network):
    with pytest.raises(ValueError):
        plane(simple_network, "values", "face", normal=_NORMAL)


def test_where_status_with_patches(four_by_four_raster):
    with pytest.raises(ValueError):
        constant(four_by_four_raster, "values", "patch", where=CORE_NODE)


def test_multiple_status_node(four_by_four_raster):
    four_by_four_raster.set_closed_boundaries_at_grid_edges(True, True, False, False)
    vals = constant(
        four_by_four_raster,
        "values",
        "node",
        where=[CORE_NODE, CLOSED_BOUNDARY],
        constant=10.,
    )
    true_array = np.array(
        [0., 0., 0., 0., 0., 10., 10., 10., 0., 10., 10., 10., 10., 10., 10., 10.]
    )
    assert_array_equal(vals, true_array)


# test all grid elements


# where is size of elements
# where is list of types
# where is list but use patches (raises value error)


# x and y of multiple types of grid elements
# ask for a plane of something without X and y.
