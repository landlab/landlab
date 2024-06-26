import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import text
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.field import GroupError
from landlab.values import constant
from landlab.values import plane
from landlab.values import random
from landlab.values import units
from landlab.values.synthetic import _plane_function
from landlab.values.synthetic import _where_to_add_values

_NORMAL = (1, 1, 1)
_POINT = (0, 0, 0)


@given(name=text(), unit_str=text())
def test_add_units_missing_field(at, name, unit_str):
    grid = RasterModelGrid((4, 4))
    units(grid, name, at=at, units=unit_str)
    assert grid[at][name] == pytest.approx(0.0)
    assert grid.field_units(at, name) == unit_str
    assert grid[at].units[name] == unit_str


def test_add_units_existing_field(at):
    grid = RasterModelGrid((4, 4))
    grid.add_empty("x", at=at, units="NONE")
    assert grid.field_units(at, "x") == "NONE"
    values = grid[at]["x"].copy()

    units(grid, "x", at=at, units="m")
    assert grid.field_units(at, "x") == "m"
    assert_array_equal(grid[at]["x"], values)


@given(name=text())
def test_add_units_without_units(at, name):
    grid = RasterModelGrid((4, 4))
    units(grid, name, at=at, units=None)
    assert grid.field_units(at, name) == "?"

    units(grid, name, at=at)
    assert grid.field_units(at, name) == "?"


def test_bad_grid_element_name(four_by_four_raster):
    with pytest.raises(KeyError):
        constant(four_by_four_raster, "some_flux", "not_a_grid_element", value=1.0)


def test_bad_distribution_name(four_by_four_raster):
    with pytest.raises(ValueError):
        random(four_by_four_raster, "values", "node", distribution="not_a_distribution")


def test_vertical_plane(four_by_four_raster):
    with pytest.raises(ValueError):
        plane(four_by_four_raster, "values", normal=(0.0, 1.0, 0.0))


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
        four_by_four_raster.xy_of_cell[:, 0],
        four_by_four_raster.xy_of_cell[:, 1],
        point=_POINT,
        normal=_NORMAL,
    )
    assert_array_equal(vals, truth)


def test_xy_link_raster(four_by_four_raster):
    vals = plane(four_by_four_raster, "values", "link", normal=_NORMAL)
    truth = _plane_function(
        four_by_four_raster.xy_of_link[:, 0],
        four_by_four_raster.xy_of_link[:, 1],
        point=_POINT,
        normal=_NORMAL,
    )
    assert_array_equal(vals, truth)


def test_xy_face_raster(four_by_four_raster):
    vals = plane(four_by_four_raster, "values", "face", normal=_NORMAL)
    truth = _plane_function(
        four_by_four_raster.xy_of_face[:, 0],
        four_by_four_raster.xy_of_face[:, 1],
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
    with pytest.raises(AttributeError):
        constant(
            four_by_four_raster,
            "values",
            "patch",
            where=four_by_four_raster.BC_NODE_IS_CORE,
        )


def test_multiple_status_node(four_by_four_raster):
    four_by_four_raster.set_closed_boundaries_at_grid_edges(True, True, False, False)
    vals = constant(
        four_by_four_raster,
        "values",
        "node",
        where=[
            four_by_four_raster.BC_NODE_IS_CORE,
            four_by_four_raster.BC_NODE_IS_CLOSED,
        ],
        value=10.0,
    )
    true_array = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            10.0,
            10.0,
            10.0,
            0.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
            10.0,
        ]
    )
    assert_array_equal(vals, true_array)


def test_where_to_add_values_where_is_none(four_by_four_raster, at):
    where = _where_to_add_values(four_by_four_raster, at, None)
    assert np.all(where == np.full(four_by_four_raster.size(at), True))


def test_where_to_add_values_where_is_everywhere(four_by_four_raster, at):
    where = np.full(four_by_four_raster.size(at), True)
    assert np.all(_where_to_add_values(four_by_four_raster, at, where) == where)


def test_where_to_add_values_where_is_nowhere(four_by_four_raster, at):
    where = np.full(four_by_four_raster.size(at), False)
    assert np.all(_where_to_add_values(four_by_four_raster, at, where) == where)


def test_where_to_add_values_where_is_somewhere(four_by_four_raster, at):
    where = np.random.randint(0, 2, size=four_by_four_raster.size(at), dtype=bool)
    assert np.all(_where_to_add_values(four_by_four_raster, at, where) == where)


def test_where_to_add_values_node_bc(four_by_four_raster, at, node_bc):
    if at == "node":
        where = _where_to_add_values(four_by_four_raster, at, node_bc)
        assert len(where) == four_by_four_raster.number_of_nodes
    elif at == "link":
        with pytest.raises(ValueError):
            _where_to_add_values(four_by_four_raster, at, node_bc)
    else:
        with pytest.raises(AttributeError):
            _where_to_add_values(four_by_four_raster, at, node_bc)


def test_where_to_add_values_link_bc(four_by_four_raster, at, link_bc):
    if at == "link":
        where = _where_to_add_values(four_by_four_raster, at, link_bc)
        assert len(where) == four_by_four_raster.number_of_links
    elif at == "node":
        with pytest.raises(ValueError):
            _where_to_add_values(four_by_four_raster, at, link_bc)
    else:
        with pytest.raises(AttributeError):
            _where_to_add_values(four_by_four_raster, at, link_bc)


def test_where_to_add_values_wrong_size(four_by_four_raster, at):
    where = np.full(four_by_four_raster.size(at) - 1, False)
    with pytest.raises(ValueError):
        _where_to_add_values(four_by_four_raster, at, where)


def test_where_to_add_values_with_bad_at(four_by_four_raster):
    with pytest.raises(GroupError):
        _where_to_add_values(four_by_four_raster, "not-a-place", None)


# test all grid elements


# where is size of elements
# where is list of types
# where is list but use patches (raises value error)


# x and y of multiple types of grid elements
# ask for a plane of something without X and y.
