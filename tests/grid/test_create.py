from collections import OrderedDict
from io import StringIO

import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import text
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import HexModelGrid
from landlab import RadialModelGrid
from landlab import RasterModelGrid
from landlab.grid.create import _parse_args_kwargs
from landlab.grid.create import add_boundary_conditions
from landlab.grid.create import add_field_from_function
from landlab.grid.create import as_list_of_tuples
from landlab.grid.create import create_grid
from landlab.grid.create import grid_from_dict
from landlab.grid.create import norm_grid_description

SIMPLE_PARAMS_STR = """
grid:
  RasterModelGrid:
    args: [[4, 5]]
    xy_spacing: [3, 4]
    fields:
      node:
        topographic__elevation:
          plane:
            - point: [1, 1, 1]
              normal: [-2, -1, 1]
    boundary_conditions:
      - set_closed_boundaries_at_grid_edges:
        - True
        - True
        - True
        - True
"""


def test_no_grid_value():
    dict_like = {"foo": "bar"}
    with pytest.raises(ValueError):
        create_grid(dict_like)


def test_bad_grid_name():
    dict_like = {"grid": "MagicModelGrid"}
    with pytest.raises(ValueError):
        create_grid(dict_like, section="grid")


def test_two_grid_types_as_dict():
    dict_like = {
        "grid": OrderedDict([("RasterModelGrid", [(4, 5)]), ("HexModelGrid", [(6, 7)])])
    }

    # dict_like = {"grid": {"RasterModelGrid": [(4, 5)], "HexModelGrid": [6, 7]}}
    # with pytest.raises(ValueError):
    grids = create_grid(dict_like, section="grid")
    assert len(grids) == 2
    assert isinstance(grids[0], RasterModelGrid)
    assert isinstance(grids[1], HexModelGrid)


def test_two_grid_types_as_list():
    dict_like = {"grid": [{"RasterModelGrid": [(4, 5)]}, {"HexModelGrid": [(6, 7)]}]}

    grids = create_grid(dict_like, section="grid")
    assert len(grids) == 2
    assert isinstance(grids[0], RasterModelGrid)
    assert isinstance(grids[1], HexModelGrid)


def test_bad_field_location():
    dict_like = {
        "grid": {
            "RasterModelGrid": {
                "args": [(4, 5)],
                "fields": {"at_bad_location": {"foo": "bar"}},
            }
        }
    }
    with pytest.raises(ValueError):
        create_grid(dict_like, section="grid")


def test_bad_field_function():
    dict_like = {
        "grid": {
            "RasterModelGrid": [
                (4, 5),
                {
                    "fields": {
                        "node": {"new_field_name": {"not_a_function": ["bar", "spam"]}}
                    }
                },
            ]
        }
    }
    with pytest.raises(ValueError):
        create_grid(dict_like, section="grid")


def test_simple_create(tmpdir):
    """Load parameters from YAML-formatted file."""
    with tmpdir.as_cwd():
        with open("params.yaml", "w") as fp:
            fp.write(SIMPLE_PARAMS_STR)
        mg = create_grid("./params.yaml", section="grid")

    assert mg.number_of_nodes == 20
    assert "topographic__elevation" in mg.at_node

    x_of_node = np.array(
        [
            0.0,
            3.0,
            6.0,
            9.0,
            12.0,
            0.0,
            3.0,
            6.0,
            9.0,
            12.0,
            0.0,
            3.0,
            6.0,
            9.0,
            12.0,
            0.0,
            3.0,
            6.0,
            9.0,
            12.0,
        ]
    )

    status_at_node = np.array(
        [4, 4, 4, 4, 4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4], dtype=np.uint8
    )
    topographic__elevation = np.array(
        [
            [-2.0, 4.0, 10.0, 16.0, 22.0],
            [2.0, 8.0, 14.0, 20.0, 26.0],
            [6.0, 12.0, 18.0, 24.0, 30.0],
            [10.0, 16.0, 22.0, 28.0, 34.0],
        ]
    )

    assert_array_equal(mg.x_of_node, x_of_node)
    assert_array_equal(status_at_node, mg.status_at_node)
    assert_array_equal(
        topographic__elevation,
        np.round(mg.at_node["topographic__elevation"].reshape(mg.shape), decimals=2),
    )


def test_esri_ascii_create(datadir):
    filename = str(datadir / "4_x_3_no_nodata_value.asc")
    dict_like = {
        "grid": {
            "RasterModelGrid": [
                (4, 3),
                {
                    "xy_spacing": 10,
                    "xy_of_lower_left": (1, 2),
                    "fields": {
                        "node": {
                            "topographic__elevation": {"read_esri_ascii": [filename]}
                        }
                    },
                },
            ]
        }
    }
    mg = create_grid(dict_like, section="grid")

    x_of_node = np.array(
        [1.0, 11.0, 21.0, 1.0, 11.0, 21.0, 1.0, 11.0, 21.0, 1.0, 11.0, 21.0]
    )
    y_of_node = np.array(
        [2.0, 2.0, 2.0, 12.0, 12.0, 12.0, 22.0, 22.0, 22.0, 32.0, 32.0, 32.0]
    )
    status_at_node = np.array([1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1], dtype=np.uint8)
    topographic__elevation = np.array(
        [9.0, 10.0, 11.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0]
    )

    assert_array_equal(mg.x_of_node, x_of_node)
    assert_array_equal(mg.y_of_node, y_of_node)
    assert_array_equal(status_at_node, mg.status_at_node)
    assert_array_equal(topographic__elevation, mg.at_node["topographic__elevation"])


def test_read_netcdf_create(datadir):
    filename = str(datadir / "test-netcdf4.nc")
    dict_like = {
        "grid": {
            "RasterModelGrid": [
                (4, 3),
                {
                    "fields": {
                        "node": {"surface__elevation": {"read_netcdf": [filename]}}
                    }
                },
            ]
        }
    }
    mg = create_grid(dict_like, section="grid")
    x_of_node = np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0])
    y_of_node = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0])
    status_at_node = np.array([1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1], dtype=np.uint8)
    surface__elevation = np.array(
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    )

    assert_array_equal(mg.x_of_node, x_of_node)
    assert_array_equal(mg.y_of_node, y_of_node)
    assert_array_equal(status_at_node, mg.status_at_node)
    assert_array_equal(surface__elevation, mg.at_node["surface__elevation"])


def test_two_boundary_condition_dicts():
    dict_like = {
        "grid": {"RasterModelGrid": [(4, 3)]},
        "boundary_conditions": [
            {
                "set_closed_boundaries_at_grid_edges": [True, True, True, True],
                "not_a_function": [True, False],
            }
        ],
    }
    with pytest.raises(ValueError):
        create_grid(dict_like)


def test_bad_boundary_condition_functions(datadir):
    filename = str(datadir / "bad_boundary.yaml")
    with pytest.raises(ValueError):
        create_grid(filename, section="grid")


def test_parse_args_kwds():
    assert _parse_args_kwargs(["foo", {"bar": "baz"}]) == (("foo",), {"bar": "baz"})


def test_parse_args_kwds_without_kwds():
    assert _parse_args_kwargs(["foo"]) == (("foo",), {})


def test_parse_args_kwds_without_args():
    assert _parse_args_kwargs([{"bar": "baz"}]) == ((), {"bar": "baz"})


def test_parse_args_kwds_with_multiple_args():
    assert _parse_args_kwargs(["foo", "bar", {"bar": "baz"}]) == (
        ("foo", "bar"),
        {"bar": "baz"},
    )


def test_parse_args_kwds_as_dict():
    assert _parse_args_kwargs({"args": "foo", "bar": "baz"}) == _parse_args_kwargs(
        ["foo", {"bar": "baz"}]
    )
    assert _parse_args_kwargs({"args": ["foo"], "bar": "baz"}) == _parse_args_kwargs(
        {"args": "foo", "bar": "baz"}
    )
    assert _parse_args_kwargs({"args": "foo"}) == (("foo",), {})
    assert _parse_args_kwargs({"bar": "baz"}) == ((), {"bar": "baz"})


@pytest.mark.parametrize("xy_of_lower_left", ((-10.0, 100.0), None))
@pytest.mark.parametrize("xy_spacing", (10.0, (100.0, 200.0), None))
def test_grid_from_dict_raster(xy_of_lower_left, xy_spacing):
    kwds = {}
    if xy_of_lower_left is not None:
        kwds["xy_of_lower_left"] = xy_of_lower_left
    if xy_spacing is not None:
        kwds["xy_spacing"] = xy_spacing

    expected = RasterModelGrid((3, 4), **kwds)
    actual = grid_from_dict("RasterModelGrid", [(3, 4), kwds])

    assert_array_almost_equal(expected.x_of_node, actual.x_of_node)
    assert_array_almost_equal(expected.y_of_node, actual.y_of_node)


@pytest.mark.parametrize("orientation", ("horizontal", "vertical", None))
@pytest.mark.parametrize("node_layout", ("hex", "rect", None))
@pytest.mark.parametrize("xy_of_lower_left", ((-10.0, 100.0), None))
def test_grid_from_dict_hex(orientation, node_layout, xy_of_lower_left):
    kwds = {"shape": (3, 4)}
    if orientation is not None:
        kwds["orientation"] = orientation
    if node_layout is not None:
        kwds["node_layout"] = node_layout
    if xy_of_lower_left is not None:
        kwds["xy_of_lower_left"] = xy_of_lower_left

    expected = HexModelGrid(**kwds)
    actual = grid_from_dict("HexModelGrid", kwds)

    assert_array_almost_equal(expected.x_of_node, actual.x_of_node)
    assert_array_almost_equal(expected.y_of_node, actual.y_of_node)


@pytest.mark.parametrize("xy_of_center", ((-10.0, 100.0), None))
def test_grid_from_dict_radial(xy_of_center):
    kwds = {}
    if xy_of_center is not None:
        kwds["xy_of_center"] = xy_of_center
    expected = RadialModelGrid(3, 12, **kwds)
    actual = grid_from_dict("RadialModelGrid", (3, 12, kwds))

    assert_array_almost_equal(expected.x_of_node, actual.x_of_node)
    assert_array_almost_equal(expected.y_of_node, actual.y_of_node)


@pytest.mark.parametrize(
    "to_list",
    (
        lambda items: ((k, v) for k, v in items),
        lambda items: tuple((k, v) for k, v in items),
        lambda items: [[k, v] for k, v in items],
        lambda items: [{k: v} for k, v in items],
        OrderedDict,
    ),
)
def test_as_list_of_tuples(to_list):
    expected = [
        ("RasterModelGrid", "raster_data"),
        ("hex_model_grid", "hex_data"),
        ("radial_model_grid", "radial_data"),
    ]
    grid_list = to_list(expected)
    assert as_list_of_tuples(grid_list) == expected


def test_as_list_of_tuples_just_one():
    assert as_list_of_tuples(("foo", "bar")) == [("foo", "bar")]
    assert as_list_of_tuples(["foo", "bar"]) == [("foo", "bar")]
    assert as_list_of_tuples({"foo": "bar"}) == [("foo", "bar")]


def test_as_list_of_tuples_mixed_types():
    expected = [
        ("RasterModelGrid", "raster_data"),
        ("hex_model_grid", "hex_data"),
        ("radial_model_grid", "radial_data"),
    ]
    grid_list = [
        ("RasterModelGrid", "raster_data"),
        OrderedDict(
            [("hex_model_grid", "hex_data"), ("radial_model_grid", "radial_data")]
        ),
    ]
    assert as_list_of_tuples(grid_list) == expected


def test_create_grid():
    contents = StringIO(
        """
RasterModelGrid:
  args: [[3, 4]]
  xy_spacing: 2.0
  xy_of_lower_left: [1.0, 2.0]
  fields:
    node:
      elevation:
        constant:
          value: 1.0
"""
    )
    expected = RasterModelGrid((3, 4), xy_spacing=2.0, xy_of_lower_left=(1, 2))
    actual = create_grid(contents)
    assert_array_almost_equal(expected.x_of_node, actual.x_of_node)
    assert_array_almost_equal(expected.y_of_node, actual.y_of_node)


def test_create_grid_multiple():
    contents = StringIO(
        """
grids:
    - RasterModelGrid:
      - [3, 4]
      - xy_spacing: 2.0
        xy_of_lower_left: [1.0, 2.0]
    - RasterModelGrid:
        args: [[3, 4]]
        xy_spacing: 2.0
        xy_of_lower_left: [1.0, 2.0]
    - - RasterModelGrid
      - args: [[3, 4]]
        xy_spacing: 2.0
        xy_of_lower_left: [1.0, 2.0]
"""
    )
    expected = RasterModelGrid((3, 4), xy_spacing=2.0, xy_of_lower_left=(1, 2))
    grids = create_grid(contents, section="grids")
    assert len(grids) == 3
    for actual in grids:
        assert_array_almost_equal(expected.x_of_node, actual.x_of_node)
        assert_array_almost_equal(expected.y_of_node, actual.y_of_node)


@pytest.mark.parametrize("to_list", (list, tuple))
@pytest.mark.parametrize("at", ("node", "link", "patch", "corner", "face", "cell"))
def test_add_field_from_function_constant(at, to_list):
    grid = RasterModelGrid((20, 30))
    function = to_list(("constant", {"value": 1.0}))
    add_field_from_function(grid, "z", function, at=at)
    assert_array_almost_equal(grid[at]["z"], 1.0)


@pytest.mark.parametrize("at", ("node", "link", "patch", "corner", "face", "cell"))
def test_add_field_from_function_constant_as_dict(at):
    grid = RasterModelGrid((20, 30))
    add_field_from_function(grid, "z", {"constant": {"value": 10.0}}, at=at)
    assert_array_almost_equal(grid[at]["z"], 10.0)


@pytest.mark.parametrize(
    "to_list",
    (
        lambda items: ((k, v) for k, v in items),
        lambda items: tuple((k, v) for k, v in items),
        lambda items: [[k, v] for k, v in items],
        lambda items: [{k: v} for k, v in items],
    ),
)
@pytest.mark.parametrize("at", ("node", "link", "patch", "corner", "face", "cell"))
def test_add_field_from_function_with_list(at, to_list):
    grid = RasterModelGrid((20, 30))
    functions = to_list(
        [
            ("constant", {"value": 1.0}),
            ("constant", {"value": 2.0}),
            ("constant", {"value": 3.0}),
            ("constant", {"value": 4.0}),
            ("constant", {"value": 5.0}),
        ]
    )
    add_field_from_function(grid, "z", functions, at=at)
    assert_array_almost_equal(grid[at]["z"], 15.0)


@pytest.mark.parametrize(
    "to_list",
    (
        lambda items: ((k, v) for k, v in items),
        lambda items: tuple((k, v) for k, v in items),
        lambda items: [[k, v] for k, v in items],
        lambda items: [{k: v} for k, v in items],
    ),
)
def test_add_boundary_conditions(to_list):
    boundary_conditions = to_list(
        [("set_closed_boundaries_at_grid_edges", (True, True, True, True))]
    )
    grid = RasterModelGrid((3, 4))
    add_boundary_conditions(grid, boundary_conditions=boundary_conditions)
    assert_array_equal(
        grid.status_at_node,
        [
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CORE,
            grid.BC_NODE_IS_CORE,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
            grid.BC_NODE_IS_CLOSED,
        ],
    )


def test_norm_grid_description():
    expected = {"args": ["arg0", "arg1"]}
    assert norm_grid_description(["arg0", "arg1"]) == expected
    assert norm_grid_description(expected) == expected


@given(units=text())
def test_field_units(units):
    params = {
        "RasterModelGrid": [
            {
                "shape": (5, 7),
                "xy_spacing": (10, 10),
            },
            {
                "fields": {
                    "node": {
                        "topographic__elevation": {
                            "constant": [{"value": 0.0}],
                            "units": {"units": units},
                        },
                    }
                },
            },
        ],
    }
    grid = create_grid(params)
    assert grid.at_node["topographic__elevation"] == pytest.approx(0.0)
    assert grid.at_node.units["topographic__elevation"] == units
    assert grid.field_units("node", "topographic__elevation") == units


@given(units=text())
def test_repeated_units(units):
    params = {
        "RasterModelGrid": [
            {
                "shape": (5, 7),
                "xy_spacing": (10, 10),
            },
            {
                "fields": {
                    "node": {
                        "topographic__elevation": [
                            ("units", {"units": units[::-1]}),
                            ("constant", [{"value": 0.0}]),
                            ("units", {"units": units}),
                        ],
                    }
                },
            },
        ],
    }
    grid = create_grid(params)
    assert grid.at_node["topographic__elevation"] == pytest.approx(0.0)
    assert grid.at_node.units["topographic__elevation"] == units
    assert grid.field_units("node", "topographic__elevation") == units
