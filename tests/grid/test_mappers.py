import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.grid.mappers import map_upwind_node_link_max_to_node


def test_map_upwind_node_link_max_to_node():
    grid = RasterModelGrid((3, 4))

    grid.at_link["grad"] = (
        [-1.0, -2.0, -1.0]
        + [0.0, 0.0, 0.0, 0.0]
        + [-1.0, -2.0, -1.0]
        + [0.0, 0.0, 0.0, 0.0]
        + [-1.0, -2.0, -1.0]
    )

    expected = [0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0]

    assert map_upwind_node_link_max_to_node(grid, "grad") == pytest.approx(expected)

    out = grid.empty(at="node")
    rtn = map_upwind_node_link_max_to_node(grid, "grad", out=out)
    assert rtn is out
    assert out == pytest.approx(expected)


def test_map_vectors_to_links_with_hex():
    C30 = np.cos(np.radians(30))

    hmg = HexModelGrid((3, 2))
    link_vecs = hmg.map_vectors_to_links(1.0, 0.0)
    expected = [1.0, -0.5, 0.5, -0.5, 0.5, 1.0, 1.0, 0.5, -0.5, 0.5, -0.5, 1.0]
    assert_array_almost_equal(link_vecs, expected)

    link_vecs = hmg.map_vectors_to_links(0.0, 1.0)
    expected = [0.0, C30, C30, C30, C30, 0.0, 0.0, C30, C30, C30, C30, 0.0]
    assert_array_almost_equal(link_vecs, expected)

    hmg = HexModelGrid((2, 3), orientation="vertical")
    link_vecs = hmg.map_vectors_to_links(1.0, 0.0)
    expected = [C30, C30, 0.0, C30, C30, 0.0, 0.0, C30, C30, 0.0, C30, C30]
    assert_array_almost_equal(link_vecs, expected)

    link_vecs = hmg.map_vectors_to_links(0.0, 1.0)
    expected = [-0.5, 0.5, 1.0, 0.5, -0.5, 1.0, 1.0, -0.5, 0.5, 1.0, 0.5, -0.5]
    assert_array_almost_equal(link_vecs, expected)


def test_map_vectors_to_links_with_raster():
    rmg = RasterModelGrid((3, 3))

    link_vecs = rmg.map_vectors_to_links(1.0, 1.0)
    expected = np.ones(rmg.number_of_links)
    assert_array_almost_equal(link_vecs, expected)

    link_vecs = rmg.map_vectors_to_links(-1.0, 1.0)
    expected[rmg.horizontal_links] = -1.0
    assert_array_almost_equal(link_vecs, expected)

    link_vecs = rmg.map_vectors_to_links(-1.0, -1.0)
    expected[:] = -1.0
    assert_array_almost_equal(link_vecs, expected)

    link_vecs = rmg.map_vectors_to_links(1.0, -1.0)
    expected[rmg.horizontal_links] = 1.0
    assert_array_almost_equal(link_vecs, expected)
