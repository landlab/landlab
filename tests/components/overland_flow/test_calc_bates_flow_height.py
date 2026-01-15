import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components.overland_flow._calc import calc_bates_flow_height


def reference_bates_flow_height(z_at_node, h_at_node, nodes_at_link, where, out):
    z_at_node = np.asarray(z_at_node).ravel()
    h_at_node = np.asarray(h_at_node).ravel()
    nodes_at_link = np.asarray(nodes_at_link)
    where = np.asarray(where).ravel()

    for link in where:
        tail = nodes_at_link[link, 0]
        head = nodes_at_link[link, 1]

        out[link] = np.maximum(
            h_at_node[tail] + z_at_node[tail], h_at_node[head] + z_at_node[head]
        ) - np.maximum(z_at_node[tail], z_at_node[head])

    return out


@pytest.mark.parametrize("z0", (-1.0, 0.0, 1.0))
def test_bates_flow_height_with_const_z(z0):
    grid = RasterModelGrid((3, 4))

    h = np.asarray(
        [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [8, 9, 10, 11],
        ],
        dtype=float,
    ).ravel()
    z = np.full_like(h, z0)

    where = np.arange(grid.number_of_links, dtype=grid.nodes_at_link.dtype)
    actual = np.empty(grid.number_of_links, dtype=z.dtype)

    calc_bates_flow_height(
        z, h, nodes_at_link=grid.nodes_at_link, where=where, out=actual
    )

    assert_array_equal(
        actual,
        [
            *[1, 2, 3],
            *[4, 5, 6, 7],
            *[5, 6, 7],
            *[8, 9, 10, 11],
            *[9, 10, 11],
        ],
    )


@pytest.mark.parametrize("h0", (-1.0, 0.0, 1.0))
def test_bates_flow_height_with_const_h(h0):
    grid = RasterModelGrid((3, 4))

    z = np.asarray(
        [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [8, 9, 10, 11],
        ],
        dtype=float,
    ).ravel()
    h = np.full_like(z, h0)

    where = np.arange(grid.number_of_links, dtype=grid.nodes_at_link.dtype)
    actual = np.empty(grid.number_of_links, dtype=z.dtype)

    calc_bates_flow_height(
        z, h, nodes_at_link=grid.nodes_at_link, where=where, out=actual
    )

    assert_array_equal(actual, h0)


@pytest.mark.parametrize(
    "here",
    (
        [
            *[0, 0, 0],
            *[0, 0, 0, 0],
            *[0, 0, 0],
            *[0, 0, 0, 0],
            *[0, 0, 0],
        ],
        [
            *[1, 1, 1],
            *[1, 1, 1, 1],
            *[1, 1, 1],
            *[1, 1, 1, 1],
            *[1, 1, 1],
        ],
        [
            *[1, 0, 0],
            *[1, 1, 1, 0],
            *[0, 0, 1],
            *[0, 0, 0, 0],
            *[1, 0, 0],
        ],
    ),
)
@pytest.mark.parametrize("values_dtype", (np.float32, np.float64))
@pytest.mark.parametrize("where_dtype", (np.int32, np.int64))
def test_bates_flow_height_where(where_dtype, values_dtype, here):
    grid = RasterModelGrid((3, 4))
    here = np.asarray(here)
    where = np.flatnonzero(here)
    actual = np.full(grid.number_of_links, -1, dtype=values_dtype)
    expected = actual.copy()

    z = np.arange(12)
    h = z[::-1] / 16

    calc_bates_flow_height(
        np.asarray(z, dtype=actual.dtype),
        np.asarray(h, dtype=actual.dtype),
        nodes_at_link=grid.nodes_at_link.astype(where_dtype),
        where=where.astype(where_dtype),
        out=actual,
    )

    assert np.all(here) or np.all(actual[~here] == -1)
    assert np.all(~here) or np.all(actual[here] >= 0)

    assert_array_equal(
        actual, reference_bates_flow_height(z, h, grid.nodes_at_link, where, expected)
    )


def test_bates_flow_height_return_value():
    grid = RasterModelGrid((3, 4))
    where = np.arange(grid.number_of_links)

    z = np.arange(grid.number_of_nodes, dtype=float)
    h = z[::-1] / 16

    actual = np.empty(grid.number_of_links, dtype=z.dtype)
    rtn = calc_bates_flow_height(
        z,
        h,
        nodes_at_link=grid.nodes_at_link,
        where=where.astype(grid.nodes_at_link.dtype),
        out=actual,
    )
    assert_array_equal(rtn, actual)
    assert rtn is actual
