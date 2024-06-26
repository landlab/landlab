import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid.ext.raster_mappers import map_max_of_link_nodes_to_link
from landlab.grid.mappers import (
    map_max_of_link_nodes_to_link as map_max_of_link_nodes_to_link_slow,
)
from landlab.grid.raster_mappers import map_max_of_inlinks_to_node
from landlab.grid.raster_mappers import map_max_of_outlinks_to_node
from landlab.grid.raster_mappers import map_mean_of_inlinks_to_node
from landlab.grid.raster_mappers import map_mean_of_outlinks_to_node
from landlab.grid.raster_mappers import map_min_of_inlinks_to_node
from landlab.grid.raster_mappers import map_min_of_outlinks_to_node
from landlab.grid.raster_mappers import map_sum_of_inlinks_to_node
from landlab.grid.raster_mappers import map_sum_of_outlinks_to_node


@pytest.mark.parametrize("n_rows", (3, 4, 8, 512))
@pytest.mark.parametrize("n_cols", (3, 4, 8, 512))
@pytest.mark.parametrize("dtype", (np.int64, np.longlong, np.float32, np.float64))
def test_boundscheck(n_rows, n_cols, dtype):
    """Test buffer over/under flow."""
    buffer_size = 1024

    n_nodes = n_rows * n_cols
    n_links = n_rows * (n_cols - 1) + (n_rows - 1) * n_cols

    value_at_node = np.ones(n_nodes, dtype=dtype)
    actual = np.full(n_links + 2 * buffer_size, -2, dtype=dtype)

    map_max_of_link_nodes_to_link(
        actual[buffer_size:-buffer_size], value_at_node, (n_rows, n_cols)
    )

    assert np.all(actual[:buffer_size] == -2)
    assert np.all(actual[-buffer_size:] == -2)
    assert np.all(actual[buffer_size:-buffer_size] == 1)


def map_max_of_link_nodes_to_link_fast(grid, value_at_node, out=None):
    if out is None:
        out = grid.empty(at="link")

    map_max_of_link_nodes_to_link(out, value_at_node, grid.shape)

    return out


@pytest.mark.benchmark(group="map_max_of_link_nodes_to_link")
@pytest.mark.parametrize(
    "func",
    (map_max_of_link_nodes_to_link_fast, map_max_of_link_nodes_to_link_slow),
    ids=("CYTHON", "PYTHON"),
)
def test_map_max_of_link_nodes_to_link_bench(benchmark, func):
    grid = RasterModelGrid((300, 4000))

    value_at_node = np.arange(grid.number_of_nodes, dtype=int)
    out = np.full(grid.number_of_links, -999, dtype=int)

    benchmark(func, grid, value_at_node, out=out)

    assert np.all(out >= 0)
    assert np.all(out < grid.number_of_nodes)


def test_map_max_of_link_nodes_to_link_cmp():
    grid = RasterModelGrid((30, 40))

    value_at_node = np.random.uniform(low=-100, high=100, size=grid.number_of_nodes)
    expected = np.empty(grid.number_of_links, dtype=float)
    actual = np.empty(grid.number_of_links, dtype=float)

    map_max_of_link_nodes_to_link_slow(grid, value_at_node, out=expected)
    map_max_of_link_nodes_to_link_fast(grid, value_at_node, out=actual)

    assert_array_equal(actual, expected)


@pytest.mark.benchmark(group="raster_mappers")
@pytest.mark.parametrize(
    "func",
    [
        map_sum_of_inlinks_to_node,
        map_mean_of_inlinks_to_node,
        map_max_of_inlinks_to_node,
        map_min_of_inlinks_to_node,
        map_sum_of_outlinks_to_node,
        map_mean_of_outlinks_to_node,
        map_max_of_outlinks_to_node,
        map_min_of_outlinks_to_node,
    ],
)
def test_mapper_bench(benchmark, func):
    grid = RasterModelGrid((300, 4000))

    value_at_link = np.arange(grid.number_of_links, dtype=int)
    out = np.full(grid.number_of_nodes, -999, dtype=int)

    benchmark(func, grid, value_at_link, out=out)

    assert np.all(out >= 0)
