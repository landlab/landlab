import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components.overland_flow._neighbors_at_link import sum_parallel_links
from landlab.graph.structured_quad.ext.at_link import fill_parallel_links_at_link


def test_sum_parallel_links_bench(benchmark):
    shape = (400, 5000)
    grid = RasterModelGrid(shape)
    value_at_link = np.arange(grid.number_of_links, dtype=int)
    actual = np.full_like(value_at_link, -999)

    benchmark(sum_parallel_links, actual, value_at_link, shape)

    sum_at_horizontal_links = actual[grid.horizontal_links].reshape((shape[0], -1))
    sum_at_vertical_links = actual[grid.vertical_links].reshape((-1, shape[1]))

    assert np.all(sum_at_horizontal_links[:, (0, -1)] == -999)
    assert np.all(sum_at_horizontal_links[:, 1:-1] != -999)
    assert np.all(sum_at_vertical_links[(0, -1), :] == -999)
    assert np.all(sum_at_vertical_links[1:-1, :] != -999)


def test_sum_parallel_links():
    shape = (4, 5)
    grid = RasterModelGrid(shape)
    value_at_link = np.arange(grid.number_of_links, dtype=int)
    actual = np.full_like(value_at_link, -999)

    sum_parallel_links(actual, value_at_link, shape)

    assert_array_equal(
        actual[grid.horizontal_links].reshape((shape[0], -1)),
        [
            [-999, 2, 4, -999],
            [-999, 20, 22, -999],
            [-999, 38, 40, -999],
            [-999, 56, 58, -999],
        ],
    )
    assert_array_equal(
        actual[grid.vertical_links].reshape((-1, shape[1])),
        [
            [-999, -999, -999, -999, -999],
            [26, 28, 30, 32, 34],
            [-999, -999, -999, -999, -999],
        ],
    )


@pytest.mark.benchmark(group="small")
def test_fill_parallel_links(benchmark):
    grid = RasterModelGrid((3, 4))

    parallel_links_at_link = np.full((grid.number_of_links, 2), -999, dtype=int)

    benchmark(fill_parallel_links_at_link, grid.shape, parallel_links_at_link)

    assert not np.any(parallel_links_at_link == -999)
    assert_array_equal(
        parallel_links_at_link[grid.horizontal_links],
        [
            [-1, 1],
            [0, 2],
            [1, -1],
            [-1, 8],
            [7, 9],
            [8, -1],
            [-1, 15],
            [14, 16],
            [15, -1],
        ],
    )
    assert_array_equal(
        parallel_links_at_link[grid.vertical_links],
        [
            [-1, 10],
            [-1, 11],
            [-1, 12],
            [-1, 13],
            [3, -1],
            [4, -1],
            [5, -1],
            [6, -1],
        ],
    )


@pytest.mark.benchmark(group="large")
def test_fill_parallel_links_speed_c(benchmark):
    grid = RasterModelGrid((300, 4000))

    parallel_links_at_link = np.full((grid.number_of_links, 2), -999, dtype=int)

    benchmark(fill_parallel_links_at_link, grid.shape, parallel_links_at_link)

    assert not np.any(parallel_links_at_link == -999)


def _parallel_links_at_link(grid, out):
    out[grid.vertical_links, 0] = grid.links_at_node[
        grid.node_at_link_tail[grid.vertical_links], 3
    ]
    out[grid.vertical_links, 1] = grid.links_at_node[
        grid.node_at_link_head[grid.vertical_links], 1
    ]
    out[grid.horizontal_links, 0] = grid.links_at_node[
        grid.node_at_link_tail[grid.horizontal_links], 2
    ]
    out[grid.horizontal_links, 1] = grid.links_at_node[
        grid.node_at_link_head[grid.horizontal_links], 0
    ]


@pytest.mark.benchmark(group="large")
def test_fill_parallel_links_python(benchmark):
    grid = RasterModelGrid((300, 4000))

    parallel_links_at_link = np.full((grid.number_of_links, 2), -999, dtype=int)

    benchmark(_parallel_links_at_link, grid, parallel_links_at_link)

    assert not np.any(parallel_links_at_link == -999)
