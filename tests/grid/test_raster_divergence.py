import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.grid.divergence import calc_flux_div_at_node
from landlab.grid.ext.raster_divergence import (
    calc_flux_div_at_node as _calc_flux_div_at_node_c,
)


def calc_flux_div_at_node_c(grid, value_at_link, out=None):
    _calc_flux_div_at_node_c(grid.shape, (grid.dx, grid.dy), value_at_link, out)
    return out


@pytest.mark.benchmark(group="small")
@pytest.mark.parametrize("func", [calc_flux_div_at_node, calc_flux_div_at_node_c])
def test_flux_div_at_node_bench(benchmark, func):
    grid = RasterModelGrid((4, 5), xy_spacing=(1.0, 2.0))

    value_at_link = np.random.uniform(size=grid.number_of_links)
    actual = grid.empty(at="node")

    actual.reshape(grid.shape)[:, (0, -1)] = 0.0
    actual.reshape(grid.shape)[(0, -1), :] = 0.0

    expected = grid.calc_flux_div_at_node(value_at_link)

    benchmark(func, grid, value_at_link, out=actual)

    assert_array_almost_equal(actual, expected)


@pytest.mark.benchmark(group="large")
@pytest.mark.parametrize("func", [calc_flux_div_at_node, calc_flux_div_at_node_c])
def test_flux_div_at_node_large_bench(benchmark, func):
    grid = RasterModelGrid((400, 5000), xy_spacing=(1.0, 2.0))

    value_at_link = np.random.uniform(size=grid.number_of_links)
    actual = grid.empty(at="node")

    actual.reshape(grid.shape)[:, (0, -1)] = 0.0
    actual.reshape(grid.shape)[(0, -1), :] = 0.0

    expected = grid.calc_flux_div_at_node(value_at_link)

    benchmark(func, grid, value_at_link, out=actual)

    assert_array_almost_equal(actual, expected)
