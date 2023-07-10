import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.grid.divergence import calc_flux_div_at_node as calc_flux_div_at_node_slow
from landlab.grid.raster_divergence import calc_flux_div_at_node


@pytest.mark.benchmark(group="small")
def test_flux_div_at_node_bench(benchmark):
    grid = RasterModelGrid((40, 50), xy_spacing=(1.0, 2.0))

    value_at_link = np.random.uniform(size=grid.number_of_links)
    actual = grid.empty(at="node")

    benchmark(calc_flux_div_at_node, grid, value_at_link, out=actual)


@pytest.mark.benchmark(group="large")
@pytest.mark.parametrize(
    "func",
    [calc_flux_div_at_node, calc_flux_div_at_node_slow],
    ids=["raster-specific", "general"],
)
def test_flux_div_at_node_large_bench(benchmark, func):
    grid = RasterModelGrid((400, 5000), xy_spacing=(1.0, 2.0))

    value_at_link = np.random.uniform(size=grid.number_of_links)
    actual = grid.empty(at="node")

    benchmark(func, grid, value_at_link, out=actual)


@pytest.mark.parametrize("shape", [(4, 5), (40, 50), (50, 40), (3, 3)])
@pytest.mark.parametrize("spacing", [(1.0, 3.0), (3.0, 1.0)])
def test_flux_div_at_node_matches(shape, spacing):
    grid = RasterModelGrid(shape, xy_spacing=spacing)

    value_at_link = np.random.uniform(size=grid.number_of_links)
    actual = grid.empty(at="node")

    expected = calc_flux_div_at_node_slow(grid, value_at_link)
    actual = calc_flux_div_at_node(grid, value_at_link)

    assert_array_almost_equal(actual, expected)
