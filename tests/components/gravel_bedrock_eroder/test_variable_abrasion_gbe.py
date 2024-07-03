import numpy as np
import pytest
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_almost_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import GravelBedrockEroder


def test_sediment_thickness():
    grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[3:] = 10.0
    sed = grid.add_zeros("soil__depth", at="node")
    sed[3:] = 10.0
    fa = FlowAccumulator(grid)
    fa.run_one_step()

    eroder = GravelBedrockEroder(
        grid,
        number_of_sediment_classes=3,
        abrasion_coefficients=[0.002, 0.0002, 0.00002],
    )

    sum_of_classes = np.sum(eroder._thickness_by_class, axis=0)

    assert_array_equal(sum_of_classes, sed)

    eroder.run_one_step(100.0)

    sum_of_classes = np.sum(eroder._thickness_by_class, axis=0)

    assert_array_almost_equal(sum_of_classes, sed)


def test_value_error():
    grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    with pytest.raises(ValueError):
        GravelBedrockEroder(
            grid, number_of_sediment_classes=3, abrasion_coefficients=[0, 0]
        )

def stats(er):
    print("z", er._elev)
    print("dzdt", er._grid.at_node["topographic__steepest_slope"])
    print("s", er._sed)
    print("sbc", er._thickness_by_class)
    try:
        print("sf", er._sediment_fraction)
    except:
        pass
    print("r", er._bedrock__elevation)
    print("dHdt", er._dHdt)
    print("dHdt-b-c", er._dHdt_by_class)
    print("ref", er._rock_exposure_fraction)
    print("rlr", er._rock_lowering_rate)

def test_response_with_one_vs_two_classes():
    """Compare changes in one step for a single class vs two identical classes."""
    grid1 = RasterModelGrid((3, 3), xy_spacing=1000.0)
    elev1 = grid1.add_zeros("topographic__elevation", at="node")
    elev1[4] = 10.0
    sed1 = grid1.add_zeros("soil__depth", at="node")
    sed1[4] = 1.0
    rock1 = grid1.add_zeros("bedrock__elevation", at="node")
    rock1[4] = elev1[4] - sed1[4]
    grid1.status_at_node[grid1.perimeter_nodes] = grid1.BC_NODE_IS_CLOSED
    grid1.status_at_node[5] = grid1.BC_NODE_IS_FIXED_VALUE
    fa1 = FlowAccumulator(grid1, runoff_rate=10.0)
    fa1.run_one_step()
    eroder1 = GravelBedrockEroder(
        grid1, sediment_porosity=0.0,
        number_of_sediment_classes=1,
        abrasion_coefficients=[0.0005],
        coarse_fractions_from_plucking=[1.0]
    )
    eroder1.run_one_step(1000.0)

    grid2 = RasterModelGrid((3, 3), xy_spacing=1000.0)
    elev2 = grid2.add_zeros("topographic__elevation", at="node")
    elev2[4] = 10.0
    sed2 = grid2.add_zeros("soil__depth", at="node")
    sed2[4] = 1.0
    rock2 = grid2.add_zeros("bedrock__elevation", at="node")
    rock2[4] = elev2[4] - sed2[4]
    grid2.status_at_node[grid2.perimeter_nodes] = grid2.BC_NODE_IS_CLOSED
    grid2.status_at_node[5] = grid2.BC_NODE_IS_FIXED_VALUE
    fa2 = FlowAccumulator(grid2, runoff_rate=10.0)
    fa2.run_one_step()
    eroder2 = GravelBedrockEroder(
        grid2, sediment_porosity=0.0,
        number_of_sediment_classes=2,
        abrasion_coefficients=[0.0005, 0.0005],
        init_fraction_per_class=[0.5, 0.5],
        coarse_fractions_from_plucking=[0.5, 0.5]
    )
    eroder2.run_one_step(1000.0)

    assert_array_equal(eroder1._elev, eroder2._elev)
    assert_array_equal(eroder1._sed, eroder2._sed)
    assert_array_equal(eroder1._sed, np.sum(eroder2._thickness_by_class, axis=0))
    assert_array_equal(eroder1._dHdt, eroder2._dHdt)
    assert_array_equal(eroder1._dHdt, np.sum(eroder2._dHdt_by_class, axis=0))
    assert_array_equal(eroder1._rock_lowering_rate, eroder2._rock_lowering_rate)

def test_equilibrium_with_two_identical_classes():
    """
    Equilibrium should be the same with 2 classes of identical sediment,
    as compared to a single class of the same sediment.
    """
    grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[4] = 1.0
    sed = grid.add_zeros("soil__depth", at="node")
    sed[4] = 1.0e4
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
    fa = FlowAccumulator(grid, runoff_rate=10.0)
    fa.run_one_step()
    eroder = GravelBedrockEroder(
        grid, sediment_porosity=0.0,
        number_of_sediment_classes=2,
        abrasion_coefficients=[0.0005, 0.0005],
        init_fraction_per_class=[0.5, 0.5],
        coarse_fractions_from_plucking=[0.5, 0.5]
    )
    rock_elev = grid.at_node["bedrock__elevation"]
    fa.run_one_step()
    dt = 10000.0
    for _ in range(200):
        rock_elev[grid.core_nodes] += 1.0e-4 * dt
        elev[grid.core_nodes] += 1.0e-4 * dt
        eroder.run_one_step(dt)
    assert_almost_equal(elev[4], 33.2, decimal=2)


if __name__ == "__main__":
    test_response_with_one_vs_two_classes()
    test_equilibrium_with_two_identical_classes()