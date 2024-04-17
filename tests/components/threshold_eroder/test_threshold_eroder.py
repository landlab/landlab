#! /usr/bin/env python
"""Unit tests for landlab.components.threshold_eroder.py"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import FlowDirectorSteepest
from landlab.components import PriorityFloodFlowRouter
from landlab.components import ThresholdEroder

try:
    PriorityFloodFlowRouter.load_richdem()
except ModuleNotFoundError:
    with_richdem = False
else:
    with_richdem = True


@pytest.mark.skipif(not with_richdem, reason="richdem is not installed")
def test_topography_rasterGrid():
    # %%
    mg = RasterModelGrid((5, 5))
    mg.set_closed_boundaries_at_grid_edges(False, False, False, False)
    z = np.array(
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 10.0, 1.0, 0.0]
        + [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    )
    _ = mg.add_field("topographic__elevation", z, at="node")

    # Instantiate Flow director (steepest slope type) and Threshold Eroder
    fdir = PriorityFloodFlowRouter(mg)
    th_ero = ThresholdEroder(mg, slope_crit=0.6)

    # Run the components for ten short timepsteps
    for _t in range(2):
        fdir.run_one_step()
        th_ero.run_one_step()

    # Check final topography

    assert_array_almost_equal(
        mg.at_node["topographic__elevation"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.6, 0.0, 0.0, 0.6, 1.2, 0.6]
            + [0.0, 0.0, 0.6, 0.6, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0]
        ),
        verbose=True,
    )


@pytest.mark.skipif(not with_richdem, reason="richdem is not installed")
def test_topo_soil_rasterGrid():
    # %%
    mg = RasterModelGrid((5, 5))
    mg.set_closed_boundaries_at_grid_edges(False, False, False, False)
    z = np.array(
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0]
        + [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    )
    topo = mg.add_zeros("topographic__elevation", at="node")
    bed = mg.add_field("bedrock__elevation", z, at="node")
    soil = mg.add_ones("soil__depth", at="node")
    soil[mg.boundary_nodes] = 0
    topo[:] = soil + bed

    # Instantiate Flow director (steepest slope type) and Threshold Eroder
    fdir = PriorityFloodFlowRouter(mg)
    th_ero = ThresholdEroder(mg, slope_crit=0.6)

    # Run the components for ten short timepsteps

    for _t in range(2):
        fdir.run_one_step()
        th_ero.run_one_step()

    # Check final topography
    assert_array_almost_equal(
        mg.at_node["bedrock__elevation"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.6, 0.0, 0.0, 0.6, 1, 0.6]
            + [0.0, 0.0, 0.6, 0.6, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0]
        ),
        verbose=True,
    )
    assert_array_almost_equal(
        mg.at_node["soil__depth"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0]
            + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
        verbose=True,
    )


# %%
def test_topography_hexGrid():
    # %%
    hmg = HexModelGrid((8, 8))

    topo = hmg.add_zeros("topographic__elevation", at="node")
    topo[hmg.core_nodes] += 100

    # Instantiate Flow director (steepest slope type) and Threshold Eroder
    fdir = FlowDirectorSteepest(hmg)
    fa = FlowAccumulator(hmg)
    th_ero = ThresholdEroder(hmg, slope_crit=0.6)

    # Run the components for ten short timepsteps
    for _t in range(5):
        fdir.run_one_step()
        fa.run_one_step()
        th_ero.run_one_step()
    assert_array_almost_equal(
        hmg.at_node["topographic__elevation"],
        np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.6, 0.6, 0.6]
            + [0.6, 0.6, 0.0, 0.0, 0.6, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 0.6, 0.0, 0.0]
            + [0.6, 1.2, 1.8, 1.8, 1.8, 1.8, 1.8, 1.2, 0.6, 0.0, 0.0, 0.6, 1.2, 1.8]
            + [1.8, 1.8, 1.8, 1.8, 1.8, 1.2, 0.6, 0.0, 0.0, 0.6, 1.2, 1.2, 1.2, 1.2]
            + [1.2, 1.2, 1.2, 0.6, 0.0, 0.0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
            + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
    )
