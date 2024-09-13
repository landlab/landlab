#!/usr/bin/env python2
"""
Created on Thu Jul 27 14:23:25 2017

@author: gtucker
"""
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import SharedStreamPower


@pytest.mark.parametrize(
    "k_bedrock,v_s", [(0.002, 1.0), (1.0, 1000.0), (0.001, 0.0001), (0.002, 0.002)]
)
@pytest.mark.parametrize("m_sp,n_sp", [(0.5, 1.0), (1.0 / 3.0, 2.0 / 3.0)])
def test_shared_stram_power_steady_state(k_bedrock, v_s, m_sp, n_sp):
    """Test steady state run."""
    grid = RasterModelGrid((5, 5))
    grid.at_node["topographic__elevation"] = 0.01 * grid.x_of_node

    fa = FlowAccumulator(grid, flow_director="FlowDirectorD8")

    ed = SharedStreamPower(
        grid,
        k_bedrock=k_bedrock,
        k_transport=k_bedrock / v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        solver="adaptive",
    )

    # run it to steady state.
    uplift = 0.001
    dt = 10.0
    for _ in range(3000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        grid.at_node["topographic__elevation"][grid.core_nodes] += uplift * dt

    sa_factor = ((1.0 + ed.v_s) * uplift / ed.K) ** (1.0 / ed.n_sp)

    actual = grid.at_node["topographic__steepest_slope"][grid.core_nodes]
    expected = sa_factor * grid.at_node["drainage_area"][grid.core_nodes] ** -(
        ed.m_sp / ed.n_sp
    )

    assert_array_almost_equal(actual, expected)


def test_erodep_slope_area_with_threshold():
    """Test steady state run with Vs = 1 and wc = 0.00001."""
    grid = RasterModelGrid((5, 5))
    grid.at_node["topographic__elevation"] = 0.01 * grid.x_of_node

    fa = FlowAccumulator(grid, flow_director="FlowDirectorD8")

    ed = SharedStreamPower(
        grid,
        k_bedrock=0.002,
        k_transport=0.002 / 1.0,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit=0.0001,
        solver="adaptive",
    )

    # run it to steady state.
    uplift = 0.001
    dt = 10.0
    for _ in range(3000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        grid.at_node["topographic__elevation"][grid.core_nodes] += uplift * dt

    sa_factor = ((1.0 + ed.v_s) * uplift + ed.sp_crit) / ed.K

    actual = grid.at_node["topographic__steepest_slope"][grid.core_nodes]
    expected = sa_factor * grid.at_node["drainage_area"][grid.core_nodes] ** -(
        ed.m_sp / ed.n_sp
    )

    assert_array_almost_equal(actual, expected)
