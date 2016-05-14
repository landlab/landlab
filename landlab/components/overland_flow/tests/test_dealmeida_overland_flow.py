# -*- coding: utf-8 -*-
"""
Unit tests for landlab.components.overland_flow.OverlandFlow

last updated: 3/14/16
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.overland_flow import OverlandFlow
from landlab.grid.structured_quad.links import left_edge_horizontal_ids

(_SHAPE, _SPACING, _ORIGIN) = ((32, 240), (25, 25), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
        from landlab import RasterModelGrid
        grid = RasterModelGrid((32, 240), spacing = 25)
        grid.add_zeros('node', 'water__depth')
        grid.add_zeros('node', 'topographic__elevation')
        grid.add_zeros('water__discharge', at='active_link')
        deAlm = OverlandFlow(grid, mannings_n = 0.01, h_init=0.001)
        globals().update({
            'deAlm': OverlandFlow(grid)})


@with_setup(setup_grid)
def test_deAlm_name():
    assert_equal(deAlm.name, 'OverlandFlow')


@with_setup(setup_grid)
def test_deAlm_input_var_names():
    assert_equal(deAlm.input_var_names,  ('water__depth',
                                          'topographic__elevation',))


@with_setup(setup_grid)
def test_deAlm_output_var_names():
    assert_equal(deAlm.output_var_names, ('water__depth', 'water__discharge',
                                         'water_surface__gradient', ))

@with_setup(setup_grid)
def test_deAlm_var_units():
    assert_equal(set(deAlm.input_var_names) |
                 set(deAlm.output_var_names),
                 set(dict(deAlm.units).keys()))

    assert_equal(deAlm.var_units('water__depth'), 'm')
    assert_equal(deAlm.var_units('water__discharge'), 'm3/s')
    assert_equal(deAlm.var_units('water_surface__gradient'),'-')
    assert_equal(deAlm.var_units('topographic__elevation'), 'm')


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(deAlm.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(deAlm.grid.number_of_node_columns, _SHAPE[1])


def test_deAlm_analytical():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((32, 240), spacing = 25)
    grid.add_zeros('node', 'water__depth')
    grid.add_zeros('node', 'topographic__elevation')
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    left_inactive_ids = left_edge_horizontal_ids(grid.shape)
    deAlm = OverlandFlow(grid, mannings_n = 0.01, h_init=0.001)
    time = 0.0

    while time < 500:
        grid['link']['water__discharge'][left_inactive_ids] =   (grid['link'][
        'water__discharge'][left_inactive_ids + 1])
        dt = deAlm.calc_time_step()
        deAlm.overland_flow(dt)
        h_boundary = (((7./3.) * (0.01**2) * (0.4**3) *
                  time) ** (3./7.))
        grid.at_node['water__depth'][grid.nodes[1: -1, 1]] = h_boundary
        time += dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = (-(7./3.) * (0.01**2) * (0.4**2) * (x - (0.4 * 500)))

    h_analytical[np.where(h_analytical > 0)] = (h_analytical[np.where(
        h_analytical > 0)] ** (3./7.))
    h_analytical[np.where(h_analytical < 0)] = 0.0

    hdeAlm = deAlm.h.reshape(grid.shape)
    hdeAlm = hdeAlm[1][1:]
    hdeAlm = np.append(hdeAlm, [0])
    np.testing.assert_almost_equal(h_analytical, hdeAlm, decimal=1)