# -*- coding: utf-8 -*-
"""
Unit tests for landlab.components.overland_flow.OverlandFlowBates

last updated: 3/14/16
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.overland_flow import OverlandFlowBates, OverlandFlow


(_SHAPE, _SPACING, _ORIGIN) = ((32, 240), (25, 25), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
        from landlab import RasterModelGrid
        grid = RasterModelGrid((32, 240), spacing = 25)
        grid.add_zeros('node', 'water__depth')
        grid.add_zeros('node', 'topographic__elevation')
        bates = OverlandFlowBates(grid, mannings_n = 0.01, h_init=0.001)
        globals().update({
            'bates': bates})


@with_setup(setup_grid)
def test_Bates_name():
    assert_equal(bates.name, 'OverlandFlowBates')


@with_setup(setup_grid)
def test_Bates_input_var_names():
    # DEJH added sets to remove reliance on ordering
    assert_equal(set(bates.input_var_names),  set(('water__depth',
                                                   'topographic__elevation')))


@with_setup(setup_grid)
def test_Bates_output_var_names():
    # DEJH added sets to remove reliance on ordering
    assert_equal(set(bates.output_var_names), set(('water__depth',
                                                   'water__discharge',
                                                   'water_surface__gradient')))

@with_setup(setup_grid)
def test_Bates_var_units():
    assert_equal(set(bates.input_var_names) |
                 set(bates.output_var_names),
                 set(dict(bates.units).keys()))

    assert_equal(bates.var_units['water__depth'], 'm')
    assert_equal(bates.var_units['water__discharge'], 'm3/s')
    assert_equal(bates.var_units['water_surface__gradient'], 'm/m')
    assert_equal(bates.var_units['topographic__elevation'], 'm')


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in bates.grid['node'].keys():
        field = bates.grid['node'][name]
        if name != 'water__depth':
            assert_true(np.all(np.isclose(field, 0.)))
        else:
            assert_true(np.all(np.isclose(field, 0.001)))
            # this remains broken, and needs JA's attention


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(bates.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(bates.grid.number_of_node_columns, _SHAPE[1])


def test_Bates_analytical():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((32, 240), spacing = 25)
    grid.add_zeros('node', 'water__depth')
    grid.add_zeros('node', 'topographic__elevation')
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    bates = OverlandFlowBates(grid, mannings_n = 0.01, h_init=0.001)
    time = 0.0
    bates.dt = 1.0
    while time < 500:
        bates.overland_flow(grid)
        h_boundary = (((7./3.) * (0.01**2) * (0.4**3) *
                  time) ** (3./7.))
        grid.at_node['water__depth'][grid.nodes[1: -1, 1]] = h_boundary
        time += bates.dt

    x = np.arange(0, ((grid.shape[1]) * grid.dx), grid.dx)
    h_analytical = (-(7./3.) * (0.01**2) * (0.4**2) * (x - (0.4 * 500)))

    h_analytical[np.where(h_analytical > 0)] = (h_analytical[np.where(
        h_analytical > 0)] ** (3./7.))
    h_analytical[np.where(h_analytical < 0)] = 0.0

    hBates = bates.h.reshape(grid.shape)
    hBates = hBates[1][1:]
    hBates = np.append(hBates, [0])
    np.testing.assert_almost_equal(h_analytical, hBates, decimal=1)