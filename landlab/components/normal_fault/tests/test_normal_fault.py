

import numpy as np
from numpy.testing import assert_array_equal

from nose.tools import assert_true

from landlab import RasterModelGrid, HexModelGrid
from landlab.components import NormalFault


def test_dx_equals_zero():
    """Test a vertical fault trace."""
    grid = RasterModelGrid((6, 6), spacing=10)

    _ = grid.add_zeros('node', 'topographic__elevation')

    param_dict = {'faulted_surface': 'topographic__elevation',
                  'fault_dip_angle': 90.0,
                  'fault_throw_rate_through_time': {'time': [0, 9, 10],
                                                    'rate': [0, 0, 0.05]},
                  'fault_trace': {'y1': 0,
                                  'x1': 30,
                                  'y2': 30,
                                  'x2': 30},
                  'include_boundaries': True}

    nf = NormalFault(grid, **param_dict)

    out = np.array([[ True,  True,  True, False, False, False],
                    [ True,  True,  True, False, False, False],
                    [ True,  True,  True, False, False, False],
                    [ True,  True,  True, False, False, False],
                    [ True,  True,  True, False, False, False],
                    [ True,  True,  True, False, False, False]], dtype=bool)

    assert_array_equal(nf.faulted_nodes.reshape(grid.shape), out)


def test_anti_aximuth_greq_2pi():
    """Test anti azimuth over 2*pi."""

    grid = RasterModelGrid((6, 6), spacing=10)

    _ = grid.add_zeros('node', 'topographic__elevation')

    param_dict = {'faulted_surface': 'topographic__elevation',
                  'fault_dip_angle': 90.0,
                  'fault_throw_rate_through_time': {'time': [0, 9, 10],
                                                    'rate': [0, 0, 0.05]},
                  'fault_trace': {'y1': 30.0,
                                  'x1': 30.0,
                                  'y2': 20.0,
                                  'x2': 0.0},
                  'include_boundaries': True}

    nf = NormalFault(grid, **param_dict)


    assert_true(nf.fault_anti_azimuth > 2.0*np.pi)

    out = np.array([[ True,  True,  True,  True,  True,  True],
                    [ True,  True,  True,  True,  True,  True],
                    [ True,  True,  True,  True,  True,  True],
                    [False, False, False, False,  True,  True],
                    [False, False, False, False, False, False],
                    [False, False, False, False, False, False]], dtype=bool)

    assert_array_equal(nf.faulted_nodes.reshape(grid.shape), out)


def test_non_raster():
    """Test a hex model grid."""

    grid = HexModelGrid(7, 3, dx=10)

    _ = grid.add_zeros('node', 'topographic__elevation')

    param_dict = {'faulted_surface': 'topographic__elevation',
                  'fault_dip_angle': 90.0,
                  'fault_throw_rate_through_time': {'time': [0, 9, 10],
                                                    'rate': [0, 0, 0.05]},
                  'fault_trace': {'y1': 30.0,
                                  'x1': 30.0,
                                  'y2': 20.0,
                                  'x2': 0.0},
                  'include_boundaries': True}

    nf = NormalFault(grid, **param_dict)

    # plotting, to test this. it works!
    #import matplotlib.pyplot as plt
    #plt.figure()
    #imshow_grid(grid, nf.faulted_nodes, color_for_background='y')
    #plt.plot(grid.x_of_node, grid.y_of_node, 'c.')
    #plt.plot([param_dict['fault_trace']['x1'], param_dict['fault_trace']['x2']],
    #         [param_dict['fault_trace']['y1'], param_dict['fault_trace']['y2']], 'r')
    #plt.show()

    out = np.array([ True,  True,  True,  True,  True,  True,  True, False,  True,
                     True,  True,  True, False, False, False, False,  True,  True,
                     False, False, False, False, False, False, False, False, False,
                     False, False, False], dtype=bool)

    assert_array_equal(nf.faulted_nodes, out)
