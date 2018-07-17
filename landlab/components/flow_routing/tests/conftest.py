import os

import pytest
import numpy as np

from landlab import RasterModelGrid
from landlab import BAD_INDEX_VALUE as XX


@pytest.fixture
def dans_grid1():
    """
    Create a 5x5 test grid.
    This is a sheet flow test.
    """
    # global fr, mg, infile
    # global z, Q_in, A_target, frcvr_target, upids_target, Q_target, \
    #     steepest_target, links2rcvr_target

    mg = RasterModelGrid((5, 5), spacing=(10., 10.))

    this_dir = os.path.abspath(os.path.dirname(__file__))
    infile = os.path.join(this_dir, 'test_fr_input.txt')

    z = mg.node_x.copy()

    Q_in = np.full(25, 2.)

    A_target = np.array([0.,  0.,  0.,  0.,  0.,
                         3.,  3.,  2.,  1.,  0.,
                         3.,  3.,  2.,  1.,  0.,
                         3.,  3.,  2.,  1.,  0.,
                         0.,  0.,  0.,  0.,  0.])*100.

    frcvr_target = np.array([0,  1,  2,  3,  4,
                             5,  5,  6,  7,  9,
                            10, 10, 11, 12, 14,
                            15, 15, 16, 17, 19,
                            20, 21, 22, 23, 24])

    upids_target = np.array([0,  1,  2,  3,  4,
                             5,  6,  7,  8,  9,
                            10, 11, 12, 13, 14,
                            15, 16, 17, 18, 19,
                            20, 21, 22, 23, 24])

    links2rcvr_target = np.full(25, XX)
    links2rcvr_target[mg.core_nodes] = np.array([ 9, 10, 11,
                                                 18, 19, 20,
                                                 27, 28, 29])

    Q_target = A_target * 2.  # only once Q_in is used

    steepest_target = np.array([0.,  0.,  0.,  0.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  0.,  0.,  0.,  0.])

    mg.add_field('node', 'topographic__elevation', z, units='-')

    class DansGrid(object):
        pass

    dans_grid1 = DansGrid()
    dans_grid1.mg = mg
    dans_grid1.z = z
    dans_grid1.infile = infile
    dans_grid1.A_target = A_target
    dans_grid1.frcvr_target = frcvr_target 
    dans_grid1.upids_target = upids_target 
    dans_grid1.Q_target = Q_target 
    dans_grid1.steepest_target = steepest_target 
    dans_grid1.links2rcvr_target = links2rcvr_target 

    return dans_grid1
