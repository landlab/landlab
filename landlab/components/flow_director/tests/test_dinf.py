
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid, VoronoiDelaunayGrid
from landlab.components import FlowDirectorDINF
from landlab.components.flow_director import flow_direction_dinf


def test_not_implemented_voroni():
    x = [0, 0.1, 0.2, 0.3, 1, 1.1, 1.2, 1.3, 2, 2.1, 2.2, 2.3]
    y = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    vmg = VoronoiDelaunayGrid(x, y)
    with pytest.raises(NotImplementedError):
        flow_direction_dinf.flow_directions_dinf(vmg)


def test_flow__distance_raster_D_infinity_low_closed_boundary_conditions():
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    z = np.array ([[0,  0,  0, 0],
                   [0, 21, 10, 0],
                   [0, 31, 20, 0],
                   [0, 32, 30, 0],
                   [0,  0,  0, 0]],dtype='float64')
    mg.add_field('node','topographic__elevation', z)
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)

    fd = FlowDirectorDINF(mg)
    fd.run_one_step()

    true_recievers = np.array([[ 0, -1],
                               [ 1, -1],
                               [ 2, -1],
                               [ 3, -1],
                               [ 4, -1],
                               [ 6, -1],
                               [ 6, -1],
                               [ 7, -1],
                               [ 8, -1],
                               [10,  6],
                               [ 6, -1],
                               [11, -1],
                               [12, -1],
                               [14, 10],
                               [10, -1],
                               [15, -1],
                               [16, -1],
                               [17, -1],
                               [18, -1],
                               [19, -1]])

    true_proportions = np.array([[  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  1.        ,   0.        ],
                                 [  1.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.06058469,   0.93941531],
                                 [  1.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   1.        ],
                                 [  1.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ],
                                 [  0.        ,   0.        ]])
    assert_array_equal(fd.receivers, true_recievers)
    assert_array_equal(np.round(fd.proportions, decimals=6),
                       np.round(true_proportions, decimals=6))
