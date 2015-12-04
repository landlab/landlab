"""Test the Voronoi compatibility of the grid."""
import os

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab import VoronoiDelaunayGrid
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder
from landlab.plot.imshow import imshow_node_grid
from pylab import show

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_voronoi():
    nnodes = 100

    np.random.seed(0)
    x = np.random.rand(nnodes)
    np.random.seed(1)
    y = np.random.rand(nnodes)
    mg = VoronoiDelaunayGrid(x, y)

    np.random.seed(2)
    z = mg.add_field('node', 'topographic__elevation',
                     np.random.rand(nnodes) / 10000., copy=False)

    fr = FlowRouter(mg)
    spe = StreamPowerEroder(mg, os.path.join(_THIS_DIR,
                                             'drive_sp_params_voronoi.txt'))

    for i in xrange(10):
        z[mg.core_nodes] += 0.01
        fr.route_flow()
        spe.erode(mg, 1.)

    z_tg = np.array([4.35994902e-05,   2.59262318e-06,   5.49662478e-05,
                     4.36094902e-05,   4.20367802e-05,   2.16664196e-03,
                     1.23484520e-02,   4.07654109e-02,   5.69974645e-02,
                     4.71086393e-02,   5.55056467e-02,   5.53655860e-02,
                     4.87984556e-02,   6.62132007e-02,   7.33132951e-02,
                     6.10003971e-02,   8.53975293e-05,   4.90709812e-02,
                     7.84728490e-02,   8.26091298e-02,   5.05246090e-05,
                     3.69289510e-02,   7.56743264e-02,   2.52850671e-02,
                     6.80517443e-02,   5.96745309e-05,   3.86335644e-02,
                     4.12526935e-02,   7.78476139e-02,   7.18165632e-02,
                     7.62359762e-02,   6.65702961e-02,   9.04684278e-02,
                     7.37203974e-02,   9.16862054e-02,   4.11168411e-02,
                     4.36498947e-02,   9.04430763e-02,   1.42123575e-02,
                     8.02734138e-02,   8.34895711e-02,   4.51564589e-02,
                     6.47124169e-02,   9.52317640e-02,   6.01917122e-05,
                     7.11667309e-02,   8.00463179e-02,   9.40668605e-02,
                     9.33247720e-02,   7.34011920e-02,   7.30924467e-02,
                     6.21387674e-02,   9.66278589e-02,   8.64266637e-02,
                     9.06100497e-02,   7.19991365e-02,   1.07162442e-02,
                     8.68685360e-02,   7.40515066e-02,   2.73264483e-02,
                     9.55406046e-02,   6.01817121e-05,   3.28663723e-02,
                     8.60706344e-02,   9.56285286e-02,   4.24741937e-02,
                     5.64139004e-03,   8.72373392e-02,   9.04304841e-02,
                     8.83708170e-02,   5.04000439e-05,   6.47845229e-02,
                     8.77304766e-02,   4.40894724e-02,   8.08029139e-02,
                     4.52732792e-02,   6.35806712e-02,   6.71004941e-02,
                     7.32682350e-02,   4.88067087e-02,   2.14920377e-02,
                     2.55509418e-06,   8.25346959e-02,   8.12610977e-02,
                     3.75778123e-02,   2.20818912e-02,   2.45940192e-02,
                     4.43365163e-02,   4.93255141e-02,   6.23517572e-02,
                     5.99594648e-02,   4.17977098e-06,   4.96450779e-02,
                     3.72952724e-02,   2.26332108e-03,   1.69925430e-02,
                     1.99835635e-02,   6.61481327e-05,   1.70477133e-05,
                     8.81652236e-05]  )

    assert_array_almost_equal(mg.at_node['topographic__elevation'], z_tg)
