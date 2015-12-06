"""Test the Voronoi compatibility of the grid."""
import os

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab import VoronoiDelaunayGrid
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import StreamPowerEroder
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

    z_tg = np.array([9.51361259e-02,   7.66610606e-02,   5.49662478e-05,
                     9.03050696e-02,   8.26053470e-02,   4.81850412e-02,
                     7.54716080e-02,   4.11015452e-02,   1.42667862e-02,
                     9.25264770e-02,   7.13743892e-02,   8.35794328e-02,
                     6.78789659e-02,   2.21271714e-02,   1.84439866e-05,
                     1.47316188e-02,   8.81752236e-05,   7.08296626e-02,
                     4.90903110e-02,   2.52943889e-02,   5.05246090e-05,
                     6.52865044e-06,   7.35158837e-02,   4.48276368e-02,
                     3.73364102e-02,   4.92939649e-02,   5.69791893e-02,
                     2.12790557e-03,   3.68422348e-02,   1.66012225e-02,
                     5.32744486e-02,   8.02589508e-02,   1.99916970e-02,
                     9.69548363e-02,   5.05236720e-05,   9.17810494e-02,
                     9.25047689e-02,   8.12319949e-02,   1.62298599e-05,
                     7.09279233e-02,   9.64551080e-05,   5.09759307e-02,
                     7.19210532e-02,   2.14645586e-02,   4.88172006e-02,
                     9.33254390e-02,   4.96294330e-02,   6.68276781e-02,
                     7.66748162e-02,   7.33502665e-02,   5.49762478e-05,
                     7.32804768e-02,   3.66342402e-05,   4.14029097e-02,
                     5.93407457e-02,   4.06372632e-02,   8.76088386e-02,
                     7.78024064e-02,   9.37008386e-02,   6.75677053e-02,
                     5.53020421e-02,   4.52790664e-02,   8.92919426e-02,
                     6.49070464e-02,   1.21519788e-02,   8.89827335e-02,
                     6.55068670e-02,   5.93375883e-02,   3.92914015e-03,
                     1.70024742e-02,   5.04000439e-05,   6.09992361e-02,
                     2.59744753e-05,   7.97259074e-02,   8.02709192e-02,
                     7.36747056e-05,   3.72440697e-02,   4.32926006e-02,
                     7.04553131e-02,   5.17799139e-02,   5.23368855e-02,
                     9.11221532e-02,   4.97713791e-02,   9.05450237e-02,
                     7.76611232e-02,   6.22424361e-02,   9.39471770e-02,
                     5.50977905e-05,   9.49390490e-02,   3.19390925e-02,
                     6.61627478e-02,   1.69771597e-02,   7.36113525e-02,
                     7.51088749e-02,   8.28367460e-02,   3.96466978e-02,
                     5.99118692e-02,   1.28542237e-02,   1.70477133e-05,
                     8.81652236e-05])

    assert_array_almost_equal(mg.at_node['topographic__elevation'], z_tg)
