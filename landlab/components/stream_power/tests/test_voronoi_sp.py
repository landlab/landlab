"""Test the Voronoi compatibility of the grid."""
import os

import numpy as np
from numpy.testing import assert_array_almost_equal
from six.moves import range

from landlab import VoronoiDelaunayGrid
from landlab.components import FlowAccumulator, StreamPowerEroder

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_voronoi():
    nnodes = 100

    np.random.seed(0)
    x = np.random.rand(nnodes)
    np.random.seed(1)
    y = np.random.rand(nnodes)
    mg = VoronoiDelaunayGrid(x, y)

    np.random.seed(2)
    z = mg.add_field(
        "node", "topographic__elevation", np.random.rand(nnodes) / 10000.0, copy=False
    )

    fr = FlowAccumulator(mg)
    spe = StreamPowerEroder(mg, os.path.join(_THIS_DIR, "drive_sp_params_voronoi.txt"))

    for i in range(10):
        z[mg.core_nodes] += 0.01
        fr.run_one_step()
        spe.erode(mg, 1.0)

    z_tg = np.array(
        [
            4.35994902e-05,
            2.59262318e-06,
            5.49662478e-05,
            6.56738615e-03,
            4.20367802e-05,
            1.21371424e-02,
            2.16596169e-02,
            4.73320898e-02,
            6.00389761e-02,
            5.22007356e-02,
            5.37507115e-02,
            5.95794752e-02,
            5.29862904e-02,
            6.76465914e-02,
            7.31720024e-02,
            6.18730861e-02,
            8.53975293e-05,
            5.32189275e-02,
            7.34302556e-02,
            8.07385044e-02,
            5.05246090e-05,
            4.08940657e-02,
            7.39971005e-02,
            3.31915602e-02,
            6.72650419e-02,
            5.96745309e-05,
            4.72752445e-02,
            3.60359567e-02,
            7.59432065e-02,
            7.24461985e-02,
            7.80305760e-02,
            4.93866869e-02,
            8.69642467e-02,
            7.21627626e-02,
            8.96368291e-02,
            4.65142080e-02,
            6.07720217e-02,
            8.83372939e-02,
            2.35887558e-02,
            7.97616193e-02,
            8.35615355e-02,
            4.61809032e-02,
            6.34634214e-02,
            9.25711770e-02,
            4.11717225e-03,
            7.24493623e-02,
            7.97908053e-02,
            9.10375623e-02,
            9.13155023e-02,
            7.10567915e-02,
            7.35271752e-02,
            6.13091341e-02,
            9.45498463e-02,
            8.48532386e-02,
            8.82702021e-02,
            7.14969941e-02,
            2.22640943e-02,
            8.53311932e-02,
            7.49161159e-02,
            3.48837223e-02,
            9.30132692e-02,
            6.01817121e-05,
            3.87455443e-02,
            8.44673586e-02,
            9.35213577e-02,
            6.76075824e-02,
            1.58614508e-02,
            8.51346837e-02,
            8.83645680e-02,
            8.69944117e-02,
            5.04000439e-05,
            5.02319084e-02,
            8.63882765e-02,
            5.00991880e-02,
            7.65156630e-02,
            5.07591983e-02,
            6.54909962e-02,
            6.91505342e-02,
            7.33358371e-02,
            5.30109890e-02,
            2.99074601e-02,
            2.55509418e-06,
            8.21523907e-02,
            8.09368483e-02,
            4.35073025e-02,
            3.04096109e-02,
            3.26298627e-02,
            4.92259177e-02,
            5.48690358e-02,
            6.44732130e-02,
            6.28133567e-02,
            4.17977098e-06,
            5.37149677e-02,
            4.32828136e-02,
            1.30559903e-02,
            2.62405261e-02,
            2.86079272e-02,
            6.61481327e-05,
            1.70477133e-05,
            8.81652236e-05,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_tg)
