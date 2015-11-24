import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab import RasterModelGrid, ModelParameterDictionary
from landlab.components.diffusion.diffusion import LinearDiffuser

"""
This tester turns over the diffuser a couple of times to ensure basic
functionality is working.
"""


def test_diffusion():
    inputs = ModelParameterDictionary('../landlab/components/diffusion/' +
                                      'tests/diffusion_params.txt')
    nrows = inputs.read_int('nrows')
    ncols = inputs.read_int('ncols')
    dx = inputs.read_float('dx')
    dt = inputs.read_float('dt')
    time_to_run = inputs.read_float('run_time')
    init_elev = inputs.read_float('init_elev')

    mg = RasterModelGrid((nrows, ncols), (dx, dx))
    uplift_rate = mg.node_y[mg.core_cells]/100000.

    # create the fields in the grid
    mg.create_node_array_zeros('topographic__elevation')
    z = mg.create_node_array_zeros() + init_elev
    np.random.seed(0)
    mg['node']['topographic__elevation'] = z + np.random.rand(len(z))/1000.

    mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

    # instantiate:
    dfn = LinearDiffuser(mg, inputs)

    # perform the loop:
    elapsed_time = 0.  # total time in simulation
    while elapsed_time < time_to_run:
        if elapsed_time+dt > time_to_run:
            dt = time_to_run - elapsed_time
        dfn.diffuse(dt)
        mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_rate*dt
        elapsed_time += dt

    z_target = np.array([5.48813504e-04,   7.15189366e-04,   6.02763376e-04,
                         5.44883183e-04,   4.23654799e-04,   6.45894113e-04,
                         4.37587211e-04,   8.91773001e-04,   9.63662761e-04,
                         3.83441519e-04,   7.91725038e-04,   9.17712569e-04,
                         1.02071232e-03,   1.10556005e-03,   1.14946096e-03,
                         1.20022436e-03,   1.12938983e-03,   1.12734463e-03,
                         1.00699946e-03,   8.70012148e-04,   9.78618342e-04,
                         1.12639399e-03,   1.41449092e-03,   2.66410106e-03,
                         2.80054364e-03,   2.82446032e-03,   2.68730913e-03,
                         2.44468915e-03,   2.03932919e-03,   4.14661940e-04,
                         2.64555612e-04,   2.14829293e-03,   2.77874345e-03,
                         3.21568734e-03,   3.45767578e-03,   4.46407381e-03,
                         4.25087300e-03,   3.82133955e-03,   3.24905386e-03,
                         6.81820299e-04,   3.59507901e-04,   3.36368061e-03,
                         4.19756389e-03,   4.81101469e-03,   5.12991298e-03,
                         5.14551834e-03,   4.82212362e-03,   5.21733341e-03,
                         4.36590155e-03,   3.63710771e-04,   5.70196770e-04,
                         4.64503306e-03,   5.67212626e-03,   6.43364397e-03,
                         6.85196170e-03,   6.84719792e-03,   6.44199511e-03,
                         5.63718526e-03,   4.53705235e-03,   2.44425592e-04,
                         1.58969584e-04,   5.85390586e-03,   7.15167925e-03,
                         8.09808225e-03,   8.58987393e-03,   8.60040722e-03,
                         8.08353544e-03,   7.11357773e-03,   5.74363061e-03,
                         9.60984079e-05,   9.76459465e-04,   6.28409752e-03,
                         7.69302638e-03,   9.77162122e-03,   1.03665461e-02,
                         1.03596559e-02,   9.77581391e-03,   8.63184402e-03,
                         7.06221620e-03,   1.18727719e-04,   3.17983179e-04,
                         7.42141185e-03,   9.16089943e-03,   1.04492602e-02,
                         1.11235269e-02,   1.21428645e-02,   1.14619409e-02,
                         1.01991149e-02,   8.52209581e-03,   9.29296198e-04,
                         3.18568952e-04,   8.66507665e-03,   1.06513243e-02,
                         1.20945739e-02,   1.28805314e-02,   1.28817803e-02,
                         1.21355585e-02,   1.16763419e-02,   9.65364442e-03,
                         4.69547619e-06,   6.77816537e-04,   9.99969455e-03,
                         1.21210512e-02,   1.37222525e-02,   1.45639747e-02,
                         1.45899214e-02,   1.37503063e-02,   1.21879324e-02,
                         1.00970967e-02,   9.52749012e-04,   4.47125379e-04,
                         1.11861039e-02,   1.35273547e-02,   1.52388381e-02,
                         1.61710069e-02,   1.61700373e-02,   1.52797248e-02,
                         1.35608698e-02,   1.12630871e-02,   6.92531590e-04,
                         7.25254280e-04,   1.14112447e-02,   1.38209264e-02,
                         1.66345746e-02,   1.75854796e-02,   1.76037603e-02,
                         1.66317633e-02,   1.48117588e-02,   1.22941311e-02,
                         2.90077607e-04,   6.18015429e-04,   1.24651035e-02,
                         1.49479969e-02,   1.67762718e-02,   1.77682020e-02,
                         1.87567995e-02,   1.77826641e-02,   1.59106401e-02,
                         1.33841504e-02,   4.31418435e-04,   8.96546596e-04,
                         1.34298871e-02,   1.58121125e-02,   1.76176213e-02,
                         1.85627061e-02,   1.85694232e-02,   1.75911811e-02,
                         1.67878955e-02,   1.43993331e-02,   9.98847007e-04,
                         1.49448305e-04,   1.40269701e-02,   1.63483655e-02,
                         1.80226762e-02,   1.89156586e-02,   1.88916976e-02,
                         1.79880403e-02,   1.62672916e-02,   1.39254090e-02,
                         6.91669955e-05,   6.97428773e-04,   1.46967049e-02,
                         1.65718319e-02,   1.79957330e-02,   1.87279193e-02,
                         1.87059832e-02,   1.79052307e-02,   1.64399258e-02,
                         1.44435378e-02,   1.71629677e-04,   5.21036606e-04,
                         1.40296771e-02,   1.54293506e-02,   1.75066749e-02,
                         1.80476077e-02,   1.79866098e-02,   1.73732550e-02,
                         1.62714602e-02,   1.47877073e-02,   3.18389295e-05,
                         1.64694156e-04,   1.41367998e-02,   1.49517470e-02,
                         1.57129817e-02,   1.59700260e-02,   1.68585204e-02,
                         1.64421520e-02,   1.58441873e-02,   1.50473253e-02,
                         3.11944995e-04,   3.98221062e-04,   2.09843749e-04,
                         1.86193006e-04,   9.44372390e-04,   7.39550795e-04,
                         4.90458809e-04,   2.27414628e-04,   2.54356482e-04,
                         5.80291603e-05,   4.34416626e-04])

    assert_array_almost_equal(mg.at_node['topographic__elevation'], z_target)
