from landlab import RasterModelGrid, HexModelGrid
from landlab.components import Space, FlowAccumulator
import numpy as np
from numpy.testing import assert_equal

def test_can_run_with_hex():
    """Test that model can run with hex model grid."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    mg = HexModelGrid(7, 7)
    z = mg.add_zeros('node', 'topographic__elevation')
    z[:] = 0.01 * mg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='FlowDirectorSteepest')

    # Parameter values for test 1
    K = 0.001
    vs = 0.0001
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    sp = Space(mg, K_sed=0.00001, K_br=0.00000000001,
                         F_f=0.5, phi=0.1, H_star=1., v_s=0.001,
                         m_sp=0.5, n_sp = 1.0, sp_crit_sed=0,
                         sp_crit_br=0, method='simple_stream_power',
                         discharge_method=None, area_field=None,
                         discharge_field=None)

    # ... and run it to steady state.
    for i in range(2000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt
