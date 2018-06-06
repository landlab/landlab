from landlab import RasterModelGrid, HexModelGrid
from landlab.components import Space, FlowAccumulator
import numpy as np
from numpy import testing

def test_matches_detachment_solution():
    """
    Test that model matches the detachment-limited analytical solution
    for slope/area relationship at steady state: S=(U/K)^(1/n)*A^(-m/n).
    """
    
    #set up a 3x3 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), 10.0)
    
    z = mg.add_zeros('node', 'topographic__elevation')
    br = mg.add_zeros('node', 'bedrock__elevation')
    soil = mg.add_zeros('node', 'soil__depth')

    mg['node']['topographic__elevation'] += mg.node_y / 10000 \
        + mg.node_x / 10000 \
        + np.random.rand(len(mg.node_y)) / 10000
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0,
                                                  mg['node']['topographic__elevation'], 
                                                  -9999.)
    br[:] = z[:] - soil[:]
        
    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='D8')
    
    #Instantiate DepressionFinderAndRouter
    
    # Parameter values for detachment-limited test
    K_br = 0.01
    U = 0.0001
    dt = 1.0
    F_f = 1.0 #all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0

    # Instantiate the Space component...
    sp = Space(mg, K_sed=0.00001, K_br=K_br,
                         F_f=F_f, phi=0.1, H_star=1., v_s=0.001,
                         m_sp=m_sp, n_sp=n_sp, sp_crit_sed=0,
                         sp_crit_br=0)

    # ... and run it to steady state (2000x1-year timesteps).
    for i in range(2000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt #m/yr
        br[mg.core_nodes] = z[mg.core_nodes] - soil[mg.core_nodes]
    
    #compare numerical and analytical slope solutions
    num_slope = mg.at_node['topographic__steepest_slope'][mg.core_nodes]
    analytical_slope = np.power(U / K_br, 1./n_sp) \
        * np.power(mg.at_node['drainage_area'][mg.core_nodes], -m_sp / n_sp)
    testing.assert_array_almost_equal(num_slope, analytical_slope, 
                                      decimal=8, 
                                      err_msg='SPACE detachment-limited test failed',
                                      verbose=True)
        
def test_can_run_with_hex():
    """Test that model can run with hex model grid."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    mg = HexModelGrid(7, 7)
    z = mg.add_zeros('node', 'topographic__elevation')
    z[:] = 0.01 * mg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='FlowDirectorSteepest')

    # Parameter values for test 1
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    sp = Space(mg, K_sed=0.00001, K_br=0.00000000001,
                         F_f=0.5, phi=0.1, H_star=1., v_s=0.001,
                         m_sp=0.5, n_sp = 1.0, sp_crit_sed=0,
                         sp_crit_br=0)

    # ... and run it to steady state.
    for i in range(2000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt
