from landlab import RasterModelGrid, HexModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator
import numpy as np
from numpy import testing
from nose.tools import assert_raises

def test_Ff_bad_vals():
    """
    Test that instantiating ErosionDeposition with a F_f value > 1 throws a 
    ValueError.
    """

    #set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), 10.0)

    mg.add_zeros('node', 'topographic__elevation')

    mg['node']['topographic__elevation'] += mg.node_y / 100000 \
        + mg.node_x / 100000 \
        + np.random.rand(len(mg.node_y)) / 10000
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0,
                                                  mg['node']['topographic__elevation'],
                                                  -9999.)

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='D8',
                         depression_finder='DepressionFinderAndRouter')


    # Instantiate the ErosionDeposition component...
    assert_raises(ValueError, ErosionDeposition, mg, K=0.01, F_f=2.0, phi=0.5, 
                  v_s=0.001, m_sp=0.5, n_sp=1.0, sp_crit=0.0, solver='basic')


def test_phi_bad_vals():
    """
    Test that instantiating ErosionDeposition with a phi value >= 1 throws a 
    ValueError.
    """

    #set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), 10.0)

    mg.add_zeros('node', 'topographic__elevation')

    mg['node']['topographic__elevation'] += mg.node_y / 100000 \
        + mg.node_x / 100000 \
        + np.random.rand(len(mg.node_y)) / 10000
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0,
                                                  mg['node']['topographic__elevation'],
                                                  -9999.)

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='D8',
                         depression_finder='DepressionFinderAndRouter')

    # Instantiate the ErosionDeposition component...
    assert_raises(ValueError, ErosionDeposition, mg, K=0.01, F_f=0.0, phi=2.0, 
                  v_s=0.001, m_sp=0.5, n_sp=1.0, sp_crit=0.0, solver='basic')


def test_q_as_field():
    """
    Test that passing in water discharge as a grid field results in self.q
    holding correct values
    """

    #set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), 10.0)

    mg.add_zeros('node', 'topographic__elevation')
    q = mg.add_zeros('node', 'user_imposed_discharge')
    q[:] += 1.0 #add 1.0 m3/yr of water

    mg['node']['topographic__elevation'] += mg.node_y / 100000 \
        + mg.node_x / 100000 \
        + np.random.rand(len(mg.node_y)) / 10000
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0,
                                                  mg['node']['topographic__elevation'],
                                                  -9999.)

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='D8',
                         depression_finder='DepressionFinderAndRouter')

    # Instantiate the ErosionDeposition component...
    ed = ErosionDeposition(mg, K=0.01, F_f=0.0, phi=0.0, v_s=0.001, m_sp=0.5, 
                           n_sp=1.0, sp_crit=0.0, 
                           discharge_field='user_imposed_discharge', 
                           solver='basic')
    
    #ensure that ed.q is everywhere equal to 1.0 m3/yr.
    testing.assert_array_equal(np.ones(mg.number_of_nodes), 
                               ed.q,
                               err_msg='E/D discharge field test failed',
                               verbose=True)

def test_q_as_array():
    """
    Test that passing in water discharge as an array results in self.q
    holding correct values
    """

    #set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), 10.0)

    mg.add_zeros('node', 'topographic__elevation')
    q = np.zeros(mg.number_of_nodes)
    q[:] += 1.0 #add 1.0 m3/yr of water

    mg['node']['topographic__elevation'] += mg.node_y / 100000 \
        + mg.node_x / 100000 \
        + np.random.rand(len(mg.node_y)) / 10000
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0,
                                                  mg['node']['topographic__elevation'],
                                                  -9999.)

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='D8',
                         depression_finder='DepressionFinderAndRouter')

    # Instantiate the ErosionDeposition component...
    ed = ErosionDeposition(mg, K=0.01, F_f=0.0, phi=0.0, v_s=0.001, m_sp=0.5, 
                           n_sp=1.0, sp_crit=0.0, 
                           discharge_field=q, 
                           solver='basic')
    
    #ensure that ed.q is everywhere equal to 1.0 m3/yr.
    testing.assert_array_equal(np.ones(mg.number_of_nodes), 
                               ed.q,
                               err_msg='E/D discharge array test failed',
                               verbose=True)

def test_sediment__flux_already_created():
    """
    Test that an existing sediment flux grid field is not changed by 
    instantiating ErosionDeposition.
    """

    #set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), 10.0)

    mg.add_zeros('node', 'topographic__elevation')
    qs = mg.add_zeros('node', 'sediment__flux')
    qs[:] += 1.0 #add 1.0 m3/yr of flux

    mg['node']['topographic__elevation'] += mg.node_y / 100000 \
        + mg.node_x / 100000 \
        + np.random.rand(len(mg.node_y)) / 10000
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0,
                                                  mg['node']['topographic__elevation'],
                                                  -9999.)

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='D8',
                         depression_finder='DepressionFinderAndRouter')

    # Instantiate the ErosionDeposition component...
    ed = ErosionDeposition(mg, K=0.01, F_f=0.0, phi=0.0, v_s=0.001, m_sp=0.5, 
                           n_sp=1.0, sp_crit=0.0, solver='basic')
    
    #ensure that 'sediment__flux' field is everywhere equal to 1.0 m3/yr.
    testing.assert_array_equal(np.ones(mg.number_of_nodes), 
                               ed.qs,
                               err_msg='E/D sediment flux field test failed',
                               verbose=True)
