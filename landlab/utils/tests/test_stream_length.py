from landlab import RasterModelGrid, FieldError, HexModelGrid, CLOSED_BOUNDARY
from landlab.components import FlowAccumulator, FastscapeEroder, FlowDirectorSteepest, FlowRouter
import numpy as np
import math
from landlab.utils.stream_length import calculate_stream_length

from nose.tools import assert_equal, assert_true, assert_false, assert_raises

def test_no_flow_recievers():
    """Test that correct error is raised when no flow recievers are on the grid."""
    # instantiate a model grid, do not run flow accumulation on it
    mg = RasterModelGrid(30, 70)
    # test that the stream length utility will fail because of a ValueError
    assert_raises(FieldError, calculate_stream_length, mg)

def test_no_upstream_array():
    """Test that correct error is raised when no flow__upstream_node_order."""
    # instantiate a model grid, do not run flow accumulation on it
    mg = RasterModelGrid(30, 70)

    fd = FlowDirectorSteepest(mg)
    fd.run_one_step()
    # test that the stream length utility will fail because of a ValueError
    assert_raises(FieldError, calculate_stream_length, mg)

def stream_length_regular_grid_d8():
    """Test to demonstrate that stream_length utility works as expected with regular grids"""
    # instantiate a model grid
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    # instantiate an elevation array
    z = np.array ([[0,0,0,0],[0,21,10,0],[0,31,20,0],[0,32,30,0],[0,0,0,0]],dtype='float64')
    # add the elevation field to the grid
    mg.add_field('node','topographic__elevation', z)
    # instantiate the expected flow_length array
    # considering flow directions calculated with D8 algorithm
    flow_length_expected = np.array ([[0,0,0,0],[0,1,0,0],[0,math.sqrt(2),1,0],[0,1+math.sqrt(2),2,0],[0,0,0,0]],dtype='float64')
    #setting boundary conditions
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, left_is_closed=True, right_is_closed=True, top_is_closed=True)
    # calculating flow directions with FlowRouter component
    fr = FlowRouter(mg, method = 'D8')
    fr.route_flow()
    # calculating flow length map
    stream__length = calculate_stream_length(mg, add_to_grid=True, noclobber=False)
    # modifying the stream_length_map because boundary and outlet nodes should not
    # have stream_length value different from 0
    flow_length = np.reshape(stream__length-stream__length[6],mg.shape)
    # test that the stream length utility works as expected
    assert_equal (flow_length_expected.all(),flow_length.all(),mg)
    
def stream_length_regular_grid_d4():
    """Test to demonstrate that stream_length utility works as expected with regular grids"""
    # instantiate a model grid
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    # instantiate an elevation array
    z = np.array ([[0,0,0,0],[0,21,10,0],[0,31,20,0],[0,32,30,0],[0,0,0,0]],dtype='float64')
    # add the elevation field to the grid
    mg.add_field('node','topographic__elevation', z)
    # instantiate the expected flow_length array
    # considering flow directions calculated with D4 algorithm
    flow_length_expected = np.array ([[0,0,0,0],[0,1,0,0],[0,2,1,0],[0,3,2,0],[0,0,0,0]],dtype='float64')
    #setting boundary conditions
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, left_is_closed=True, right_is_closed=True, top_is_closed=True)
    # calculating flow directions with FlowRouter component
    fr = FlowRouter(mg, method = 'D4')
    fr.route_flow()
    # calculating flow length map
    stream__length = calculate_stream_length(mg, add_to_grid=True, noclobber=False)
    # modifying the stream_length_map because boundary and outlet nodes should not
    # have stream_length value different from 0
    flow_length = np.reshape(stream__length-stream__length[6],mg.shape)
    # test that the stream length utility works as expected
    assert_equal (flow_length_expected.all(),flow_length.all(),mg)
    
def stream_length_irregular_grid_d8():
    """Test to demonstrate that stream_length utility works as expected with irregular grids"""
    # instantiate a model grid
    dx=(2./(3.**0.5))**0.5
    hmg = HexModelGrid(5,3, dx)
    # instantiate and add the elevation field
    _ = hmg.add_field('topographic__elevation', hmg.node_x + np.round(hmg.node_y), at = 'node')
    # calculating flow directions with FlowRouter component: D8 algorithm
    fr = FlowRouter(hmg, method = 'D8')
    fr.route_flow()
    # calculating flow length map
    stream__length = calculate_stream_length(hmg, add_to_grid=True, noclobber=False)

def stream_length_irregular_grid_d4():
    """Test to demonstrate that stream_length utility works as expected with irregular grids"""
    # instantiate a model grid
    dx=(2./(3.**0.5))**0.5
    hmg = HexModelGrid(5,3, dx)
    # instantiate and add the elevation field
    _ = hmg.add_field('topographic__elevation', hmg.node_x + np.round(hmg.node_y), at = 'node')
    # calculating flow directions with FlowRouter component: D4 algorithm
    fr = FlowRouter(hmg, method = 'D4')
    fr.route_flow()
    # calculating flow length map
    stream__length = calculate_stream_length(hmg, add_to_grid=True, noclobber=False)