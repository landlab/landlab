#Importing the requested modules

import numpy as np
import math
from nose.tools import assert_raises
from numpy.testing import assert_array_equal
from landlab import RasterModelGrid, FieldError, HexModelGrid
from landlab.components import FlowAccumulator, FlowDirectorSteepest
from landlab.utils.stream_length import calculate_stream_length

def test_no_flow_recievers():
    """Test that correct error is raised when no flow recievers are 
    on the grid."""
    
    # instantiate a model grid, do not run flow accumulation on it
    
    mg = RasterModelGrid(30, 70)
    
    # test that the stream length utility will fail because of a ValueError
    
    assert_raises(FieldError, calculate_stream_length, mg)


def test_no_upstream_array():
    """Test that correct error is raised when no flow__upstream_node_order."""
    
    # instantiate a model grid, do not run flow accumulation on it
    
    mg = RasterModelGrid(30, 70)
    
    #Add a field called topographic__elevation to mg
    
    z = mg.add_ones('node','topographic__elevation')
    
    #Run the FlowDirectorSteepest component
    
    fd = FlowDirectorSteepest(mg)   
    fd.run_one_step()
    
    # test that the stream length utility will fail because of a ValueError
    
    assert_raises(FieldError, calculate_stream_length, mg)


def test_stream_length_regular_grid_d8():
    """Test to demonstrate that stream_length utility works as expected with 
    regular grids"""
    
    # instantiate a model grid
    
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    
    # instantiate an elevation array
    
    z = np.array([[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0],
                  [0, 0, 0, 0]], dtype='float64')
    
    # add the elevation field to the grid
    
    mg.add_field('node', 'topographic__elevation', z)
    
    # instantiate the expected flow_length array
    # considering flow directions calculated with D8 algorithm
    
    flow_length_expected = np.array([[0, 0, 0, 0], [0, 1, 0, 0], 
                                     [0, math.sqrt(2), 1, 0], 
                                     [0, 1+math.sqrt(2), 2, 0], [0, 0, 0, 0]], 
                                    dtype='float64')    
    flow_length_expected = np.reshape(flow_length_expected,
                                      mg.number_of_node_rows * 
                                      mg.number_of_node_columns)

    #setting boundary conditions
    
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, 
                                           left_is_closed=True, 
                                           right_is_closed=True, 
                                           top_is_closed=True)
    
    # calculating flow directions with FlowAccumulator component
    
    fr = FlowAccumulator(mg, flow_director='D8')
    fr.run_one_step()
    
    # calculating flow length map
    
    flow_length = calculate_stream_length(mg, add_to_grid=True,
                                          noclobber=False)
    flow_length = np.reshape(flow_length, mg.number_of_node_rows * 
                             mg.number_of_node_columns)
    
    # modifying the flow_length_map because boundary and outlet nodes should 
    # not have stream_length value different from 0
    
    flow_length[mg.boundary_nodes] = 0
    outlet_id = 6
    flow_length[outlet_id] = 0
    
    #Add the flow length to the grid
    
    mg.add_field('node', 'flow_length', flow_length)
    
    # test that the stream length utility works as expected
    
    assert_array_equal(flow_length_expected, flow_length)


def test_stream_length_regular_grid_d4():
    """Test to demonstrate that stream_length utility works as expected with
    regular grids"""
    
    # instantiate a model grid
    
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    
    # instantiate an elevation array
    
    z = np.array([[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0],
                  [0, 0, 0, 0]], dtype='float64')
        
    # add the elevation field to the grid
    
    mg.add_field('node', 'topographic__elevation', z)
    
    # instantiate the expected flow_length array
    # considering flow directions calculated with D4 algorithm
    
    flow_length_expected = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 2, 1, 0], 
                                     [0, 3, 2, 0], [0, 0, 0, 0]], 
                                    dtype='float64')
    flow_length_expected = np.reshape(flow_length_expected, 
                                      mg.number_of_node_rows * 
                                      mg.number_of_node_columns)
    
    #setting boundary conditions
    
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, 
                                           left_is_closed=True, 
                                           right_is_closed=True, 
                                           top_is_closed=True)
    
    # calculating flow directions with FlowAccumulator component
    
    fr = FlowAccumulator(mg, flow_director='D4')
    fr.run_one_step()
    
    # calculating flow length map
    
    flow_length = calculate_stream_length(mg, add_to_grid=True, 
                                          noclobber=False)
    flow_length = np.reshape(flow_length,mg.number_of_node_rows * 
                             mg.number_of_node_columns)
    
    # modifying the stream_length_map because boundary and outlet nodes 
    # should not have stream_length value different from 0
    
    flow_length[mg.boundary_nodes] = 0
    outlet_id = 6
    flow_length[outlet_id] = 0
    
    #Add the flow length to the grid
    
    mg.add_field('node', 'flow_length', flow_length)
    
    # test that the stream length utility works as expected
    
    assert_array_equal(flow_length_expected, flow_length)


def test_stream_length_irregular_grid_d4():
    """Test to demonstrate that stream_length utility works as expected with 
    irregular grids"""
    
    # instantiate a model grid
    
    dx = (2./(3.**0.5))**0.5
    hmg = HexModelGrid(5, 3, dx)
    
    # instantiate and add the elevation field
    
    hmg.add_field('topographic__elevation', hmg.node_x + np.round(hmg.node_y),
                  at='node')
    
    # instantiate the expected flow_length array
    
    flow_length_expected = np.array([0, 0, 0, 0, 0, dx, 0, 0, dx, dx, 2*dx, 0,
                                     0, 2*dx, 2*dx, 0, 0, 0, 0])
        
    #setting boundary conditions   
    
    hmg.set_closed_nodes(hmg.boundary_nodes)
    
    # calculating flow directions with FlowAccumulator component: D4 algorithm
    
    fr = FlowAccumulator(hmg, flow_director = 'D4')
    fr.run_one_step()
    
    # calculating flow length map
    
    flow_length = calculate_stream_length(hmg, add_to_grid=True, 
                                          noclobber=False)
    
   #Add the flow length to the grid
    
    hmg.add_field('node', 'flow_length', flow_length)
    
    # test that the stream length utility works as expected
    
    assert_array_equal(flow_length_expected, flow_length)


def test_stream_length_raster_MFD_diagonals_true():
    """Test of stream length utility with a raster grid and MFD."""
    
    # instantiate a model grid
    
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    
    # instantiate an elevation array
    
    z = np.array([[0 ,0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0],
                  [0, 0, 0, 0]], dtype='float64')
        
    # add the elevation field to the grid
    
    mg.add_field('node', 'topographic__elevation', z)
    
    # instantiate the expected flow_length array
    # considering flow directions calculated with MFD algorithm
    
    flow_length_expected = np.array([[0, 0, 0, 0], [0, 1, 0, 0], 
                                     [0, math.sqrt(2), 1, 0], 
                                     [0, 1+math.sqrt(2), 2, 0], [0, 0, 0, 0]], 
                                    dtype='float64')
    flow_length_expected = np.reshape(flow_length_expected, 
                                      mg.number_of_node_rows * 
                                      mg.number_of_node_columns)
    
    #setting boundary conditions
    
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True, 
                                           right_is_closed=True, 
                                           top_is_closed=True)
    
    # calculating flow directions with FlowAccumulator component
    
    fa = FlowAccumulator(mg, 'topographic__elevation', flow_director='MFD',
                         diagonals=True)
    fa.run_one_step()
    
    # calculating flow length map
    
    flow_length = calculate_stream_length(mg, add_to_grid=True, 
                                          noclobber=False)
    
    # test that the stream length utility works as expected
    
    assert_array_equal(flow_length_expected, flow_length)

def test_stream_length_raster_MFD_diagonals_false():
    """Test of stream length utility with a raster grid and MFD."""
    
    # instantiate a model grid
    
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    
    # instantiate an elevation array
    
    z = np.array([[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0],
                  [0, 0, 0, 0]], dtype='float64')
        
    # add the elevation field to the grid
    
    mg.add_field('node', 'topographic__elevation', z)
    
    # instantiate the expected flow_length array
    # considering flow directions calculated with MFD algorithm
    
    flow_length_expected = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [0, 2, 1, 0],
                                     [0, 3, 2, 0], [0, 0, 0, 0]], 
                                    dtype='float64')
    flow_length_expected = np.reshape(flow_length_expected, 
                                      mg.number_of_node_rows * 
                                      mg.number_of_node_columns)
    
    #setting boundary conditions
    
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, 
                                           left_is_closed=True, 
                                           right_is_closed=True, 
                                           top_is_closed=True)
    
    # calculating flow directions with FlowAccumulator component
    
    fa = FlowAccumulator(mg, 'topographic__elevation', flow_director='MFD',
                         diagonals=False)
    fa.run_one_step()
    
    # calculating flow length map
    
    flow_length = calculate_stream_length(mg, add_to_grid=True, 
                                          noclobber=False)
    
    # test that the stream length utility works as expected
    
    assert_array_equal(flow_length_expected, flow_length)
    
def test_stream_length_raster_D_infinity():
    """Test of stream length utility with a raster grid and D infinity."""
    
    # instantiate a model grid
    
    mg = RasterModelGrid((5, 4), spacing=(1, 1))
    
    # instantiate an elevation array
    
    z = np.array([[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0],
                  [0, 0, 0, 0]], dtype='float64')
        
    # add the elevation field to the grid
    
    mg.add_field('node', 'topographic__elevation', z)
    
    #setting boundary conditions
    
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, 
                                           left_is_closed=True, 
                                           right_is_closed=True, 
                                           top_is_closed=True)
    
    # calculating flow directions with FlowAccumulator component
    
    fa = FlowAccumulator(mg, 'topographic__elevation', flow_director='DINF')
    fa.run_one_step()
    
    # calculating flow length map
    
    flow_length = calculate_stream_length(mg, add_to_grid=True, 
                                          noclobber=False)
