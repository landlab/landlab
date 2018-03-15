from landlab import RasterModelGrid, FieldError
from landlab.components import FlowAccumulator, FastscapeEroder, FlowDirectorSteepest
import numpy as np
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
    z = mg.add_zeros('topographic__elevation', at='node')
    fd = FlowDirectorSteepest(mg)
    fd.run_one_step()
    # test that the stream length utility will fail because of a ValueError
    assert_raises(FieldError, calculate_stream_length, mg)
