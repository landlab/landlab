"""Test the flow director components.

@author: krb
"""
# Created on Thurs Nov 12, 2015
import os

from six.moves import range

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from nose.tools import (with_setup, assert_true, assert_false, assert_raises,
                        assert_almost_equal, assert_equal)
try:
    from nose.tools import (assert_is, assert_set_equal, assert_dict_equal)
except ImportError:
    from landlab.testing.tools import (assert_is, assert_set_equal,
                                       assert_dict_equal)

import landlab
from landlab import RasterModelGrid, HexModelGrid, RadialModelGrid, FieldError
from landlab.components.flow_director import(FlowDirectorD8, 
                                             FlowDirectorDINF, 
                                             FlowDirectorMFD, 
                                             FlowDirectorSteepest)

from landlab.components.flow_director.flow_director import _FlowDirector
from landlab.components.flow_director.flow_director_to_many import _FlowDirectorToMany
from landlab.components.flow_director.flow_director_to_one import _FlowDirectorToOne

from landlab import CLOSED_BOUNDARY
from landlab import BAD_INDEX_VALUE as XX


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))

def test_not_implemented():
    """Test that private run_one_step is not implemented"""
    
    mg = RasterModelGrid((10,10), spacing=(1, 1))
    z = mg.add_field('topographic__elevation', mg.node_x**2 + mg.node_y**2, at = 'node')
    
    fd0 = _FlowDirector(mg, 'topographic__elevation')
    assert_raises(NotImplementedError, fd0.run_one_step)
    
    fd1 = _FlowDirectorToMany(mg, 'topographic__elevation')
    assert_raises(NotImplementedError, fd1.run_one_step)

    fd2 = _FlowDirectorToOne(mg, 'topographic__elevation')
    assert_raises(NotImplementedError, fd2.run_one_step)
    

def test_grid_type_testing():
    """Test that only the right grids can be implemented."""
    dx=(2./(3.**0.5))**0.5
    hmg = HexModelGrid(9,5, dx)
    z = hmg.add_field('topographic__elevation', hmg.node_x + np.round(hmg.node_y), at = 'node')
        
     # D8 is ONLY RASTER
    assert_raises(NotImplementedError, FlowDirectorD8, hmg)
    
    # DINF IS ONLY RASTER RASTER
    assert_raises(NotImplementedError, FlowDirectorDINF, hmg)
    
def test_check_fields():
    """Check to make sure the right fields have been created."""
    
    mg = RasterModelGrid((10,10), spacing=(1, 1))
    z = mg.add_field('topographic__elevation', mg.node_x**2 + mg.node_y**2, at = 'node')
    
    fd0 = _FlowDirector(mg, 'topographic__elevation')
    mg.at_node.keys()
    

    fd1 = _FlowDirectorToMany(mg, 'topographic__elevation')

    fd2 = _FlowDirectorToOne(mg, 'topographic__elevation')
    
    fd3 = FlowDirectorMFD(mg, 'topographic__elevation')

    fd4 = FlowDirectorDINF(mg, 'topographic__elevation')
    
    fd5 = FlowDirectorSteepest(mg, 'topographic__elevation')

    fd6 = FlowDirectorD8(mg, 'topographic__elevation')
    