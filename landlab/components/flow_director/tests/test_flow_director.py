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
    """Check to make sure the right fields have been created.

    Check that the sizes of at least one are also correct (they are at node
    so if one is the right size, then they all should be)
    """

    mg0 = RasterModelGrid((10,10), spacing=(1, 1))
    z0 = mg0.add_field('topographic__elevation', mg0.node_x**2 + mg0.node_y**2, at = 'node')
    fd0 = _FlowDirector(mg0, 'topographic__elevation')
    assert_equal(list(mg0.at_node.keys()), ['topographic__elevation'])
    assert_equal(np.size(mg0.at_node['topographic__elevation']), mg0.number_of_nodes)

    mg1 = RasterModelGrid((10,10), spacing=(1, 1))
    z1 = mg1.add_field('topographic__elevation', mg1.node_x**2 + mg1.node_y**2, at = 'node')
    fd1 = _FlowDirectorToMany(mg1, 'topographic__elevation')
    assert_equal(sorted(list(mg1.at_node.keys())), ['flow__link_to_receiver_node',
                                                     'flow__receiver_node',
                                                     'flow__sink_flag',
                                                     'topographic__elevation',
                                                     'topographic__steepest_slope'])
    assert_equal(np.size(mg1.at_node['topographic__elevation']), mg1.number_of_nodes)

    mg2 = RasterModelGrid((10,10), spacing=(1, 1))
    z2 = mg2.add_field('topographic__elevation', mg2.node_x**2 + mg2.node_y**2, at = 'node')
    fd2 = _FlowDirectorToOne(mg2, 'topographic__elevation')
    assert_equal(sorted(list(mg2.at_node.keys())), ['flow__link_to_receiver_node',
                                                    'flow__receiver_node',
                                                    'flow__sink_flag',
                                                    'topographic__elevation',
                                                    'topographic__steepest_slope'])
    assert_equal(np.size(mg2.at_node['topographic__elevation']), mg2.number_of_nodes)


    mg3 = RasterModelGrid((10,10), spacing=(1, 1))
    z3 = mg3.add_field('topographic__elevation', mg3.node_x**2 + mg3.node_y**2, at = 'node')
    fd3 = FlowDirectorMFD(mg3, 'topographic__elevation')
    assert_equal(sorted(list(mg3.at_node.keys())), ['flow__link_to_receiver_node',
                                                     'flow__links_to_receiver_nodes',
                                                     'flow__receiver_node',
                                                     'flow__receiver_nodes',
                                                     'flow__receiver_proportions',
                                                     'flow__sink_flag',
                                                     'topographic__elevation',
                                                     'topographic__steepest_slope'])
    assert_equal(np.size(mg3.at_node['topographic__elevation']), mg3.number_of_nodes)



    mg4 = RasterModelGrid((10,10), spacing=(1, 1))
    z4 = mg4.add_field('topographic__elevation', mg4.node_x**2 + mg4.node_y**2, at = 'node')
    fd4 = FlowDirectorDINF(mg4, 'topographic__elevation')
    assert_equal(sorted(list(mg4.at_node.keys())), ['flow__link_to_receiver_node',
                                                     'flow__links_to_receiver_nodes',
                                                     'flow__receiver_node',
                                                     'flow__receiver_nodes',
                                                     'flow__receiver_proportions',
                                                     'flow__sink_flag',
                                                     'topographic__elevation',
                                                     'topographic__steepest_slope'])
    assert_equal(np.size(mg4.at_node['topographic__elevation']), mg4.number_of_nodes)

    mg5 = RasterModelGrid((10,10), spacing=(1, 1))
    z5 = mg5.add_field('topographic__elevation', mg5.node_x**2 + mg5.node_y**2, at = 'node')
    fd5 = FlowDirectorSteepest(mg5, 'topographic__elevation')
    assert_equal(sorted(list(mg5.at_node.keys())), ['flow__link_to_receiver_node',
                                                    'flow__receiver_node',
                                                    'flow__sink_flag',
                                                    'topographic__elevation',
                                                    'topographic__steepest_slope'])
    assert_equal(np.size(mg5.at_node['topographic__elevation']), mg5.number_of_nodes)


    mg6 = RasterModelGrid((10,10), spacing=(1, 1))
    z6 = mg6.add_field('topographic__elevation', mg6.node_x**2 + mg6.node_y**2, at = 'node')
    fd6 = FlowDirectorD8(mg6, 'topographic__elevation')
    assert_equal(sorted(list(mg6.at_node.keys())), ['flow__link_to_receiver_node',
                                                    'flow__receiver_node',
                                                    'flow__sink_flag',
                                                    'topographic__elevation',
                                                    'topographic__steepest_slope'])
    assert_equal(np.size(mg6.at_node['topographic__elevation']), mg6.number_of_nodes)
