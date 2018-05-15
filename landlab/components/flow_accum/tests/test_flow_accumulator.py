"""Test the flow accumulator component.

@author: krb
"""
# Created on Thurs Nov 12, 2015
import os

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
from landlab.components.flow_accum import FlowAccumulator

from landlab.components.flow_routing.lake_mapper import DepressionFinderAndRouter

from landlab.components.flow_director import(FlowDirectorD8,
                                             FlowDirectorDINF,
                                             FlowDirectorMFD,
                                             FlowDirectorSteepest)

from landlab import CLOSED_BOUNDARY
from landlab import BAD_INDEX_VALUE as XX

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_check_fields():
    """Check to make sure the right fields have been created."""

    mg = RasterModelGrid((10,10), spacing=(1, 1))
    z = mg.add_field('topographic__elevation', mg.node_x**2 + mg.node_y**2, at = 'node')

    fa = FlowAccumulator(mg)
    assert_array_equal(z, mg.at_node['topographic__elevation'])
    assert_array_equal(np.zeros(100), mg.at_node['drainage_area'])
    assert_array_equal(np.ones(100), mg.at_node['water__unit_flux_in'])

    fa = FlowAccumulator(mg, runoff_rate=2.)
    assert_array_equal(np.full(100, 2.), mg.at_node['water__unit_flux_in'])


def test_director_adding_methods_are_equivalent_Steepest():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10,10), spacing=(1, 1))
    z0 = mg0.add_field('topographic__elevation', mg0.node_x**2 + mg0.node_y**2, at = 'node')
    fa0 = FlowAccumulator(mg0, flow_director='D4')
    fa0.run_one_step()

    mg1 = RasterModelGrid((10,10), spacing=(1, 1))
    z1 = mg1.add_field('topographic__elevation', mg1.node_x**2 + mg1.node_y**2, at = 'node')
    fa1 = FlowAccumulator(mg1, flow_director='Steepest')
    fa1.run_one_step()

    mg2 = RasterModelGrid((10,10), spacing=(1, 1))
    z2 = mg2.add_field('topographic__elevation', mg2.node_x**2 + mg2.node_y**2, at = 'node')
    fa2 = FlowAccumulator(mg2, flow_director=FlowDirectorSteepest)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10,10), spacing=(1, 1))
    z3 = mg3.add_field('topographic__elevation', mg3.node_x**2 + mg3.node_y**2, at = 'node')
    fd = FlowDirectorSteepest(mg3)
    fa3 = FlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key],
                           mg1.at_node[key])

        assert_array_equal(mg1.at_node[key],
                           mg2.at_node[key])

        assert_array_equal(mg2.at_node[key],
                           mg3.at_node[key])


def test_director_adding_methods_are_equivalent_D8():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10,10), spacing=(1, 1))
    z0 = mg0.add_field('topographic__elevation', mg0.node_x**2 + mg0.node_y**2, at = 'node')
    fa0 = FlowAccumulator(mg0, flow_director='D8')
    fa0.run_one_step()

    mg1 = RasterModelGrid((10,10), spacing=(1, 1))
    z1 = mg1.add_field('topographic__elevation', mg1.node_x**2 + mg1.node_y**2, at = 'node')
    fa1 = FlowAccumulator(mg1, flow_director='FlowDirectorD8')
    fa1.run_one_step()

    mg2 = RasterModelGrid((10,10), spacing=(1, 1))
    z2 = mg2.add_field('topographic__elevation', mg2.node_x**2 + mg2.node_y**2, at = 'node')
    fa2 = FlowAccumulator(mg2, flow_director=FlowDirectorD8)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10,10), spacing=(1, 1))
    z3 = mg3.add_field('topographic__elevation', mg3.node_x**2 + mg3.node_y**2, at = 'node')
    fd = FlowDirectorD8(mg3)
    fa3 = FlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key],
                           mg1.at_node[key])

        assert_array_equal(mg1.at_node[key],
                           mg2.at_node[key])

        assert_array_equal(mg2.at_node[key],
                           mg3.at_node[key])


def test_director_adding_methods_are_equivalent_Dinf():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10,10), spacing=(1, 1))
    z0 = mg0.add_field('topographic__elevation', mg0.node_x**2 + mg0.node_y**2, at = 'node')
    fa0 = FlowAccumulator(mg0, flow_director='DINF')
    fa0.run_one_step()

    mg1 = RasterModelGrid((10,10), spacing=(1, 1))
    z1 = mg1.add_field('topographic__elevation', mg1.node_x**2 + mg1.node_y**2, at = 'node')
    fa1 = FlowAccumulator(mg1, flow_director='FlowDirectorDINF')
    fa1.run_one_step()

    mg2 = RasterModelGrid((10,10), spacing=(1, 1))
    z2 = mg2.add_field('topographic__elevation', mg2.node_x**2 + mg2.node_y**2, at = 'node')
    fa2 = FlowAccumulator(mg2, flow_director=FlowDirectorDINF)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10,10), spacing=(1, 1))
    z3 = mg3.add_field('topographic__elevation', mg3.node_x**2 + mg3.node_y**2, at = 'node')
    fd = FlowDirectorDINF(mg3)
    fa3 = FlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key],
                           mg1.at_node[key])

        assert_array_equal(mg1.at_node[key],
                           mg2.at_node[key])

        assert_array_equal(mg2.at_node[key],
                           mg3.at_node[key])


def test_director_adding_methods_are_equivalent_MFD():
    """Check that different methods to specifying the director are the same."""

    mg0 = RasterModelGrid((10,10), spacing=(1, 1))
    z0 = mg0.add_field('topographic__elevation', mg0.node_x**2 + mg0.node_y**2, at = 'node')
    fa0 = FlowAccumulator(mg0, flow_director='MFD')
    fa0.run_one_step()

    mg1 = RasterModelGrid((10,10), spacing=(1, 1))
    z1 = mg1.add_field('topographic__elevation', mg1.node_x**2 + mg1.node_y**2, at = 'node')
    fa1 = FlowAccumulator(mg1, flow_director='FlowDirectorMFD')
    fa1.run_one_step()

    mg2 = RasterModelGrid((10,10), spacing=(1, 1))
    z2 = mg2.add_field('topographic__elevation', mg2.node_x**2 + mg2.node_y**2, at = 'node')
    fa2 = FlowAccumulator(mg2, flow_director=FlowDirectorMFD)
    fa2.run_one_step()

    mg3 = RasterModelGrid((10,10), spacing=(1, 1))
    z3 = mg3.add_field('topographic__elevation', mg3.node_x**2 + mg3.node_y**2, at = 'node')
    fd = FlowDirectorMFD(mg3)
    fa3 = FlowAccumulator(mg3, flow_director=fd)
    fa3.run_one_step()

    for key in mg0.at_node.keys():
        assert_array_equal(mg0.at_node[key],
                           mg1.at_node[key])

        assert_array_equal(mg1.at_node[key],
                           mg2.at_node[key])

        assert_array_equal(mg2.at_node[key],
                           mg3.at_node[key])


def test_passing_a_bad_component():
    """Check that a random component can't be a director."""
    from landlab.components import ChiFinder

    mg = RasterModelGrid((10,10), spacing=(1, 1))
    _ = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')
    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    assert_raises(ValueError, FlowAccumulator, mg, 'topographic__elevation', flow_director=ChiFinder)


def test_error_for_to_many_with_depression():
    """Check that an error is thrown when to_many methods started DF."""

    mg0 = RasterModelGrid((10,10), spacing=(1, 1))
    z0 = mg0.add_field('topographic__elevation', mg0.node_x**2 + mg0.node_y**2, at = 'node')

    mg1 = RasterModelGrid((10,10), spacing=(1, 1))
    z1 = mg1.add_field('topographic__elevation', mg1.node_x**2 + mg1.node_y**2, at = 'node')


    assert_raises(ValueError, FlowAccumulator, mg0, flow_director='MFD', depression_finder='DepressionFinderAndRouter')
    assert_raises(ValueError, FlowAccumulator, mg0, flow_director='DINF', depression_finder='DepressionFinderAndRouter')

    fa0 = FlowAccumulator(mg0, flow_director='MFD')
    fa0.run_one_step()
    assert_raises(ValueError, DepressionFinderAndRouter, mg0)

    fa1 = FlowAccumulator(mg1, flow_director='DINF')
    fa1.run_one_step()
    assert_raises(ValueError, DepressionFinderAndRouter, mg1)


def test_fields():
    """Check to make sure the right fields have been created.

    Check that the sizes are also correct.
    """
    mg = RasterModelGrid((10,10), spacing=(1, 1))
    _ = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')
    fa = FlowAccumulator(mg)
    fa.run_one_step()

    assert_equal(sorted(list(mg.at_node.keys())), ['drainage_area',
                                                   'flow__data_structure_delta',
                                                   'flow__link_to_receiver_node',
                                                   'flow__receiver_node',
                                                   'flow__sink_flag',
                                                   'flow__upstream_node_order',
                                                   'surface_water__discharge',
                                                   'topographic__elevation',
                                                   'topographic__steepest_slope',
                                                   'water__unit_flux_in'])
    assert_equal(sorted(list(mg.at_link.keys())), ['flow__data_structure_D'])

    mg2 = RasterModelGrid((10,10), spacing=(1, 1))
    _ = mg2.add_field('topographic__elevation', mg2.node_x + mg2.node_y, at = 'node')
    fa2 = FlowAccumulator(mg2, flow_director='MFD')
    fa2.run_one_step()
    assert_equal(sorted(list(mg2.at_node.keys())), ['drainage_area',
                                                    'flow__data_structure_delta',
                                                    'flow__link_to_receiver_node',
                                                    'flow__links_to_receiver_nodes',
                                                    'flow__receiver_node',
                                                    'flow__receiver_nodes',
                                                    'flow__receiver_proportions',
                                                    'flow__sink_flag',
                                                    'flow__upstream_node_order',
                                                    'surface_water__discharge',
                                                    'topographic__elevation',
                                                    'topographic__steepest_slope',
                                                    'water__unit_flux_in'])

    assert_equal(sorted(list(mg2.at_link.keys())), ['flow__data_structure_D'])


def test_accumulated_area_closes():
    """Check that accumulated area is area of core nodes."""

    fds= ['Steepest','D8','MFD','DINF']

    for fd in fds:
        mg = RasterModelGrid((10,10), spacing=(1, 1))
        z = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')
        fa = FlowAccumulator(mg)
        fa.run_one_step()

        drainage_area = mg.at_node['drainage_area']
        drained_area = np.sum(drainage_area[mg.boundary_nodes])
        core_area = np.sum(mg.cell_area_at_node[mg.core_nodes])
        assert_equal(drained_area, core_area)

# def test_passing_unnecessary_kwarg():
#     """Test that passing a bad kwarg raises a ValueError."""
#     mg = RasterModelGrid((10,10), spacing=(1, 1))
#     z = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')
#     assert_raises(ValueError, FlowAccumulator, mg, bad_kwarg='woo')

def test_specifying_routing_method_wrong():
    """Test specifying incorrect method for routing compatability with DepressionFinderAndRouter."""
    mg = RasterModelGrid((10,10), spacing=(1, 1))
    z = mg.add_field('topographic__elevation', mg.node_x + mg.node_y, at = 'node')

    assert_raises(ValueError, FlowAccumulator, mg,
                  flow_director='D4',
                  depression_finder='DepressionFinderAndRouter')

    df = DepressionFinderAndRouter(mg)
    assert_raises(ValueError, FlowAccumulator, mg,
                  flow_director='D4',
                  depression_finder=df)


def test_field_name_array_float_case1():
    """Topography as field, runoff rate as float"""
    mg = RasterModelGrid((5,4), spacing=(1, 1))
    topographic__elevation = np.array([0.,  0.,  0., 0.,
                                       0., 21., 10., 0.,
                                       0., 31., 20., 0.,
                                       0., 32., 30., 0.,
                                       0.,  0.,  0., 0.])
    _ = mg.add_field('node', 'topographic__elevation', topographic__elevation)
    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = FlowAccumulator(mg, 'topographic__elevation', runoff_rate=10.)
    assert_array_equal(mg.at_node['water__unit_flux_in'], 10.*np.ones(mg.size('node')))

    fa.run_one_step()
    reciever = np.array([ 0,  1,  2,  3,
                          4,  1,  2,  7,
                          8, 10,  6, 11,
                         12, 14, 10, 15,
                         16, 17, 18, 19])

    da = np.array([ 0.,  1.,  5.,  0.,
                    0.,  1.,  5.,  0.,
                    0.,  1.,  4.,  0.,
                    0.,  1.,  2.,  0.,
                    0.,  0.,  0.,  0.])

    q = np.array([ 0.,  10.,  50.,  0.,
                   0.,  10.,  50.,  0.,
                   0.,  10.,  40.,  0.,
                   0.,  10.,  20.,  0.,
                   0.,   0.,   0.,  0.])

    assert_array_equal(mg.at_node['flow__receiver_node'], reciever)
    assert_array_equal(mg.at_node['drainage_area'], da)
    assert_array_equal(mg.at_node['surface_water__discharge'], q)


def test_field_name_array_float_case2():
    """Topography as field, runoff rate as field name"""
    mg = RasterModelGrid((5,4), spacing=(1, 1))
    topographic__elevation = np.array([0.,  0.,  0., 0.,
                                       0., 21., 10., 0.,
                                       0., 31., 20., 0.,
                                       0., 32., 30., 0.,
                                       0.,  0.,  0., 0.])

    runoff_rate = ([1.,  1.,  1., 1.,
                    2.,  2.,  2., 2.,
                    3.,  3.,  3., 3.,
                    4.,  4.,  4., 4.,
                    5.,  5.,  5., 5.])

    _ = mg.add_field('node', 'topographic__elevation', topographic__elevation)
    _ = mg.add_field('node', 'runoff_rate', runoff_rate)

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = FlowAccumulator(mg, 'topographic__elevation', runoff_rate='runoff_rate')

    fa.run_one_step()
    reciever = np.array([ 0,  1,  2,  3,
                          4,  1,  2,  7,
                          8, 10,  6, 11,
                         12, 14, 10, 15,
                         16, 17, 18, 19])

    da = np.array([ 0.,  1.,  5.,  0.,
                    0.,  1.,  5.,  0.,
                    0.,  1.,  4.,  0.,
                    0.,  1.,  2.,  0.,
                    0.,  0.,  0.,  0.])

    q = np.array([  0.,   2.,  16.,   0.,   # KRB double checked these numbers by hand 5/15/18 - OK
                    0.,   2.,  16.,   0.,
                    0.,   3.,  14.,   0.,
                    0.,   4.,   8.,   0.,
                    0.,   0.,   0.,   0.])

    assert_array_equal(mg.at_node['flow__receiver_node'], reciever)
    assert_array_equal(mg.at_node['drainage_area'], da)
    assert_array_equal(mg.at_node['surface_water__discharge'], q)


def test_field_name_array_float_case3():
    """Topography as field, runoff rate as array"""
    mg = RasterModelGrid((5,4), spacing=(1, 1))
    topographic__elevation = np.array([0.,  0.,  0., 0.,
                                       0., 21., 10., 0.,
                                       0., 31., 20., 0.,
                                       0., 32., 30., 0.,
                                       0.,  0.,  0., 0.])

    runoff_rate = ([1.,  1.,  1., 1.,
                    2.,  2.,  2., 2.,
                    3.,  3.,  3., 3.,
                    4.,  4.,  4., 4.,
                    5.,  5.,  5., 5.])

    _ = mg.add_field('node', 'topographic__elevation', topographic__elevation)

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = FlowAccumulator(mg, 'topographic__elevation', runoff_rate=runoff_rate)

    fa.run_one_step()
    reciever = np.array([ 0,  1,  2,  3,
                          4,  1,  2,  7,
                          8, 10,  6, 11,
                         12, 14, 10, 15,
                         16, 17, 18, 19])

    da = np.array([ 0.,  1.,  5.,  0.,
                    0.,  1.,  5.,  0.,
                    0.,  1.,  4.,  0.,
                    0.,  1.,  2.,  0.,
                    0.,  0.,  0.,  0.])

    q = np.array([  0.,   2.,  16.,   0.,   # KRB double checked these numbers by hand 5/15/18 - OK
                    0.,   2.,  16.,   0.,
                    0.,   3.,  14.,   0.,
                    0.,   4.,   8.,   0.,
                    0.,   0.,   0.,   0.])

    assert_array_equal(mg.at_node['flow__receiver_node'], reciever)
    assert_array_equal(mg.at_node['drainage_area'], da)
    assert_array_equal(mg.at_node['surface_water__discharge'], q)


def test_field_name_array_float_case4():
    """Topography as array, runoff rate as float"""
    mg = RasterModelGrid((5,4), spacing=(1, 1))
    topographic__elevation = np.array([0.,  0.,  0., 0.,
                                       0., 21., 10., 0.,
                                       0., 31., 20., 0.,
                                       0., 32., 30., 0.,
                                       0.,  0.,  0., 0.])
    _ = mg.add_field('node', 'topographic__elevation', topographic__elevation)
    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = FlowAccumulator(mg, topographic__elevation, runoff_rate=10.)
    assert_array_equal(mg.at_node['water__unit_flux_in'], 10.*np.ones(mg.size('node')))

    fa.run_one_step()
    reciever = np.array([ 0,  1,  2,  3,
                          4,  1,  2,  7,
                          8, 10,  6, 11,
                         12, 14, 10, 15,
                         16, 17, 18, 19])

    da = np.array([ 0.,  1.,  5.,  0.,
                    0.,  1.,  5.,  0.,
                    0.,  1.,  4.,  0.,
                    0.,  1.,  2.,  0.,
                    0.,  0.,  0.,  0.])

    q = np.array([ 0.,  10.,  50.,  0.,
                   0.,  10.,  50.,  0.,
                   0.,  10.,  40.,  0.,
                   0.,  10.,  20.,  0.,
                   0.,   0.,   0.,  0.])

    assert_array_equal(mg.at_node['flow__receiver_node'], reciever)
    assert_array_equal(mg.at_node['drainage_area'], da)
    assert_array_equal(mg.at_node['surface_water__discharge'], q)


def test_field_name_array_float_case5():
    """Topography as array, runoff rate as field name"""
    mg = RasterModelGrid((5,4), spacing=(1, 1))
    topographic__elevation = np.array([0.,  0.,  0., 0.,
                                       0., 21., 10., 0.,
                                       0., 31., 20., 0.,
                                       0., 32., 30., 0.,
                                       0.,  0.,  0., 0.])

    runoff_rate = ([1.,  1.,  1., 1.,
                    2.,  2.,  2., 2.,
                    3.,  3.,  3., 3.,
                    4.,  4.,  4., 4.,
                    5.,  5.,  5., 5.])

    _ = mg.add_field('node', 'topographic__elevation', topographic__elevation)
    _ = mg.add_field('node', 'runoff_rate', runoff_rate)

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = FlowAccumulator(mg, 'topographic__elevation', runoff_rate='runoff_rate')

    fa.run_one_step()
    reciever = np.array([ 0,  1,  2,  3,
                          4,  1,  2,  7,
                          8, 10,  6, 11,
                         12, 14, 10, 15,
                         16, 17, 18, 19])

    da = np.array([ 0.,  1.,  5.,  0.,
                    0.,  1.,  5.,  0.,
                    0.,  1.,  4.,  0.,
                    0.,  1.,  2.,  0.,
                    0.,  0.,  0.,  0.])

    q = np.array([  0.,   2.,  16.,   0.,   # KRB double checked these numbers by hand 5/15/18 - OK
                    0.,   2.,  16.,   0.,
                    0.,   3.,  14.,   0.,
                    0.,   4.,   8.,   0.,
                    0.,   0.,   0.,   0.])

    assert_array_equal(mg.at_node['flow__receiver_node'], reciever)
    assert_array_equal(mg.at_node['drainage_area'], da)
    assert_array_equal(mg.at_node['surface_water__discharge'], q)


def test_field_name_array_float_case6():
    """Topography as array, runoff rate as array"""
    mg = RasterModelGrid((5,4), spacing=(1, 1))
    topographic__elevation = np.array([0.,  0.,  0., 0.,
                                       0., 21., 10., 0.,
                                       0., 31., 20., 0.,
                                       0., 32., 30., 0.,
                                       0.,  0.,  0., 0.])

    runoff_rate = ([1.,  1.,  1., 1.,
                    2.,  2.,  2., 2.,
                    3.,  3.,  3., 3.,
                    4.,  4.,  4., 4.,
                    5.,  5.,  5., 5.])

    _ = mg.add_field('node', 'topographic__elevation', topographic__elevation)

    mg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    fa = FlowAccumulator(mg, topographic__elevation, runoff_rate=runoff_rate)

    fa.run_one_step()
    reciever = np.array([ 0,  1,  2,  3,
                          4,  1,  2,  7,
                          8, 10,  6, 11,
                         12, 14, 10, 15,
                         16, 17, 18, 19])

    da = np.array([ 0.,  1.,  5.,  0.,
                    0.,  1.,  5.,  0.,
                    0.,  1.,  4.,  0.,
                    0.,  1.,  2.,  0.,
                    0.,  0.,  0.,  0.])

    q = np.array([  0.,   2.,  16.,   0.,   # KRB double checked these numbers by hand 5/15/18 - OK
                    0.,   2.,  16.,   0.,
                    0.,   3.,  14.,   0.,
                    0.,   4.,   8.,   0.,
                    0.,   0.,   0.,   0.])

    assert_array_equal(mg.at_node['flow__receiver_node'], reciever)
    assert_array_equal(mg.at_node['drainage_area'], da)
    assert_array_equal(mg.at_node['surface_water__discharge'], q)
