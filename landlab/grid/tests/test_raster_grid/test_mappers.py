import numpy as np
from numpy.testing import assert_array_equal
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab import RasterModelGrid
import landlab.grid.mappers as maps


class TestLinkEndsToLink():

    def test_max(self):
        rmg = RasterModelGrid(4, 5)
        rmg.add_empty('node', 'values')
        node_values = rmg.at_node['values']
        node_values[:] = np.arange(rmg.number_of_nodes)

        values_at_links = maps.map_max_of_link_nodes_to_link(rmg, 'values')

        assert_array_equal(
            values_at_links,
            np.array([1,  2,  3,  4,
                      5,  6,  7,  8,  9,
                      6,  7,  8,  9,
                      10, 11, 12, 13, 14,
                      11, 12, 13, 14,
                      15, 16, 17, 18, 19,
                      16, 17, 18, 19]))

    def test_min(self):
        rmg = RasterModelGrid(4, 5)
        rmg.add_empty('node', 'values')
        node_values = rmg.at_node['values']
        node_values[:] = np.arange(rmg.number_of_nodes)

        values_at_links = maps.map_min_of_link_nodes_to_link(rmg, 'values')

        assert_array_equal(
            values_at_links,
            np.array([0,  1,  2,  3,
                      0,  1,  2,  3, 4,
                      5,  6,  7,  8,
                      5,  6,  7,  8, 9,
                      10, 11, 12, 13,
                      10, 11, 12, 13, 14,
                      15, 16, 17, 18]))


class TestNodeToLinkMappers():

    def test_to_node(self):
        rmg = RasterModelGrid((4, 5), spacing=(1., 1.))
        rmg.add_empty('node', 'values')

        node_values = rmg.at_node['values']
        node_values[:] = np.arange(rmg.number_of_nodes)

        link_values = maps.map_link_tail_node_to_link(rmg, 'values')

        assert_array_equal(link_values,
                           np.array([ 0,  1,  2,  3,  0,
                                      1,  2,  3,  4,  5,
                                      6,  7,  8,  5,  6,
                                      7,  8,  9, 10, 11,
                                     12, 13, 10, 11, 12,
                                     13, 14, 15, 16, 17,
                                     18]))
        out = np.empty_like(link_values)
        rtn = maps.map_link_tail_node_to_link(rmg, 'values', out=out)
        assert_array_equal(out, link_values)
        assert_is(rtn, out)

    def test_from_node(self):
        rmg = RasterModelGrid(4, 5, 1.)
        rmg.add_empty('node', 'values')

        node_values = rmg.at_node['values']
        node_values[:] = np.arange(rmg.number_of_nodes)

        link_values = maps.map_link_head_node_to_link(rmg, 'values')

        #link_values = rmg.at_link['values']
        assert_array_equal(link_values,
                           np.array([1,  2,  3,  4,
                                     5,  6,  7,  8,  9,
                                     6,  7,  8,  9, 
                                     10, 11, 12, 13, 14, 
                                     11, 12, 13, 14, 
                                     15, 16, 17, 18, 19,
                                     16, 17, 18, 19]))

    def test_mean_node(self):
        rmg = RasterModelGrid(4, 5, 1.)
        rmg.add_empty('node', 'values')

        node_values = rmg.at_node['values']
        node_values[:] = np.arange(rmg.number_of_nodes)

        link_values = maps.map_mean_of_link_nodes_to_link(rmg, 'values')

        assert_array_equal(link_values,
                           np.array([0.5,  1.5,  2.5,  3.5,  2.5,
                                     3.5,  4.5,  5.5,  6.5,  5.5,
                                     6.5,  7.5,  8.5,  7.5,  8.5,
                                     9.5, 10.5, 11.5, 10.5,
                                     11.5, 12.5, 13.5, 12.5,
                                     13.5, 14.5, 15.5, 16.5,
                                     15.5, 16.5, 17.5, 18.5]))

    def test_cell(self):
        rmg = RasterModelGrid(4, 5, 1.)
        rmg.add_empty('node', 'values')

        node_values = rmg.at_node['values']
        node_values[:] = np.arange(rmg.number_of_nodes)

        cell_values = maps.map_node_to_cell(rmg, 'values')

        assert_array_equal(np.array([6., 7., 8., 11., 12., 13.]), cell_values)
