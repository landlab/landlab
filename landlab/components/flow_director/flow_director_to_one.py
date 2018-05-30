#! /usr/env/python

"""
flow_director_to_one.py provides a private class to help create FlowDirectors.

Provides the _FlowDirectorToOne component which makes sure all model grid
fields are set up correctly.
"""

from landlab import FieldError
from landlab.components.flow_director.flow_director import _FlowDirector
import numpy
from landlab import BAD_INDEX_VALUE


class _FlowDirectorToOne(_FlowDirector):

    """
    Private class for creating components to calculate flow directions.

    This class is not meant to be used directly in modeling efforts. It
    inherits from the _FlowDirector class and builds on it to provide the
    functionality that all flow direction calculators need if they direct flow
    only to one nodes, as in steepest descent or D8 direction finding. It
    exists in contrast to the other intermediate flow director class
    _FlowDirectorToMany which provides equivalent functionality for flow
    direction algorithms such as D infinity or D trig that route flow from one
    cell to multiple nodes. As the primary difference between these two methods
    is the names of the fields they create and use, the primary function of
    this class is to create model grid fields.

    Specifically, it stores as ModelGrid fields:

    -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
       there is no receiver: *'flow__receiver_node'*
    -  Node array of steepest downhill slopes:
       *'topographic__steepest_slope'*
    -  Node array containing ID of link that leads from each node to its
       receiver, or BAD_INDEX_VALUE if no link:
       *'flow__link_to_receiver_node'*
    -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class, :func:`run_one_step` is not implemented.

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to direct flow across.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_director.flow_director_to_one import(
    ... _FlowDirectorToOne)
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y,
    ...                  at = 'node')
    >>> fd = _FlowDirectorToOne(mg, 'topographic__elevation')
    >>> fd.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> sorted(list(mg.at_node.keys()))
    ['flow__link_to_receiver_node',
     'flow__receiver_node',
     'flow__sink_flag',
     'topographic__elevation',
     'topographic__steepest_slope']
    """

    _name = 'FlowDirectorToOne'

    _input_var_names = ('topographic__elevation',
                        )

    _output_var_names = ('flow__receiver_node',
                         'topographic__steepest_slope',
                         'flow__link_to_receiver_node',
                         'flow__sink_flag',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'flow__receiver_node': '-',
                  'topographic__steepest_slope': '-',
                  'flow__link_to_receiver_node': '-',
                  'flow__sink_flag': '-',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_node': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow__link_to_receiver_node': 'node',
                    'flow__sink_flag': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__sink_flag': 'Boolean array, True at local lows',
    }

    def __init__(self, grid, surface):
        """Initialize the _FlowDirectorTo_One class."""
        # run init for the inherited class
        super(_FlowDirectorToOne, self).__init__(grid, surface)
        self.to_n_receivers = 'one'
        # initialize new fields
        try:
            self.receiver = grid.add_field('flow__receiver_node',
                                           BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                           at='node', dtype=int)
        except FieldError:
            self.receiver = grid.at_node['flow__receiver_node']

        try:
            self.steepest_slope = grid.add_zeros(
                'topographic__steepest_slope', at='node', dtype=float)
        except FieldError:
            self.steepest_slope = grid.at_node['topographic__steepest_slope']

        try:
            self.links_to_receiver = grid.add_field('flow__link_to_receiver_node',
                                                    BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                                    at='node', dtype=int)

        except FieldError:
            self.links_to_receiver = grid.at_node[
                'flow__link_to_receiver_node']

        grid.add_zeros('flow__sink_flag', at='node', dtype=numpy.int8,
                       noclobber=False)

    def run_one_step(self):
        """run_one_step is not implemented for this component."""
        raise NotImplementedError('run_one_step()')

    def _determine_link_directions(self):
        """Determine link directions and set flow__link_direction field.

        This routine is slightly different between the route-to-one and
        route-to-many methods.

        An additional issue may be if this also needs to be updated by
        depression finding.
        """
        # start by re-setting all links to zero.
        self.flow__link_direction[:] = 0

        # grid.link_dirs_at_node: 0 = NA, -1 = out of (this node is HEAD)
        # -1 = into (this node is TAIL)

        # grid.links_at_node = Link ID at node

        # grid.node_at_link_head = node ID at link HEAD
        # grid.node_at_link_tail = node ID at link tail

        # links that are used are in `flow__link_to_receiver_node`.
        # if `flow__link_to_receiver_node` is leading to the link HEAD node,
        # then -1, otherwise 1 (double check this is not the reverse).

        # identify where flow is active on links
        is_active_flow_link = self.links_to_receiver != BAD_INDEX_VALUE

        # make an array that says which link ID is active
        active_flow_links = self.links_to_receiver[is_active_flow_link]

        # for each of those links, the position is the upstream node
        upstream_node_of_active_flow_link = numpy.where(is_active_flow_link)[0]

        # get the head node
        head_node_at_active_flow_link = self._grid.node_at_link_head[active_flow_links]

        # if head node is upstream node = -1, else 1
        self.flow__link_direction[active_flow_links[head_node_at_active_flow_link == upstream_node_of_active_flow_link]] = -1
        self.flow__link_direction[active_flow_links[head_node_at_active_flow_link != upstream_node_of_active_flow_link]] = 1


    # set properties. These are the same for all DirectToOne Directors
    @property
    def node_receiving_flow(self):
        """Return the node id of the node receiving flow."""
        return self._grid['node']['flow__receiver_node']

    @property
    def node_steepest_slope(self):
        """Return the steepest link slope at a node."""
        return self._grid['node']['topographic__steepest_slope']

    @property
    def link_to_flow_receiving_node(self):
        """Return the link id along the link transporting flow."""
        return self._grid['node']['flow__link_to_receiver_node']

    @property
    def sink_flag(self):
        """Return the array with sink flags."""
        return self._grid['node']['flow__sink_flag']


if __name__ == '__main__':
    import doctest
    doctest.testmod()
