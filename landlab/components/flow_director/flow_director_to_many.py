#! /usr/env/python

"""
flow_director_to_many.py provides a private class to help create FlowDirectors.

Provides the _FlowDirectorToMany component which makes sure all model grid
fields are set up correctly.
"""

from landlab import FieldError
from landlab.components.flow_director.flow_director import _FlowDirector
import numpy
from landlab import BAD_INDEX_VALUE


class _FlowDirectorToMany(_FlowDirector):

    """
    Private class for creating components to calculate flow directions.

    This class is not meant to be used directly in modeling efforts. It
    inherits from the _FlowDirector class and builds on it to provide the
    functionality that all flow direction calculators need if they direct flow
    only to multiple nodes, as in  D infinity or MFD direction finding. It 
    exists in contrast to the other intermediate flow director class 
    _FlowDirectorToOne which provides equivalent functionality for flow 
    direction algorithms such as D8 or steepest descent which directs flow only
    to one other node. As the primary difference between these two methods is 
    the names of the fields they create and use, the primary function of this 
    class is to create model grid fields.

    Specifically, it stores as ModelGrid fields:

    -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
       there is no receiver: *'flow__receiver_nodes'*
    -  Node array of flow proportion: *'flow__receiver_proportions'*
    -  Node array of steepest downhill slope from each reciever:
       *'topographic__steepest_slope'*
    -  Node array containing ID of steepest link that leads from each node to a
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
    ['flow__link_to_receiver_node', 'flow__receiver_node',
           'flow__sink_flag', 'topographic__elevation',
           'topographic__steepest_slope']
    """

    _name = 'FlowDirectorToMany'

    _input_var_names = ('topographic__elevation',
                        )

    _output_var_names = ('flow__receiver_nodes',
                         'flow__receiver_proportions'
                         'topographic__steepest_slope',
                         'flow__link_to_receiver_node',
                         'flow__sink_flag',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'flow__receiver_nodes': '-',
                  'flow__receiver_proportions': '-',
                  'topographic__steepest_slope': '-',
                  'flow__link_to_receiver_node': '-',
                  'flow__sink_flag': '-',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_nodes': 'node',
                    'flow__receiver_proportions': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow__link_to_receiver_node': 'node',
                    'flow__sink_flag': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'flow__receiver_nodes':
            'Node array of receivers (nodes that receives flow from current '
            'node). This array is of dimension (number of nodes x max number '
            'of receivers.',
         'flow__receiver_proportions':
            'Node array of proportion of flow sent from current node to '
            'downstream nodes.',
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__sink_flag': 'Boolean array, True at local lows',
    }

    def __init__(self, grid, surface):
        """Initialize the _FlowDirectorTo_One class."""
        # run init for the inherited class
        super(_FlowDirectorToMany, self).__init__(grid, surface)
        self.to_n_receivers = 'many'
        # initialize new fields
        try:
            self.receivers = grid.add_field('flow__receiver_nodes',
                                            BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                            at='node', dtype=int)
        except FieldError:
            self.receiver = grid.at_node['flow__receiver_nodes']
            
        try:
            self.proportions = grid.add_field('flow__receiver_proportions',
                                            BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                            at='node', dtype=int)
        except FieldError:
            self.proportions = grid.at_node['flow__receiver_proportions']

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

    # set properties. These are the same for all DirectToMany Directors
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
