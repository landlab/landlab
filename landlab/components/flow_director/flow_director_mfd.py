#! /usr/env/python

"""
flow_director_mfd.py: provides the component FlowDirectorMFD.

This components finds the steepest single-path steepest descent flow
directions. It is equivalent to D4 method in the special case of a raster grid
in that it does not consider diagonal links between nodes. For that capability,
use FlowDirectorD8.
"""

from landlab.components.flow_director.flow_director_to_many import(
_FlowDirectorToMany)
from landlab.components.flow_director import flow_direction_mfd
from landlab import VoronoiDelaunayGrid
from landlab import (FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY,
                     BAD_INDEX_VALUE)
import numpy


class FlowDirectorMFD(_FlowDirectorToMany):

    """Single-path (steepest direction) flow direction without diagonals.

    This components finds the steepest single-path steepest descent flow
    directions. It is equivalent to D4 method in the special case of a raster
    grid in that it does not consider diagonal links between nodes. For that
    capability, use FlowDirectorD8.

    Specifically, it stores as ModelGrid fields:

    -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
       there is no receiver: *'flow__receiver_nodes'*. This array is 2D, and is
       of dimension (number of nodes x max number of receivers).
    -  Node array of flow proportions: *'flow__receiver_proportions'*. This
       array is 2D, and is of dimension (number of nodes x max number of
       receivers).
    -  Node array of links carrying flow:  *'flow__links_to_receiver_nodes'*.
       This array is 2D, and is of dimension (number of nodes x max number of
       receivers).
    -  Node array of the steepest downhill receiver. *'flow__receiver_nodes'*
    -  Node array of steepest downhill slope from each receiver:
       *'topographic__steepest_slope'*
    -  Node array containing ID of steepest link that leads from each node to a
       receiver, or BAD_INDEX_VALUE if no link:
       *'flow__link_to_receiver_node'*
    -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class is :func:`run_one_step`.

    Examples
    --------

    This method works for both raster and irregular grids. First we will look
    at a raster example, and then an irregular example.

    >>> import numpy as numpy
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorMFD
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y,
    ...                  at = 'node')

    The MFD flow director can be uses for raster and irregular grids. For
    raster grids, use of diagonal links is specified with the keyword
    *diagonals* (default is False).

    >>> fd = FlowDirectorMFD(mg, 'topographic__elevation', diagonals = True)
    >>> fd.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> fd.run_one_step()

    Unlike flow directors that only direct flow to one node, FlowDirectorMFD
    directs flow to all downstream nodes. It stores the receiver information
    is a (number of nodes x maximum number or receivers) shape field at nodes.

    >>> mg.at_node['flow__receiver_nodes']
    array([[ 0, -1, -1, -1, -1, -1, -1, -1],
           [ 1, -1, -1, -1, -1, -1, -1, -1],
           [ 2, -1, -1, -1, -1, -1, -1, -1],
           [ 3, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1,  1, -1, -1,  0, -1],
           [ 5, -1, -1, -1, -1, -1, -1, -1],
           [ 6, -1, -1, -1, -1, -1, -1, -1],
           [ 7, -1, -1, -1, -1, -1, -1, -1],
           [ 8, -1, -1, -1, -1, -1, -1, -1]])

    It also stores the proportions of flow going to each receiver and the link
    on which the flow moves in at node arrays.

    >>> mg.at_node['flow__receiver_proportions'] # doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ],
           [ 0.        ,  0.        ,  0.        ,  0.41421356,  0.        ,
             0.        ,  0.58578644,  0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ,  0.        ,  0.        ]])
    >>> mg.at_node['flow__links_to_receiver_nodes']
    array([[-1, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1,  3, -1, -1, 12, -1],
           [-1, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1, -1]])

    Like flow directors that only direct flow to one downstream node,
    FlowDirectorMFD identifies and stores the steepest slope leading downhill
    from each node, the link carrying that flow, and the receiver receiving
    flow on that link.

    >>> mg.at_node['topographic__steepest_slope'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  1.41421356,
            0.        ,  0.        ,  0.        ,  0.        ])
    >>> mg.at_node['flow__link_to_receiver_node']
    array([-1, -1, -1, -1, 12, -1, -1, -1, -1])
    >>> mg.at_node['flow__receiver_node']
    array([0, 1, 2, 3, 0, 5, 6, 7, 8])

    Finally, FlowDirectorMFD identifies sinks, or local lows.

    >>> mg.at_node['flow__sink_flag']
    array([1, 1, 1, 1, 0, 1, 1, 1, 1], dtype=int8)

    The flow directors also have the ability to return the flow receiver nodes.
    For this example, we will turn the diagonals off. This is the default
    value.

    >>> fd = FlowDirectorMFD(mg, 'topographic__elevation')
    >>> fd.run_one_step()
    >>> receivers, proportions = fd.direct_flow()
    >>> receivers
    array([[ 0, -1, -1, -1],
           [ 1, -1, -1, -1],
           [ 2, -1, -1, -1],
           [ 3, -1, -1, -1],
           [-1, -1, -1,  1],
           [ 5, -1, -1, -1],
           [ 6, -1, -1, -1],
           [ 7, -1, -1, -1],
           [ 8, -1, -1, -1]])
    >>> proportions # doctest: +NORMALIZE_WHITESPACE
    array([[ 1.,  0.,  0.,  0.],
           [ 1.,  0.,  0.,  0.],
           [ 1.,  0.,  0.,  0.],
           [ 1.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  1.],
           [ 1.,  0.,  0.,  0.],
           [ 1.,  0.,  0.,  0.],
           [ 1.,  0.,  0.,  0.],
           [ 1.,  0.,  0.,  0.]])

    For each donor node (represented by each row) the proportions should sum to
    one.

    >>> proportions.sum(axis=1)
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

    For the second example we will use a Hexagonal Model Grid, a special type
    of Voroni Grid that has regularly spaced hexagonal cells. FlowDirectorMFD
    has multiple ways to partition flow based on slope. The default method is
    based on the slope angle. A secondary methods is to partion based on the
    square root of slope. This represents the solution to a steady kinematic
    wave.

    >>> from landlab import HexModelGrid
    >>> mg = HexModelGrid(5,3)
    >>> _ = mg.add_field('topographic__elevation',
    ...                  mg.node_x + numpy.round(mg.node_y),
    ...                  at = 'node')
    >>> fd = FlowDirectorMFD(mg, 'topographic__elevation',
    ...                      partition_method='square_root_of_slope')
    >>> fd.surface_values # doctest: +NORMALIZE_WHITESPACE
    array([ 0. ,  1. ,  2. ,
            0.5,  1.5,  2.5,  3.5,
            1. ,  2. ,  3. ,  4. ,  5. ,
            2.5,  3.5,  4.5,  5.5,
            3. ,  4. ,  5. ])
    >>> fd.run_one_step()
    >>> mg.at_node['flow__receiver_nodes']
    array([[ 0, -1, -1, -1, -1, -1],
           [ 1, -1, -1, -1, -1, -1],
           [ 2, -1, -1, -1, -1, -1],
           [ 3, -1, -1, -1, -1, -1],
           [-1, -1, -1,  3,  0,  1],
           [-1, -1, -1,  4,  1,  2],
           [ 6, -1, -1, -1, -1, -1],
           [ 7, -1, -1, -1, -1, -1],
           [-1, -1, -1,  7,  3,  4],
           [-1, -1, -1,  8,  4,  5],
           [-1, -1, -1,  9,  5,  6],
           [11, -1, -1, -1, -1, -1],
           [12, -1, -1, -1, -1, -1],
           [-1, -1, 16, 12,  8,  9],
           [-1, -1, 17, 13,  9, 10],
           [15, -1, -1, -1, -1, -1],
           [16, -1, -1, -1, -1, -1],
           [17, -1, -1, -1, -1, -1],
           [18, -1, -1, -1, -1, -1]])
    >>> mg.at_node['flow__receiver_proportions'] # doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 0.        ,  0.        ,  0.        ,  0.34108138,  0.41773767,
             0.24118095],
           [ 0.        ,  0.        ,  0.        ,  0.34108138,  0.41773767,
             0.24118095],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 0.        ,  0.        ,  0.        ,  0.34108138,  0.41773767,
             0.24118095],
           [ 0.        ,  0.        ,  0.        ,  0.34108138,  0.41773767,
             0.24118095],
           [ 0.        ,  0.        ,  0.        ,  0.34108138,  0.41773767,
             0.24118095],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 0.        ,  0.        ,  0.19431571,  0.27480391,  0.33656468,
             0.19431571],
           [ 0.        ,  0.        ,  0.19431571,  0.27480391,  0.33656468,
             0.19431571],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ],
           [ 1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
             0.        ]])
    >>> mg.at_node['flow__links_to_receiver_nodes']
    array([[-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1,  8,  3,  4],
           [-1, -1, -1,  9,  5,  6],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, 19, 12, 13],
           [-1, -1, -1, 20, 14, 15],
           [-1, -1, -1, 21, 16, 17],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, 35, 31, 25, 26],
           [-1, -1, 37, 32, 27, 28],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1]])
    >>> mg.at_node['topographic__steepest_slope'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0. ,  0. ,  0. ,  0. ,  1.5,  1.5,  0. ,  0. ,  1.5,  1.5,  1.5,
            0. ,  0. ,  1.5,  1.5,  0. ,  0. ,  0. ,  0. ])
    >>> mg.at_node['flow__link_to_receiver_node']
    array([-1, -1, -1,
           -1,  3,  5, 10,
           -1, 12, 14, 16, 22,
           24, 25, 27, 29,
           -1, 36, 38])
    >>> mg.at_node['flow__sink_flag']
    array([1, 1, 1,
           1, 0, 0, 1,
           1, 0, 0, 0, 1,
           1, 0, 0, 1,
           1, 1, 1], dtype=int8)
    """

    _name = 'FlowDirectorMFD'

    def __init__(self, grid, surface='topographic__elevation', **kwargs):
        """
        Parameters
        ----------
        grid : ModelGrid
            A grid.
        surface : field name at node or array of length node, optional
            The surface to direct flow across, default is field at node:
            topographic__self.surface_valuesation.
        partition_method: string, optional
            Method for partitioning flow. Options include 'slope' (default) and
            'square_root_of_slope'.
        """
        # unpack kwargs:
        try:
            partition_method = kwargs.pop('partition_method')
        except:
            partition_method = 'slope'
        try:
            diagonals = kwargs.pop('diagonals')
        except:
            diagonals = False

        self.method = 'MFD'
        super(FlowDirectorMFD, self).__init__(grid, surface)
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        if self._is_Voroni:
            diagonals = False
        self.updated_boundary_conditions()
        self.partition_method = partition_method
        self.diagonals = diagonals

        if self._is_Voroni == False and diagonals == False:
            self.max_receivers = 4
        if self._is_Voroni == False and diagonals == True:
            self.max_receivers = 8
        else:
            self.max_receivers = self._grid.adjacent_nodes_at_node.shape[1]

        # set the number of recievers, proportions, and receiver links with the
        # right size.
        self.receivers = grid.add_field('flow__receiver_nodes',
                                        BAD_INDEX_VALUE*numpy.ones((self._grid.number_of_nodes,
                                                                    self.max_receivers),
                                                                   dtype=int),
                                        at='node',
                                        dtype=int,
                                        noclobber=False)

        self.receiver_links = grid.add_field('flow__links_to_receiver_nodes',
                                             BAD_INDEX_VALUE*numpy.ones((self._grid.number_of_nodes,
                                                                         self.max_receivers),
                                                                        dtype=int),
                                        at='node',
                                        dtype=int,
                                        noclobber=False)

        self.proportions = grid.add_field('flow__receiver_proportions',
                                        BAD_INDEX_VALUE*numpy.ones((self._grid.number_of_nodes,
                                                                    self.max_receivers),
                                                                   dtype=float),
                                        at='node',
                                        dtype=int,
                                        noclobber=False)

    def updated_boundary_conditions(self):
        """
        Method to update FlowDirectorMFD when boundary conditions change.

        Call this if boundary conditions on the grid are updated after the
        component is instantiated.
        """
        self._active_links = self.grid.active_links
        self._activelink_tail = self.grid.node_at_link_tail[self.grid.active_links]
        self._activelink_head = self.grid.node_at_link_head[self.grid.active_links]

    def run_one_step(self):
        """
        Find flow directions and save to the model grid.

        run_one_step() checks for updated boundary conditions, calculates
        slopes on links, finds basself.surface_valuesel nodes based on the
        status at node, calculates flow directions, and saves results to the
        grid.

        An alternative to run_one_step() is direct_flow() which does the same
        things but also returns the receiver nodes not return values.
        """
        self.direct_flow()

    def direct_flow(self):
        """
        Find flow directions, save to the model grid, and return receivers.

        direct_flow() checks for updated boundary conditions, calculates
        slopes on links, finds basself.surface_valuesel nodes based on the status at node,
        calculates flow directions, saves results to the grid, and returns a
        at-node array  of receiver nodes. This array is stored in the grid at:
        grid['node']['flow__receiver_nodes']

        An alternative to direct_flow() is run_one_step() which does the same
        things but also returns a at-node array  of receiver nodes. This array
        is stored in the grid at:
        grid['node']['flow__receiver_nodes']
        """
        # step 0. Check and update BCs
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code

        # step 1. Required inumpyuts for flow_directions_MFD
        # this is where diagonals are or are not included in
        # flow direction calculations

        # Option for no diagonals (default)
        if self.diagonals == False:
            neighbors_at_node = self.grid.adjacent_nodes_at_node
            links_at_node = self.grid.links_at_node
            active_link_dir_at_node = self.grid.active_link_dirs_at_node

            # this needs to be the gradient
            link_slope = self.grid.calc_grad_at_link(self.surface_values)

        # Option with diagonals.
        else:

            # need to create a list of diagonal links since it doesn't exist.
            diag_links = numpy.sort(numpy.unique(self.grid.d8s_at_node[:, 4:]))
            diag_links = diag_links[diag_links > 0]

            # get diagonal active links (though this actually includes ALL
            # active links)
            dal = self.grid.active_d8

            # calculate graidents across diagonals
            diag_grads = numpy.zeros(diag_links.shape)
            where_active_diag = dal >= diag_links.min()
            active_diags_inds = dal[where_active_diag]-diag_links.min()
            active_diag_grads = self.grid._calculate_gradients_at_d8_active_links(self.surface_values)
            diag_grads[active_diags_inds] = active_diag_grads[where_active_diag]

            # calculate gradients on orthogonal links
            ortho_grads = self.grid.calc_grad_at_link(self.surface_values)

            # concatenate the diagonal and orthogonal grid elements
            neighbors_at_node = numpy.hstack((self.grid.adjacent_nodes_at_node,
                                              self.grid.diagonal_adjacent_nodes_at_node))
            active_link_dir_at_node = numpy.hstack((self.grid.active_link_dirs_at_node,
                                                    self.grid.active_diagonal_dirs_at_node))
            link_slope = numpy.hstack((ortho_grads,
                                       diag_grads))

            links_at_node = self.grid.d8s_at_node

        # Step 2. Find and save base level nodes.
        (baselevel_nodes, ) = numpy.where(
            numpy.logical_or(self._grid.status_at_node == FIXED_VALUE_BOUNDARY,
                             self._grid.status_at_node == FIXED_GRADIENT_BOUNDARY))

        # Calculate flow directions
        (self.receivers, self.proportions, steepest_slope,
        steepest_receiver, sink,
        receiver_links, steepest_link) = \
        flow_direction_mfd.flow_directions_mfd(self.surface_values,
                                               neighbors_at_node,
                                               links_at_node,
                                               active_link_dir_at_node,
                                               link_slope,
                                               baselevel_nodes=baselevel_nodes,
                                               partition_method=self.partition_method)

        # Save the four ouputs of this component.
        self._grid['node']['flow__receiver_nodes'][:] = self.receivers
        self._grid['node']['flow__receiver_node'][:] = steepest_receiver
        self._grid['node']['flow__receiver_proportions'][:] = self.proportions
        self._grid['node']['topographic__steepest_slope'][:] = steepest_slope
        self._grid['node']['flow__link_to_receiver_node'][:] = steepest_link
        self._grid['node']['flow__links_to_receiver_nodes'][:] = receiver_links
        self._grid['node']['flow__sink_flag'][:] = numpy.zeros_like(steepest_link,
                                                                    dtype=bool)
        self._grid['node']['flow__sink_flag'][sink] = True

        return (self.receivers, self.proportions)

    @property
    def nodes_receiving_flow(self):
        """Return the node ids of the nodes receiving flow."""
        return self._grid['node']['flow__receiver_nodes']

    @property
    def proportions_of_flow(self):
        """Return the proportion of flow going to receivers."""
        return self._grid['node']['flow__receiver_proportions']

if __name__ == '__main__':
    import doctest
    doctest.testmod()
