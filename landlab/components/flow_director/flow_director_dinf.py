#! /usr/env/python

"""
flow_director_dinf.py: provides the component FlowDirectorDINF.

Directs flow on raster grids only using the Dinfinity algorithm of
Tarboton 1997.
"""

import numpy

from landlab import (
    BAD_INDEX_VALUE,
    FIXED_GRADIENT_BOUNDARY,
    FIXED_VALUE_BOUNDARY,
    VoronoiDelaunayGrid,
)
from landlab.components.flow_director import flow_direction_dinf
from landlab.components.flow_director.flow_director_to_many import _FlowDirectorToMany


class FlowDirectorDINF(_FlowDirectorToMany):

    """Flow direction on a raster grid by the D infinity method.

    Directs flow by the D infinity method (Tarboton, 1997). Each node is
    assigned two flow directions, toward the two neighboring nodes that are on
    the steepest subtriangle. Partitioning of flow is done based on the aspect
    of the subtriangle.

    Specifically, it stores as ModelGrid fields:

    -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
       there is no receiver: *'flow__receiver_node'*. This array is 2D, and is
       of dimension (number of nodes x max number of receivers).
    -  Node array of flow proportions: *'flow__receiver_proportions'*. This
       array is 2D, and is of dimension (number of nodes x max number of
       receivers).
    -  Node array of links carrying flow:  *'flow__link_to_receiver_node'*.
       This array is 2D, and is of dimension (number of nodes x max number of
       receivers).
    -  Node array of downhill slopes corresponding to each receiver.:
       *'topographic__steepest_slope'* This array is 2D, and is
       of dimension (number of nodes x max number of receivers). If the slope is
       uphill or flat, the value is assigned zero.
    -  Boolean node array of all local lows: *'flow__sink_flag'*
    -  Link array identifing if flow goes with (1) or against (-1) the link
       direction: *'flow_link_direction'*

    The primary method of this class is :func:`run_one_step`.

    Examples
    --------

    This method works for both raster and irregular grids. First we will look
    at a raster example, and then an irregular example.

    >>> import numpy as numpy
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorDINF
    >>> mg = RasterModelGrid((4,4), xy_spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field(
    ...     'topographic__elevation',
    ...     mg.node_x**2 + mg.node_y**2,
    ...     at = 'node'
    ... )

    The DINF flow director can be uses for raster grids only.

    >>> fd = FlowDirectorDINF(mg, 'topographic__elevation')
    >>> fd.surface_values.reshape(mg.shape)
    array([[  0.,   1.,   4.,   9.],
           [  1.,   2.,   5.,  10.],
           [  4.,   5.,   8.,  13.],
           [  9.,  10.,  13.,  18.]])
    >>> fd.run_one_step()

    Unlike flow directors that only direct flow to one node or to all
    downstream nodes, FlowDirectorDINF directs flow two nodes only. It stores
    the receiver information is a (number of nodes x 2) shape field at nodes.

    >>> mg.at_node['flow__receiver_node']
    array([[ 0, -1],
           [ 1, -1],
           [ 2, -1],
           [ 3, -1],
           [ 4, -1],
           [ 0,  -1],
           [ 5,  1],
           [ 7, -1],
           [ 8, -1],
           [ 5, -1],
           [ 5, -1],
           [11, -1],
           [12, -1],
           [13, -1],
           [14, -1],
           [15, -1]])

    It also stores the proportions of flow going to each receiver, the link on
    which the flow moves in at node arrays, and the slope of each link.

    >>> mg.at_node['flow__receiver_proportions'] # doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 0.59033447,  0.40966553],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ]])
    >>> mg.at_node['flow__link_to_receiver_node']
    array([[-1, -1],
           [-1, -1],
           [-1, -1],
           [-1, -1],
           [ 3, 25],
           [24,  4],
           [ 8, 26],
           [ 9, 28],
           [14, 31],
           [11, 33],
           [32, 12],
           [16, 34],
           [21, 37],
           [18, 39],
           [19, 38],
           [20, 40]])
    >>> mg.at_node['topographic__steepest_slope'] # doctest: +NORMALIZE_WHITESPACE
    array([[ -1.00000000e+00,  -1.41421356e+00],
           [  1.00000000e+00,  -7.12763635e+02],
           [  3.00000000e+00,   1.41421356e+00],
           [  5.00000000e+00,   2.82842712e+00],
           [  1.00900000e+03,   7.12763635e+02],
           [  1.41421356e+00,   1.00000000e+00],
           [  3.00000000e+00,   2.82842712e+00],
           [  1.00400000e+03,   7.10642315e+02],
           [  1.00400000e+03,   7.12056529e+02],
           [  3.00000000e+00,   0.00000000e+00],
           [  4.24264069e+00,   3.00000000e+00],
           [  1.00100000e+03,   7.09935208e+02],
           [ -0.00000000e+00,   7.09935208e+02],
           [  1.00400000e+03,   7.07813888e+02],
           [  1.00100000e+03,   7.09935208e+02],
           [  0.00000000e+00,   7.07813888e+02]])

    Finally, FlowDirectorDINF identifies sinks, or local lows.

    >>> mg.at_node['flow__sink_flag']
    array([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=int8)

    The flow directors also have the ability to return the flow receiver nodes
    through a function called direct_flow()

    >>> fd = FlowDirectorDINF(mg, 'topographic__elevation')
    >>> fd.run_one_step()
    >>> receivers, proportions = fd.direct_flow()
    >>> receivers
    array([[ 0, -1],
           [ 1, -1],
           [ 2, -1],
           [ 3, -1],
           [ 4, -1],
           [ 0, -1],
           [ 5,  1],
           [ 7, -1],
           [ 8, -1],
           [ 5, -1],
           [ 5, -1],
           [11, -1],
           [12, -1],
           [13, -1],
           [14, -1],
           [15, -1]])
    >>> proportions # doctest: +NORMALIZE_WHITESPACE
    array([[ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 0.59033447,  0.40966553],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ],
           [ 1.        ,  0.        ]])

    For each donor node (represented by each row) the proportions should sum to
    one.

    >>> proportions.sum(axis=1)
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
            1.,  1.,  1.])
    """

    _name = "FlowDirectorDINF"

    def __init__(self, grid, surface="topographic__elevation"):
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

        self.method = "DINF"
        self.max_receivers = 2
        super(FlowDirectorDINF, self).__init__(grid, surface)
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        if self._is_Voroni:
            raise NotImplementedError(
                "FlowDirectorDINF is not implemented" " for irregular grids."
            )

        self.updated_boundary_conditions()

        # set the number of recievers, proportions, and receiver links with the
        # right size.
        self.receivers = grid.add_field(
            "flow__receiver_node",
            BAD_INDEX_VALUE
            * numpy.ones((self._grid.number_of_nodes, self.max_receivers), dtype=int),
            at="node",
            dtype=int,
            noclobber=False,
        )

        self.steepest_slope = grid.add_field(
            "topographic__steepest_slope",
            BAD_INDEX_VALUE
            * numpy.ones((self._grid.number_of_nodes, self.max_receivers), dtype=float),
            at="node",
            dtype=float,
            noclobber=False,
        )

        self.receiver_links = grid.add_field(
            "flow__link_to_receiver_node",
            BAD_INDEX_VALUE
            * numpy.ones((self._grid.number_of_nodes, self.max_receivers), dtype=int),
            at="node",
            dtype=int,
            noclobber=False,
        )

        self.proportions = grid.add_field(
            "flow__receiver_proportions",
            BAD_INDEX_VALUE
            * numpy.ones((self._grid.number_of_nodes, self.max_receivers), dtype=float),
            at="node",
            dtype=int,
            noclobber=False,
        )

    def updated_boundary_conditions(self):
        """Method to update FlowDirectorDINF when boundary conditions change.

        Call this if boundary conditions on the grid are updated after
        the component is instantiated.
        """
        self._active_links = self.grid.active_links
        self._activelink_tail = self.grid.node_at_link_tail[self.grid.active_links]
        self._activelink_head = self.grid.node_at_link_head[self.grid.active_links]

    def run_one_step(self):
        """Find flow directions and save to the model grid.

        run_one_step() checks for updated boundary conditions, calculates
        slopes on links, finds basself.surface_valuesel nodes based on the status at node,
        calculates flow directions, and saves results to the grid.

        An alternative to direct_flow() is direct_flow() which does the same
        things but also returns the receiver nodes not return values.
        """
        self.direct_flow()

    def direct_flow(self):
        """Find flow directions, save to the model grid, and return receivers.

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
        self._check_updated_bc()

        # Step 1. Find and save base level nodes.
        (baselevel_nodes,) = numpy.where(
            numpy.logical_or(
                self._grid.status_at_node == FIXED_VALUE_BOUNDARY,
                self._grid.status_at_node == FIXED_GRADIENT_BOUNDARY,
            )
        )

        # Calculate flow directions
        (
            self.receivers,
            self.proportions,
            slopes_to_receivers,
            steepest_slope,
            steepest_receiver,
            sink,
            receiver_links,
            steepest_link,
        ) = flow_direction_dinf.flow_directions_dinf(
            self._grid, self.surface_values, baselevel_nodes=baselevel_nodes
        )

        # Save the four ouputs of this component.
        self._grid["node"]["flow__receiver_node"][:] = self.receivers
        self._grid["node"]["flow__receiver_proportions"][:] = self.proportions
        self._grid["node"]["topographic__steepest_slope"][:] = slopes_to_receivers
        self._grid["node"]["flow__link_to_receiver_node"][:] = receiver_links
        self._grid["node"]["flow__sink_flag"][:] = False
        self._grid["node"]["flow__sink_flag"][sink] = True

        return (self.receivers, self.proportions)


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
