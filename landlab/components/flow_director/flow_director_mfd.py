#! /usr/env/python

"""
flow_director_steepest.py: provides the component FlowDirectorsSteepest.

This components finds the steepest single-path steepest descent flow
directions. It is equivalent to D4 method in the special case of a raster grid
in that it does not consider diagonal links between nodes. For that capability,
use FlowDirectorD8.
"""

from landlab.components.flow_director.flow_director_to_many import(
_FlowDirectorToMany)
from landlab.components.flow_director import flow_direction_mfd
from landlab import VoronoiDelaunayGrid
from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY
import numpy


class FlowDirectorMFD(_FlowDirectorToMany):

    """
    Single-path (steepest direction) flow direction without diagonals.

    This components finds the steepest single-path steepest descent flow
    directions. It is equivalent to D4 method in the special case of a raster
    grid in that it does not consider diagonal links between nodes. For that
    capability, use FlowDirectorD8.

     Stores as ModelGrid fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_node'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Node array containing ID of link that leads from each node to its
           receiver, or BAD_INDEX_VALUE if no link:
           *'flow__link_to_receiver_node'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class is :func:`run_one_step`.

    Construction::

        FlowDirectorMFD(grid, surface='topographic__self.surface_valuesation')

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

    Examples
    --------

    This method works for both raster and irregular grids. First we will look
    at a raster example, and then an irregular example.

    >>> import numpy as numpy
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorMFD
    >>> mg = RasterModelGrid((3,3), spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field('topographic__self.surface_valuesation',
    ...                  mg.node_x + mg.node_y,
    ...                  at = 'node')
    >>> fd = FlowDirectorMFD(mg, 'topographic__self.surface_valuesation')
    >>> fd.surface_values
    array([ 0.,  1.,  2.,  1.,  2.,  3.,  2.,  3.,  4.])
    >>> fd.run_one_step()
    >>> mg.at_node['flow__receiver_nodes']
    array([0, 1, 2, 3, 1, 5, 6, 7, 8])
    >>> mg.at_node['topographic__steepest_slope']
    array([ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.])
    >>> mg.at_node['flow__link_to_receiver_node']
    array([-1, -1, -1, -1,  3, -1, -1, -1, -1])
    >>> mg.at_node['flow__sink_flag']
    array([1, 1, 1, 1, 0, 1, 1, 1, 1], dtype=int8)
    >>> mg_2 = RasterModelGrid((5, 4), spacing=(1, 1))
    >>> topographic__self.surface_valuesation = numpy.array([0.,  0.,  0., 0.,
    ...                                    0., 21., 10., 0.,
    ...                                    0., 31., 20., 0.,
    ...                                    0., 32., 30., 0.,
    ...                                    0.,  0.,  0., 0.])
    >>> _ = mg_2.add_field('node',
    ...                    'topographic__self.surface_valuesation',
    ...                    topographic__self.surface_valuesation)
    >>> mg_2.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fd_2 = FlowDirectorSteepest(mg_2)
    >>> fd_2.run_one_step()
    >>> mg_2.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  1,  2,  3,
            4,  1,  2,  7,
            8, 10,  6, 11,
           12, 14, 10, 15,
           16, 17, 18, 19])

    The flow directors also have the ability to return the flow receiver nodes

    >>> receiver = fd.direct_flow()
    >>> receiver
    array([0, 1, 2,
           3, 1, 5,
           6, 7, 8])

    For the second example we will use a Hexagonal Model Grid, a special type
    of Voroni Grid that has regularly spaced hexagonal cells.

    >>> from landlab import HexModelGrid
    >>> mg = HexModelGrid(5,3)
    >>> _ = mg.add_field('topographic__self.surface_valuesation',
    ...                  mg.node_x + numpy.round(mg.node_y),
    ...                  at = 'node')
    >>> fd = FlowDirectorSteepest(mg, 'topographic__self.surface_valuesation')
    >>> fd.surface_values
    array([ 0. ,  1. ,  2. ,
        0.5,  1.5,  2.5,  3.5,
      1. ,  2. ,  3. ,  4. , 5. ,
        2.5,  3.5,  4.5,  5.5,
            3. ,  4. ,  5. ])
    >>> fd.run_one_step()
    >>> mg.at_node['flow__receiver_node']
    array([ 0,  1,  2,
          3,  0,  1,  6,
        7,  3,  4,  5,  11,
          12,  8,  9, 15,
            16, 17, 18])
    >>> mg.at_node['topographic__steepest_slope']
    array([ 0. ,  0. ,  0. ,
        0. ,  1.5,  1.5,   0. ,
      0. ,  1.5,  1.5,  1.5,  0. ,
        0. ,  1.5,  1.5,  0. ,
            0. ,  0. ,  0. ])
    >>> mg.at_node['flow__link_to_receiver_node']
    array([-1, -1, -1,
         -1,  3,  5, -1,
       -1, 12, 14, 16, -1,
         -1, 25, 27, -1,
           -1, -1, -1])
    >>> mg.at_node['flow__sink_flag']
    array([1, 1, 1,
          1, 0, 0, 1,
         1, 0, 0, 0, 1,
          1, 0, 0, 1,
            1, 1, 1], dtype=int8)
    >>> receiver = fd.direct_flow()
    >>> receiver
    array([ 0,  1,  2,
          3,  0,  1,  6,
        7,  3,  4,  5, 11,
         12,  8,  9, 15,
          16, 17, 18])
    """

    _name = 'FlowDirectorMFD'

    def __init__(self, grid, surface='topographic__self.surface_valuesation',
                 partition_method='slope', diagonals=False):
        """
        Initialize FlowDirectorMFD
        """
        self.method = 'MFD'
        super(FlowDirectorMFD, self).__init__(grid, surface)
        self._is_Voroni = isinstance(self._grid, VoronoiDelaunayGrid)
        self.updated_boundary_conditions()
        self.partition_method =  partition_method
        self.diagonals = diagonals
        
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
        slopes on links, finds basself.surface_valuesel nodes based on the status at node,
        calculates flow directions, and saves results to the grid.

        An alternative to direct_flow() is direct_flow() which does the same
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
        if self.diagonals == False:
            neighbors_at_node = self.grid.neighbors_at_node
            links_at_node = self.grid.links_at_node
            active_link_dir_at_node = self.grid.active_link_dirs_at_node
            
            # this needs to change from the gradient to the slope. 
            link_slope = numpy.arctan(self.grid.calc_grad_at_link(self.surface_values))
            
    
        else: # if diagonals are true
            self.grid._create_diag_links_at_node()
            dal, d8t, d8h = self.grid._d8_active_links()
                    
            neighbors_at_node = numpy.hstack((self.grid.neighbors_at_node, 
                                              self.grid._diagonal_neighbors_at_node))
            links_at_node = numpy.hstack((self.grid.links_at_node,
                                          self.grid._diagonal_links_at_node))
            active_link_dir_at_node = numpy.hstack((self.grid.active_link_dirs_at_node,
                                                 self.grid._diag__active_link_dirs_at_node))
            
            # need to create a list of diagonal links since it doesn't exist. 
            diag_links = numpy.sort(numpy.unique(self.grid._diag_links_at_node))
            diag_links = diag_links[diag_links>0]
            
            # calculate graidents across diagonals
            diag_grads = numpy.zeros(diag_links.shape)
            where_active_diag = dal>=diag_links.min()
            active_diags_inds = dal[where_active_diag]-diag_links.min()
        
            active_diag_grads = self.grid._calculate_gradients_at_d8_active_links(self.surface_values)
            
            diag_grads[active_diags_inds] = active_diag_grads[where_active_diag]
            ortho_grads = self.grid.calc_grad_at_link(self.surface_values)
                                 
            link_slope = numpy.hstack((numpy.arctan(ortho_grads),
                                                     numpy.arctan(diag_grads)))
        # Step 2. Find and save base level nodes.
        (basself.surface_valuesel_nodes, ) = numpy.where(
            numpy.logical_or(self._grid.status_at_node == FIXED_VALUE_BOUNDARY,
                             self._grid.status_at_node == FIXED_GRADIENT_BOUNDARY))

        # Calculate flow directions
        (receivers, proportions, steepest_slope, 
        steepest_receiver, sink, 
        receiver_links, steepest_link)= \
        flow_direction_mfd.flow_directions_mfd(self.surface_values, 
                                               neighbors_at_node,
                                               links_at_node,
                                               active_link_dir_at_node,
                                               link_slope, 
                                               basself.surface_valuesel_nodes=None,
                                               partition_method='slope')

        # Save the four ouputs of this component.
        self._grid['node']['flow__receiver_nodes'][:] = receivers
        self._grid['node']['flow__receiver_proportions'][:] = proportions
        self._grid['node']['topographic__steepest_slope'][:] = steepest_slope
        self._grid['node']['flow__link_to_receiver_node'][:] = steepest_link
        self._grid['node']['flow__sink_flag'][:] = numpy.zeros_like(receiver,
                                                                    dtype=bool)
        self._grid['node']['flow__sink_flag'][sink] = True

        return (receivers, proportions)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
