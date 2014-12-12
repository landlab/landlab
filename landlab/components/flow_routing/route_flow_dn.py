#! /usr/env/python

"""
route_flow_dn.py:
    
Given a ModelGrid, calculates single-path (steepest direction) flow directions,
drainage area, and (optionally) discharge.

The "dn" in the name means that this is a generalization of the D8 algorithm,
for a grid in which a node has N neighbors (N might happen to be 8, or not).

Created GT Nov 2013
Modified to save data to grid directly, DEJH March 2014
"""

import landlab
#from landlab import RasterModelGrid
from landlab.components.flow_routing import flow_direction_DN
#reload(flow_direction_DN)
from landlab.components.flow_accum import flow_accum_bw
import numpy
from scipy import weave
#from scipy.weave.build_tools import CompileError

#output_suppression_flag = True

class FlowRouter():
    """
    This class implements single-path (steepest direction) flow routing, and 
    calculates flow directions, drainage area, and (optionally) discharge. 
    
    It initializes with a reference to a ModelGrid of any kind. Optionally, it
    can also take *value_field*, the string which is the name of the elvation 
    field in the model grid to use to route flow. If *value_field* is not
    provided, it defaults to 'topographic_elevation'.
    
    The primary method of this class is :func:`route_flow`.
    """
    
    def __init__(self, model_grid, value_field='topographic_elevation'):
        
        # We keep a local reference to the grid
        self._grid = model_grid
        self.value_field = value_field
        
        # We'll also keep track of the active links; if raster, then these are
        # the "D8" links; otherwise, it's just activelinks
        if type(model_grid) is landlab.grid.raster.RasterModelGrid:
            dal, d8f, d8t = model_grid.d8_active_links()
            self._active_links = dal
            self._activelink_from = d8f
            self._activelink_to = d8t
        else:
            self._active_links = model_grid.active_links
            self._activelink_from = model_grid.activelink_fromnode
            self._activelink_to = model_grid.activelink_tonode
            
        # Keep track of the following variables:
        #   - drainage area at each node
        #   - receiver of each node
        self.drainage_area = model_grid.add_zeros('node', 'drainage_area')
        self.receiver = model_grid.create_node_array_zeros('flow_receiver')
        self.steepest_slope = model_grid.create_node_array_zeros('steepest_slope')
        self.discharges = model_grid.create_node_array_zeros('water_discharges')
        self.upstream_ordered_nodes = model_grid.create_node_array_zeros('upstream_ID_order')
        self.links_to_receiver = model_grid.create_node_array_zeros('links_to_flow_receiver')
        
        self.weave_flag = model_grid.weave_flag
        
        
    def route_flow(self, elevs=None, grid=None, runoff_rate=1.0,
                   boundary_nodes=None):
        """
        Routes surface-water flow by (1) assigning to each node a single 
        drainage direction, and then (2) adding up the number of nodes that
        contribute flow to each node on the grid (including the node itself).
        If a scalar is specified for cell_area, computes the total surface 
        contributing area by assuming that each cell has the same surface area.
        If an array is given (with length equal to the number of nodes), these
        areas are used for each cell. Likewise, runoff_rate, in length per
        time (volume per area per time) may be given either as a scalar
        (identical for each cell) or as an array whose length is the number of
        nodes in the grid.
        
        Takes:
            - Either *elevs*, an array of node elevations, or *grid*, a 
              reference to a ModelGrid.
              
        Takes as optional inputs:
            - *runoff_rate*, a float (for constant rainfall) or array (for
              spatially variable rainfall) of runoff rates, such that drainage
              area is in volume, rather than number of upstream cells.
             
            - Note that this module **no longer** takes *node_cell_area* as an
              optional input. Node cell area will always be read from the
              grid supplied.
        
        Stores as ModelGrid fields, or returns, if *elevs* was provided rather
        than *grid*:
            - Node array of receivers (nodes that receive flow): 
              *'flow_receiver'*
            - Node array of drainage areas: *'drainage_area'*
            - Node array of discharges: *'water_discharges'*
            - Node array of steepest downhill slopes: *'steepest_slope'*
            - Node array containing downstream-to-upstream ordered list of node
              IDs: *'upstream_ID_order'*
            - Node array containing ID of link that leads from each node to its
              receiver (or ITS OWN ID if there is no receiver):
              *'links_to_flow_receiver'*
        
        Returns, if *grid* was provided:
            - the modified grid object
        
        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing.route_flow_dn import FlowRouter
        >>> mg = RasterModelGrid(5, 4, 1.0)
        >>> elev = np.array([0.,  0.,  0., 0.,
        ...                  0., 21., 10., 0.,
        ...                  0., 31., 20., 0.,
        ...                  0., 32., 30., 0.,
        ...                  0.,  0.,  0., 0.])
        >>> mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
        >>> fr = FlowRouter(mg)
        >>> r, a, q, ss, s, rl = fr.route_flow(elevs=elev)
        >>> r
        array([ 0,  1,  2,  3,  4,  1,  2,  7,  8,  6,  6, 11, 12, 10, 10, 15, 16,
               17, 18, 19])
        >>> a
        array([ 1.,  2.,  6.,  1.,  1.,  1.,  5.,  1.,  1.,  1.,  3.,  1.,  1.,
                1.,  1.,  1.,  1.,  1.,  1.,  1.])

        Now let's change the cell area and the runoff rates:
        
        >>> mg = RasterModelGrid(5, 4, 10.) #so cell area==100.
        >>> mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
        >>> fr = FlowRouter(mg)
        >>> runoff_rate = np.arange(mg.number_of_nodes)
        >>> r, a, q, ss, s, rl = fr.route_flow(elevs=elev, runoff_rate=runoff_rate)
        >>> q
        array([    0.,   600.,  5400.,   300.,   400.,   500.,  5200.,   700.,
                 800.,   900.,  3700.,  1100.,  1200.,  1300.,  1400.,  1500.,
                1600.,  1700.,  1800.,  1900.])
        
        """
        
        #if elevs is not provided, default to stored grid values, which must be provided as grid
        if elevs is None:
            if grid is not None:
                self._grid = grid
                elevs = grid['node'][self.value_field]
            else:
                raise ValueError('Either an elevation array or a copy of the grid must be provided!')
        
        node_cell_area = self._grid.forced_cell_areas
            
        
        # Calculate the downhill-positive slopes at the d8 active links
        #TODO: generalize to use EITHER D8, if raster, or just active links,
        # otherwise.
        link_slope = -self._grid.calculate_gradients_at_d8_active_links(elevs)
        # Find the baselevel nodes
        (baselevel_nodes, ) = numpy.where(numpy.logical_or(self._grid.node_status==1, self._grid.node_status==2))

        # Calculate flow directions
        receiver, steepest_slope, sink, recvr_link  = \
            flow_direction_DN.flow_directions(elevs, self._active_links, 
                                         self._activelink_from,
                                         self._activelink_to, link_slope, 
                                         grid=grid,
                                         baselevel_nodes=baselevel_nodes, 
                                         use_weave=self.weave_flag)
#############grid=None???
        
        # TODO: either need a way to calculate and return the *length* of the
        # flow links, OR the caller has to handle the raster / non-raster case.
        
        #print 'sinks:', sink

        # Calculate drainage area, discharge, and ...
        a, q, s = flow_accum_bw.flow_accumulation(receiver, sink,
                                                  node_cell_area, runoff_rate,
                                                  boundary_nodes, self.weave_flag)
                                                  
        #added DEJH March 2014:
        #store the generated data in the grid
        self._grid['node']['drainage_area'] = a
        self._grid['node']['flow_receiver'] = receiver
        self._grid['node']['steepest_slope'] = steepest_slope
        self._grid['node']['water_discharges'] = q
        self._grid['node']['upstream_ID_order'] = s
        self._grid['node']['links_to_flow_receiver'] = recvr_link
        self._grid['node']['flow_sinks'] = numpy.zeros_like(receiver, dtype=bool)
        self._grid['node']['flow_sinks'][sink] = True
        

        if grid:
            return self._grid
        else:
            return receiver, a, q, steepest_slope, s, recvr_link

    @property
    def node_drainage_area(self):
        return self._grid['node']['drainage_area']

    @property
    def node_receiving_flow(self):
        return self._grid['node']['flow_receiver']

    @property
    def node_steepest_slope(self):
        return self._grid['node']['steepest_slope']

    @property
    def node_water_discharge(self):
        return self._grid['node']['water_discharges']

    @property
    def node_order_upstream(self):
        return self._grid['node']['upstream_ID_order']

    @property
    def link_to_flow_receiving_node(self):
        return self._grid['node']['links_to_flow_receiver']


if __name__ == '__main__':
    import doctest
    doctest.testmod()
