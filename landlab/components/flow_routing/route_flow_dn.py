#! /usr/env/python

"""
route_flow_dn.py:
    
Given a ModelGrid, calculates single-path (steepest direction) flow directions,
drainage area, and (optionally) discharge.

The "dn" in the name means that this is a generalization of the D8 algorithm,
for a grid in which a node has N neighbors (N might happen to be 8, or not).

Created GT Nov 2013
"""

import landlab
from landlab import RasterModelGrid
from landlab.components.flow_routing import flow_direction_DN
reload(flow_direction_DN)
from landlab.components.flow_accum import flow_accum_bw
import numpy


class FlowRouter():
    
    def __init__(self, model_grid):
        
        # We keep a local reference to the grid
        self._grid = model_grid
        
        # We'll also keep track of the active links; if raster, then these are
        # the "D8" links; otherwise, it's just activelinks
        if type(model_grid) is landlab.grid.raster.RasterModelGrid:
            dal, d8f, d8t = model_grid.d8_active_links()
            self._activelink_from = d8f
            self._activelink_to = d8t
        else:
            self._activelink_from = model_grid.activelink_fromnode
            self._activelink_to = model_grid.activelink_tonode
            
        # Keep track of the following variables:
        #   - drainage area at each node
        #   - receiver of each node
        self.drainage_area = model_grid.add_zeros('node', 'drainage_area')
        self.receiver = model_grid.create_node_array_zeros('receiver')
        self.steepest_slope = model_grid.create_node_array_zeros('steepest_slope')
        
    def route_flow(self, elevs, node_cell_area=1.0, runoff_rate=1.0):
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
        
        Returns:
            - Node array of receivers (nodes that receive flow)
            - Node array of drainage areas
            - Node array of discharges
            - Node array of steepest downhill slopes
            - Node array containing downstream-to-upstream ordered list of node
              IDs
        
        Example:
            >>> mg = RasterModelGrid(5, 4, 1.0)
            >>> elev = numpy.array([0.,  0.,  0., 0., \
                                 0., 21., 10., 0., \
                                 0., 31., 20., 0., \
                                 0., 32., 30., 0., \
                                 0., 0., 0., 0.])
            >>> mg.set_inactive_boundaries(False, True, True, True)
            >>> fr = FlowRouter(mg)
            >>> r, a, q, ss, s = fr.route_flow(elev)
            >>> r
            array([ 0,  1,  2,  3,  4,  1,  2,  7,  8,  6,  6, 11, 12, 10, 10, 15, 16,
                   17, 18, 19])
            >>> a
            array([ 1.,  2.,  6.,  1.,  1.,  1.,  5.,  1.,  1.,  1.,  3.,  1.,  1.,
                    1.,  1.,  1.,  1.,  1.,  1.,  1.])
            >>> r, a, q, ss, s = fr.route_flow(elev, 10.0)
            >>> a
            array([ 10.,  20.,  60.,  10.,  10.,  10.,  50.,  10.,  10.,  10.,  30.,
                    10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.])
            >>> cell_areas = 10.0 + numpy.arange(mg.number_of_nodes)
            >>> r, a, q, ss, s = fr.route_flow(elev, cell_areas)
            >>> a
            array([  10.,   26.,  114.,   13.,   14.,   15.,  102.,   17.,   18.,
                     19.,   67.,   21.,   22.,   23.,   24.,   25.,   26.,   27.,
                     28.,   29.])
            >>> runoff_rate = numpy.arange(mg.number_of_nodes)
            >>> r, a, q, ss, s = fr.route_flow(elev, 100.0, runoff_rate)
            >>> q
            array([    0.,   600.,  5400.,   300.,   400.,   500.,  5200.,   700.,
                     800.,   900.,  3700.,  1100.,  1200.,  1300.,  1400.,  1500.,
                    1600.,  1700.,  1800.,  1900.])
            >>> r, a, q, ss, s = fr.route_flow(elev, cell_areas, runoff_rate)
            >>> q
            array([    0.,    86.,  1126.,    39.,    56.,    75.,  1102.,   119.,
                     144.,   171.,   835.,   231.,   264.,   299.,   336.,   375.,
                     416.,   459.,   504.,   551.])
        """
        
        # Calculate the downhill-positive slopes at the d8 active links
        link_slope = -self._grid.calculate_gradients_at_d8_active_links(elevs)
        
        # Find the baselevel nodes
        (baselevel_nodes, ) = numpy.where( self._grid.node_status==1 )

        # Calculate flow directions
        receiver, steepest_slope, sink = flow_direction_DN.flow_directions(
                                         elevs, self._activelink_from,
                                         self._activelink_to, link_slope, 
                                         baselevel_nodes)

        # Calculate drainage area, discharge, and ...
        a, q, s = flow_accum_bw.flow_accumulation(receiver, sink,
                                                  node_cell_area, runoff_rate)

        return receiver, a, q, steepest_slope, s


if __name__ == '__main__':
    import doctest
    doctest.testmod()
