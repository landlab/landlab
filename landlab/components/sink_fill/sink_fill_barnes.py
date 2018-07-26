#!/usr/env/python

"""
fill_sinks_barnes.py

Fill sinks in a landscape to the brim, following the Barnes et al. (2014) algos.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component, BAD_INDEX_VALUE
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.components import LakeMapperBarnes
from landlab.utils.return_array import return_array_at_node
from landlab.core.messages import warning_message
from landlab.components.lake_fill.lake_fill_barnes import StablePriorityQueue
from landlab.components import LakeMapperBarnes

from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY
from landlab import CLOSED_BOUNDARY
from collections import deque
import six
import numpy as np
import heapq
import itertools

LOCAL_BAD_INDEX_VALUE = BAD_INDEX_VALUE


class SinkFillerBarnes(LakeMapperBarnes):
    """
    Uses the Barnes et al (2014) algorithms to replace pits with flats, or
    optionally to very shallow gradient surfaces to allow continued draining.

    This component is NOT intended for use iteratively as a model runs;
    rather, it is to fill in an initial topography. If you want to fill pits
    as a landscape develops, you are after the LakeMapperBarnes component.

    The locations and depths etc. of the fills will be tracked, and properties
    are provided to access this information.

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to fill.
    method : {'steepest', 'd8'}
        Whether or not to recognise diagonals as valid flow paths, if a raster.
        Otherwise, no effect.
    fill_flat : bool
        If True, pits will be filled to perfectly horizontal. If False, the new
        surface will be slightly inclined to give steepest descent flow paths
        to the outlet.
    ignore_overfill : bool
        If True, suppresses the Error that would normally be raised during
        creation of a gentle incline on a fill surface (i.e., if not
        fill_flat). Typically this would happen on a synthetic DEM where more
        than one outlet is possible at the same elevation. If True, the
        was_there_overfill property can still be used to see if this has
        occurred.
    """
    def __init__(self, grid, surface='topographic__elevation',
                 method='d8', fill_flat=True,
                 ignore_overfill=False):
        # Most of the functionality of this component is directly inherited
        # from LakeMapperBarnes, so
        super().__init__(grid, surface=surface,
                     method=method, fill_flat=fill_flat,
                     fill_surface=surface,
                     route_flow_steepest_descent=False,
                     calc_slopes=False, ignore_overfill=ignore_overfill,
                     track_lakes=True)
        # note we will always track the fills, since we're only doing this
        # once... Likewise, no need for flow routing; this is not going to
        # get used dynamically.

    def run_one_step(self):
        super().run_one_step()

    @property
    def fill_dict(self):
        """
        Return a dictionary where the keys are the outlet nodes of each filled
        area, and the values are deques of nodes within each. Items are not
        returned in ID order.
        """
        return super().fill_dict

    @property
    def fill_outlets(self):
        """
        Returns the outlet for each filled area, not necessarily in ID order.
        """
        return super().lake_outlets

    @property
    def number_of_fills(self):
        """
        Return the number of individual filled areas.
        """
        return super().number_of_lakes

    @property
    def fill_map(self):
        """
        Return an array of ints, where each filled node is labelled
        with its outlet node ID.
        Nodes not in a filled area are labelled with LOCAL_BAD_INDEX_VALUE
        (default -1).
        """
        return super().lake_map

    @property
    def fill_at_node(self):
        """
        Return a boolean array, True if the node is filled, False otherwise.
        """
        return super().lake_at_node

    @property
    def fill_depth(self):
        """Return the change in surface elevation at each node this step.
        """
        return super().lake_depth

    @property
    def fill_areas(self):
        """
        A nlakes-long array of the area of each fill. The order is the same as
        that of the keys in fill_dict, and of fill_outlets.
        """
        return super().lake_areas

    @property
    def fill_volumes(self):
        """
        A nlakes-long array of the volume of each fill. The order is the same
        as that of the keys in fill_dict, and of fill_outlets.
        """
        return super().lake_volumes






#         """
#         Examples
#         --------
#         >>> import numpy as np
#         >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
#         >>> # from landlab.components import LakeMapperBarnes
#         >>> mg = RasterModelGrid((5, 6), 1.)
#         >>> for edge in ('left', 'top', 'bottom'):
#         ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
#         >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
#         >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
#         >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
#         >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
#         >>> z_init = z.copy()
#         >>> lmb = LakeMapperBarnes(mg, method='steepest')
#         """
#         self._grid = grid
#         self._open = StablePriorityQueue()
#         self._pit = []
#         self._closed = self.grid.zeros('node', dtype=bool)
#         self._gridclosednodes = self.grid.status_at_node == CLOSED_BOUNDARY
#         gridopennodes = np.logical_not(self._gridclosednodes)
#         # close up the CLOSED_BOUNDARY permanently:
#         self._closed[self._gridclosednodes] = True
# 
#         # this component maintains its own internal count of how often it has
#         # been called. This is to enable "cheap" data access of the various
#         # available data structures without needless recalculation
#         self._runcounter = itertools.count()
#         self._runcount = -1  # not yet run
#         self._lastcountforlakemap = -1  # lake_map has not yet been called
#         self._trackingcounter = itertools.count()  # sequ counter if tracking
#         self._trackercount = -1  # tracking calls not yet run
#         self._lastcountfortracker = -1
#         self._ignore_overfill = ignore_overfill
# 
#         # get the neighbour call set up:
#         assert method in {'steepest', 'd8'}
#         if method == 'd8':
#             try:
#                 self._allneighbors = self.grid.d8s_at_node
#             except AttributeError:
#                 self._allneighbors = self.grid.adjacent_nodes_at_node
#         else:
#             self._allneighbors = self.grid.adjacent_nodes_at_node
# 
#         # A key difference from the "pure" Barnes algorithm for LL is that
#         # we must'n flood from all the edges. Instead, we can only flood from
#         # a true grid edge, i.e., only the FIXED boundary types. (Both
#         # CLOSED and LOOPED assume flow either can't get out there, or at
#         # least, there's more land in that direction that will be handled
#         # otherwise.) Note we'l add a test that there is at least some kind
#         # of outlet!!
#         self._edges = np.where(np.logical_or(
#             self.grid.status_at_node == FIXED_VALUE_BOUNDARY,
#             self.grid.status_at_node == FIXED_GRADIENT_BOUNDARY))[0]
#         if self._edges.size == 0:
#             raise ValueError(
#                 'No valid outlets found on the grid! You cannot run the ' +
#                 'filling algorithms!')
#         # and finally, close these up permanently as well (edges will always
#         # be edges...)
#         self._closed[self._edges] = True
# 
#         # now, work out what constitutes a "surface" under various input opts:
#         self._dontreroute = not route_flow_steepest_descent
# # NOTE: check the overlap with the rerouting method here.
# # do we have to still implicitly do the filling??
#         # check if we are modifying in place or not:
#         self._inplace = surface is fill_surface
#         # then
#         self._surface = return_array_at_node(grid, surface)
#         self._fill_surface = return_array_at_node(grid, surface)
#         if fill_flat:
#             self._fill_one_node = self._fill_one_node_to_flat
#         else:
#             self._fill_one_node = self._fill_one_node_to_slant
#         if track_filled_nodes:
#             self._tracking_fill = True
#             self._orig_surface = self._surface.copy()
#             self._lakefiller = LakeMapperBarnes(
#                 grid, surface=self._orig_surface, method=method, fill=True,
#                 fill_flat=fill_flat, fill_surface=surface,
#                 route_flow_steepest_descent=route_flow_steepest_descent,
#                 calc_slopes=calc_slopes, ignore_overfill=ignore_overfill)
#         else:
#             self._tracking_fill = False
# 
#     def run_one_step(self):
#         "Fills the surface to remove all pits."
#         if self._tracking_fill:
#             self._lakefiller.run_one_step()
#         else:
#             # increment the run count
#             self._runcount = next(self._runcounter)
#             # First get _fill_surface in order.
#             self._fill_surface[:] = self._surface  # surfaces begin identical
#             # note this is nice & efficent if _fill_surface is _surface
#             # now, return _closed to its initial cond, w only the
#             # CLOSED_BOUNDARY and grid draining nodes pre-closed:
#             closedq = self._closed.copy()
#             if self._dontreroute:  # i.e. algos 2 & 3
#                 for edgenode in self._edges:
#                     self._open.add_task(
#                         edgenode, priority=self._surface[edgenode])
#                 self._closed[self._edges] = True
#                 while True:
#                     try:
#                         self._fill_one_node(self.grid, )
#                     except KeyError:  # run out of nodes to fill...
#                         break

    # @property
    # def fill_dict(self):
    #     """
    #     Return a dictionary where the keys are the outlet nodes of each filled
    #     area, and the values are deques of nodes within each. Items are not
    #     returned in ID order.
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.lake_dict
    # 
    # @property
    # def fill_outlets(self):
    #     """
    #     Returns the outlet for each filled area, not necessarily in ID order.
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.lake_outlets
    # 
    # @property
    # def number_of_fills(self):
    #     """
    #     Return the number of individual filled areas.
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.number_of_lakes
    # 
    # @property
    # def fill_map(self):
    #     """
    #     Return an array of ints, where each filled node is labelled
    #     with its outlet node ID.
    #     Nodes not in a filled area are labelled with LOCAL_BAD_INDEX_VALUE
    #     (default -1).
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.lake_map
    # 
    # @property
    # def fill_at_node(self):
    #     """
    #     Return a boolean array, True if the node is filled, False otherwise.
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.lake_at_node
    # 
    # @property
    # def fill_depth(self):
    #     """Return the change in surface elevation at each node this step.
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.fill_depth
    # 
    # @property
    # def fill_areas(self):
    #     """
    #     A nlakes-long array of the area of each fill. The order is the same as
    #     that of the keys in fill_dict, and of fill_outlets.
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.lake_areas
    # 
    # @property
    # def fill_volumes(self):
    #     """
    #     A nlakes-long array of the volume of each fill. The order is the same
    #     as that of the keys in fill_dict, and of fill_outlets.
    #     """
    #     assert self._tracking_fill
    #     return self._lakefiller.lake_volumes
