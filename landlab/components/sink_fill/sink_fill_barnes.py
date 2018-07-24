#!/usr/env/python

"""
fill_sinks_barnes.py

Fill sinks in a landscape to the brim, following the Barnes et al. (2014) algos.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component, BAD_INDEX_VALUE
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
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


def _fill_one_node_to_flat(grid, fill_surface, all_neighbors,
                           pitq, openq, closedq):
    """
    Implements the Barnes et al. algorithm for a simple fill. Assumes the
    _open and _closed lists have already been updated per Barnes algos 2&3,
    lns 1-7.

    Parameters
    ----------
    grid : ModelGrid
        In this function, behaves as a dummy parameter.
    fill_surface : 1-D array of length nnodes
        The surface to fill in LL node order. Modified in place.
    all_neighbors : (nnodes, max_nneighbours) array
        Adjacent nodes at each node.
    pitq : heap queue (i.e., a structured list)
        Current nodes known to be in a lake, if already identified.
    openq : StablePriorityQueue object
        Ordered queue of nodes remaining to be checked out by the algorithm
        that are known not to be in a lake.
    closedq : 1-D boolean array of length nnodes
        Nodes already or not to be explored by the algorithm.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
    >>> mg = RasterModelGrid((5, 6), 1.)
    >>> for edge in ('left', 'top', 'bottom'):
    ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
    >>> z = mg.zeros('node', dtype=float)
    >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
    >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
    >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
    >>> zw = z.copy()
    >>> openq = StablePriorityQueue()
    >>> pitq = []
    >>> closedq = mg.zeros('node', dtype=bool)
    >>> closedq[mg.status_at_node == CLOSED_BOUNDARY] = True
    >>> edges = np.array([11, 17, 23])
    >>> for edgenode in edges:
    ...     openq.add_task(edgenode, priority=z[edgenode])
    >>> closedq[edges] = True
    >>> while True:
    ...     try:
    ...         _fill_one_node_to_flat(mg, zw, mg.adjacent_nodes_at_node,
    ...                                pitq, openq, closedq)
    ...     except KeyError:
    ...         break
    ...     else:
    ...         print(np.sort(openq.nodes_currently_in_queue()), pitq) ### REMOVE POST DEBUG

    Now check the values make sense.

    >>> lake = np.array([False, False, False, False, False, False,
    ...                  False, False,  True,  True, False, False,
    ...                  False, False,  True,  True, False, False,
    ...                  False, False,  True,  True, False, False,
    ...                  False, False, False, False, False, False])
    >>> np.allclose(zw[lake], z[16])
    True
    >>> np.all(np.greater(zw[lake], z[lake]))
    True
    >>> np.allclose(zw[np.logical_not(lake)], z[np.logical_not(lake)])
    True
    """
    try:
        c = heapq.heappop(pitq)
    except IndexError:
        c = openq.pop_task()
        # this will raise a KeyError once it's exhausted both queues
    cneighbors = all_neighbors[c]
    openneighbors = cneighbors[
        np.logical_not(closedq[cneighbors])]  # for efficiency
    closedq[openneighbors] = True
    for n in openneighbors:
        if fill_surface[n] <= fill_surface[c]:
            fill_surface[n] = fill_surface[c]
            heapq.heappush(pitq, n)
        else:
            openq.add_task(n, priority=fill_surface[n])


class LakeFillerBarnes(Component):
    """
    Uses the Barnes et al (2014) algorithms to replace pits with flats, or
    optionally to very shallow gradient surfaces to allow continued draining.

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
    route_flow_steepest_descent : bool
        If True, the component outputs the 'flow__receiver_node' and
        'flow__link_to_receiver_node' fields.
    calc_slopes : bool
        For direct comparison with the flow_director components, if True, the
        component outputs the 'topographic__steepest_slope' field.
    track_filled_nodes : bool
        If True, permits use of a suite of component properties to query
        which nodes were filled, and by how much. Set to False to enhance
        component speed.
    ignore_overfill : bool
        If True, suppresses the Error that would normally be raised during
        creation of a gentle incline on a fill surface (i.e., if
        fill_flat=False). Typically this would happen on a synthetic DEM
        where more than one outlet is possible at the same elevation.
    """
    def __init__(self, grid, surface='topographic__elevation',
                 method='d8', fill_flat=True,
                 route_flow_steepest_descent=True,
                 calc_slopes=True, track_filled_nodes=True,
                 ignore_overfill=False):
        """
        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> # from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        """
        self._grid = grid
        self._open = StablePriorityQueue()
        self._pit = []
        self._closed = self.grid.zeros('node', dtype=bool)
        self._gridclosednodes = self.grid.status_at_node == CLOSED_BOUNDARY
        gridopennodes = np.logical_not(self._gridclosednodes)
        # close up the CLOSED_BOUNDARY permanently:
        self._closed[self._gridclosednodes] = True

        # this component maintains its own internal count of how often it has
        # been called. This is to enable "cheap" data access of the various
        # available data structures without needless recalculation
        self._runcounter = itertools.count()
        self._runcount = -1  # not yet run
        self._lastcountforlakemap = -1  # lake_map has not yet been called
        self._trackingcounter = itertools.count()  # sequ counter if tracking
        self._trackercount = -1  # tracking calls not yet run
        self._lastcountfortracker = -1
        self._ignore_overfill = ignore_overfill

        # get the neighbour call set up:
        assert method in {'steepest', 'd8'}
        if method == 'd8':
            try:
                self._allneighbors = self.grid.d8s_at_node
            except AttributeError:
                self._allneighbors = self.grid.adjacent_nodes_at_node
        else:
            self._allneighbors = self.grid.adjacent_nodes_at_node

        # A key difference from the "pure" Barnes algorithm for LL is that
        # we must'n flood from all the edges. Instead, we can only flood from
        # a true grid edge, i.e., only the FIXED boundary types. (Both
        # CLOSED and LOOPED assume flow either can't get out there, or at
        # least, there's more land in that direction that will be handled
        # otherwise.) Note we'l add a test that there is at least some kind
        # of outlet!!
        self._edges = np.where(np.logical_or(
            self.grid.status_at_node == FIXED_VALUE_BOUNDARY,
            self.grid.status_at_node == FIXED_GRADIENT_BOUNDARY))[0]
        if self._edges.size == 0:
            raise ValueError(
                'No valid outlets found on the grid! You cannot run the ' +
                'filling algorithms!')
        # and finally, close these up permanently as well (edges will always
        # be edges...)
        self._closed[self._edges] = True

        # now, work out what constitutes a "surface" under various input opts:
        self._dontreroute = not route_flow_steepest_descent
# NOTE: check the overlap with the rerouting method here.
# do we have to still implicitly do the filling??
        # check if we are modifying in place or not:
        self._inplace = surface is fill_surface
        # then
        self._surface = return_array_at_node(grid, surface)
        self._fill_surface = return_array_at_node(grid, surface)
        if fill_flat:
            self._fill_one_node = self._fill_one_node_to_flat
        else:
            self._fill_one_node = self._fill_one_node_to_slant
        if track_filled_nodes:
            self._tracking_fill = True
            self._orig_surface = self._surface.copy()
            self._lakefiller = LakeMapperBarnes(
                grid, surface=self._orig_surface, method=method, fill=True,
                fill_flat=fill_flat, fill_surface=surface,
                route_flow_steepest_descent=route_flow_steepest_descent,
                calc_slopes=calc_slopes, ignore_overfill=ignore_overfill)
        else:
            self._tracking_fill = False

    def run_one_step(self):
        "Fills the surface to remove all pits."
        if self._tracking_fill:
            self._lakefiller.run_one_step()
        else:
            # increment the run count
            self._runcount = next(self._runcounter)
            # First get _fill_surface in order.
            self._fill_surface[:] = self._surface  # surfaces begin identical
            # note this is nice & efficent if _fill_surface is _surface
            # now, return _closed to its initial cond, w only the
            # CLOSED_BOUNDARY and grid draining nodes pre-closed:
            closedq = self._closed.copy()
            if self._dontreroute:  # i.e. algos 2 & 3
                for edgenode in self._edges:
                    self._open.add_task(
                        edgenode, priority=self._surface[edgenode])
                self._closed[self._edges] = True
                while True:
                    try:
                        self._fill_one_node(self.grid, )
                    except KeyError:  # run out of nodes to fill...
                        break

    @property
    def fill_dict(self):
        """
        Return a dictionary where the keys are the outlet nodes of each filled
        area, and the values are deques of nodes within each. Items are not
        returned in ID order.
        """
        assert self._tracking_fill
        return self._lakefiller.lake_dict

# NOTE: need additional counter on the tracking mechanism to ensure that's getting run, else _lakemappings could be wrong or out-of-date

    @property
    def fill_outlets(self):
        """
        Returns the outlet for each filled area, not necessarily in ID order.
        """
        assert self._tracking_fill
        return self._lakefiller.lake_outlets

    @property
    def number_of_fills(self):
        """
        Return the number of individual filled areas.
        """
        assert self._tracking_fill
        return self._lakefiller.number_of_lakes

    @property
    def fill_map(self):
        """
        Return an array of ints, where each filled node is labelled
        with its outlet node ID.
        Nodes not in a filled area are labelled with LOCAL_BAD_INDEX_VALUE
        (default -1).
        """
        assert self._tracking_fill
        return self._lakefiller.lake_map

    @property
    def fill_at_node(self):
        """
        Return a boolean array, True if the node is filled, False otherwise.
        """
        assert self._tracking_fill
        return self._lakefiller.lake_at_node

    @property
    def fill_depth(self):
        """Return the change in surface elevation at each node this step.
        """
        assert self._tracking_fill
        return self._lakefiller.fill_depth

    @property
    def fill_areas(self):
        """
        A nlakes-long array of the area of each fill. The order is the same as
        that of the keys in fill_dict, and of fill_outlets.
        """
        assert self._tracking_fill
        return self._lakefiller.lake_areas

    @property
    def fill_volumes(self):
        """
        A nlakes-long array of the volume of each fill. The order is the same
        as that of the keys in fill_dict, and of fill_outlets.
        """
        assert self._tracking_fill
        return self._lakefiller.lake_volumes


class LakeEvaporator(Component):
    """
    This component will permit the lowering of a lake surface under
    evaporation, following the "slab removal" method used in the
    LakeFlooderBarnes. This is needed in case evaporation can change the
    surface area significantly.
    It can also be used by the LakeFlooderBarnes to permit the "special case"
    there.
    """


class LakeFlooderBarnes(LakeMapperBarnes):
    """
    This algorithm honours the actual discharges supplying water to
    depressions, and so honours water mass balance. It requires both a surface
    and a water_surface to already exist, and flow routing to have already
    occured over the surface (i.e., ignoring that water_surface).
    """
    # mechanism: expand outwards from the outflow, as in the normal LMB. Once
    # the expansion meets a pit, ingest it.
    # Continue to sum ingested pits as we go, and sum how many.
    # Once single lake is filled, check we had enough water to do this.
    # If not, flag as problematic.

    # If we're OK, will still need to reduce outflow, but do this after?

    # Problematic lakes need 2nd pass.

    # Easier if only one pit in lake (special case). Take list of all elevs
    # in the lake & sort.
    # (If Voronoi, also need to track area at each.)
    # This represents increments we can drop the surface by w/o changing
    # area yet. So, given we know the water shortfall, take of water slabs til
    # we get to the right vol.

    # Else, we must iterate. Easiest to work outwards from the pits, per the
    # approach in the old LakeMapper. i.e., take the slabs, but work upwards.
