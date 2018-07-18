#!/usr/env/python

"""
fill_sinks_barnes.py

Fill sinks in a landscape to the brim, following the Barnes et al. (2014) algos.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.return_array import return_array_at_node
from landlab.core.messages import warning_message

from landlab import BAD_INDEX_VALUE
from collections import deque
import six
import numpy as np
import heapq
import itertools
import deque


class StablePriorityQueue():
    """
    Implements a stable priority queue, that tracks insertion order; i.e., this
    is used to break ties.

    See https://docs.python.org/2/library/heapq.html#priority-queue-implementation-notes
    & https://www.sciencedirect.com/science/article/pii/S0098300413001337
    """
    _pq = []                          # list of entries arranged in a heap
    _entry_finder = {}                # mapping of tasks to entries
    _REMOVED = '<removed-task>'       # placeholder for a removed task
    _counter = itertools.count()      # unique sequence count
    _nodes_ever_in_queue = deque([])  # tracks what has ever been added

    def add_task(self, task, priority=0):
        "Add a new task or update the priority of an existing task"
        if task in self._entry_finder:
            self.remove_task(task)
        count = next(self._counter)
        entry = [priority, count, task]
        self._entry_finder[task] = entry
        heapq.heappush(self._pq, entry)
        _nodes_ever_in_queue.append(task)

    def remove_task(self, task):
        "Mark an existing task as _REMOVED.  Raise KeyError if not found."
        entry = self._entry_finder.pop(task)
        entry[-1] = self._REMOVED

    def pop_task(self):
        "Remove and return the lowest priority task. Raise KeyError if empty."
        while self._pq:
            priority, count, task = heapq.heappop(self._pq)
            if task is not self._REMOVED:
                del self._entry_finder[task]
                return task
        raise KeyError('pop from an empty priority queue')

    def peek_at_task(self):
        """
        Return the lowest priority task without removal. Raise KeyError if
        empty.
        """
        while self._pq:
            priority, count, task = self._pq[0]
            if task is not self._REMOVED:
                return task
        raise KeyError('peeked at an empty priority queue')

    def nodes_ever_in_queue(self):
        """
        Return array of all nodes ever added to this queue object. Repeats
        are permitted.
        """
        return np.array(self._nodes_ever_in_queue)


class LakeMapperBarnes(Component):
    """
    Note this will also flag our lake nodes, since we can access
    nodes_ever_in_queue for each spill point.

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to direct flow across.
    method : {'steepest', 'd8'}
        Whether or not to recognise diagonals as valid flow paths, if a raster.
        Otherwise, no effect.
    fill : bool
        If True, surface will be modified to convert pits to flats or slightly
        inclined surfaces, as determined by fill_flat.
    fill_flat : bool
        If True, pits will be filled to perfectly horizontal. If False, the new
        surface will be slightly inclined to give steepest descent flow paths
        to the outlet. Only has effect when fill == True.
    fill_surface : bool
        If fill == True, sets the field or array to fill. If fill_surface is
        surface, this operation occurs in place, and is faster.
        Note that the component will overwrite fill_surface if it exists; to
        supply an existing water level to it, supply it as surface, not
        fill_surface.
    route_flow_steepest_descent : bool
        If True, the component outputs the 'flow__receiver_node' and
        'flow__link_to_receiver_node' fields.
    calc_slopes : bool
        For direct comparison with the flow_director components, if True, the
        component outputs the 'topographic__steepest_slope' field.
    """
    def __init__(self, grid, surface='topographic__elevation',
                 method='d8', fill=False, fill_flat=True,
                 fill_surface='topographic__elevation',
                 route_flow_steepest_descent=True,
                 calc_slopes=True):
        self._grid = grid
        self._open = StablePriorityQueue()
        self._pit = heapq([])
        self._closed = self.grid.zeros('node', dtype=bool)

        # get the neighbour call set up:
        assert method in {'steepest', 'd8'}
        if method == 'd8':
            try:
                self._allneighbors = self.grid.d8s_at_node
            except AttributeError:
                self._allneighbors = self.grid.adjacent_nodes_at_node
        else:
            self._allneighbors = self.grid.adjacent_nodes_at_node
        # now a special case; iterate in once to find the "actual" (non-closed)
        # edge, and assume this does not change.
        # note that any closed "islands" in the grid will get operated on
        # every time; this just trims the edges.
        # this gives us self._edges.
        closededges = []
        self._edges = set()
        tested = set()
        for node in self.grid.boundary_nodes:
            if self.grid.status_at_node != CLOSED_BOUNDARY:
                self._edges.add(node)
            else:
                closededges.append(node)
            tested.add(node)
        queuetotest = heapq.heapify(closededges)
        for curr_edge in queuetotest:
            neighbors = set(self._allneighbors[curr_edge].flat)
            for ngb in neighbors.difference(tested):
                if self.grid.status_at_node != CLOSED_BOUNDARY:
                    self._edges.add(node)
                else:
                    queuetotest.heappush(ngb)
                tested.add(ngb)
        self._edges = np.array(self._edges)

        # now, work out what constitutes a "surface" under various input opts:
        self._dontreroute = not route_flow_steepest_descent
# NOTE: check the overlap with the rerouting method here.
# do we have to still implicitly do the filling??
        if fill:
            self._surface = return_array_at_node(grid, surface)
            self._fill_surface = return_array_at_node(grid, fill_surface)
            if fill_flat:
                self._fill_one_node = self._fill_one_node_to_flat
            else:
                self._fill_one_node = self._fill_one_node_to_slant

    def _fill_one_node_to_flat(self):
        """
        Implements the Barnes et al. algorithm for a simple fill. Assumes the
        _open and _closed lists have already been updated per Barnes algos 2&3,
        lns 1-7.
        """
        try:
            c = self._pit.heappop()
        except KeyError:
            c = self._open.pop_task()
            # this will raise a KeyError once it's exhausted both queues
        cneighbors = self._allneighbors[c]
        openneighbors = cneighbors[
            np.logical_not(self._closed[cneighbors])]  # for efficiency
        self._closed[openneighbors] = True
        for n in openneighbors:
            if self._fill_surface[n] <= self._fill_surface[c]:
                self._fill_surface[n] = self._fill_surface[c]
                self._pit.heappush(n)
            else:
                self._open.add_task(n, priority=self._fill_surface[n])

    def _fill_one_node_to_slant(self):
        """
        Implements the Barnes et al. algorithm to obtain a naturally draining
        surface. Assumes the _open and _closed lists have already been updated
        per Barnes algos 2&3, lns 1-7.
        """
        try:
            topopen = self._open.peek_at_task()
        except KeyError:
            pass_on = True
        try:
            toppit = self._pit[0]
        except IndexError:
            pass_on = True
        else:
            pass_on = not np.isclose(topopen, toppit)  # careful; does this invalidate the incrementing part??
        if not pass_on:
                c = self._open.pop_task()
                PitTop = None
        else:
            try:
                c = self._pit.heappop()
                try:
                    if PitTop is None:  ######## not clear exactly when PitTop gets first defined; check this.
                        PitTop = self._fill_surface[c]
                except NameError:
                    pass
            except KeyError:
                c = self._open.pop_task()  # again, returns KeyError if empty
                PitTop = None

        cneighbors = self._allneighbors[c]
        openneighbors = cneighbors[
            np.logical_not(self._closed[cneighbors])]  # for efficiency
        self._closed[openneighbors] = True
        for n in openneighbors:
            nextval = np.nextafter(
                self._fill_surface[c], 9999999999.)
            if self.grid.status_at_node[n] == CLOSED_BOUNDARY:
                self._pit.heappush(n)
            elif self._fill_surface[n] <= nextval:
                if (PitTop < self._fill_surface[n] and
                        nextval >= self._fill_surface[n]):
                    raise ValueError(
                        "Pit is overfilled due to very low gradients in DEM"
                    )
                self._fill_surface[n] = nextval
                self._pit.heappush(n)
            else:
                self._open.add_task(n, priority=self._fill_surface[n])

    def run_one_step(self):
        "Fills the surface to remove all pits.
        "
        # First get _fill_surface in order.
        self._fill_surface[:] = self._surface  # surfaces begin identical
        # note this is nice & efficent if _fill_surface is _surface
        if self._dontreroute:  # i.e. algos 2 & 3
            for edgenode in self._edges:
                self._open.add_task(edgenode, priority=self._surface[edgenode])
            self._closed[self._edges] = True
            while True:
                try:
                    self._fill_one_node()
                except KeyError:  # run out of nodes to fill...
                    break
        # NOTE: we need to track what goes into _pit, not _open!!


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
