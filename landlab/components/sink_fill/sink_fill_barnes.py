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

from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY
from landlab import CLOSED_BOUNDARY
from collections import deque
import six
import numpy as np
import heapq
import itertools


class StablePriorityQueue():
    """
    Implements a stable priority queue, that tracks insertion order; i.e., this
    is used to break ties.

    See https://docs.python.org/2/library/heapq.html#priority-queue-implementation-notes
    & https://www.sciencedirect.com/science/article/pii/S0098300413001337
    """
    def __init__(self):
        self._pq = []                          # list of entries arranged in a heap
        self._entry_finder = {}                # mapping of tasks to entries
        self._REMOVED = BAD_INDEX_VALUE        # placeholder for a removed task
        self._counter = itertools.count()      # unique sequence count
        self._nodes_ever_in_queue = deque([])  # tracks what has ever been added

    def add_task(self, task, priority=0):
        "Add a new task or update the priority of an existing task"
        if task in self._entry_finder:
            self.remove_task(task)
        count = next(self._counter)
        entry = [priority, count, task]
        self._entry_finder[task] = entry
        heapq.heappush(self._pq, entry)
        self._nodes_ever_in_queue.append(task)

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

    def nodes_currently_in_queue(self):
        "Return array of nodes currently in the queue."
        mynodes = [task for (priority, count, task) in self._pq]
        return np.array(mynodes)

    def nodes_ever_in_queue(self):
        """
        Return array of all nodes ever added to this queue object. Repeats
        are permitted.
        """
        return np.array(self._nodes_ever_in_queue)


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
            return n
        else:
            openq.add_task(n, priority=fill_surface[n])


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
        if fill:
            self._surface = return_array_at_node(grid, surface)
            self._fill_surface = return_array_at_node(grid, fill_surface)
            if fill_flat:
                self._fill_one_node = self._fill_one_node_to_flat
            else:
                self._fill_one_node = self._fill_one_node_to_slant

        # HACKY TWO LINES TO REMOVE POST DEBUG
        self._surface = return_array_at_node(grid, surface)
        self._fill_surface = return_array_at_node(grid, fill_surface)

    def _fill_one_node_to_slant(self, grid, fill_surface, all_neighbors,
                                pitq, openq, closedq):
        """
        Implements the Barnes et al. algorithm to obtain a naturally draining
        surface. Assumes the _open and _closed lists have already been updated
        per Barnes algos 2&3, lns 1-7.
        """
        try:
            topopen = openq.peek_at_task()
        except KeyError:
            pass_on = True
        try:
            toppit = pitq[0]
        except IndexError:
            pass_on = True
        else:
            pass_on = not np.isclose(topopen, toppit)  # careful; does this invalidate the incrementing part??
        if not pass_on:
                c = openq.pop_task()
                PitTop = None
        else:
            try:
                c = heapq.heappop(pitq)
                try:
                    if PitTop is None:  ######## not clear exactly when PitTop gets first defined; check this.
                        PitTop = fill_surface[c]
                except NameError:
                    pass
            except IndexError:
                c = openq.pop_task()  # again, returns KeyError if empty
                PitTop = None

        cneighbors = all_neighbors[c]
        openneighbors = cneighbors[
            np.logical_not(closedq[cneighbors])]  # for efficiency
        closedq[openneighbors] = True
        for n in openneighbors:
            nextval = np.nextafter(
                fill_surface[c], 9999999999.)
            if self._gridclosednodes[n]:
                heapq.heappush(pitq, n)
            elif fill_surface[n] <= nextval:
                if (PitTop < fill_surface[n] and
                        nextval >= fill_surface[n]):
                    raise ValueError(
                        "Pit is overfilled due to very low gradients in DEM"
                    )
                fill_surface[n] = nextval
                heapq.heappush(pitq, n)
            else:
                openq.add_task(n, priority=fill_surface[n])


    def _fill_to_flat_with_tracking(self, fill_surface, all_neighbors,
                                    pitq, openq, closedq):
        """
        Implements the Barnes et al. algorithm for a simple fill. Assumes the
        _open and _closed lists have already been updated per Barnes algos 2&3,
        lns 1-7.

        This version runs a little more slowly to enable tracking of which nodes
        are linked to which outlets.
        
        Parameters
        ----------
        fill_surface : 1-D array
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

        Returns
        -------
        lakemappings : {outlet_ID : [nodes draining to outlet]}
            Dict with outlet nodes of individual lakes as keys, and lists of
            each node inundated (i.e., depth > 0.) by that lake. Note
            len(keys) is the number of individually mapped lakes.
        
        """
        lakemappings = dict()
        outlet_ID = BAD_INDEX_VALUE
        while True:
            try:
                c = heapq.heappop(pitq)
            except IndexError:
                try:
                    c = openq.pop_task()
                    outlet_ID = c
                except KeyError:
                    break
            else:
                try:
                    lakemappings[outlet_ID].append(c)  # add this node to lake
                except KeyError:  # this is the first node of a new lake
                    lakeappings[outlet_ID] = deque([c, ])

            cneighbors = all_neighbors[c]
            openneighbors = cneighbors[
                np.logical_not(closedq[cneighbors])]  # for efficiency
            closedq[openneighbors] = True
            for n in openneighbors:
                if fill_surface[n] <= fill_surface[c]:
                    fill_surface[n] = fill_surface[c]
                    heapq.heappush(pitq, n)
                    return n
                else:
                    openq.add_task(n, priority=fill_surface[n])
            print(np.sort(openq.nodes_currently_in_queue()), pitq)
        return lakemappings


    def run_one_step(self):
        "Fills the surface to remove all pits."
        # First get _fill_surface in order.
        self._fill_surface[:] = self._surface  # surfaces begin identical
        # note this is nice & efficent if _fill_surface is _surface
        # now, return _closed to its initial cond, w only the CLOSED_BOUNDARY
        # and grid draining nodes pre-closed:
        closedq = self._closed.copy()
        if self._dontreroute:  # i.e. algos 2 & 3
            for edgenode in self._edges:
                self._open.add_task(edgenode, priority=self._surface[edgenode])
            self._closed[self._edges] = True
            while True:
                try:
                    self._fill_one_node(self.grid, )
                except KeyError:  # run out of nodes to fill...
                    break
    
    def run_one_step_to_flat_with_tracking(self):
        """
        Fills the surface, but also creates a property that is a dict of lake
        outlets and inundated nodes.

        Examples
        --------
        
        Elementary test:

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
        >>> lmb.run_one_step_to_flat_with_tracking()
        """
        
        # do the prep:
        # First get _fill_surface in order.
        self._fill_surface[:] = self._surface  # surfaces begin identical
        # note this is nice & efficent if _fill_surface is _surface
        # now, return _closed to its initial cond, w only the CLOSED_BOUNDARY
        # and grid draining nodes pre-closed:
        closedq = self._closed.copy()
        for edgenode in self._edges:
            self._open.add_task(edgenode, priority=self._surface[edgenode])
        self._closed[self._edges] = True
        for edgenode in self._edges:
            self._open.add_task(edgenode, priority=self._surface[edgenode])
        self._closed[self._edges] = True
        self._lakemappings = self._fill_to_flat_with_tracking(
            self._fill_surface, self._allneighbors, self._pit, self._open, 
            self._closed)


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
