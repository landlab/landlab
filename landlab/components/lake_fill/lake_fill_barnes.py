#!/usr/env/python

"""
lake_fill_barnes.py

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

LOCAL_BAD_INDEX_VALUE = BAD_INDEX_VALUE
LARGE_ELEV = 9999999999.

# TODO: Needs to have rerouting functionality...


class StablePriorityQueue():
    """
    Implements a stable priority queue, that tracks insertion order; i.e., this
    is used to break ties.

    See https://docs.python.org/2/library/heapq.html#priority-queue-implementation-notes
    & https://www.sciencedirect.com/science/article/pii/S0098300413001337
    """
    def __init__(self):
        self._pq = []                          # list of entries as a heap
        self._entry_finder = {}                # mapping of tasks to entries
        self._REMOVED = BAD_INDEX_VALUE        # placeholder for a removed task
        self._counter = itertools.count()      # unique sequence count
        self._nodes_ever_in_queue = deque([])
        # last one tracks all nodes that have ever been added

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


def _fill_one_node_to_flat(fill_surface, all_neighbors,
                           pitq, openq, closedq, dummy):
    """
    Implements the Barnes et al. algorithm for a simple fill. Assumes the
    _open and _closed lists have already been updated per Barnes algos 2&3,
    lns 1-7.

    Parameters
    ----------
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
    dummy : any Python object
        Necessary for direct comparison with _fill_one_node_to_slant.

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
    ...         _fill_one_node_to_flat(zw, mg.adjacent_nodes_at_node,
    ...                                pitq, openq, closedq, None)
    ...     except KeyError:
    ...         break

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
    fill_flat : bool
        If True, pits will be filled to perfectly horizontal. If False, the new
        surface will be slightly inclined to give steepest descent flow paths
        to the outlet.
    fill_surface : bool
        Sets the field or array to fill. If fill_surface is surface, this
        operation occurs in place, and is faster.
        Note that the component will overwrite fill_surface if it exists; to
        supply an existing water level to it, supply that water level field as
        surface, not fill_surface.
    route_flow_steepest_descent : bool
        If True, the component outputs the 'flow__receiver_node' and
        'flow__link_to_receiver_node' fields.
    calc_slopes : bool
        For direct comparison with the flow_director components, if True, the
        component outputs the 'topographic__steepest_slope' field.
    ignore_overfill : bool
        If True, suppresses the Error that would normally be raised during
        creation of a gentle incline on a fill surface (i.e., if not
        fill_flat). Typically this would happen on a synthetic DEM where more
        than one outlet is possible at the same elevation. If True, the
        was_there_overfill property can still be used to see if this has
        occurred.
    track_lakes : bool
        If True, the component permits a slight hit to performance in order to
        explicitly track which nodes have been filled, and to enable queries
        on that data in retrospect. Set to False to simply fill the surface
        and be done with it.
    """
    def __init__(self, grid, surface='topographic__elevation',
                 method='d8', fill_flat=True,
                 fill_surface='topographic__elevation',
                 route_flow_steepest_descent=True,
                 calc_slopes=True, ignore_overfill=False, track_lakes=True):
        """
        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
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
        self._PitTop = LARGE_ELEV  # variable to not overfill slanted surfaces
        self._ignore_overfill = ignore_overfill
        self._overfill_flag = False
        self._track_lakes = track_lakes

        # get the neighbour call set up:
        assert method in {'steepest', 'd8', 'D8'}
        if method in ('d8', 'D8'):
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

        # NOTE: this component can't yet do this rerouting, so...
        if route_flow_steepest_descent:
            raise ValueError, "Component can't yet do rerouting, sorry..."
        
        # check if we are modifying in place or not. This gets used to check
        # it makes sense to calculate various properties.
        self._inplace = surface is fill_surface
        # then
        self._surface = return_array_at_node(grid, surface)
        self._fill_surface = return_array_at_node(grid, fill_surface)
        
        # NOTE: buggy functionality of return_array_at_node here means
        # component can't yet handle arrays as opposed to fields...
        # This will be resolved by a modification to return_array_at_node


        self._fill_flat = fill_flat
        if fill_flat:
            self._fill_one_node = _fill_one_node_to_flat
        else:
            self._fill_one_node = self._fill_one_node_to_slant

    def _fill_one_node_to_slant(self, fill_surface, all_neighbors,
                                pitq, openq, closedq, ignore_overfill):
        """
        Implements the Barnes et al. algorithm to obtain a naturally draining
        surface. Assumes the _open and _closed lists have already been updated
        per Barnes algos 2&3, lns 1-7.

        Parameters
        ----------
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
        ignore_overfill : bool
            If False, method will raise a ValueError if adding an increment
            to the node's elevation would fundamentally alter the resulting
            drainage pattern (e.g., it would create a new outlet somewhere).
            If True, the elevation of the node will be changed regardless.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([11, 17, 23])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> first_nodes_checked = []

        >>> for i in range(3):  # run a couple of steps
        ...     lmb._fill_one_node_to_slant(
        ...         z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...         lmb._closed, False)
        ...     print(lmb._open.peek_at_task())
        ...     assert lmb._pit == []  # these steps don't find pits
        17
        23
        16

        >>> lmb._fill_one_node_to_slant(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False)
        >>> lmb._pit == [15, ]  # Popping 16 off "open" puts 15 in "pit"
        True
        >>> np.isclose(z[15], z[16])  # 15 has also been filled in this step
        True
        >>> z[15] > z[16]  # ...but 15 is incrementally greater than 16
        True

        >>> lmb._fill_one_node_to_slant(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False)
        >>> lmb._pit == [9, 21, 14]  # 15 pops of pit, these neighbors go in
        True
        >>> np.allclose(z[15], [z[9], z[21], z[14]])  # now filled
        True
        >>> np.all([z[9] == z[21], z[21] == z[14]])  # these perfectly level
        True
        >>> z[9] > z[15]  # ...but incrementally above 15
        True

        >>> lmb._fill_one_node_to_slant(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False)
        >>> lmb._pit == [8, 21, 14]  # 9 popped off pit, 8 went in. And so on.

        Test a failing example. This behaviour exists to prevent the
        application of the gradient from fundamentally altering the drainage
        pattern that "should" result.

        >>> mg = RasterModelGrid((3, 7), 1.)
        >>> for edge in ('top', 'right', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1,
        ...                                 1.0000000000000004, 1.5]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([7, ])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> while True:
        ...     try:
        ...         lmb._fill_one_node_to_slant(
        ...             z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...             lmb._closed, False)
        ...     except KeyError:
        ...         break
        ValueError: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.
        """
        try:
            topopen = openq.peek_at_task()
        except KeyError:
            noopen = True
        else:
            noopen = False
        try:
            toppit = pitq[0]
        except IndexError:
            nopit = True
        else:
            nopit = False
        if not (nopit or noopen):
            if topopen == toppit:  # intentionally tight comparison
                c = openq.pop_task()
                self._PitTop = LARGE_ELEV
        if not nopit:
            c = heapq.heappop(pitq)
            if np.isclose(self._PitTop, LARGE_ELEV):
                self._PitTop = fill_surface[c]
        else:
            c = openq.pop_task()  # again, returns KeyError if empty
            self._PitTop = LARGE_ELEV

        for n in all_neighbors[c]:
            if closedq[n]:
                continue
            else:
                closedq[n] = True
            nextval = np.nextafter(
                fill_surface[c], LARGE_ELEV)
            if self._gridclosednodes[n]:
                heapq.heappush(pitq, n)
            elif fill_surface[n] <= nextval:
                if (self._PitTop < fill_surface[n] and
                        nextval >= fill_surface[n]):
                    if ignore_overfill:
                        self._overfill_flag = True
                    else:
                        raise ValueError(
                            "Pit is overfilled due to creation of two " +
                            "outlets as the minimum gradient gets applied. " +
                            "Suppress this Error with the ignore_overfill " +
                            "flag at component instantiation."
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

        This version runs a little more slowly to enable tracking of which
        nodes are linked to which outlets.

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

        Examples
        --------
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([11, 17, 23])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> lmb._fill_to_flat_with_tracking(z, mg.adjacent_nodes_at_node,
        ...                                 lmb._pit, lmb._open, lmb._closed)
        {8: deque([7]), 16: deque([15, 9, 14, 22])}
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
                    lakemappings[outlet_ID] = deque([c, ])

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
            # print(np.sort(openq.nodes_currently_in_queue()), pitq)
        return lakemappings

    def _fill_to_slant_with_optional_tracking(self, fill_surface,
                                              all_neighbors, pitq, openq,
                                              closedq, ignore_overfill,
                                              track_lakes):
        """
        Parameters
        ----------
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
        ignore_overfill : bool
            If False, method will raise a ValueError if adding an increment
            to the node's elevation would fundamentally alter the resulting
            drainage pattern (e.g., it would create a new outlet somewhere).
            If True, the elevation of the node will be changed regardless.
        track_lakes : bool
            If True, returns a dict with data on the lakes created. If false,
            returns an empty dict.

        Returns
        -------
        lakemappings : dict
            If track_lakes, {outlet_ID : [nodes draining to outlet]}. This is
            a dict with outlet nodes of individual lakes as keys, and lists
            (strictly, deques) of each node inundated (i.e., depth > 0.) by
            that lake. Note len(keys) is the number of individually mapped
            lakes.
            If not track_lakes, an empty dict.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([11, 17, 23])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> lmb._fill_to_slant_with_optional_tracking(
        >>>         z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        >>>         lmb._closed, False, True)
        {16: deque([15, 9, 8, 14, 20, 21])}

        Test two pits:

        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([11, 17, 23])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> lmb._fill_to_slant_with_optional_tracking(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False, True)
        {8: deque([7]), 16: deque([15, 9, 14, 22])}
        >>> fr = FlowRouter(mg, method='D4')
        >>> fr.route_flow()
        >>> np.all(mg.at_node['flow__sink_flag'][mg.core_nodes] == 0)
        True
        >>> drainage_area = np.array([  0.,   0.,   0.,   0.,   0.,   0.,
        ...                             0.,   1.,   2.,   3.,   1.,   1.,
        ...                             0.,   1.,   4.,   9.,  10.,  10.,
        ...                             0.,   1.,   2.,   1.,   1.,   1.,
        ...                             0.,   0.,   0.,   0.,   0.,   0.])
        >>> np.allclose(mg.at_node['drainage_area'], drainage_area)
        True

        With track_lakes == False, fill still works just fine, but the dict
        returned is empty:

        >>> z[:] = mg.node_x.max() - mg.node_x  # all this as above
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([11, 17, 23])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> lmb._fill_to_slant_with_optional_tracking(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False, False)  # empty dict now
        {}

        >>> fr.route_flow()  # drains fine still, as above
        >>> np.allclose(mg.at_node['drainage_area'], drainage_area)
        True

        Test a failing example:

        >>> mg = RasterModelGrid((3, 7), 1.)
        >>> for edge in ('top', 'right', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1,
        ...                                 1.0000000000000004, 1.5]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([7, ])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> lmb._fill_to_slant_with_optional_tracking(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False, True)
        ValueError: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.
        """
        lakemappings = dict()
        outlet_ID = BAD_INDEX_VALUE
        while True:
            try:
                topopen = openq.peek_at_task()
            except KeyError:
                noopen = True
                topopen = None
            else:
                noopen = False
            try:
                toppit = pitq[0]
            except IndexError:
                nopit = True
                toppit = None
            else:
                nopit = False
            if (not (nopit or noopen)) and (topopen == toppit):
                # intentionally tight comparison
                c = openq.pop_task()
                outlet_ID = c
                self._PitTop = LARGE_ELEV
            elif not nopit:
                c = heapq.heappop(pitq)
                if np.isclose(self._PitTop, LARGE_ELEV):
                    self._PitTop = fill_surface[c]
                if track_lakes:
                    try:
                        lakemappings[outlet_ID].append(c)
                        # ^add this node to lake
                    except KeyError:
                        # ^this is the first node of a new lake
                        lakemappings[outlet_ID] = deque([c, ])
            else:
                try:
                    c = openq.pop_task()
                    # ^again, returns KeyError if empty
                except KeyError:
                    break
                outlet_ID = c
                self._PitTop = LARGE_ELEV

            for n in all_neighbors[c]:
                if closedq[n]:
                    continue
                else:
                    closedq[n] = True
                nextval = np.nextafter(
                    fill_surface[c], LARGE_ELEV)
                if self._gridclosednodes[n]:
                    heapq.heappush(pitq, n)
                elif fill_surface[n] <= nextval:
                    if (self._PitTop < fill_surface[n] and
                            nextval >= fill_surface[n]):
                        if ignore_overfill:
                            self._overfill_flag = True
                        else:
                            raise ValueError(
                                "Pit is overfilled due to creation of two " +
                                "outlets as the minimum gradient gets " +
                                "applied. Suppress this Error with the " +
                                "ignore_overfill flag at component " +
                                "instantiation."
                                )
                    fill_surface[n] = nextval
                    heapq.heappush(pitq, n)
                else:
                    openq.add_task(n, priority=fill_surface[n])
        return lakemappings

    def run_one_step(self):
        """
        Fills the surface to fill all pits. Note that a second run on a
        surface that has already been filled will *not* "see" any existing
        lakes correctly - it will see lakes, but with zero depths. In
        particular, if fill_flat is False, an attempt to fill a
        surface a second time will raise a ValueError unless ignore_overfill.
        (In this case, setting ignore_overfill is True will give the expected
        behaviour.)

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        # >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        # ...                        fill_surface=z, fill_flat=False,
        # ...                        route_flow_steepest_descent=False,
        # ...                        calc_slopes=False, track_lakes=False)
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=False)
        >>> lmb.run_one_step()
        >>> z_out = np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        ...                    0. ,  2.1,  1.5,  1.5,  1.6,  0. ,
        ...                    0. ,  2. ,  1.5,  1.5,  1.5,  0. ,
        ...                    0. ,  2.2,  1.5,  1.5,  1.7,  0. ,
        ...                    0. ,  0. ,  0. ,  0. ,  0. ,  0. ])
        >>> np.allclose(z, z_out)
        True
        >>> np.all(np.equal(z, z_out))  # those 1.5's are actually a bit > 1.5
        False
        >>> lmb.lake_map  # not created, as we aren't tracking
        AssertionError: Enable tracking to access information about lakes
        >>> lmb.was_there_overfill  # everything fine with slope adding
        False

        >>> fr = FlowRouter(mg, method='D4')  # routing will work fine now
        >>> fr.route_flow()
        >>> np.all(mg.at_node['flow__sink_flag'][mg.core_nodes] == 0)
        True
        >>> drainage_area = np.array([  0.,   0.,   0.,   0.,   0.,   0.,
        ...                             0.,   1.,   2.,   3.,   1.,   1.,
        ...                             0.,   1.,   4.,   9.,  10.,  10.,
        ...                             0.,   1.,   2.,   1.,   1.,   1.,
        ...                             0.,   0.,   0.,   0.,   0.,   0.])
        >>> np.allclose(mg.at_node['drainage_area'], drainage_area)
        True

        Test two pits:

        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=True,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict
        {8: deque([7]), 16: deque([15, 9, 14, 22])}
        >>> lmb.number_of_lakes
        2
        >>> lmb.lake_depths  # z was both surface and 'fill_surface'
        AssertionError: surface and fill_surface must be different fields or arrays to enable the property fill_depth!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        fill_flat=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict
        {8: deque([7]), 16: deque([15, 9, 14, 22])}
        >>> np.allclose(lmb.lake_areas, np.array([ 4.,  1.]))
        True
        >>> lmb.run_one_step()  # z already filled, so...
        ValueError: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.

        Suppress this behaviour with ignore_overfill:

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        fill_flat=False, track_lakes=True,
        ...                        ignore_overfill=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict
        {8: deque([7]), 16: deque([15, 9, 14, 22])}
        >>> lmb.run_one_step()
        >>> np.allclose(lmb.lake_areas, np.array([ 4.,  1.]))  # found them!
        True
        """
        # do the prep:
        # increment the run counter
        self._runcount = next(self._runcounter)
        # First get _fill_surface in order.
        self._fill_surface[:] = self._surface  # surfaces begin identical
        # note this is nice & efficent if _fill_surface is _surface
        # now, return _closed to its initial cond, w only the CLOSED_BOUNDARY
        # and grid draining nodes pre-closed:
        closedq = self._closed.copy()
# NOTE: several of the below were formerly self._closed...????
# Are these updates even necessary???
        if self._track_lakes:
            for edgenode in self._edges:
                self._open.add_task(edgenode, priority=self._surface[edgenode])
            closedq[self._edges] = True
            if self._fill_flat:
                self._lakemappings = self._fill_to_flat_with_tracking(
                    self._fill_surface, self._allneighbors, self._pit,
                    self._open, closedq)
            else:
                self._lakemappings = (
                    self._fill_to_slant_with_optional_tracking(
                        self._fill_surface, self._allneighbors, self._pit,
                        self._open, closedq,
                        ignore_overfill=self._ignore_overfill,
                        track_lakes=True))

        else:  # not tracked
            if self._dontreroute:  # i.e. algos 2 & 3
                for edgenode in self._edges:
                    self._open.add_task(edgenode,
                                        priority=self._surface[edgenode])
                closedq[self._edges] = True
                while True:
                    try:
                        self._fill_one_node(
                            self._surface, self._allneighbors, self._pit,
                            self._open, closedq, self._ignore_overfill)
                    except KeyError:  # run out of nodes to fill...
                        break

    @property
    def lake_dict(self):
        """
        Return a dictionary where the keys are the outlet nodes of each lake,
        and the values are deques of nodes within each lake. Items are not
        returned in ID order. The outlet nodes are NOT part of the lake.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=False)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict
        AssertionError: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict
        {16: deque([15, 9, 8, 14, 20, 21])}
        """
        assert self._track_lakes, \
            "Enable tracking to access information about lakes"
        return self._lakemappings

# NOTE: need additional counter on the tracking mechanism to ensure that's getting run, else _lakemappings could be wrong or out-of-date

    @property
    def lake_outlets(self):
        """
        Returns the outlet for each lake, not necessarily in ID order.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=False)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets
        AssertionError: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets == [16, ]
        True
        """
        assert self._track_lakes, \
            "Enable tracking to access information about lakes"
        return self._lakemappings.keys()

    @property
    def number_of_lakes(self):
        """
        Return the number of individual lakes. Lakes sharing outlet nodes are
        considered part of the same lake.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=False)
        >>> lmb.run_one_step()
        >>> lmb.number_of_lakes
        AssertionError: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.number_of_lakes
        2
        """
        assert self._track_lakes, \
            "Enable tracking to access information about lakes"
        return len(self._lakemappings)

    @property
    def lake_map(self):
        """
        Return an array of ints, where each node within a lake is labelled
        with its outlet node ID. The outlet nodes are NOT part of the lakes.
        Nodes not in a lake are labelled with LOCAL_BAD_INDEX_VALUE
        (default -1).

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=False)
        >>> lmb.run_one_step()
        >>> lmb.lake_map
        AssertionError: Enable tracking to access information about lakes

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lake_map = np.array([-1, -1, -1, -1, -1, -1,
        ...                      -1,  8, -1, 16, -1, -1,
        ...                      -1, -1, 16, 16, -1, -1,
        ...                      -1, -1, -1, -1, 16, -1,
        ...                      -1, -1, -1, -1, -1, -1])
        >>> np.all(np.equal(lmb.lake_map, lake_map))
        True

        Note that the outlet node is NOT part of the lake.

        Updating the elevations works fine:

        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> lmb.run_one_step()
        >>> new_lake_map = np.array([-1, -1, -1, -1, -1, -1,
        ...                          -1, -1, 16, 16, -1, -1,
        ...                          -1, -1, 16, 16, -1, -1,
        ...                          -1, -1, 16, 16, -1, -1,
        ...                          -1, -1, -1, -1, -1, -1])
        >>> np.all(np.equal(lmb.lake_map, new_lake_map))
        True
        """
        if self._runcount > self._lastcountforlakemap:
            # things have changed since last call to lake_map
            self._lake_map = np.full(mg.number_of_nodes, LOCAL_BAD_INDEX_VALUE,
                                     dtype=int)
            for (outlet, lakenodes) in self.lake_dict.iteritems():
                self._lake_map[lakenodes] = outlet
        else:
            pass  # old map is fine
        self._lastcountforlakemap = self._runcount
        return self._lake_map

    @property
    def lake_at_node(self):
        """
        Return a boolean array, True if the node is flooded, False otherwise.
        The outlet nodes are NOT part of the lakes.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lake_at_node = np.array([False, False, False, False, False, False,
        ...                          False,  True, False,  True, False, False,
        ...                          False, False,  True,  True, False, False,
        ...                          False, False, False, False,  True, False,
        ...                          False, False, False, False, False, False],
        ...                          dtype=bool)
        >>> np.all(np.equal(lmb.lake_at_node, lake_at_node))
        True
        """
        return self.lake_map != LOCAL_BAD_INDEX_VALUE

    @property
    def lake_depths(self):
        """Return the change in surface elevation at each node this step.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_depths  # won't work as surface & fill_surface are both z
        AssertionError: surface and fill_surface must be different fields or arrays to enable the property lake_depths!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
                                   surface=z_init,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lake_depths = np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        ...                          0. ,  1. ,  0. ,  0.5,  0. ,  0. ,
        ...                          0. ,  0. ,  0.4,  0.7,  0. ,  0. ,
        ...                          0. ,  0. ,  0. ,  0. ,  0.1,  0. ,
        ...                          0. ,  0. ,  0. ,  0. ,  0. ,  0. ])
        >>> np.all(np.equal(lmb.lake_depths,
        ...                 lake_depths))  # slope applied, so...
        False
        >>> np.allclose(lmb.lke_depths, lake_depths)
        True
        """
        assert not self._inplace, \
            "surface and fill_surface must be different fields or arrays " + \
            "to enable the property lake_depths!"
        return self._fill_surface - self._surface

    @property
    def lake_areas(self):
        """
        A nlakes-long array of the area of each lake. The order is the same as
        that of the keys in lake_dict, and of lake_outlets. Note that outlet
        nodes are not parts of the lakes.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=False)
        >>> lmb.run_one_step()
        >>> lmb.lake_areas  # note track_lakes=False
        AssertionError: Enable tracking to access information about lakes

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets
        [16, 8]
        >>> np.allclose(lmb.lake_areas, np.array([ 4.,  1.]))
        True
        """
        lakeareas = np.empty(self.number_of_lakes, dtype=float)
        for (i, (outlet, lakenodes)) in enumerate(
             self.lake_dict.iteritems()):
            lakeareas[i] = self.grid.cell_area_at_node[lakenodes].sum()
        return lakeareas

    @property
    def lake_volumes(self):
        """
        A nlakes-long array of the volume of each lake. The order is the same
        as that of the keys in lake_dict, and of lake_outlets.
        Note that this calculation is performed relative to the initial
        surface, so is only a true lake volume if the initial surface was the
        rock suface (not an earlier water level).

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowRouter
        >>> mg = RasterModelGrid((5, 6), 1.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_volumes  # won't work as surface & fill_surface are both z
        AssertionError: surface and fill_surface must be different fields or arrays to enable the property lake_volumes!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
                                   surface=z_init,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets
        [16, 8]
        >>> np.allclose(lmb.lake_volumes, np.array([ 1.7,  1. ]))
        True
        """
        lake_vols = np.empty(self.number_of_lakes, dtype=float)
        col_vols = self.grid.cell_area_at_node * self.lake_depths
        for (i, (outlet, lakenodes)) in enumerate(
             self.lake_dict.iteritems()):
            lake_vols[i] = col_vols[lakenodes].sum()
        return lake_vols

    @property
    def was_there_overfill(self):
        """
        If the ignore_overfill flag was set to True at instantiation, this
        property indicates if any depression in the grid has, at any point,
        been overfilled.

        Examples
        --------
        >>> mg = RasterModelGrid((3, 7), 1.)
        >>> for edge in ('top', 'right', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1,
        ...                                 1.0000000000000004, 1.5]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=True,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, ignore_overfill=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.was_there_overfill
        AssertionError: was_there_overfill is only defined if filling to an inclined surface!

        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, ignore_overfill=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        ValueError: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.

        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, ignore_overfill=True,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.was_there_overfill
        True

        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1, 1., 1.5]
        >>> lmb.run_one_step()
        >>> lmb.was_there_overfill  # still true as was in the previous example
        True

        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1, 1., 1.5]
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        route_flow_steepest_descent=False,
        ...                        calc_slopes=False, ignore_overfill=True,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.was_there_overfill  # Now reset
        False

        >>> lmb.run_one_step()  # 2nd run on same fill_surface creates overfill
        >>> lmb.was_there_overfill
        True

        Note however that in this last example, values have NOT changed.
        """
        assert self._fill_flat is False, \
            "was_there_overfill is only defined if filling to an " + \
            "inclined surface!"
        return self._overfill_flag


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
