#!/usr/env/python

"""
lake_fill_barnes.py

Fill sinks in a landscape to the brim, following the Barnes et al. (2014)
algorithms.
"""

from __future__ import print_function

from warnings import warn

from landlab import FieldError, Component, BAD_INDEX_VALUE
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.return_array import return_array_at_node
from landlab.core.messages import warning_message

from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY
from landlab import CLOSED_BOUNDARY, CORE_NODE
from landlab.components import FlowDirectorSteepest, FlowAccumulator
# ^ this simply in case Katy updates to add more fields, that we would also
# need to update...
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
    redirect_flow_steepest_descent : bool
        If True, the component outputs modified versions of the
        'flow__receiver_node', 'flow__link_to_receiver_node',
        'flow__sink_flag', and 'topographic__steepest_slope' fields. These
        are the fields output by the FlowDirector components, so set to
        True if you wish to pass this LakeFiller to the FlowAccumulator,
        or if you wish to work directly with the new, correct flow directions
        and slopes without rerunning these components on your new surface.
        Ensure the necessary fields already exist, and have already been
        calculated by a FlowDirector! This also means you need to instantiate
        your FlowDirector **before** you instantiate the LakeFillerBarnes.
        Note that the new topographic__steepest_slope will always be set to
        zero, even if fill_flat=False (i.e., there is actually a miniscule
        gradient on the new topography, which gets ignored).
    reaccumulate_flow : bool
        If True, and redirect_flow_steepest_descent is True, the run method
        will (re-)accumulate the flow after redirecting the flow. This means
        the 'drainage_area' and 'surface_water__discharge' fields will now
        reflect the new drainage patterns without having to manually
        reaccumulate the discharge. If True but redirect_flow_steepest_descent
        is False, raises an AssertionError.
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
                 redirect_flow_steepest_descent=False,
                 reaccumulate_flow=False,
                 ignore_overfill=False, track_lakes=True):
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
        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        redirect_flow_steepest_descent=False)
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
                self._allneighbors = np.concatenate(
                    (self.grid.adjacent_nodes_at_node,
                     self.grid.diagonal_adjacent_nodes_at_node), axis=1)
            except AttributeError:  # not a raster
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
        # ...note there's a slight of hand here, Because of the ordering of LL
        # grids, the last node will always be a boundary node, even for very
        # odd Voronois. This enables us to treat out -1s in the neighbour
        # arrays as always True. But, just in case...
        assert self._closed[-1] == True

        # check if we are modifying in place or not. This gets used to check
        # it makes sense to calculate various properties.
        self._inplace = surface is fill_surface
        # then
        self._surface = return_array_at_node(grid, surface)
        self._fill_surface = return_array_at_node(grid, fill_surface)

        # NOTE: buggy functionality of return_array_at_node here means
        # component can't yet handle arrays as opposed to fields...
        # This will be resolved by a modification to return_array_at_node

        # now, work out what constitutes a "surface" under various input opts:
        self._dontredirect = not redirect_flow_steepest_descent

        if redirect_flow_steepest_descent:
            # this routine only permitted if we store the lake info, so
            assert track_lakes, "You must track_lakes to redirect the flow!"
            # Check we have the necessary fields already existing.
            # These will raise FieldErrors if they don't.
            # This will cause a bunch of our tests to break, so users will
            # never see this.
            assert len(FlowDirectorSteepest.output_var_names) == 4
            self._receivers = self.grid.at_node['flow__receiver_node']
            assert len(self._receivers.shape) == 1, \
                'The LakeFillerBarnes does not yet work with one-to-many ' + \
                'flow directing schemes!'
            self._receiverlinks = self.grid.at_node[
                'flow__link_to_receiver_node']
            self._steepestslopes = self.grid.at_node[
                'topographic__steepest_slope']
            # if raster, do the neighbors & diagonals separate when rerouting
            # so we'll need to pull these separately:
            try:
                self._neighbor_arrays = (
                    self.grid.adjacent_nodes_at_node,
                    self.grid.diagonal_adjacent_nodes_at_node)
                self._link_arrays = (
                    self.grid.links_at_node,
                    self.grid.d8s_at_node[:, 4:])
                self._neighbor_lengths = self.grid.length_of_d8
            except AttributeError:  # this wasn't a raster
                self._neighbor_arrays = (self.grid.adjacent_nodes_at_node, )
                self._link_arrays = (self.grid.links_at_node, )
                self._neighbor_lengths = self.grid.length_of_link

        if reaccumulate_flow:
            assert redirect_flow_steepest_descent, (
                "You must also redirect_flow_steepest_descent if you want " +
                "to reaccumulate_flow!")
            self._reaccumulate = True
            self._fa = FlowAccumulator(self.grid)
        else:
            self._reaccumulate = False

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
        True

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
        ...     except ValueError:
        ...         print('ValueError was raised: here we overfilled')
        ValueError was raised: here we overfilled
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
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
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
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([11, 17, 23])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> out = lmb._fill_to_flat_with_tracking(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open, lmb._closed)
        >>> out == {8: deque([7]), 16: deque([15, 9, 14, 22])}
        True
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
        >>> out = lmb._fill_to_slant_with_optional_tracking(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False, True)
        >>> out == {16: deque([15, 9, 8, 14, 20, 21])}
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
        >>> lmb = LakeMapperBarnes(mg, method='steepest')
        >>> lmb._closed = mg.zeros('node', dtype=bool)
        >>> lmb._closed[mg.status_at_node == CLOSED_BOUNDARY] = True
        >>> edges = np.array([11, 17, 23])
        >>> for edgenode in edges:
        ...     lmb._open.add_task(edgenode, priority=z[edgenode])
        >>> lmb._closed[edges] = True
        >>> out = lmb._fill_to_slant_with_optional_tracking(
        ...     z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...     lmb._closed, False, True)
        >>> out == {8: deque([7]), 16: deque([15, 9, 14, 22])}
        True
        >>> fr = FlowRouter(mg, method='D4')
        >>> fr.run_one_step()
        >>> np.all(mg.at_node['flow__sink_flag'][mg.core_nodes] == 0)
        True
        >>> drainage_area = np.array([  0.,   0.,   0.,   0.,   0.,   0.,
        ...                             0.,   1.,   2.,   3.,   1.,   1.,
        ...                             0.,   1.,   4.,   9.,  11.,  11.,
        ...                             0.,   1.,   2.,   1.,   1.,   0.,
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

        >>> fr.run_one_step()  # drains fine still, as above
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
        >>> try:
        ...     lmb._fill_to_slant_with_optional_tracking(
        ...         z, mg.adjacent_nodes_at_node, lmb._pit, lmb._open,
        ...         lmb._closed, False, True)
        ... except ValueError:
        ...     print('ValueError was raised: Pit is overfilled due to ' +
        ...           'creation of two outlets as the minimum gradient ' +
        ...           'gets applied. Suppress this Error with the ' +
        ...           'ignore_overfill flag at component instantiation.')
        ValueError was raised: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.
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

    def _track_original_surface(self):
        """
        This helper method ensures that if flow is to be redircted, the
        _redirect_flowdirs() method can still get access to this information
        when it needs it. The idea here is that the operation is essentially
        free when surface and fill_surface were different to start with,
        which should make us faster.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6), 2.)
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z_new = mg.add_zeros('node', 'topographic__fill', dtype=float)
        >>> lmb = LakeMapperBarnes(mg, surface='topographic__elevation',
        ...                        fill_surface='topographic__fill',
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> orig_surf = lmb._track_original_surface()
        >>> z is orig_surf
        True
        >>> lmb = LakeMapperBarnes(mg, surface='topographic__elevation',
        ...                        fill_surface='topographic__elevation',
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> orig_surf = lmb._track_original_surface()
        >>> z is orig_surf
        False
        """
        if self._inplace:
            orig_surf = self._surface.copy()
        else:
            orig_surf = self._surface
        return orig_surf

    def _redirect_flowdirs(self, surface, lake_dict):
        """
        For nodes within lakes that have already been defined, modifies
        existing FlowDirector fields to account for the lake filling, viz.
        'flow__receiver_node', 'flow__link_to_receiver_node',
        'flow__sink_flag', and 'topographic__steepest_slope'.

        Note that the topographic__steepest_slope of a lake node will always
        be exactly 0., even if fill_flat is False.

        Examples
        --------
        >>> import numpy as np
        >>> from collections import deque
        >>> from landlab import RasterModelGrid
        >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
        >>> from landlab.components import FlowDirectorSteepest
        >>> from landlab.components import FlowAccumulator
        >>> mg = RasterModelGrid((5, 6), 2.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[23] = 1.3
        >>> z[15] = -2.  # this deep pit causes the outlet to first drain *in*
        >>> z[10] = 1.3  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own, if d8
        >>> z[9] = -1.
        >>> z[14] = 0.6  # [9, 14, 15] is a lake in both methods
        >>> z[16] = 1.2
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16 if d8
        >>> z_init = z.copy()
        >>> fd = FlowDirectorSteepest(mg)
        >>> fa = FlowAccumulator(mg)
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=True,
        ...                        redirect_flow_steepest_descent=True,
        ...                        track_lakes=True)

        In this test, we won't run the lmb. Instead, directly specify the
        correct answer:

        >>> lake_dict = {8: deque([7]), 16: deque([15, 9, 14, 22])}
        >>> fd.run_one_step()  # fill the director fields
        >>> fa.run_one_step()  # get a drainage_area
        >>> np.alltrue(mg.at_node['flow__sink_flag'][[7, 15, 22]])  # sinks
        True
        >>> nodes_in_lakes = np.array([7, 8, 9, 14, 15, 16, 22])
        >>> nodes_not_in_lakes = np.setdiff1d(mg.nodes.flat, nodes_in_lakes)

        Note we're here defining the outlets as inside the lakes, which isn't
        actually the behaviour of the component, but helps us demonstrate
        what changes, below.

        Now save the info we already have on the Flow fields:

        >>> receivers_init = mg.at_node['flow__receiver_node'].copy()
        >>> rec_links_init = mg.at_node['flow__link_to_receiver_node'].copy()
        >>> steepest_init = mg.at_node['topographic__steepest_slope'].copy()
        >>> drainage_area = mg.at_node['drainage_area'].copy()
        >>> orig_surf = lmb._track_original_surface()

        Note flow doesn't make it to the outlets:

        >>> outlets = np.where(mg.status_at_node == FIXED_VALUE_BOUNDARY)
        >>> drainage_area[outlets].sum() == mg.cell_area_at_node[
        ...     mg.core_nodes].sum()
        False

        Now, run the method:

        >>> lmb._redirect_flowdirs(orig_surf, lake_dict)

        Now, we haven't touched the nodes outside the lakes:

        >>> np.all(np.equal(receivers_init[nodes_not_in_lakes],
        ...                 mg.at_node['flow__receiver_node'][
        ...                     nodes_not_in_lakes]))
        True
        >>> np.all(np.equal(rec_links_init[nodes_not_in_lakes],
        ...     mg.at_node['flow__link_to_receiver_node'][
        ...         nodes_not_in_lakes]))
        True
        >>> np.all(np.equal(steepest_init[nodes_not_in_lakes],
        ...     mg.at_node['topographic__steepest_slope'][nodes_not_in_lakes]))
        True

        ...but the ones inside have been rewired to new patterns:

        >>> np.all(np.equal(receivers_init[nodes_in_lakes],
        ...                 mg.at_node['flow__receiver_node'][nodes_in_lakes]))
        False
        >>> np.all(np.equal(rec_links_init[nodes_in_lakes],
        ...     mg.at_node['flow__link_to_receiver_node'][nodes_in_lakes]))
        False
        >>> np.all(np.equal(steepest_init[nodes_in_lakes],
        ...     mg.at_node['topographic__steepest_slope'][nodes_in_lakes]))
        False

        There are now no pits:

        >>> np.any(mg.at_node['flow__sink_flag'][[7, 15, 22]])
        False

        The lake nodes now flow out:

        >>> mg.at_node['flow__receiver_node'][lake_dict[16]]
        array([16, 16, 15, 16])
        >>> mg.at_node['flow__receiver_node'][lake_dict[8]]
        array([8])

        ...and any outlet nodes that used to drain into the lake now drain
        out. Note that 8 has also rewired itself, since it now drains directly
        into a lake. However, in this case, node 8 already drained outwards;
        it now just drains to a slightly different lake node.

        >>> receivers_init[8]
        9
        >>> mg.at_node['flow__receiver_node'][8]
        15
        >>> receivers_init[16]
        15
        >>> mg.at_node['flow__receiver_node'][16]
        17
        >>> np.is_close(mg.at_node['topographic__steepest_slope'][16], 0.6)
        True
        

        If we reaccumulate the flow, we'll now see that the boundary nodes do
        now accumulate the total available discharge:

        >>> area, discharge = fa.accumulate_flow(update_flow_director=False)
        >>> mg.at_node['drainage_area'][outlets].sum() == (
        ...     mg.cell_area_at_node[mg.core_nodes].sum())
        True
        """
        openq = self._open
        closedq = self.grid.ones('node', dtype=int)
        # Using a slightly different approach. We recognise three types: lake
        # (0), lake margin (1), and closed (2). This lets us work the
        # perimeter too. Open each lake as needed.
        # close the known boundary nodes:
        closedq[self.grid.status_at_node != CORE_NODE] = 2

        # now the actual loop. Work forward lake by lake to avoid unnecessary
        # processing (nodes outside lakes are already correct, by definition).
        for (outlet, lakenodes) in lake_dict.items():
            # open the lake:
            closedq[lakenodes] = 0
            # make a deque for liminal nodes:
            liminal_nodes = deque([])
            openq.add_task(outlet, priority=surface[outlet])

            # it's possible the outlet used to drain *into* the lake,
            # so it needs separate consideration. Likewise, the gradients
            # of the perimeter nodes, although going the right way, are likely
            # to be wrong.
            if self.grid.status_at_node[outlet] != CORE_NODE:
                # don't do anything if the outlet happens to be a boundary
                pass
            else:
                out_elev = LARGE_ELEV
                for neighbor_set, link_set in zip(
                        self._neighbor_arrays, self._link_arrays):
                    not_lake_neighbors = np.not_equal(
                        closedq[neighbor_set[outlet]], 0)
                    minusones = np.equal(neighbor_set[outlet], -1)
                    not_lake_neighbors[minusones] = False
                    closednodes = np.equal(
                        self.grid.status_at_node[neighbor_set[outlet]],
                        CLOSED_BOUNDARY)  # closed BCs can't count
                    not_lake_neighbors[closednodes] = False
                    try:
                        min_val = np.amin(
                            surface[neighbor_set[outlet][not_lake_neighbors]])
                    except ValueError:
                        continue
                    if min_val < out_elev:
                        viable_nodes = neighbor_set[outlet][not_lake_neighbors]
                        min_neighbor_byTrue = np.argmin(
                            surface[viable_nodes])
                        min_neighbor = viable_nodes[min_neighbor_byTrue]
                        min_link = link_set[outlet][not_lake_neighbors][
                            min_neighbor_byTrue]
                        out_elev = min_val
                self._receivers[outlet] = min_neighbor
                self._receiverlinks[outlet] = min_link
                self._steepestslopes[outlet] = (
                    (surface[outlet] - surface[min_neighbor]) /
                    self._neighbor_lengths[min_link])

            print(closedq.reshape(self.grid.shape))

            while True:
                try:
                    c = openq.pop_task()
                except KeyError:
                    break
                else:
                    closedq[c] = 2  # close it
                    # if raster, do the neighbors & diagonals separate...
                    for neighbor_set, link_set in zip(
                            self._neighbor_arrays, self._link_arrays):
                        for (n, l) in zip(neighbor_set[c, :], link_set[c, :]):
                            if closedq[n] == 2:  # fully closed
                                continue
                            elif n == -1:
                                continue
                            elif self.grid.status_at_node[n] != CORE_NODE:
                                closedq[n] = 2
                                continue
                            else:
                                print(n)
                                if closedq[n] == 0:
                                    self._receivers[n] = c
                                    self._receiverlinks[n] = l
                                    self._steepestslopes[n] = 0.
                                    closedq[n] = 2  # close it
                                    openq.add_task(n, priority=surface[n])
                                else:  # it's liminal (1); grads likely wrong
                                    # ...but it's not if set by the outlet...
                                    if c == outlet:
                                        # still need these nodes to be explored
                                        # by other lake nodes as needed, so
                                        # don't close either
                                        pass
                                    else:  # liminal to actual lake
                                        closedq[n] = 2
                                        liminal_nodes.append(n)
                                        # ...& don't add to the queue
            print(liminal_nodes)

            # now know which nodes we need to reassess. So:
            for liminal in liminal_nodes:
                min_elev = LARGE_ELEV
                min_link = -1
                for neighbor_set, link_set in zip(
                        self._neighbor_arrays, self._link_arrays):
                    neighbors = neighbor_set[liminal]
                    neighbors_valid = np.not_equal(neighbors, -1)
                    closednodes = np.equal(
                        self.grid.status_at_node[neighbors],
                        CLOSED_BOUNDARY)  # closed BCs can't count
                    neighbors_valid[closednodes] = False
                    neighbors_to_check = neighbors[neighbors_valid]
                    if len(neighbors_to_check) == 0:
                        continue
                    else:
                        min_neighbor_now = np.amin(
                            self._fill_surface[neighbors_to_check])
                        if min_neighbor_now < min_elev:
                            min_elev = min_neighbor_now
                            links_available = link_set[liminal][neighbors_valid]
                            min_link_of_valid = np.argmin(
                                self._fill_surface[neighbors_to_check])
                            min_receiver = neighbors_to_check[min_link_of_valid]
                            min_link = links_available[min_link_of_valid]
                            max_grad = ((
                                self._fill_surface[liminal] - min_elev) /
                                self._neighbor_lengths[min_link])
                        else:
                            pass
                assert min_link != -1, neighbors_valid  # link successfully found
                self._receivers[liminal] = min_receiver
                self._receiverlinks[liminal] = min_link
                self._steepestslopes[liminal] = max_grad

            # by the time we get here, we've removed all the pits! So...
            print('lakenodes ', lakenodes)
            self.grid.at_node['flow__sink_flag'][lakenodes] = 0
            # reclose the lake:
            closedq[outlet] = 1
            closedq[lakenodes] = 1
            closedq[liminal_nodes] = 1

# this all an absolute mess
# 
#             # reopen the lake for now
#             closedq[lakenodes] = 0
#             print(closedq.reshape(self.grid.shape))
#             # redo the grads of anything in the deque:
#             for liminal in liminal_nodes:
#                 #out_elev = LARGE_ELEV
#                 out_elev_in = self._fill_surface[self.grid.at_node[
#                     'flow__receiver_node'][liminal]]
#                 out_elev = out_elev_in
#                 # the existing lowest elev
#                 min_val = LARGE_ELEV
#                 for neighbor_set, link_set in zip(
#                         self._neighbor_arrays, self._link_arrays):
#                     lake_neighbors = np.equal(
#                         closedq[neighbor_set[liminal]], 0)
#                     #lake_neighbors = np.ones_like(neighbor_set[liminal], dtype=bool)
#                     minusones = np.equal(neighbor_set[liminal], -1)
#                     lake_neighbors[minusones] = False
#                     closednodes = np.equal(
#                         self.grid.status_at_node[neighbor_set[liminal]],
#                         CLOSED_BOUNDARY)  # closed BCs can't count
#                     lake_neighbors[closednodes] = False
#                     try:
#                         min_val = np.amin(
#                             self._fill_surface[neighbor_set[liminal][
#                                 lake_neighbors]])
#                     except ValueError:  # if lake is only in diags, or whatever
#                         continue
#                     print('lim, out_elev, min_val', liminal, out_elev, min_val)
#                     if min_val <= out_elev:
#                         viable_nodes = neighbor_set[liminal][
#                             lake_neighbors]
#                         min_neighbor_byTrue = np.argmin(
#                             self._fill_surface[viable_nodes])
#                         min_neighbor = viable_nodes[min_neighbor_byTrue]
#                         min_link = link_set[liminal][lake_neighbors][
#                             min_neighbor_byTrue]
#                         out_elev = min_val
#                 # if out_elev == out_elev_in:
#                 #     continue
# #####ALL MESSED UP
# # run last run_one_step eg
#                     # if closedq[min_neighbor] == 3:
#                     #     self._receivers[min_neighbor] = liminal
#                     #     self._receiverlinks[min_neighbor] = min_link
#                     #     self._steepestslopes[min_neighbor] = (
#                     #         (self._fill_surface[min_neighbor] -
#                     #          self._fill_surface[liminal]) /
#                     #         self._neighbor_lengths[min_link])
#                     # elif closedq[min_neighbor] == 0:
# 
#                 # now we have the lowest node inside the lake. But we only
#                 # want to overprint if this is actually deeper than what
#                 # we have already. So:
#                 if min_val < out_elev_in:
#                     lake_slope = ((self._fill_surface[liminal] -
#                                    self._fill_surface[min_neighbor]) /
#                                   self._neighbor_lengths[min_link])
#                     print(lake_slope)
#                     if lake_slope > self._steepestslopes[liminal]:
#                         self._receivers[liminal] = min_neighbor
#                         self._receiverlinks[liminal] = min_link
#                         self._steepestslopes[liminal] = lake_slope
#                     # else:  # these cases are found nr the outlet
#                     #     pass


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
        >>> from landlab.components import LakeMapperBarnes, FlowAccumulator
        >>> from landlab.components import FlowDirectorSteepest
        >>> mg = RasterModelGrid((5, 6), 2.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)

        TODO: once return_array_at_node is fixed, this example should also
        take fill_surface...

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
        >>> try:
        ...     lmb.lake_map  # not created, as we aren't tracking
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        AssertionError was raised: Enable tracking to access information about lakes
        >>> lmb.was_there_overfill  # everything fine with slope adding
        False

        >>> fd = FlowDirectorSteepest(mg)
        >>> fa = FlowAccumulator(mg)  # routing will work fine now
        >>> fd.run_one_step()
        >>> fa.run_one_step()
        >>> np.all(mg.at_node['flow__sink_flag'][mg.core_nodes] == 0)
        True
        >>> drainage_area = np.array([  0.,   0.,   0.,   0.,   0.,   0.,
        ...                             0.,   4.,   8.,  12.,   4.,   4.,
        ...                             0.,   4.,  16.,  36.,  40.,  40.,
        ...                             0.,   4.,   8.,   4.,   4.,   4.,
        ...                             0.,   0.,   0.,   0.,   0.,   0.])
        >>> np.allclose(mg.at_node['drainage_area'], drainage_area)
        True

        Test two pits:

        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[23] = 1.3
        >>> z[15] = 0.3
        >>> z[10] = 1.3  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own, if d8
        >>> z[9] = 0.5
        >>> z[14] = 0.6  # [9, 14, 15] is a lake in both methods
        >>> z[16] = 1.2
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16 if d8
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='d8', fill_flat=True,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()  # note the d8 routing now
        >>> lmb.lake_dict == {22: deque([15, 9, 14])}
        True
        >>> lmb.number_of_lakes
        1
        >>> try:
        ...     lmb.lake_depths  # z was both surface and 'fill_surface'
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'surface and fill_surface must be different fields ' +
        ...           'or arrays to enable the property fill_depth!')
        AssertionError was raised: surface and fill_surface must be different fields or arrays to enable the property fill_depth!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        fill_flat=False, track_lakes=True)
        >>> lmb.run_one_step()  # compare to the method='d8' lakes, above...
        >>> lmb.lake_dict == {8: deque([7]), 16: deque([15, 9, 14, 22])}
        True
        >>> lmb.number_of_lakes
        2
        >>> np.allclose(lmb.lake_areas, np.array([ 16.,  4.]))
        True
        >>> try:
        ...     lmb.run_one_step()  # z already filled, so...
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Pit is overfilled due to creation of two outlets as ' +
        ...           'the minimum gradient gets applied. Suppress this ' +
        ...           'Error with the ignore_overfill flag at component ' +
        ...           'instantiation.')
        ValueError was raised: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.

        Suppress this behaviour with ignore_overfill:

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        fill_flat=False, track_lakes=True,
        ...                        ignore_overfill=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict == {8: deque([7]), 16: deque([15, 9, 14, 22])}
        True
        >>> lmb.run_one_step()
        >>> np.allclose(lmb.lake_areas, np.array([ 16.,  4.]))  # found them!
        True

        The component can redirect flow to account for the fills that have
        been carried out (all necessary fields get updated):

        >>> z[:] = z_init
        >>> fd.run_one_step()
        >>> init_flowdirs = mg.at_node['flow__receiver_node'].copy()
        >>> fa.run_one_step()
        >>> init_areas = mg.at_node['drainage_area'].copy()
        >>> init_qw = mg.at_node['surface_water__discharge'].copy()

        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        fill_flat=False, track_lakes=True,
        ...                        redirect_flow_steepest_descent=False,
        ...                        ignore_overfill=True)
        >>> lmb.run_one_step()
        >>> np.all(mg.at_node['flow__receiver_node'] == init_flowdirs)
        True

        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        fill_flat=False, track_lakes=True,
        ...                        redirect_flow_steepest_descent=True,
        ...                        ignore_overfill=True)
        >>> lmb.run_one_step()
        >>> np.all(mg.at_node['flow__receiver_node'] == init_flowdirs)
        False

        However, note that unless the reaccumulate_flow argument is also
        set, the 'drainage_area' and 'surface_water__discharge' fields
        *won't* also get updated:

        >>> np.all(mg.at_node['drainage_area'] == init_areas)
        True
        >>> np.all(mg.at_node['surface_water__discharge'] == init_qw)
        True

        >>> lmb = LakeMapperBarnes(mg, method='steepest',
        ...                        fill_flat=False, track_lakes=True,
        ...                        redirect_flow_steepest_descent=True,
        ...                        reaccumulate_flow=True,
        ...                        ignore_overfill=True)
        >>> lmb.run_one_step()
        >>> np.all(mg.at_node['drainage_area'] == init_areas)
        False
        >>> np.all(mg.at_node['surface_water__discharge'] == init_qw)
        False

        Be sure to set both redirect_flow_steepest_descent and
        reaccumulate_flow to True if you want to reaccumulate flow...

        >>> try:
        ...     lmb = LakeMapperBarnes(mg, method='steepest',
        ...                            fill_flat=False, track_lakes=True,
        ...                            redirect_flow_steepest_descent=False,
        ...                            reaccumulate_flow=True,
        ...                            ignore_overfill=True)
        ... except AssertionError:
        ...     print('Oops!')
        Oops!

        The component is completely happy with irregular grids:

        >>> from landlab import HexModelGrid, FieldError
        >>> hmg = HexModelGrid(5, 4, dx=2.)
        >>> z_hex = hmg.add_zeros('node', 'topographic__elevation')
        >>> z_hex[:] = hmg.node_x
        >>> z_hex[11] = -3.
        >>> z_hex[12] = -1.  # middle nodes become a pit
        >>> z_hex_init = z_hex.copy()
        >>> lmb = LakeMapperBarnes(hmg, fill_flat=True, track_lakes=False)
        >>> lmb.run_one_step()
        >>> np.allclose(z_hex[10:13], 0.)
        True
        >>> z_hex[:] = z_hex_init
        >>> try:
        ...     lmb = LakeMapperBarnes(hmg, fill_flat=False,
        ...                            surface=z_hex_init,
        ...                            redirect_flow_steepest_descent=True,
        ...                            track_lakes=True)
        ... except FieldError:
        ...     print("Oops!")  # flowdir field must already exist!
        Oops!
        >>> fd = FlowDirectorSteepest(hmg)
        >>> lmb = LakeMapperBarnes(hmg, fill_flat=False, surface=z_hex_init,
        ...                        redirect_flow_steepest_descent=True,
        ...                        track_lakes=True)
        >>> fd.run_one_step()
        >>> lmb.run_one_step()
        >>> np.allclose(z_hex[10:13], 0.)
        True
        >>> z_hex[11] > z_hex[10]
        True
        >>> z_hex[12] > z_hex[11]
        True
        >>> np.allclose(lmb.lake_depths[10:14], np.array([ 0.,  3.,  1.,  0.]))
        True
        >>> np.round(lmb.lake_volumes, 4)
        array([ 13.8564])

        Together, all this means that we can now run a topographic growth
        model that permits flooding as it runs:

        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowAccumulator
        >>> from landlab.components import FlowDirectorSteepest
        >>> from landlab.components import StreamPowerEroder
        >>> mg = RasterModelGrid((6, 8))
        >>> for edge in ('right', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x
        >>> z[11] = 1.5
        >>> z[19] = 0.5
        >>> z[34] = 1.1
        >>> z_init = z.copy()
        >>> fd = FlowDirectorSteepest(mg)
        >>> fa = FlowAccumulator(mg)
        >>> lmb = LakeMapperBarnes(mg, fill_flat=True, surface=z_init,
        ...                        redirect_flow_steepest_descent=True,
        ...                        reaccumulate_flow=True,
        ...                        track_lakes=True)
        >>> sp = StreamPowerEroder(mg, K_sp=1., m_sp=1., n_sp=1.)
        >>> fd.run_one_step()
        >>> fa.run_one_step()
        >>> np.isclose(mg.at_node['topographic__steepest_slope'][11], 2.)
        True
        >>> np.allclose(mg.at_node['drainage_area'].reshape(mg.shape)[1, :],
        ...             np.array([ 2.,  2.,  1.,  1.,  4.,  2.,  1.,  0.]))
        True
        >>> lmb.run_one_step()
        >>> np.isclose(mg.at_node['topographic__steepest_slope'][11], 1.)
        True
        >>> np.allclose(mg.at_node['drainage_area'].reshape(mg.shape)[1, :],
        ...             np.array([ 6.,  6.,  5.,  4.,  3.,  2.,  1.,  0.]))
        True
        >>> sp.run_one_step(1.)
        >>> z
        
        """
        # do the prep:
        # increment the run counter
        self._runcount = next(self._runcounter)
        # First get _fill_surface in order.
        self._fill_surface[:] = self._surface  # surfaces begin identical
        # note this is nice & efficent if _fill_surface is _surface
        # if we're doing a redirect, we're going to need to preserve this
        # initial topo, so let's do that:
        if not self._dontredirect:
            orig_topo = self._track_original_surface()
        # now, return _closed to its initial cond, w only the CLOSED_BOUNDARY
        # and grid draining nodes pre-closed:
        closedq = self._closed.copy()
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
            if not self._dontredirect:
                self._redirect_flowdirs(orig_topo, self._lakemappings)
                if self._reaccumulate:
                    _, _ = self._fa.accumulate_flow(update_flow_director=False)

        else:  # not tracked
            # note we've already checked _dontredirect is True in setup,
            # so we don't need to worry about these cases.
            for edgenode in self._edges:
                self._open.add_task(edgenode,
                                    priority=self._surface[edgenode])
            closedq[self._edges] = True
            while True:
                try:
                    self._fill_one_node(
                        self._fill_surface, self._allneighbors, self._pit,
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_dict
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        AssertionError was raised: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict == {16: deque([15, 9, 8, 14, 20, 21])}
        True
        """
        assert self._track_lakes, \
            "Enable tracking to access information about lakes"
        return self._lakemappings

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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_outlets
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        AssertionError was raised: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.number_of_lakes
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        AssertionError was raised: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_map
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        AssertionError was raised: Enable tracking to access information about lakes

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
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
            self._lake_map = np.full(self.grid.number_of_nodes,
                                     LOCAL_BAD_INDEX_VALUE, dtype=int)
            for (outlet, lakenodes) in self.lake_dict.items():
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> try:  # won't work as surface & fill_surface are both z
        ...     lmb.lake_depths
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'surface and fill_surface must be different fields ' +
        ...           'or arrays to enable the property lake_depths!')
        AssertionError was raised: surface and fill_surface must be different fields or arrays to enable the property lake_depths!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        surface=z_init,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lake_depths = np.array([ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        ...                          0. ,  1. ,  0. ,  0.5,  0. ,  0. ,
        ...                          0. ,  0. ,  0.4,  0.7,  0. ,  0. ,
        ...                          0. ,  0. ,  0. ,  0. ,  0.1,  0. ,
        ...                          0. ,  0. ,  0. ,  0. ,  0. ,  0. ])
        >>> np.all(np.equal(lmb.lake_depths,
        ...                 lake_depths))  # slope applied, so...
        False
        >>> np.allclose(lmb.lake_depths, lake_depths)
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_areas  # note track_lakes=False
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        AssertionError was raised: Enable tracking to access information about lakes

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets
        [16, 8]
        >>> np.allclose(lmb.lake_areas, np.array([ 4.,  1.]))
        True
        """
        lakeareas = np.empty(self.number_of_lakes, dtype=float)
        for (i, (outlet, lakenodes)) in enumerate(
             self.lake_dict.items()):
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> try:  # won't work as surface & fill_surface are both z
        ...     lmb.lake_volumes
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'surface and fill_surface must be different fields ' +
        ...           'or arrays to enable the property lake_volumes!')
        AssertionError was raised: surface and fill_surface must be different fields or arrays to enable the property lake_volumes!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        surface=z_init,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets
        [16, 8]
        >>> np.allclose(lmb.lake_volumes, np.array([ 1.7,  1. ]))
        True
        """
        lake_vols = np.empty(self.number_of_lakes, dtype=float)
        col_vols = self.grid.cell_area_at_node * self.lake_depths
        for (i, (outlet, lakenodes)) in enumerate(
             self.lake_dict.items()):
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
        ...                        redirect_flow_steepest_descent=False,
        ...                        ignore_overfill=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.was_there_overfill
        ... except AssertionError:
        ...     print('AssertionError was raised: ' +
        ...           'was_there_overfill is only defined if filling to an ' +
        ...           'inclined surface!')
        AssertionError was raised: was_there_overfill is only defined if filling to an inclined surface!

        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        ignore_overfill=False, track_lakes=True)
        >>> try:
        ...     lmb.run_one_step()
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Pit is overfilled due to creation of two outlets ' +
        ...           'as the minimum gradient gets applied. Suppress this ' +
        ...           'Error with the ignore_overfill flag at component ' +
        ...           'instantiation.')
        ValueError was raised: Pit is overfilled due to creation of two outlets as the minimum gradient gets applied. Suppress this Error with the ignore_overfill flag at component instantiation.

        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        ignore_overfill=True, track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.was_there_overfill
        True

        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1, 1., 1.5]
        >>> lmb.run_one_step()
        >>> lmb.was_there_overfill  # still true as was in the previous example
        True

        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1, 1., 1.5]
        >>> lmb = LakeMapperBarnes(mg, method='steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        ignore_overfill=True, track_lakes=True)
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
