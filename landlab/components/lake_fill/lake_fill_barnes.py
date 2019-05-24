#!/usr/env/python

"""
lake_fill_barnes.py

Fill sinks in a landscape to the brim, following the Barnes et al. (2014)
algorithms.
"""

from __future__ import print_function

import heapq
import itertools

# ^ this simply in case Katy updates to add more fields, that we would also
# need to update...
from collections import deque

import numpy as np
from six import iteritems

from landlab import (
    BAD_INDEX_VALUE,
    CLOSED_BOUNDARY,
    CORE_NODE,
    FIXED_GRADIENT_BOUNDARY,
    FIXED_VALUE_BOUNDARY,
    Component,
    RasterModelGrid,
)
from landlab.components import FlowAccumulator, FlowDirectorSteepest
from landlab.utils import StablePriorityQueue
from landlab.utils.return_array import return_array_at_node

LOCAL_BAD_INDEX_VALUE = BAD_INDEX_VALUE
LARGE_ELEV = 9999999999.0

# TODO: Needs to have rerouting functionality...


def _fill_one_node_to_flat(fill_surface, all_neighbors, pitq, openq, closedq, dummy):
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
    >>> mg = RasterModelGrid((5, 6))
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
    openneighbors = cneighbors[np.logical_not(closedq[cneighbors])]  # for efficiency
    closedq[openneighbors] = True
    for n in openneighbors:
        if fill_surface[n] <= fill_surface[c]:
            fill_surface[n] = fill_surface[c]
            heapq.heappush(pitq, n)
        else:
            openq.add_task(n, priority=fill_surface[n])


class LakeMapperBarnes(Component):
    """
    A Landlab implementation of the Barnes et al. (2014) lake filling & lake
    routing algorithms, lightly modified and adapted for Landlab by DEJH. This
    component is designed as a direct replacement for the LakeMapper as
    existing pre-Aug 2018, and provides a suite of properties to access
    information about the lakes created each time it is run. Only significant
    difference is the way the lakes are coded: this component uses the
    (unique) ID of the outlet node, whereas DepressionFinderAndRouter uses
    one of the pit node IDs. Note also this component does not offer the
    `lake_codes` or `display_depression_map` options, for essentially this
    reason. Use `lake_map` instead for both. It also uses a much more
    Landlabbian `run_one_step()` method as its driver, superceding
    DepressionFinderAndRouter's `map_depressions()`.

    A variety of options is provided. Flow routing is route-to-one in this
    implementation, but can be either D4 ("steepest") or D8 on a raster.
    The surface can be filled to either flat or a very slight downward
    incline, such that subsequent flow routing will run over the lake surface.
    This incline is applied at machine precision to minimise the chances of
    creating false outlets by overfill, and note that the gradient as
    calculated on such surfaces may still appear to be zero.
    The filling can either be performed in place, or on a new (water) surface
    distinct from the original (rock) surface. For efficiency, data structures
    describing the lakes and their properties are only created, and existing
    flow direction and flow accumulation fields modified, if those flags are
    set at instantiation.

    With care, this component can be used to create a dynamically flooding
    surface in a fluvial landscape (interacting with, e.g., the
    StreamPowerEroder). See the run_one_step docstring for an example.

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    surface : field name at node or array of length node
        The surface to direct flow across.
    method : {'Steepest', 'D8'}
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
        the 'drainage_area', 'surface_water__discharge',
        'flow__upstream_node_order', and the other various flow accumulation
        fields (see output field names) will now reflect the new drainage
        patterns without having to manually reaccumulate the discharge. If
        True but redirect_flow_steepest_descent is False, raises an
        ValueError.
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

    _name = "LakeMapperBarnes"

    _cite_as = """@article{BARNES2014117,
        title = "Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models",
        journal = "Computers & Geosciences",
        volume = "62",
        pages = "117 - 127",
        year = "2014",
        issn = "0098-3004",
        doi = "https://doi.org/10.1016/j.cageo.2013.04.024",
        url = "http://www.sciencedirect.com/science/article/pii/S0098300413001337",
        author = "Richard Barnes and Clarence Lehman and David Mulla",
        keywords = "Pit filling, Terrain analysis, Hydrology, Drainage network, Modeling, GIS"
        }"""

    _input_var_names = (
        "topographic__elevation",
        "drainage_area",
        "surface_water__discharge",
        "flow__link_to_receiver_node",
        "flow__upstream_node_order",
        "flow__data_structure_delta",
        "flow__data_structure_D",
        "flow__receiver_node",
        "flow__sink_flag",
    )

    _output_var_names = (
        "topographic__elevation",
        "drainage_area",
        "surface_water__discharge",
        "flow__link_to_receiver_node",
        "flow__upstream_node_order",
        "flow__data_structure_delta",
        "flow__data_structure_D",
        "flow__receiver_node",
        "flow__sink_flag",
    )

    _var_units = {
        "topographic__elevation": "m",
        "drainage_area": "m**2",
        "surface_water__discharge": "m**3/s",
        "flow__link_to_receiver_node": "-",
        "flow__upstream_node_order": "-",
        "flow__data_structure_delta": "-",
        "flow__data_structure_D": "-",
        "flow__receiver_node": "-",
        "flow__sink_flag": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "drainage_area": "node",
        "surface_water__discharge": "node",
        "flow__link_to_receiver_node": "node",
        "flow__upstream_node_order": "node",
        "flow__data_structure_delta": "node",
        "flow__data_structure_D": "link",
        "flow__receiver_node": "node",
        "flow__sink_flag": "node",
    }

    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "drainage_area": "Upstream accumulated surface area contributing to the node's "
        "discharge",
        "surface_water__discharge": "Discharge of water through each node",
        "flow__link_to_receiver_node": "ID of link downstream of each node, which carries the discharge",
        "flow__upstream_node_order": "Node array containing downstream-to-upstream ordered list of "
        "node IDs",
        "flow__data_structure_delta": "Node array containing the elements delta[1:] of the data "
        'structure "delta" used for construction of the downstream-to-'
        "upstream node array",
        "flow__data_structure_D": "Link array containing the data structure D used for construction"
        "of the downstream-to-upstream node array",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current "
        "node)",
        "flow__sink_flag": "Boolean array, True at local lows",
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        method="Steepest",
        fill_flat=True,
        fill_surface="topographic__elevation",
        redirect_flow_steepest_descent=False,
        reaccumulate_flow=False,
        ignore_overfill=False,
        track_lakes=True,
    ):
        """
        Initialize the component.
        """
        if "flow__receiver_node" in grid.at_node:
            if grid.at_node["flow__receiver_node"].size != grid.size("node"):
                msg = (
                    "A route-to-multiple flow director has been "
                    "run on this grid. The landlab development team has not "
                    "verified that LakeMapperBarnes is compatible with "
                    "route-to-multiple methods. Please open a GitHub Issue "
                    "to start this process."
                )
                raise NotImplementedError(msg)
        self._grid = grid
        self._open = StablePriorityQueue()
        self._pit = []
        self._closed = self.grid.zeros("node", dtype=bool)
        self._gridclosednodes = self.grid.status_at_node == CLOSED_BOUNDARY
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
        if method not in {"Steepest", "D8"}:
            raise ValueError(
                "{method}: method must be 'Steepest' or 'D8'".format(method=method)
            )
        if method == "D8":
            if isinstance(grid, RasterModelGrid):
                self._allneighbors = np.concatenate(
                    (
                        self.grid.adjacent_nodes_at_node,
                        self.grid.diagonal_adjacent_nodes_at_node,
                    ),
                    axis=1,
                )
            else:  # not a raster
                raise ValueError(
                    (
                        "D8 is not a valid value for method if grid type is "
                        + "{gridtype}!"
                    ).format(gridtype=type(grid))
                )
        else:
            self._allneighbors = self.grid.adjacent_nodes_at_node

        # A key difference from the "pure" Barnes algorithm for LL is that
        # we must'n flood from all the edges. Instead, we can only flood from
        # a true grid edge, i.e., only the FIXED boundary types. (Both
        # CLOSED and LOOPED assume flow either can't get out there, or at
        # least, there's more land in that direction that will be handled
        # otherwise.) Note we add a test that there is at least some kind
        # of outlet!!
        self._edges = np.where(
            np.logical_or(
                self.grid.status_at_node == FIXED_VALUE_BOUNDARY,
                self.grid.status_at_node == FIXED_GRADIENT_BOUNDARY,
            )
        )[0]
        if self._edges.size == 0:
            raise ValueError(
                "No valid outlets found on the grid! You cannot run the "
                + "filling algorithms!"
            )
        # and finally, close these up permanently as well (edges will always
        # be edges...)
        self._closed[self._edges] = True
        # ...note there's a slight of hand here, Because of the ordering of LL
        # grids, the last node will always be a boundary node, even for very
        # odd Voronois. This enables us to treat out -1s in the neighbour
        # arrays as always True. But, just in case...
        assert self._closed[-1]

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
            if not track_lakes:
                raise ValueError("You must track_lakes to redirect the flow!")
            # Check we have the necessary fields already existing.
            # These will raise FieldErrors if they don't.
            # This will cause a bunch of our tests to break, so users will
            # never see this.
            assert len(FlowDirectorSteepest.output_var_names) == 4
            self._receivers = self.grid.at_node["flow__receiver_node"]
            self._receiverlinks = self.grid.at_node["flow__link_to_receiver_node"]
            self._steepestslopes = self.grid.at_node["topographic__steepest_slope"]
            # if raster, do the neighbors & diagonals separate when rerouting
            # so we'll need to pull these separately:
            if method == "D8":  # Raster test unnecessary given tests above
                self._neighbor_arrays = (
                    self.grid.adjacent_nodes_at_node,
                    self.grid.diagonal_adjacent_nodes_at_node,
                )
                self._link_arrays = (
                    self.grid.links_at_node,
                    self.grid.d8s_at_node[:, 4:],
                )
                self._neighbor_lengths = self.grid.length_of_d8
            else:
                self._neighbor_arrays = (self.grid.adjacent_nodes_at_node,)
                self._link_arrays = (self.grid.links_at_node,)
                self._neighbor_lengths = self.grid.length_of_link

        if reaccumulate_flow:
            if not redirect_flow_steepest_descent:
                raise ValueError(
                    "You must also redirect_flow_steepest_descent if you "
                    + "want to reaccumulate_flow!"
                )
            self._reaccumulate = True
            self._fa = FlowAccumulator(self.grid, flow_director=method)
        else:
            self._reaccumulate = False

        self._fill_flat = fill_flat
        if fill_flat:
            self._fill_one_node = _fill_one_node_to_flat
        else:
            self._fill_one_node = self._fill_one_node_to_slant

    def _fill_one_node_to_slant(
        self, fill_surface, all_neighbors, pitq, openq, closedq, ignore_overfill
    ):
        """
        Implements the Barnes et al. algorithm to obtain a naturally draining
        surface, updating a single node. Assumes the _open and _closed lists
        have already been updated per Barnes algos 2&3, lns 1-7.

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
        >>> mg = RasterModelGrid((5, 6))
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest')
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

        >>> mg = RasterModelGrid((3, 7))
        >>> for edge in ('top', 'right', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1,
        ...                                 1.0000000000000004, 1.5]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest')
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
            # not clear how this occurs, but present in Barnes ->
            # DEJH suspects this should be an elevation comparison given the
            # text description. Regardless, this is only to ensure
            # repeatability, so it's not vital even if these cases don't
            # trigger
            if topopen == toppit:  # intentionally tight comparison
                # print('yessssss')
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
            nextval = np.nextafter(fill_surface[c], LARGE_ELEV)
            # DEJH believes that in the LL use cases this is impossible,
            # but retained as comments since present in Barnes algorithm
            # if self._gridclosednodes[n]:
            #     heapq.heappush(pitq, n)
            if fill_surface[n] <= nextval:  # former elif
                if self._PitTop < fill_surface[n] and nextval >= fill_surface[n]:
                    if ignore_overfill:
                        self._overfill_flag = True
                    else:
                        raise ValueError(
                            "Pit is overfilled due to creation of two "
                            + "outlets as the minimum gradient gets applied. "
                            + "Suppress this Error with the ignore_overfill "
                            + "flag at component instantiation."
                        )
                fill_surface[n] = nextval
                heapq.heappush(pitq, n)
            else:
                openq.add_task(n, priority=fill_surface[n])

    def _fill_to_flat_with_tracking(
        self, fill_surface, all_neighbors, pitq, openq, closedq
    ):
        """
        Implements the Barnes et al. algorithm for a simple fill over the
        grid. Assumes the _open and _closed lists have already been updated
        per Barnes algos 2&3, lns 1-7.

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
        >>> mg = RasterModelGrid((5, 6))
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest')
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
                    lakemappings[outlet_ID] = deque([c])

            cneighbors = all_neighbors[c]
            openneighbors = cneighbors[
                np.logical_not(closedq[cneighbors])
            ]  # for efficiency
            closedq[openneighbors] = True
            for n in openneighbors:
                if fill_surface[n] <= fill_surface[c]:
                    fill_surface[n] = fill_surface[c]
                    heapq.heappush(pitq, n)
                else:
                    openq.add_task(n, priority=fill_surface[n])
            # print(np.sort(openq.tasks_currently_in_queue()), pitq)
        return lakemappings

    def _fill_to_slant_with_optional_tracking(
        self,
        fill_surface,
        all_neighbors,
        pitq,
        openq,
        closedq,
        ignore_overfill,
        track_lakes,
    ):
        """
        Implements the Barnes et al. algorithm to obtain a naturally draining
        surface over the grid. Assumes the _open and _closed lists have
        already been updated per Barnes algos 2&3, lns 1-7.

        This version runs a little more slowly to enable tracking of which
        nodes are linked to which outlets.

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
        >>> from landlab.components import LakeMapperBarnes, FlowAccumulator
        >>> mg = RasterModelGrid((5, 6))
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest')
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest')
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
        >>> fr = FlowAccumulator(mg, flow_director='D4')
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest')
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

        >>> mg = RasterModelGrid((3, 7))
        >>> for edge in ('top', 'right', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1,
        ...                                 1.0000000000000004, 1.5]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest')
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
            # as above, DEJH is unclear how this clause triggers, so
            # retained but untested ->
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
                        lakemappings[outlet_ID] = deque([c])
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
                nextval = np.nextafter(fill_surface[c], LARGE_ELEV)
                # as in non-tracker, DEJH believes this is redundant in LL
                # if self._gridclosednodes[n]:
                #     heapq.heappush(pitq, n)
                if fill_surface[n] <= nextval:  # formerly elif
                    if self._PitTop < fill_surface[n] and nextval >= fill_surface[n]:
                        if ignore_overfill:
                            self._overfill_flag = True
                        else:
                            raise ValueError(
                                "Pit is overfilled due to creation of two "
                                + "outlets as the minimum gradient gets "
                                + "applied. Suppress this Error with the "
                                + "ignore_overfill flag at component "
                                + "instantiation."
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
        >>> mg = RasterModelGrid((5, 6), xy_spacing=2.)
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z_new = mg.add_zeros('node', 'topographic__fill', dtype=float)
        >>> lmb = LakeMapperBarnes(mg, method='D8',
        ...                        surface='topographic__elevation',
        ...                        fill_surface='topographic__fill',
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> orig_surf = lmb._track_original_surface()
        >>> z is orig_surf
        True
        >>> lmb = LakeMapperBarnes(mg, method='D8',
        ...                        surface='topographic__elevation',
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
        be exactly 0., even if fill_flat is False. This is because we are
        adding an increment to elevation at machine precision.

        Examples
        --------
        >>> import numpy as np
        >>> from collections import deque
        >>> from landlab import RasterModelGrid
        >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
        >>> from landlab.components import FlowDirectorSteepest
        >>> from landlab.components import FlowAccumulator
        >>> mg = RasterModelGrid((5, 6), xy_spacing=2.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[23] = 1.3
        >>> z[15] = -2.  # this deep pit causes the outlet to first drain *in*
        >>> z[10] = 1.3  # raise "guard" exit nodes
        >>> z[7] = 2.  # is a lake on its own, if D8
        >>> z[9] = -1.
        >>> z[14] = 0.6  # [9, 14, 15] is a lake in both methods
        >>> z[16] = 1.2
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16 if D8
        >>> z_init = z.copy()
        >>> fd = FlowDirectorSteepest(mg)
        >>> fa = FlowAccumulator(mg)
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=True,
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

        Now the flow directions all ignore the pits:

        >>> np.all(np.equal(mg.at_node['flow__receiver_node'],
        ...                 np.array([ 0,  1,  2,  3,  4,  5,
        ...                            6,  8,  9, 15,  9, 11,
        ...                           12, 14, 15, 16, 17, 17,
        ...                           18, 20, 14, 15, 16, 23,
        ...                           24, 25, 26, 27, 28, 29])))
        True

        (Note the filling of the pits might redirect the occasional node
        not in a lake, but on its perimeter - if the node used to drain
        into the lake, but now has a steeper descent path elsewhere.)

        There are now no pits:

        >>> np.any(mg.at_node['flow__sink_flag'][[7, 15, 22]])
        False

        The lake nodes now flow out:

        >>> mg.at_node['flow__receiver_node'][lake_dict[16]]
        array([16, 15, 15, 16])
        >>> mg.at_node['flow__receiver_node'][lake_dict[8]]
        array([8])

        ...and any outlet nodes that used to drain into the lake now drain
        out.

        >>> receivers_init[16]
        15
        >>> mg.at_node['flow__receiver_node'][16]
        17
        >>> np.isclose(mg.at_node['topographic__steepest_slope'][16], 0.6)
        True

        If we reaccumulate the flow, we'll now see that the boundary nodes do
        now accumulate the total available discharge:

        >>> area, discharge = fa.accumulate_flow(update_flow_director=False)
        >>> mg.at_node['drainage_area'][outlets].sum() == (
        ...     mg.cell_area_at_node[mg.core_nodes].sum())
        True
        """
        openq = self._open
        closedq = self.grid.ones("node", dtype=int)
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
            # of the perimeter nodes are likely to be wrong.
            if self.grid.status_at_node[outlet] != CORE_NODE:
                # don't do anything if the outlet happens to be a boundary
                pass
            else:
                out_elev = LARGE_ELEV
                for neighbor_set, link_set in zip(
                    self._neighbor_arrays, self._link_arrays
                ):
                    not_lake_neighbors = np.not_equal(closedq[neighbor_set[outlet]], 0)
                    minusones = np.equal(neighbor_set[outlet], -1)
                    not_lake_neighbors[minusones] = False
                    closednodes = np.equal(
                        self.grid.status_at_node[neighbor_set[outlet]], CLOSED_BOUNDARY
                    )  # closed BCs can't count
                    not_lake_neighbors[closednodes] = False
                    try:
                        min_val = np.amin(
                            surface[neighbor_set[outlet][not_lake_neighbors]]
                        )
                    except ValueError:
                        continue
                    if min_val < out_elev:
                        viable_nodes = neighbor_set[outlet][not_lake_neighbors]
                        min_neighbor_byTrue = np.argmin(surface[viable_nodes])
                        min_neighbor = viable_nodes[min_neighbor_byTrue]
                        min_link = link_set[outlet][not_lake_neighbors][
                            min_neighbor_byTrue
                        ]
                        out_elev = min_val
                self._receivers[outlet] = min_neighbor
                self._receiverlinks[outlet] = min_link
                self._steepestslopes[outlet] = (
                    surface[outlet] - surface[min_neighbor]
                ) / self._neighbor_lengths[min_link]

            while True:
                try:
                    c = openq.pop_task()
                except KeyError:
                    break
                else:
                    closedq[c] = 2  # close it
                    # if raster, do the neighbors & diagonals separate...
                    for neighbor_set, link_set in zip(
                        self._neighbor_arrays, self._link_arrays
                    ):
                        for (n, l) in zip(neighbor_set[c, :], link_set[c, :]):
                            if closedq[n] == 2:  # fully closed
                                continue
                            elif n == -1:
                                continue
                            elif self.grid.status_at_node[n] != CORE_NODE:
                                closedq[n] = 2
                                continue
                            else:
                                if closedq[n] == 0:
                                    self._receivers[n] = c
                                    self._receiverlinks[n] = l
                                    self._steepestslopes[n] = 0.0
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

            # TODO: obvious case for Cython accel here
            # Now know which nodes we need to reassess. So:
            for liminal in liminal_nodes:
                min_elev = LARGE_ELEV
                min_link = -1
                for neighbor_set, link_set in zip(
                    self._neighbor_arrays, self._link_arrays
                ):
                    neighbors = neighbor_set[liminal]
                    neighbors_valid = np.not_equal(neighbors, -1)
                    closednodes = np.equal(
                        self.grid.status_at_node[neighbors], CLOSED_BOUNDARY
                    )  # closed BCs can't count
                    neighbors_valid[closednodes] = False
                    neighbors_to_check = neighbors[neighbors_valid]
                    if len(neighbors_to_check) == 0:
                        continue
                    else:
                        min_neighbor_now = np.amin(
                            self._fill_surface[neighbors_to_check]
                        )
                        if min_neighbor_now < min_elev:
                            min_elev = min_neighbor_now
                            links_available = link_set[liminal][neighbors_valid]
                            min_link_of_valid = np.argmin(
                                self._fill_surface[neighbors_to_check]
                            )
                            min_receiver = neighbors_to_check[min_link_of_valid]
                            min_link = links_available[min_link_of_valid]
                            max_grad = (
                                self._fill_surface[liminal] - min_elev
                            ) / self._neighbor_lengths[min_link]
                        else:
                            pass
                assert min_link != -1, neighbors_valid
                # ^link successfully found
                self._receivers[liminal] = min_receiver
                self._receiverlinks[liminal] = min_link
                self._steepestslopes[liminal] = max_grad

            # by the time we get here, we've removed all the pits! So...
            self.grid.at_node["flow__sink_flag"][lakenodes] = 0
            # reclose the lake:
            closedq[outlet] = 1
            closedq[lakenodes] = 1
            closedq[liminal_nodes] = 1

    def run_one_step(self):
        """
        Fills the surface to fill all pits. Note that a second run on a
        surface that has already been filled will *not* "see" any existing
        lakes correctly - it will see lakes, but with zero depths. In
        particular, if fill_flat is False, an attempt to fill a
        surface a second time will raise a ValueError unless ignore_overfill.
        (In this case, setting ignore_overfill is True will give the expected
        behaviour.)

        If reaccumulate_flow was True at instantiation, run_one_step also
        updates all the various flow fields produced by the FlowDirector and
        FlowAccumulator components.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowAccumulator
        >>> from landlab.components import FlowDirectorSteepest
        >>> mg = RasterModelGrid((5, 6), xy_spacing=2.)
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init,
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
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        ValueError was raised: Enable tracking to access information about lakes
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
        >>> z[7] = 2.  # is a lake on its own, if D8
        >>> z[9] = 0.5
        >>> z[14] = 0.6  # [9, 14, 15] is a lake in both methods
        >>> z[16] = 1.2
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16 if D8
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='D8', fill_flat=True,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()  # note the D8 routing now
        >>> lmb.lake_dict == {22: deque([15, 9, 14])}
        True
        >>> lmb.number_of_lakes
        1
        >>> try:
        ...     lmb.lake_depths  # z was both surface and 'fill_surface'
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'surface and fill_surface must be different fields ' +
        ...           'or arrays to enable the property fill_depth!')
        ValueError was raised: surface and fill_surface must be different fields or arrays to enable the property fill_depth!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='Steepest',
        ...                        fill_flat=False, track_lakes=True)
        >>> lmb.run_one_step()  # compare to the method='D8' lakes, above...
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest',
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

        >>> lmb = LakeMapperBarnes(mg, method='Steepest',
        ...                        fill_flat=False, track_lakes=True,
        ...                        redirect_flow_steepest_descent=False,
        ...                        ignore_overfill=True)
        >>> lmb.run_one_step()
        >>> np.all(mg.at_node['flow__receiver_node'] == init_flowdirs)
        True

        >>> lmb = LakeMapperBarnes(mg, method='Steepest',
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

        >>> lmb = LakeMapperBarnes(mg, method='Steepest',
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
        ...     lmb = LakeMapperBarnes(mg, method='Steepest',
        ...                            fill_flat=False, track_lakes=True,
        ...                            redirect_flow_steepest_descent=False,
        ...                            reaccumulate_flow=True,
        ...                            ignore_overfill=True)
        ... except ValueError:
        ...     print('Oops!')
        Oops!

        The component is completely happy with irregular grids:

        >>> from landlab import HexModelGrid, FieldError
        >>> hmg = HexModelGrid(5, 4, dx=2.)
        >>> z_hex = hmg.add_zeros('node', 'topographic__elevation')
        >>> z_hex[:] = hmg.node_x
        >>> z_hex[11] = -3.
        >>> z_hex[12] = -1.
        >>> z_hex_init = z_hex.copy()
        >>> z_hex
        array([   2.,   4.,   6.,   8.,
               1.,   3.,   5.,   7.,   9.,
            0.,   2.,   -3.,  -1.,   8.,  10.,
               1.,   3.,   5.,   7.,   9.,
                  2.,   4.,  6.,   8.])

        As you can see, nodes 11 and 12 are now a pit. If they were to fill
        they would fill to the level of 2, the lowest downstream value.

        >>> lmb = LakeMapperBarnes(
        ...     hmg,
        ...     method='Steepest',
        ...     fill_flat=True,
        ...     track_lakes=False)
        >>> lmb.run_one_step()
        >>> np.allclose(z_hex[10:13], 2.)
        True
        >>> z_hex[:] = z_hex_init
        >>> try:
        ...     lmb = LakeMapperBarnes(hmg, method='Steepest',
        ...                            fill_flat=False,
        ...                            surface=z_hex_init,
        ...                            redirect_flow_steepest_descent=True,
        ...                            track_lakes=True)
        ... except FieldError:
        ...     print("Oops!")  # flowdir field must already exist!
        Oops!
        >>> fd = FlowDirectorSteepest(hmg)
        >>> lmb = LakeMapperBarnes(hmg, method='Steepest',
        ...                        fill_flat=False, surface=z_hex_init,
        ...                        redirect_flow_steepest_descent=True,
        ...                        track_lakes=True)
        >>> fd.run_one_step()
        >>> lmb.run_one_step()
        >>> np.allclose(z_hex[10:13], 2.)
        True
        >>> z_hex[11] > z_hex[10]
        True
        >>> z_hex[12] > z_hex[11]
        True
        >>> np.allclose(lmb.lake_depths[10:14], np.array([ 0.,  5.,  3.,  0.]))
        True
        >>> np.round(lmb.lake_volumes, 4)
        array([ 27.7128])

        Together, all this means that we can now run a topographic growth
        model that permits flooding as it runs:

        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes, FlowAccumulator
        >>> from landlab.components import FlowDirectorSteepest
        >>> from landlab.components import FastscapeEroder
        >>> mg = RasterModelGrid((6, 8))
        >>> for edge in ('right', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY

        Because it is what we want the FastscapeEroder to see and work on,
        it's actually the water surface that needs to go in as
        'topographic__elevation'. We'll also need to keep track of the bed
        elevation though, since the LakeMapper will need it. We start them
        equal (i.e., topo starts dry).

        >>> z_water = mg.add_zeros(
        ...     'node', 'topographic__elevation', dtype=float)
        >>> z_water[:] = mg.node_x
        >>> z_water[11] = 1.5
        >>> z_water[19] = 0.5
        >>> z_water[34] = 1.1
        >>> z_bed = mg.add_zeros(
        ...     'node', 'bedrock__elevation', dtype=float)
        >>> z_bed[:] = z_water  # topo starts dry

        Let's just take a look:

        >>> np.all(np.equal(
        ...     np.round(z_water, 2),
        ...     np.array([0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 1.5, 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 0.5, 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 1.1, 3. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ])))
        True

        >>> fd = FlowDirectorSteepest(mg)
        >>> fa = FlowAccumulator(mg)
        >>> lmb = LakeMapperBarnes(mg, method='D8', fill_flat=True,
        ...                        surface='bedrock__elevation',
        ...                        fill_surface='topographic__elevation',
        ...                        redirect_flow_steepest_descent=True,
        ...                        reaccumulate_flow=True,
        ...                        track_lakes=True)
        >>> sp = FastscapeEroder(mg, K_sp=1., m_sp=0., n_sp=1.)
        >>> fd.run_one_step()
        >>> fa.run_one_step()  # node 18 is draining into the pit...
        >>> np.isclose(mg.at_node['topographic__steepest_slope'][18], 1.5)
        True
        >>> np.allclose(mg.at_node['drainage_area'],
        ...             np.array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        ...                        2.,  2.,  1.,  4.,  3.,  2.,  1.,  0.,
        ...                        1.,  1.,  1., 13.,  3.,  2.,  1.,  0.,
        ...                        2.,  2.,  1.,  4.,  3.,  2.,  1.,  0.,
        ...                        6.,  6.,  5.,  4.,  3.,  2.,  1.,  0.,
        ...                        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]))
        True
        >>> lmb.run_one_step()  # now node 18 drains correctly, outward ->
        >>> np.isclose(mg.at_node['topographic__steepest_slope'][18], 1.)
        True
        >>> np.allclose(mg.at_node['drainage_area'],
        ...             np.array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        ...                       13., 13., 12.,  4.,  3.,  2.,  1.,  0.,
        ...                        2.,  2.,  1.,  7.,  3.,  2.,  1.,  0.,
        ...                        2.,  2.,  1.,  1.,  3.,  2.,  1.,  0.,
        ...                        7.,  7.,  6.,  4.,  3.,  2.,  1.,  0.,
        ...                        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]))
        True
        >>> np.all(np.equal(
        ...     np.round(z_water, 2),
        ...     np.array([0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 2. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 2. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 1.1, 3. , 4. , 5. , 6. , 7. ,
        ...               0. , 1. , 2. , 3. , 4. , 5. , 6. , 7. ])))
        True

        >>> sp.run_one_step(0.05)  # note m=0 to illustrate effect of slopes
        >>> np.all(np.equal(
        ...     np.round(z_water, 2),
        ...     np.array([0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ,
        ...               0.  , 0.95, 1.95, 2.  , 3.9 , 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.95, 2.  , 3.9 , 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.95, 2.93, 3.93, 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.09, 2.91, 3.95, 4.95, 5.95, 7.  ,
        ...               0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ])))
        True

        If we want to keep this going honouring the depths of the lakes try
        this next in your loop:

        >>> z_bed[:] = np.minimum(z_water, z_bed)
        >>> np.all(np.equal(
        ...     np.round(z_bed, 2),
        ...     np.array([0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ,
        ...               0.  , 0.95, 1.95, 1.5 , 3.9 , 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.95, 0.5 , 3.9 , 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.95, 2.93, 3.93, 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.09, 2.91, 3.95, 4.95, 5.95, 7.  ,
        ...               0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ])))
        True
        >>> fd.run_one_step()
        >>> fa.run_one_step()
        >>> lmb.run_one_step()

        Lake node depths are now updated in lmb:

        >>> np.round(
        ...     [lmb.lake_depths[lake] for lake in lmb.lake_dict.values()], 2)
        array([[ 0.45,  1.45]])

        ...and the "topography" (i.e., water surface) at the flooded nodes
        has lowered itself as the lip of the outlet was eroded in the last
        step:

        >>> np.all(np.equal(
        ...     np.round(z_water, 2),
        ...     np.array([0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ,
        ...               0.  , 0.95, 1.95, 1.95 , 3.9 , 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.95, 1.95 , 3.9 , 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.95, 2.93, 3.93, 4.95, 5.95, 7.  ,
        ...               0.  , 0.95, 1.09, 2.91, 3.95, 4.95, 5.95, 7.  ,
        ...               0.  , 1.  , 2.  , 3.  , 4.  , 5.  , 6.  , 7.  ])))
        True

        >>> sp.run_one_step(0.05)

        ...and so on.

        Note that this approach, without passing `flooded_nodes` to the
        FastscapeEroder run method, is both more "Landlabbic" and also
        ensures the information about the lake and the water surface
        topography are all updated cleanly and correctly.
        """
        if "flow__receiver_node" in self._grid.at_node:
            if self._grid.at_node["flow__receiver_node"].size != self._grid.size(
                "node"
            ):
                msg = (
                    "A route-to-multiple flow director has been "
                    "run on this grid. The landlab development team has not "
                    "verified that LakeMapperBarnes is compatible with "
                    "route-to-multiple methods. Please open a GitHub Issue "
                    "to start this process."
                )
                raise NotImplementedError(msg)
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
                    self._fill_surface,
                    self._allneighbors,
                    self._pit,
                    self._open,
                    closedq,
                )
            else:
                self._lakemappings = self._fill_to_slant_with_optional_tracking(
                    self._fill_surface,
                    self._allneighbors,
                    self._pit,
                    self._open,
                    closedq,
                    ignore_overfill=self._ignore_overfill,
                    track_lakes=True,
                )
            if not self._dontredirect:
                self._redirect_flowdirs(orig_topo, self._lakemappings)
                if self._reaccumulate:
                    _, _ = self._fa.accumulate_flow(update_flow_director=False)

        else:  # not tracked
            # note we've already checked _dontredirect is True in setup,
            # so we don't need to worry about these cases.
            for edgenode in self._edges:
                self._open.add_task(edgenode, priority=self._surface[edgenode])
            closedq[self._edges] = True
            while True:
                try:
                    self._fill_one_node(
                        self._fill_surface,
                        self._allneighbors,
                        self._pit,
                        self._open,
                        closedq,
                        self._ignore_overfill,
                    )
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
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_dict
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        ValueError was raised: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_dict == {16: deque([15, 9, 8, 14, 20, 21])}
        True
        """
        if not self._track_lakes:
            raise ValueError("Enable tracking to access information about lakes")
        return self._lakemappings

    @property
    def lake_outlets(self):
        """
        Returns the outlet for each lake, not necessarily in ID order.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
        >>> for edge in ('left', 'top', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[2, 1:-1] = [2., 1., 0.5, 1.5]
        >>> z.reshape(mg.shape)[1, 1:-1] = [2.1, 1.1, 0.6, 1.6]
        >>> z.reshape(mg.shape)[3, 1:-1] = [2.2, 1.2, 0.7, 1.7]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_outlets
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        ValueError was raised: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets == [16, ]
        True
        """
        if not self._track_lakes:
            raise ValueError("Enable tracking to access information about lakes")
        return list(self._lakemappings.keys())

    @property
    def number_of_lakes(self):
        """
        Return the number of individual lakes. Lakes sharing outlet nodes are
        considered part of the same lake.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.number_of_lakes
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        ValueError was raised: Enable tracking to access information about lakes

        >>> lmb = LakeMapperBarnes(mg, method='Steepest', surface=z_init,
        ...                        fill_surface=z, fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.number_of_lakes
        2
        """
        if not self._track_lakes:
            raise ValueError("Enable tracking to access information about lakes")
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
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_map
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        ValueError was raised: Enable tracking to access information about lakes

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
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
            self._lake_map = np.full(
                self.grid.number_of_nodes, LOCAL_BAD_INDEX_VALUE, dtype=int
            )
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
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
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
        """
        Return the change in surface elevation at each node this step.
        Requires that surface and fill_surface were not the same array at
        instantiation.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> try:  # won't work as surface & fill_surface are both z
        ...     lmb.lake_depths
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'surface and fill_surface must be different fields ' +
        ...           'or arrays to enable the property lake_depths!')
        ValueError was raised: surface and fill_surface must be different fields or arrays to enable the property lake_depths!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
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
        if self._inplace:
            raise ValueError(
                "surface and fill_surface must be different fields or "
                + "arrays to enable the property lake_depths!"
            )
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
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=False)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.lake_areas  # note track_lakes=False
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'Enable tracking to access information about lakes')
        ValueError was raised: Enable tracking to access information about lakes

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> lmb.lake_outlets
        [16, 8]
        >>> np.allclose(lmb.lake_areas, np.array([ 4.,  1.]))
        True
        """
        lakeareas = np.empty(self.number_of_lakes, dtype=float)
        for (i, (outlet, lakenodes)) in enumerate(iteritems(self.lake_dict)):
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

        Requires that surface and fill_surface were not the same array at
        instantiation.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import LakeMapperBarnes
        >>> mg = RasterModelGrid((5, 6))
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
        ...                        redirect_flow_steepest_descent=False,
        ...                        track_lakes=True)
        >>> lmb.run_one_step()
        >>> try:  # won't work as surface & fill_surface are both z
        ...     lmb.lake_volumes
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'surface and fill_surface must be different fields ' +
        ...           'or arrays to enable the property lake_volumes!')
        ValueError was raised: surface and fill_surface must be different fields or arrays to enable the property lake_volumes!

        >>> z[:] = z_init
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
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
        for (i, (outlet, lakenodes)) in enumerate(iteritems(self.lake_dict)):
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
        >>> mg = RasterModelGrid((3, 7))
        >>> for edge in ('top', 'right', 'bottom'):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
        >>> z = mg.add_zeros('node', 'topographic__elevation', dtype=float)
        >>> z.reshape(mg.shape)[1, 1:-1] = [1., 0.2, 0.1,
        ...                                 1.0000000000000004, 1.5]
        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=True,
        ...                        redirect_flow_steepest_descent=False,
        ...                        ignore_overfill=False, track_lakes=True)
        >>> lmb.run_one_step()
        >>> try:
        ...     lmb.was_there_overfill
        ... except ValueError:
        ...     print('ValueError was raised: ' +
        ...           'was_there_overfill is only defined if filling to an ' +
        ...           'inclined surface!')
        ValueError was raised: was_there_overfill is only defined if filling to an inclined surface!

        >>> z_init = z.copy()
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
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
        >>> lmb = LakeMapperBarnes(mg, method='Steepest', fill_flat=False,
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
        if self._fill_flat is True:
            raise ValueError(
                "was_there_overfill is only defined if filling to an "
                + "inclined surface!"
            )
        return self._overfill_flag
