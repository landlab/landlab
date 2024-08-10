cimport cython

from cython.parallel import prange

from libc.math cimport INFINITY
from libc.stdint cimport uint8_t


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef find_lowest_node_on_lake_perimeter_c(
    const cython.integral [:, :] neighbor_nodes_at_node,
    cython.integral [:] flood_status_at_node,
    const cython.floating [:] value_at_node,
    cython.integral [:] nodes_in_pit,
    long n_pit_nodes,
):
    """Locate the lowest node on the margin of the "lake".

    Parameters
    ----------
    neighbor_nodes_at_node : (nnodes, 4) or (nnodes, 8) array of int
        The node neighbors, as stored by a DepressionFinderAndRouter
        component
    flood_status_at_node : nnodes array of int
        The node flooded status at the point of the function call, as stored
        by a DepressionFinderAndRouter component
    value_at_node : nnodes array of float
        The value of each node in the grid
    nodes_in_pit : nnodes array of int
        Nodes that form a pit, followed by padding values to make up a nnodes-
        long array. This should be passed in with the first value as the
        pit node, then padding values, but it will be updated in place to
        reflect the nodes in the lake.
    n_pit_nodes : int
        The number of nodes currently in the lake.

    Returns
    -------
    (int, int)
        Lowest node on the perimeter of a depression, updated n_pit_nodes

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, DepressionFinderAndRouter
    >>> from landlab.components.flow_routing.cfuncs import (
    ...     find_lowest_node_on_lake_perimeter_c
    ... )
    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field("topographic__elevation", mg.node_x.copy(), at="node")
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node["topographic__elevation"].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> fr = FlowAccumulator(mg, flow_director='D8')
    >>> fr.run_one_step()  # the flow "gets stuck" in the hole
    >>> df = DepressionFinderAndRouter(mg)

    >>> node_nbrs = df._node_nbrs
    >>> flood_status = df.flood_status
    >>> elev = df._elev
    >>> nodes_this_depression = mg.zeros('node', dtype=int)
    >>> nodes_this_depression[0] = 16
    >>> pit_count = 1

    >>> find_lowest_node_on_lake_perimeter_c(
    ...     node_nbrs, flood_status, elev, nodes_this_depression, pit_count,
    ... )
    (23, 1)
    >>> nodes_this_depression[1] = 8
    >>> pit_count = 2
    >>> find_lowest_node_on_lake_perimeter_c(
    ...     node_nbrs, flood_status, elev, nodes_this_depression, pit_count,
    ... )
    (0, 2)
    """
    # Start with the first node on the list, and an arbitrarily large value
    cdef long lowest_node = nodes_in_pit[0]
    cdef double lowest_value = INFINITY

    cdef long i

    # Codes for depression status
    cdef int UNFLOODED = 0
    cdef int PIT = 1
    cdef int CURRENT_LAKE = 2
    cdef int FLOODED = 3

    cdef int n_neighbors = neighbor_nodes_at_node.shape[1]
    cdef long node
    cdef long neighbor_node
    cdef long neighbor

    with nogil:
        i = 0
        while i < n_pit_nodes:
            node = nodes_in_pit[i]
            for neighbor in range(n_neighbors):
                neighbor_node = neighbor_nodes_at_node[node, neighbor]
                if neighbor_node != -1:
                    if flood_status_at_node[neighbor_node] == UNFLOODED:
                        if value_at_node[neighbor_node] < lowest_value:
                            lowest_node = neighbor_node
                            lowest_value = value_at_node[neighbor_node]
                    elif (
                        flood_status_at_node[neighbor_node] == PIT
                        or flood_status_at_node[neighbor_node] == FLOODED
                    ):
                        nodes_in_pit[n_pit_nodes] = neighbor_node
                        flood_status_at_node[neighbor_node] = CURRENT_LAKE
                        n_pit_nodes += 1
            i += 1

    return lowest_node, n_pit_nodes


@cython.boundscheck(False)
@cython.wraparound(False)
def find_pits(
    const cython.floating [:] value_at_node,
    const cython.integral [:, :] nodes_at_link,
    const cython.integral [:] links,
    const uint8_t [:] is_open_node,
    uint8_t [:] is_pit_node,
):
    """
    A node is defined as being a pit if and only if:

    1. All neighboring core nodes have equal or greater elevation, and
    2. Any neighboring open boundary nodes have a greater elevation.

    The algorithm starts off assuming that all core nodes are pits. We then
    loop through all active links. For each link, if one node is higher
    than the other, the higher one cannot be a pit, so we flag it False.
    We also look at cases in which an active link's nodes have equal
    elevations. If one is an open boundary, then the other must be a core
    node, and we declare the latter not to be a pit (via rule 2 above).
    """
    cdef long i
    cdef long link
    cdef long node_at_head
    cdef long node_at_tail
    cdef long n_links = len(links)

    for i in prange(n_links, nogil=True, schedule="static"):
        link = links[i]

        node_at_head = nodes_at_link[link, 0]
        node_at_tail = nodes_at_link[link, 1]

        if value_at_node[node_at_head] > value_at_node[node_at_tail]:
            is_pit_node[node_at_head] = False
        elif value_at_node[node_at_tail] > value_at_node[node_at_head]:
            is_pit_node[node_at_tail] = False
        else:
            if is_open_node[node_at_head]:
                is_pit_node[node_at_tail] = False
            if is_open_node[node_at_tail]:
                is_pit_node[node_at_head] = False
