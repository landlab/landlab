cimport cython
from libc.math cimport INFINITY

from landlab.core.messages import warning_message

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
cpdef find_lowest_node_on_lake_perimeter_c(
    const id_t [:, :] neighbor_nodes_at_node,
    id_t [:] flood_status_at_node,
    const cython.floating [:] value_at_node,
    id_t [:] nodes_in_pit,
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
    cdef int _UNFLOODED = 0
    cdef int _PIT = 1
    cdef int _CURRENT_LAKE = 2
    cdef int _FLOODED = 3

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
                    if flood_status_at_node[neighbor_node] == _UNFLOODED:
                        if value_at_node[neighbor_node] < lowest_value:
                            lowest_node = neighbor_node
                            lowest_value = value_at_node[neighbor_node]
                    elif (
                        flood_status_at_node[neighbor_node] == _PIT
                        or flood_status_at_node[neighbor_node] == _FLOODED
                    ):
                        nodes_in_pit[n_pit_nodes] = neighbor_node
                        flood_status_at_node[neighbor_node] = _CURRENT_LAKE
                        n_pit_nodes += 1
            i += 1

    return lowest_node, n_pit_nodes
