import numpy as np
cimport numpy as np
cimport cython

from landlab.core.messages import warning_message


DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t


@cython.boundscheck(False)
cpdef find_lowest_node_on_lake_perimeter_c(
        np.ndarray[DTYPE_INT_t, ndim=2] node_nbrs,
        np.ndarray[DTYPE_INT_t, ndim=1] flood_status,
        np.ndarray[DTYPE_FLOAT_t, ndim=1] elev,
        np.ndarray[DTYPE_INT_t, ndim=1] nodes_this_depression,
        DTYPE_INT_t pit_count,
        DTYPE_FLOAT_t BIG_ELEV
    ):
    """Locate the lowest node on the margin of the "lake".

    Parameters
    ----------
    node_nbrs : (nnodes, 4) or (nnodes, 8) array of int
        The node neighbors, as stored by a DepressionFinderAndRouter
        component
    flood_status : nnodes array of int
        The node flooded status at the point of the function call, as stored
        by a DepressionFinderAndRouter component
    elev : nnodes array of float
        The elevations of each node in the grid
    nodes_this_depression : nnodes array of int
        Nodes that form a pit, followed by padding values to make up a nnodes-
        long array. This should be passed in with the first value as the
        pit node, then padding values, but it will be updated in place to
        reflect the nodes in the lake.
    pit_count : int
        The number of nodes currently in the lake.

    Returns
    -------
    (int, int)
        (Lowest node on the perimeter of a depression, updated pit_count)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, DepressionFinderAndRouter
    >>> from landlab.components.flow_routing.cfuncs import find_lowest_node_on_lake_perimeter_c
    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node['topographic__elevation'].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> fr = FlowAccumulator(mg, flow_director='D8')
    >>> fr.run_one_step()  # the flow "gets stuck" in the hole
    >>> df = DepressionFinderAndRouter(mg)

    >>> node_nbrs = df._node_nbrs
    >>> flood_status = df.flood_status
    >>> elev = df._elev
    >>> BIG_ELEV = df._BIG_ELEV
    >>> nodes_this_depression = mg.zeros('node', dtype=int)
    >>> nodes_this_depression[0] = 16
    >>> pit_count = 1

    >>> find_lowest_node_on_lake_perimeter_c(
    ...     node_nbrs, flood_status, elev, nodes_this_depression, pit_count,
    ...     BIG_ELEV
    ... )
    (23, 1)
    >>> nodes_this_depression[1] = 8
    >>> pit_count = 2
    >>> find_lowest_node_on_lake_perimeter_c(
    ...     node_nbrs, flood_status, elev, nodes_this_depression, pit_count,
    ...     BIG_ELEV
    ... )
    (0, 2)
    """
    # Start with the first node on the list, and an arbitrarily large elev
    cdef int lowest_node = nodes_this_depression[0]
    cdef DTYPE_FLOAT_t lowest_elev = BIG_ELEV

    # set up a worst-case scanario array for the pits, and a counter to pull
    # the good entries later:
    cdef int current_iter = 0

    # Codes for depression status
    cdef int _UNFLOODED = 0
    cdef int _PIT = 1
    cdef int _CURRENT_LAKE = 2
    cdef int _FLOODED = 3

    cdef int n
    cdef int nbr
    cdef int i

    while current_iter < pit_count:
        n = nodes_this_depression[current_iter]
        for nbr in node_nbrs[n]:
            if nbr != -1:
                if flood_status[nbr] == _UNFLOODED:
                    if elev[nbr] < lowest_elev:
                        lowest_node = nbr
                        lowest_elev = elev[nbr]
                elif (
                    flood_status[nbr] == _PIT
                    or flood_status[nbr] == _FLOODED
                ):
                    nodes_this_depression[pit_count] = nbr
                    pit_count += 1
                    flood_status[nbr] = _CURRENT_LAKE
        current_iter += 1
    if lowest_elev == BIG_ELEV:
        print("Unable to find drainage outlet for a lake.")
        print("In lake with " + str(len(nodes_this_depression)), "nodes:")
        print(str(nodes_this_depression))

        for i in nodes_this_depression:
            print("Node ID: ", i)
            print("Node Elevation: ", elev[i])
            print("Node Flood Status: ", flood_status[i])
            print("Node Neigbors: ", node_nbrs[i])
            print("Neighbor Elevations: ", elev[node_nbrs[i]])
            print("Neigbor Flood Status: ", flood_status[node_nbrs[i]])
        warning_message(
            """If you see no data values in any of the elevation terms
            this may because you have disconnected open nodes (which
            sometimes occurs during raster clipping.

            Consider running
            set_open_nodes_disconnected_from_watershed_to_closed
            which will remove isolated open nodes."""
        )
    return lowest_node, pit_count
