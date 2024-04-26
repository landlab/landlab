# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
# distutils: extra_compile_args = -std=c++11
# distutils: extra_link_args = -std=c++11
""" Contains the cython functions for the component method
Flow_router.run_directions(). Flow directions and depression handling are done
adapting Barnes et al., 2014 algorithm #4."""

import numpy as np

cimport cython
cimport numpy as cnp
from libcpp cimport bool
from libcpp.pair cimport pair


# 1-3. Instantiate the Queues (Steps noted #1-3, #5-6 in Barnes, 2014's algorithm #4
####################################################################################
cdef extern from "_priority_queue.hpp" nogil:
    cdef cppclass _priority_queue:
        _priority_queue(...) except +
        void push(pair[cnp.int_t, cnp.float_t])
        pair[cnp.int_t, cnp.float_t] top() except +
        void pop()
        bool empty()
        cnp.int_t size()


cdef bool _compare_second(pair[int, double] a, pair[int, double] b) nogil:
    return a.second > b.second


@cython.boundscheck(False)  # turn off bounds-checking for entire function
cdef void _init_flow_direction_queues(
    const cnp.int_t [:] base_level_nodes,
    const cnp.int_t [:] closed_nodes,
    cnp.float_t [:] z,
    _priority_queue& to_do,
    cnp.int_t [:] receivers,
    cnp.int_t [:] outlet_nodes,
    cnp.int_t [:] done,
    cnp.int_t* done_n_ptr,
) nogil:
    """
    Add the base-level nodes to the queue and update receivers for base-level and
    closed nodes. Updates to_do, receivers, outlet_nodes, done and the value pointed
    by done_n_ptr.

    Implementation remarks: It's important to pass the priority_queue to_do as a
    reference, otherwise it won't be modified. I didn't manage to use the bint type
    for the variable done, so here 0 is for False and 1 for True.

    Params
    ------
    base_level_nodes: memoryview(long).
        Base-level nodes.
    closed_nodes: memoryview(long).
        Closed nodes.
    z: memoryview(double).
        Elevation of eachn
    to_do: _priority_queue&
        Queue which handles the nodes where to direct flow
    receivers: memoryview(long).
        Receiver of each node (ordered by node id).
    outlet_nodes: memoryview(long).
        Base-level outlet of each node (ordered by node id).
    done: memoryview(long).
        Flag for each node (ordered by node id). If 1, the node has been done,
        otherwise 0.
    done_n_ptr: long*
        Pointer to the number of done nodes. Necessary for multithreading
        (future evolution).
    """

    cdef:
        cnp.int_t node_id, i, n = len(base_level_nodes), m = len(closed_nodes)
        pair[cnp.int_t, cnp.float_t] node_pair

    for i in range(n):
        # NB: for node_i in open_boundary raises a compiling error with nogil.
        node_id = base_level_nodes[i]
        node_pair = pair[cnp.int_t, cnp.float_t](node_id, z[node_id])
        to_do.push(node_pair)
        receivers[node_id] = node_id
        outlet_nodes[node_id] = node_id
        done[node_id] = 1
        done_n_ptr[0] += 1

    for i in range(m):
        node_id = closed_nodes[i]
        receivers[node_id] = node_id
        outlet_nodes[node_id] = node_id
        done[node_id] = 1
        done_n_ptr[0] += 1
        # NB: done_n_ptr[0]: method to dereference the pointer. The
        # cython.operator.dereference doesn't seem to work on this pointer to a
        # long-type variable.

# 4. Functions necessary for the flow direction processing (Steps #11 - 20)
###########################################################################


@cython.boundscheck(False)
cdef void _set_flooded_and_outlet(
    cnp.int_t donor_id,
    cnp.float_t [:] z,
    cnp.int_t [:] receivers,
    cnp.int_t [:] outlet_nodes,
    cnp.int_t [:] depression_outlet_nodes,
    cnp.int_t [:] flooded_nodes,
    cnp.float_t [:] depression_depths,
    cnp.float_t [:] depression_free_elevations,
    cnp.int_t flooded_status,
    cnp.int_t bad_index,
    cnp.float_t min_elevation_relative_diff,
) nogil:
    """ Updates the base-level outlet nodes (outlet_nodes), the depression outlet
    nodes (depression_outlet_nodes), the flooded status (flooded_nodes), and the
    depths of the depressions (depression_depths) for the node donor_id depending
    on its surface z value and the one of its receiver outlet and neighbors.

    Parameters
    ----------
    donor_id: long.
        Id of the node giving the flow.
    z: memoryview(double).
        Values of the surface where the flow is directed, for each node.
    receivers: memoryview(long).
        Ids of the receiver nodes for each node considered as a donor, ordered by
        donor id.
    outlet_nodes: memoryview(long).
        Ids of the base-level outlet node, for each donor node, ordered by the donor
        node ids.
    depression_outlet_nodes: memoryview(long).
        Outlet nodes of the depression.
    flooded_nodes: memoryview(bool).
        Flooded status for each node.
    depression_depths: memoryview(double).
        Depths of the depression (if existing) below the level of the depression outlet
        for each node, ordered by node id.
    depression_free_elevations: memoryview(double).
        Elevation of the surface corrected from depressions.
    flooded_status: long.
        Constant for flooded status.
    bad_index: long.
        Constant for bad index.
    min_elevation_relative_diff: double
        Minimum relative difference in elevation for the depression_free_elevations
        surface.
    """
    cdef:
        cnp.int_t receiver_id = receivers[donor_id]
        cnp.int_t receiver_depression_outlet = (
            receiver_id if (
                depression_outlet_nodes[receiver_id] == bad_index
            ) else depression_outlet_nodes[receiver_id]
        )
        cnp.int_t receiver_outlet = (
            receiver_id if (
                outlet_nodes[receiver_id] == bad_index
            ) else outlet_nodes[receiver_id]
        )
    outlet_nodes[donor_id] = receiver_outlet
    if z[donor_id] < z[receiver_depression_outlet]:
        depression_outlet_nodes[donor_id] = receiver_depression_outlet
        flooded_nodes[donor_id] = flooded_status
        depression_depths[donor_id] = z[receiver_depression_outlet] - z[donor_id]
        depression_free_elevations[donor_id] = (
            (1 + min_elevation_relative_diff) * depression_free_elevations[receiver_id]
        )


@cython.boundscheck(False)
cdef void _set_receiver(
    cnp.int_t donor_id,
    cnp.int_t receiver_id,
    cnp.int_t [:] receivers,
    cnp.int_t [:] done,
    cnp.int_t* done_n_ptr,
) nogil:
    """ Updates the receiver (receivers) and the process statuses (done) of the donor
    node donor_id.

    Parameters
    ----------
    donor_id: long
        Id of the donor node.
    receiver_id: long
        Id of the receiver node.
    receivers: memoryview(long)
        Ids of the receivers for all donor nodes, ordered by the id of the donor nodes.
    done: memoryview(bool)
        Process statuses for all nodes. 1 for done. 0 otherwise.
    """
    receivers[donor_id] = receiver_id
    done[donor_id] = 1  # 1 For True
    done_n_ptr[0] += 1


@cython.boundscheck(False)
cdef void _set_donor_properties(
    cnp.int_t donor_id,
    cnp.int_t receiver_id,
    cnp.int_t [:] sorted_pseudo_tails,
    const cnp.int_t [:, :] head_start_end_indexes,
    const cnp.int_t [:] sorted_dupli_links,
    cnp.float_t [:] sorted_dupli_gradients,
    cnp.float_t [:] z,
    cnp.float_t [:] steepest_slopes,
    cnp.int_t [:] links_to_receivers,
) nogil:
    """ Updates the steepest_slopes and the links_to_receivers of a donor in function
    of the slopes with its neighbors and the head-tail links of the grid. Steepest
    slope is set to 0. if the donor is in a depression.

    Parameters
    ----------
    donor_id: long
        Id of the donor node.
    receiver_id: long.
        Id of the receiver node.
    sorted_pseudo_tails: memoryview(long).
        All tails and heads for all existing links in the grid, ordered by
        increasing head and tail ids in sorted_pseudo_heads (see FlowRouter.py file).
    head_start_end_indexes: memoryview(longxlong).
        Start (array[0, :] and end (array[1, :]) index of all possible value nodes in
        the sorted_pseudo_heads (ordered by increasing head and tail ids, see
        FlowRouter.py file).
    sorted_dupli_links: memoryview(long).
        All existing links in the grid, duplicated, ordered by
        increasing head and tail ids in sorted_pseudo_heads (see FlowRouter.py file).
    sorted_dupli_gradients: memoryview(double).
        Gradients of the surface for all links duplicated, ordered by  increasing
        head and tail ids in sorted_pseudo_heads.
    z: memoryview(double).
        Values of the surface where the flow is directed, for each node.
    steepest_slopes: memoryview(double).
        Steepest slope for each node.
    links_to_receivers: memoryview(long).
        Ids of the links between a donor and a receiver, ordered by the ids of the
        donor nodes.
    """

    # range of indexes where donor_id is founded in sorted_pseudo_heads
    cdef:
        cnp.int_t idx1 = head_start_end_indexes[0, donor_id]
        cnp.int_t idx2 = head_start_end_indexes[1, donor_id] + 1

        cnp.int_t [:] s = sorted_pseudo_tails[idx1:idx2]
        cnp.int_t n = len(s), c = -1, i

    # loop to bypass the impossibility to use
    # sorted_pseudo_tails[idx1:idx2] == receiver_id with memoryviews
    for i in range(n):
        if s[i] == receiver_id:
            c = idx1 + i
            break
    steepest_slopes[donor_id] = (
        sorted_dupli_gradients[c] if (
            z[receiver_id] <= z[donor_id]
        ) else 0
    )
    links_to_receivers[donor_id] = sorted_dupli_links[c]

# ######################################################################################
# Main functions to direct flow


@cython.boundscheck(False)
cdef void _direct_flow_c(
    cnp.int_t nodes_n,
    const cnp.int_t[:] base_level_nodes,
    const cnp.int_t[:] closed_nodes,
    cnp.int_t[:] sorted_pseudo_tails,
    cnp.float_t[:] sorted_dupli_gradients,
    const cnp.int_t[:] sorted_dupli_links,
    const cnp.int_t[:, :] head_start_end_indexes,
    cnp.int_t [:] outlet_nodes,
    cnp.int_t [:] depression_outlet_nodes,
    cnp.int_t[:] flooded_nodes,
    cnp.float_t[:] depression_depths,
    cnp.float_t[:] depression_free_elevations,
    cnp.int_t[:] links_to_receivers,
    cnp.int_t[:] receivers,
    cnp.float_t[:] steepest_slopes,
    cnp.float_t[:] z,
    cnp.int_t flooded_status,
    cnp.int_t bad_index,
    cnp.int_t neighbors_max_number,
    cnp.float_t min_elevation_relative_diff,
):
    """
    Main function implementing the flow directing through breaching depressions.
    Updates outlet_nodes, depression_outlet_nodes, flooded_nodes, links_to_receivers,
    receivers, steepest_slopes.

    Params
    ------
    nodes_n: long.
        Number of nodes to handle.
    base_level_nodes: memoryview(long).
        Base-level nodes.
    closed_nodes: memoryview(long).
        Closed nodes.
    sorted_pseudo_tails: memoryview(long).
        All tails and heads for all existing links in the grid, ordered by
        increasing head and tail ids in sorted_pseudo_heads (see FlowRouter.py file).
    sorted_dupli_gradients: memoryview(double).
        Gradients of the surface for all links duplicated, ordered by increasing
        head and tail ids in sorted_pseudo_heads.
    sorted_dupli_links: memoryview(long).
        All existing links in the grid, duplicated, ordered by
        increasing head and tail ids in sorted_pseudo_heads (see FlowRouter.py file).
    head_start_end_indexes: memoryview(longxlong).
        Start (array[0, :] and end (array[1, :]) index of all possible value nodes in
        the sorted_pseudo_heads (ordered by increasing head and tail ids, see
        FlowRouter.py file).
    outlet_nodes: memoryview(long).
        Ids of the base-level outlet node, for each donor node, ordered by the donor
        node ids.
    depression_outlet_nodes: memoryview(long).
        Outlet nodes of the depression, ordered by the donor node ids.
    flooded_nodes: memoryview(bool).
        Flooded status for each node, ordered by node ids.
    depression_depths: memoryview(double).
        Depths of the depression (if existing) below the level of the depression outlet
        for each node, ordered by node ids.
    depression_free_elevations: memoryview(double).
        Elevation of the surface corrected from depressions.
    links_to_receivers: memoryview(long).
        Ids of the links between a donor and a receiver, ordered by the ids of the
        donor nodes.
    receivers: memoryview(long)
        Ids of the receivers for all donor nodes, ordered by the id of the donor nodes.
    steepest_slopes: memoryview(double).
        Steepest slope for each node, ordered by node id.
    z: memoryview(double).
        Values of the surface where the flow is directed, for each node, ordered by
        node id.
    flooded_status: long.
        Constant for flooded status.
    bad_index: long.
        Constant for bad index.
    neighbors_max_number: long
        Maximum number of neighbors in the grid.
    min_elevation_relative_diff: double
        Minimum relative difference in elevation for the depression_free_elevations
        surface.
    """
    cdef:
        cnp.int_t [:] done
        # cnp.int_t [:] tmp_neighbors
        # cnp.int_t [:] neighbors_to_do
        _priority_queue to_do = _priority_queue(_compare_second)
        cnp.int_t receiver_id, donor_id, i, j, done_n
        # cnp.int_t [:] neighbors
        pair[cnp.int_t, cnp.float_t] node_pair

    done = np.full(nodes_n, 0, dtype=int)
    # tmp_neighbors = np.full(neighbors_max_number, 0, dtype=int)
    # neighbors_to_do = np.array([], dtype=int)

    # done_n is input only for MULTITHREADING, a future evolution, and
    # is not checked in this function.
    done_n = 0

    _init_flow_direction_queues(
        base_level_nodes, closed_nodes, z, to_do, receivers, outlet_nodes, done, &done_n
    )

    for j in range(nodes_n):
        # a while loop is possible here, but prefer for loop, with future
        # multithreading evolution.
        if to_do.empty():
            break
        receiver_id = to_do.top().first
        to_do.pop()
        done[0] = 1

        # Get the neighbors to handle.
        idx1 = head_start_end_indexes[0, receiver_id]
        idx2 = head_start_end_indexes[1, receiver_id] + 1

        # Handle each neighbor.
        for i in range(idx1, idx2):
            donor_id = sorted_pseudo_tails[i]
            if done[donor_id] == 1:
                continue  # 0 for False.
            _set_receiver(donor_id, receiver_id, receivers, done, &done_n)
            _set_flooded_and_outlet(
                donor_id,
                z,
                receivers,
                outlet_nodes,
                depression_outlet_nodes,
                flooded_nodes,
                depression_depths,
                depression_free_elevations,
                flooded_status,
                bad_index,
                min_elevation_relative_diff,
            )
            _set_donor_properties(
                donor_id,
                receiver_id,
                sorted_pseudo_tails,
                head_start_end_indexes,
                sorted_dupli_links,
                sorted_dupli_gradients,
                z,
                steepest_slopes,
                links_to_receivers,
            )
            node_pair = pair[cnp.int_t, cnp.float_t](donor_id, z[donor_id])
            to_do.push(node_pair)


def _direct_flow(
    cnp.int_t nodes_n,
    const cnp.int_t[:] base_level_nodes,
    const cnp.int_t[:] closed_nodes,
    cnp.int_t[:] sorted_pseudo_tails,
    cnp.float_t[:] sorted_dupli_gradients,
    const cnp.int_t[:] sorted_dupli_links,
    const cnp.int_t[:, :] head_start_end_indexes,
    cnp.int_t [:] outlet_nodes,
    cnp.int_t [:] depression_outlet_nodes,
    cnp.int_t[:] flooded_nodes,
    cnp.float_t[:] depression_depths,
    cnp.float_t[:] depression_free_elevations,
    cnp.int_t[:] links_to_receivers,
    cnp.int_t[:] receivers,
    cnp.float_t[:] steepest_slopes,
    cnp.float_t[:] z,
    cnp.int_t flooded_status,
    cnp.int_t bad_index,
    cnp.int_t neighbors_max_number,
    cnp.float_t min_elevation_relative_diff,
):
    """
    Main function calling the function that implements flow directing through
    breaching depressions. Updates outlet_nodes, depression_outlet_nodes,
    flooded_nodes, links_to_receivers, receivers, steepest_slopes.

    Params
    ------
    nodes_n: long.
        Number of nodes to handle.
    base_level_nodes: memoryview(long).
        Base-level nodes.
    closed_nodes: memoryview(long).
        Closed nodes.
    sorted_pseudo_tails: memoryview(long).
        All tails and heads for all existing links in the grid, ordered by
        increasing head and tail ids in sorted_pseudo_heads (see FlowRouter.py file).
    sorted_dupli_gradients: memoryview(double).
        Gradients of the surface for all links duplicated, ordered by increasing
        head and tail ids in sorted_pseudo_heads.
    sorted_dupli_links: memoryview(long).
        All existing links in the grid, duplicated, ordered by
        increasing head and tail ids in sorted_pseudo_heads (see FlowRouter.py file).
    head_start_end_indexes: memoryview(longxlong).
        Start (array[0, :] and end (array[1, :]) index of all possible value nodes in
        the sorted_pseudo_heads (ordered by increasing head and tail ids, see
        FlowRouter.py file).
    outlet_nodes: memoryview(long).
        Ids of the base-level outlet node, for each donor node, ordered by the donor
        node ids.
    depression_outlet_nodes: memoryview(long).
        Outlet nodes of the depression, ordered by the donor node ids.
    flooded_nodes: memoryview(bool).
        Flooded status for each node, ordered by node ids.
    depression_depths: memoryview(double).
        Depths of the depression (if existing) below the level of the depression outlet
        for each node, ordered by node ids.
    depression_free_elevations: memoryview(double).
        Elevation of the surface corrected from depressions.
    links_to_receivers: memoryview(long).
        Ids of the links between a donor and a receiver, ordered by the ids of the
        donor nodes.
    receivers: memoryview(long)
        Ids of the receivers for all donor nodes, ordered by the id of the donor nodes.
    steepest_slopes: memoryview(double).
        Steepest slope for each node, ordered by node id.
    z: memoryview(double).
        Values of the surface where the flow is directed, for each node, ordered by
        node id.
    flooded_status: long.
        Constant for flooded status.
    bad_index: long.
        Constant for bad index.
    neighbors_max_number: long
        Maximum number of neighbors in the grid.
    min_elevation_relative_diff: double
        Minimum relative difference in elevation for the depression_free_elevations
        surface.
    """
    _direct_flow_c(
        nodes_n,
        base_level_nodes,
        closed_nodes,
        sorted_pseudo_tails,
        sorted_dupli_gradients,
        sorted_dupli_links,
        head_start_end_indexes,
        outlet_nodes,
        depression_outlet_nodes,
        flooded_nodes,
        depression_depths,
        depression_free_elevations,
        links_to_receivers,
        receivers,
        steepest_slopes,
        z,
        flooded_status,
        bad_index,
        neighbors_max_number,
        min_elevation_relative_diff,
    )
