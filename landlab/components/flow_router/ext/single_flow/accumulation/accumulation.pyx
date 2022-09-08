#distutils: language = c++
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

""" Contains the cython functions for the component method 
Flow_router.run_accumulations(). Flow accumulation and upstream node ordering are done
adapting Braun and Willett, 2013 algorithm."""

import numpy as np
cimport numpy as cnp
cimport cython
from libcpp cimport bool


cdef cnp.int_t _add_to_upstream_ordered_nodes(cnp.int_t receiver_id,
    cnp.int_t node_idx_in_stack,
    cnp.int_t [:] upstream_ordered_nodes, cnp.int_t [:] donors_start_indexes,
    cnp.int_t [:] donors):
    """
    Adapted version of the recursive algorithm of Braun and Willett, 2013.
    This function updates node_idx_in_stack and upstream_ordered_nodes.
    
    NB for future evolutions: this recursive function might be not parallelizable as
        such.

    Params
    ------
    receiver_id: long
        Id of the flow receiver node. (noted l in Braun and Willett, 2013).
    node_idx_in_stack: long.
        Current index in the stack to be used to add a new receiver id (noted j).
    upstream_ordered_nodes: memoryview(long).
        Stack of nodes ordered downstream to upstream, built from the recursive call
        of _add_to_upstream_ordered_nodes().
    donors_start_indexes: memoryview(long).
        Starting indexes in donors for each receiver. E.g. 3 for receiver n and 4 for
        receiver n+1 means that the donors for receiver n are in index 3 in donors
        (noted delta_i).
    donors: memoryview(long).
        Donors ordered by the id of their receivers (noted D).
        
    Returns
    -------
    node_idx_in_stack: long.
        Current index in the stack to be used to add a new receiver id (noted j).
    """
    cdef cnp.int_t donor_index
    cdef cnp.int_t donor_id
    cdef cnp.int_t idx1 = donors_start_indexes[receiver_id]
    cdef cnp.int_t idx2 = donors_start_indexes[receiver_id + 1]
    cdef cnp.int_t _return_node_idx_in_stack = 0
    upstream_ordered_nodes[node_idx_in_stack] = receiver_id
    node_idx_in_stack += 1

    for donor_index in range(idx1, idx2): 
        donor_id = donors[donor_index]      
        if donor_id != receiver_id: 
            node_idx_in_stack = _add_to_upstream_ordered_nodes(donor_id,
                                    node_idx_in_stack, upstream_ordered_nodes,
                                    donors_start_indexes, donors)
    return node_idx_in_stack


cdef void _calc_upstream_order_for_nodes_c(cnp.int_t[:] base_level_and_closed_nodes,
    cnp.int_t[:] upstream_ordered_nodes, cnp.int_t[:] donors_start_indexes,
    cnp.int_t[:] donors):
    cdef cnp.int_t node_idx_in_stack = 0
    cdef cnp.int_t n = len(base_level_and_closed_nodes)
    cdef cnp.int_t receiver_id
    cdef cnp.int_t i
    """
    Orders nodes downstream to upstream for each base-level and add closed nodes.
    
    Params
    ------
    base_level_and_closed_nodes: memoryview(long).
        Base level and closed nodes
    upstream_ordered_nodes: memoryview(long).
        Stack of nodes ordered downstream to upstream, built from the recursive call
        of _add_to_upstream_ordered_nodes().
    donors_start_indexes: memoryview(long).
        Starting indexes in donors for each receiver. E.g. 3 for receiver n and 4 for
        receiver n+1 means that the donors for receiver n are in index 3 in donors
        (noted delta_i).
    donors: memoryview(long).
        Donors ordered by the id of their receivers (noted D).
    
    Return
    ------
    void
    """

    for i in range(n):
        receiver_id = base_level_and_closed_nodes[i]
        node_idx_in_stack = _add_to_upstream_ordered_nodes(receiver_id,
                            node_idx_in_stack,
                            upstream_ordered_nodes,
                            donors_start_indexes, donors)


def _calc_upstream_order_for_nodes(cnp.int_t[:] base_level_and_closed_nodes,
    cnp.int_t[:] upstream_ordered_nodes,
    cnp.int_t[:] donors_start_indexes, cnp.int_t[:] donors):
    """
    Orders nodes downstream to upstream for each base-level and add closed nodes.
    Updates upstream_ordered_nodes.
    
    Params
    ------
    base_level_and_closed_nodes: memoryview(long).
        Base level and closed nodes.
    upstream_ordered_nodes: memoryview(long).
        Stack of nodes ordered downstream to upstream, built from the recursive call
        of _add_to_upstream_ordered_nodes().
    donors_start_indexes: memoryview(long).
        Starting indexes in donors for each receiver. E.g. 3 for receiver n and 4 for
        receiver n+1 means that the donors for receiver n are in index 3 in donors
        (noted delta_i).
    donors: memoryview(long).
        Donors ordered by the id of their receivers (noted D).
    
    Return
    ------
    void
    """
    _calc_upstream_order_for_nodes_c(base_level_and_closed_nodes, 
        upstream_ordered_nodes, donors_start_indexes, donors)


def _calc_drainage_areas(cnp.int_t [:] downstream_ordered_nodes,
                         cnp.int_t [:] receivers, cnp.float_t [:] drainage_areas):
    """
    Calculates drainage areas for each node.
    Updates drainage_areas.
    Params
    ------
    downstream_ordered_nodes: memoryview(long).
        Nodes ordered upstream to downstream.
    receivers: memoryview(long).
        receivers for each donor, this latter ordered by node id.
    drainage_areas: memoryview(long).
        drainage areas of each node (ordered by node id)
    """
    cdef cnp.int_t n = len(downstream_ordered_nodes), i, donor_id, receiver_id

    for i in range(n):
        donor_id = downstream_ordered_nodes[i]
        receiver_id = receivers[donor_id];
        if receiver_id != donor_id:
            drainage_areas[receiver_id] += drainage_areas[donor_id]

