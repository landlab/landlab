# distutils: language = c++
#cython: language_level=3

import numpy as np
cimport numpy as cnp
cimport cython
from libcpp cimport bool
from libcpp.pair cimport pair

# 1-3. Instantiate the Queues (Steps #1-3 and #5-6 of Barnes, 2014 see Python FlowRouter)
#########################################################################################
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


@cython.boundscheck(False) # turn off bounds-checking for entire function    
cdef void _init_flow_direction_queues(cnp.int_t nodes_n, const cnp.int_t [:] base_level_nodes, const cnp.int_t [:] closed_nodes, 
    cnp.float_t [:] z, _priority_queue& to_do, cnp.int_t [:] receivers, 
    cnp.int_t [:] outlet_nodes, cnp.int_t [:] done, cnp.int_t* done_n_ptr) nogil: 
        # important to pass the priority_queue to_do as a reference, otherwise it won't be modified
    # &: to_do is passed as a reference
    # I didn't manage to use the bint type for the variable done, so here 0 is for False and 1 for True
    cdef:
        cnp.int_t node_id, i, n = len(base_level_nodes), m = len(closed_nodes)
        pair[cnp.int_t, cnp.float_t] node_pair
            
    for i in range(n): # for node_i in open_boundary: raise a compiling error with nogil
        node_id = base_level_nodes[i]
        node_pair = pair[cnp.int_t, cnp.float_t](node_id, z[node_id])                   
        to_do.push(node_pair)
        receivers[node_id] = node_id
        outlet_nodes[node_id] = node_id
        done[node_id] = 1; done_n_ptr[0] +=1
        
    for i in range(m):
        node_id = closed_nodes[i]
        receivers[node_id] = node_id
        outlet_nodes[node_id] = node_id
        done[node_id] = 1; done_n_ptr[0] +=1 # done_n_ptr[0]: method to dereference the pointer. 
        #cython.operator.dereference doesn't seem to work on this pointer to a long-type variable

# 4. Functions necessary for the flow direction processing (Steps #11 - 20)
###########################################################################

@cython.boundscheck(False)
cdef void _set_flooded_and_outlet(cnp.int_t donor_id, cnp.float_t [:] z, cnp.int_t [:] receivers, 
    cnp.int_t [:] outlet_nodes, cnp.int_t [:] depression_outlet_nodes, cnp.int_t [:] flooded_nodes, 
    cnp.float_t [:] depression_depths, cnp.int_t flooded_status, cnp.int_t bad_index) nogil:
    """ Modify the flooded status (flooded_nodes), the outlet nodes and depths of the depressions (outlet_nodes,
    depression_depths) for the node donor_id depending on its surface z value and the one of 
    its receiver outlet and neighbors.
    Parameters:
    -----------
    donor_id: int
        Id of the node giving the flow
    z: nd.array(float)
        Surface where the flow is directed
    receivers: nd.array(int)
        Ids of the receiver nodes ordered by their donor id
    outlet_nodes: nd.Array(int)
        Id of the outlet node, ordered by the donor id, if the donor is in a depression
    flooded_nodes: nd.array(bool)
        Flooded status for all nodes
    depression_depths: nd.array(float)
        Depths of the depression (if existing) below flow level for all nodes
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    """
    cdef: 
        cnp.int_t receiver_id = receivers[donor_id]
        cnp.int_t receiver_depression_outlet = receiver_id if depression_outlet_nodes[receiver_id] == bad_index \
            else depression_outlet_nodes[receiver_id]
        cnp.int_t receiver_outlet = receiver_id if outlet_nodes[receiver_id] == bad_index \
            else outlet_nodes[receiver_id]
    outlet_nodes[donor_id] = receiver_outlet
    if z[donor_id] < z[receiver_depression_outlet]:
        depression_outlet_nodes[donor_id] = receiver_depression_outlet
        flooded_nodes[donor_id] = flooded_status        
        depression_depths[donor_id] = z[receiver_depression_outlet] - z[donor_id]                

@cython.boundscheck(False)
cdef void _set_receiver(cnp.int_t donor_id, cnp.int_t receiver_id, cnp.int_t [:] receivers, cnp.int_t [:] done, cnp.int_t* done_n_ptr) nogil:
    """ Modify the receiver (receivers) and the process statuses (done) of the donor node donor_id
    Parameters:
    -----------
    donor_id: int
        Id of the donor node
    receiver_id: int
        Id of the receiver node
    receivers: nd.array(int)
        Ids of the receivers for all donor nodes, ordered by the id of the donor nodes
    done: nd.array(bool)
        Process statuses for all nodes
    """
    cdef cnp.int_t done_n
    
    receivers[donor_id] = receiver_id
    done[donor_id] = 1 # 1 For True
    done_n_ptr[0] += 1

@cython.boundscheck(False)
cdef void _set_donor_properties(cnp.int_t donor_id, cnp.int_t receiver_id, cnp.int_t [:] sorted_pseudo_tails, 
    const cnp.int_t [:,:] head_start_end_indexes, const cnp.int_t [:] sorted_dupli_links, 
    cnp.float_t [:] sorted_dupli_gradients, cnp.float_t [:] z, cnp.float_t [:] steepest_slopes, 
    cnp.int_t [:] links_to_receivers) nogil:
    """ Update the steepest_slopes and the links_to_receivers of a donor in function of the
    slopes with its neighbors and the head-tail links of the grid. Two tests are done whether
    donor is head and receiver is tail, or the reverse.
    Parameters:
    -----------
    donor_id: int
        Id of the donor node
    receiver_id: int
        Id of the receiver node
    sorted_pseudo_tails: nd.array(int)
        All tails and heads for all existing links in the grid, ordered by 
        increasing head and tail ids in sorted_pseudo_heads
    head_start_end_indexes: nd.array(intxint)
        Start (array[0, :] and end (array[1, :]) inde of all possible value nodes in the sorted_pseudo_heads 
            (ordered by increasing head and tail ids
    sorted_dupli_links: nd.array(int)
        All existing links in the grid, duplicated, ordered by 
        increasing head and tail ids in sorted_pseudo_heads
    sorted_dupli_gradients: nd.array(float)
        Gradients of the surface for all links duplicated, ordered by 
        increasing head and tail ids in sorted_pseudo_heads
    z: nd.array(int)
        Surface where the flow is directed             
    steepest_slopes: nd.array(float)
        Steepest slope for all nodes
    links_to_receivers: nd.array(int)
        Ids of the links between a donor and a receiver, ordered by the ids of the donor nodes
    """

    # range of indexes where donor_id is founded in sorted_pseudo_heads
    cdef:
        cnp.int_t idx1 = head_start_end_indexes[0, donor_id]
        cnp.int_t idx2 = head_start_end_indexes[1, donor_id] + 1 
    
    # loop to bypass the impossibility to use sorted_pseudo_tails[idx1:idx2] == receiver_id with memoryviews
        cnp.int_t [:] s = sorted_pseudo_tails[idx1:idx2]
        cnp.int_t n = len(s), c = -1, i
    
    for i in range(n):
        if s[i] == receiver_id: 
            c = idx1 + i
            break
    steepest_slopes[donor_id] = sorted_dupli_gradients[c] if z[receiver_id] <= z[donor_id] else 0
    links_to_receivers[donor_id] = sorted_dupli_links[c]
    
#################################################################################################################

# Main functions to direct flow
###############################

@cython.boundscheck(False)
cdef void _direct_flow_c(cnp.int_t nodes_n, const cnp.int_t[:] base_level_nodes, const cnp.int_t[:] closed_nodes, 
                cnp.int_t[:] sorted_pseudo_tails, cnp.float_t[:] sorted_dupli_gradients, 
                const cnp.int_t[:] sorted_dupli_links, const cnp.int_t[:, :] head_start_end_indexes, 
                cnp.int_t [:] outlet_nodes, cnp.int_t [:] depression_outlet_nodes, 
                cnp.int_t[:] flooded_nodes, cnp.float_t[:] depression_depths, cnp.int_t[:] links_to_receivers, 
                cnp.int_t[:] receivers, 
                cnp.float_t[:] steepest_slopes, cnp.float_t[:] z, cnp.int_t flooded_status, cnp.int_t bad_index,
                cnp.int_t neighbors_max_number):
    cdef:
        cnp.int_t [:] done
        cnp.int_t [:] tmp_neighbors
        cnp.int_t [:] neighbors_to_do
        _priority_queue to_do = _priority_queue(_compare_second)
        cnp.int_t receiver_id, donor_id, n, i, j, done_n
        cnp.int_t [:] neighbors
        pair[cnp.int_t, cnp.float_t] node_pair

    done = np.full(nodes_n, 0, dtype=int)
    tmp_neighbors = np.full(neighbors_max_number, 0, dtype=int) 
    neighbors_to_do = np.array([], dtype=int)
    done_n = 0 # done_n is input only for MULTITHREADING and is not checked in this function
    
    _init_flow_direction_queues(nodes_n, base_level_nodes, closed_nodes, z, to_do, receivers, 
                            outlet_nodes, done, &done_n)

    for j in range(nodes_n):
        if to_do.empty() == True:
            break          
        receiver_id = to_do.top().first; to_do.pop(); done[0] = 1 
        
        #neighbors = _get_neighbors(receiver_id, sorted_pseudo_tails, head_start_end_indexes)    
       # neighbors_to_do = _get_neighbors_to_do(neighbors, done, tmp_neighbors)
       # n = len(neighbors_to_do)
        idx1 = head_start_end_indexes[0, receiver_id]; idx2 = head_start_end_indexes[1, receiver_id] + 1
        for i in range(idx1, idx2):
            donor_id = sorted_pseudo_tails[i]; 
            if done[donor_id] == 1: continue # 0 for False
        #for i in range(n): #, nogil=True): # len(neighbors_to_do) within the prange instruction triggers the stop of the kernel when execution
            #donor_id = neighbors_to_do[i]
            _set_receiver(donor_id, receiver_id, receivers, done, &done_n)     
            _set_flooded_and_outlet(donor_id, z, receivers, outlet_nodes, \
                depression_outlet_nodes, flooded_nodes, depression_depths, flooded_status,
                bad_index)
            _set_donor_properties(donor_id, receiver_id, sorted_pseudo_tails, head_start_end_indexes,
                        sorted_dupli_links, sorted_dupli_gradients, z, steepest_slopes, links_to_receivers)
            node_pair = pair[cnp.int_t, cnp.float_t](donor_id, z[donor_id])
            to_do.push(node_pair)

cpdef void _direct_flow(cnp.int_t nodes_n, const cnp.int_t[:] base_level_nodes, const cnp.int_t[:] closed_nodes, 
                cnp.int_t[:] sorted_pseudo_tails, cnp.float_t[:] sorted_dupli_gradients, 
                const cnp.int_t[:] sorted_dupli_links, const cnp.int_t[:, :] head_start_end_indexes, 
                cnp.int_t[:] outlet_nodes, cnp.int_t[:] depression_outlet_nodes, cnp.int_t[:] flooded_nodes, 
                cnp.float_t[:] depression_depths, cnp.int_t[:] links_to_receivers, cnp.int_t[:] receivers, 
                cnp.float_t[:] steepest_slopes, cnp.float_t[:] z, cnp.int_t flooded_status, cnp.int_t bad_index,
                cnp.int_t neighbors_max_number):
    _direct_flow_c(nodes_n, base_level_nodes, closed_nodes, 
                    sorted_pseudo_tails, sorted_dupli_gradients, 
                    sorted_dupli_links, head_start_end_indexes, 
                    outlet_nodes, depression_outlet_nodes,
                    flooded_nodes, depression_depths, links_to_receivers, receivers, 
                    steepest_slopes, z, flooded_status, bad_index,
                    neighbors_max_number)