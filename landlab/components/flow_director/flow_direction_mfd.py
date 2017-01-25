#! /usr/env/python

"""
flow_direction_mfd.py: calculate multiple-flow-direction flow directions.

Works on both a regular or irregular grid. Also calculates flow proportions.

KRB Jan 2017
"""

import numpy as np

from landlab import BAD_INDEX_VALUE
from landlab.core.utils import as_id_array

UNDEFINED_INDEX = BAD_INDEX_VALUE


def flow_directions_mfd(elev, 
                        neighbors_at_node,
                        links_at_node,
                        active_link_dir_at_node,
                        tail_node, 
                        head_node, 
                        link_slope, 
                        baselevel_nodes=None,
                        partition_method='slope'):

    """
    Find multiple-flow-direction flow directions on a grid.

    Finds and returns flow directions and proportions for a given elevation 
    grid. Each node is assigned multiple flow directions, toward all of the N 
    neighboring nodes that are lower than it. If none of the neighboring nodes
    are lower, it is assigned to itself. Flow proportions can be calculated as
    proportional to slope (default) or proportional to the square root of 
    slope, which is the solution to a steady kinematic wave. 

    Parameters
    ----------
    elev : array_like
        Elevations at nodes.
    neighbors_at_node : array_like (num nodes, max neighbors at node)
        For each node, the link IDs of active links.
    links_at_node : array_like (num nodes, max neighbors at node)
    
    link_dir_at_node: array_like (num nodes, max neighbors at node)
    
    tail_node : array_like
        IDs of the tail node for each link.
    head_node : array_like
        IDs of the head node for each link.
    link_slope : array_like
        slope of each link, defined POSITIVE DOWNHILL (i.e., a negative value
        means the link runs uphill from the fromnode to the tonode).        
    baselevel_nodes : array_like, optional
        IDs of open boundary (baselevel) nodes.
    partition_method: string, optional
        Method for partitioning flow. Options include 'slope' (default) and 
        'square_root_of_slope'. 

    Returns
    -------
    receivers : ndarray of size (num nodes, max neighbors at node)
        For each node, the IDs of the nodes that receive its flow. For nodes
        that do not direct flow to all neighbors, BAD_INDEX_VALUE is given as
        a placeholder. The ID of the node itself is given if no other receiver
        is assigned.
    proportions : ndarray of size (num nodes, max neighbors at node)
        For each receiver, the proportion of flow (between 0 and 1) is given.
        A proportion of zero indicates that the link does not have flow along 
        it. 
    steepest_slope : ndarray
        The slope value (positive downhill) in the direction of flow.
    steepest_receiver : ndarray
        For each node, the node ID of the node connected by the steepest link. 
        BAD_INDEX_VALUE is given if no flow emmanates from the node. 
    sink : ndarray
        IDs of nodes that are flow sinks (they are their own receivers)
    receiver_links : ndarray of size (num nodes, max neighbors at node)
        ID of links that leads from each node to its receiver, or
        UNDEFINED_INDEX if no flow occurs on this link.
    steepest_link : ndarray
        For each node, the link ID of the steepest link. 
        BAD_INDEX_VALUE is given if no flow emmanates from the node.

    Examples
    --------
    The example below assigns elevations to the 10-node example network in
    Braun and Willett (2012), so that their original flow pattern should be
    re-created.

    >>> import numpy as np
    >>> from landlab.components.flow_routing import flow_directions
    >>> z = np.array([2.4, 1.0, 2.2, 3.0, 0.0, 1.1, 2.0, 2.3, 3.1, 3.2])
    >>> fn = np.array([1,4,4,0,1,2,5,1,5,6,7,7,8,6,3,3,2,0])
    >>> tn = np.array([4,5,7,1,2,5,6,5,7,7,8,9,9,8,8,6,3,3])
    >>> s = z[fn] - z[tn]  # slope with unit link length, positive downhill
    >>> active_links = np.arange(len(fn))
    >>> r, ss, snk, rl = flow_directions(z, active_links, fn, tn, s)
    >>> r
    array([1, 4, 1, 6, 4, 4, 5, 4, 6, 7])
    >>> ss
    array([ 1.4,  1. ,  1.2,  1. ,  0. ,  1.1,  0.9,  2.3,  1.1,  0.9])
    >>> snk
    array([4])
    >>> rl[3:8]
    array([15, -1,  1,  6,  2])

    """
    # Calculate the number of nodes.
    num_nodes = len(elev)

    # Create a node array
    node_id = np.arange(num_nodes)
    
    # Calculate the maximum number of neighbors at node. 
    max_number_of_neighbors = neighbors_at_node.shape[1]
    
    # Make a copy of neighbors_at_node so we can change it into the receiver
    # array. 
    receivers = neighbors_at_node.copy()  
    
    # Construct the array of slope to neighbors at node. This also will adjust
    # for the slope convention based on the direction of the link. 
    slopes_to_neighbors_at_node = -link_slope[links_at_node]*active_link_dir_at_node
    
    # Make a copy so this can be changed based on where no flow occurs. 
    receiver_links = links_at_node.copy()
    
    # some of these potential recievers may have already been assigned as 
    # UNDEFINED_INDEX because the link was inactive. Make a mask of these for
    # future use. Also find the close nodes. 
    inactive_link_to_neighbor = active_link_dir_at_node == 0 
    closed_nodes = np.sum(np.abs(active_link_dir_at_node), 1) == 0
    # Now calculate where flow occurs. 
    # First, make an elevation array of potential receivers. 
    potential_receiver_elev = elev[neighbors_at_node]
    
    # now make an array of the same shape (for direct comparison) of the source
    # node elevation.
    source_node_elev = elev[np.tile(node_id, (max_number_of_neighbors,1)).T]
    
    # find where flow does not occur (source is lower that receiver)
    flow_does_not_occur = source_node_elev<=potential_receiver_elev
    
    # Where the source is lower, set receivers to UNDEFINED_INDEX
    receivers[flow_does_not_occur] = UNDEFINED_INDEX
     
    # Where the link is not active, set receivers to UNDEFINED_INDEX
    receivers[inactive_link_to_neighbor] = UNDEFINED_INDEX
    
    # Next, find where a node drains to itself 
    drains_to_self = receivers.sum(1) == -1*max_number_of_neighbors
    
    # Where this occurs, set the receiver ID in the first column of receivers
    # to the node ID. 
    receivers[drains_to_self, 0] = node_id[drains_to_self]
    
    # Finally, set the first element of the closed nodes to themselves. 
    receivers[closed_nodes, 0] = node_id[closed_nodes]
    
    # Next, calculate flow proportions.     
    # Copy slope array and mask by where flow is not occuring and where the 
    # link is inactive. 
    flow_slopes = slopes_to_neighbors_at_node.copy()
    flow_slopes[flow_does_not_occur] = 0.
    flow_slopes[inactive_link_to_neighbor] = 0.          
             
    if partition_method == 'square_root_of_slope':
        values_for_partitioning = flow_slopes**0.5
    elif partition_method == 'slope':
        values_for_partitioning = flow_slopes
    else:
        raise ValueError ('Keyword argument to partition_method invalid.')
    
    # Calculate proportions by normalizing by rowsums. 
    denom = np.tile(values_for_partitioning.sum(1), (max_number_of_neighbors,1)).T
    denom[denom<=0] = 1  # to prevent runtime errors
    proportions = values_for_partitioning/denom                      
    proportions[drains_to_self, 0] = 1
    proportions[drains_to_self, 1:] = 0       
                                  
    # Might need to sort by proportions and rearrange to follow expectations 
    # of no UNDEFINED_INDEX value in first column. KRB NOT SURE
    
    # mask the receiver_links by where flow doesn't occur to return
    receiver_links[flow_does_not_occur] = UNDEFINED_INDEX
    receiver_links[inactive_link_to_neighbor] = UNDEFINED_INDEX
    
    # identify the steepest link so that the steepest receiver, link, and slope
    # can be returned. 
    slope_sort = np.argsort(flow_slopes, 1) == max_number_of_neighbors-11
    steepest_slope = flow_slopes[slope_sort]
    
    ## identify the steepest link and steepest receiever. 
    steepest_link = receiver_links[slope_sort]
    steepest_receiver = receivers[slope_sort]
    
    # Optionally, handle baselevel nodes: they are their own receivers
    if baselevel_nodes is not None:
        receivers[baselevel_nodes,:] = node_id[baselevel_nodes,]
        receiver_links[baselevel_nodes,:] = UNDEFINED_INDEX
        steepest_slope[baselevel_nodes,:] = 0.
    
    # The sink nodes are those that are their own receivers (this will normally
    # include boundary nodes as well as interior ones; "pits" would be sink
    # nodes that are also interior nodes).
    (sink, ) = np.where(node_id==receivers[:,0])
    sink = as_id_array(sink)
    
    return (receivers, proportions, steepest_slope, steepest_receiver, sink, 
            receiver_links, steepest_link)

if __name__ == '__main__':
    import doctest
    doctest.testmod()