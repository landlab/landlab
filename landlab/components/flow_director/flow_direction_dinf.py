# def direct_dinf(grid, elevs='topographic_elevation', baselevel_nodes=None):


"""
flow_direction_mfd.py: calculate multiple-flow-direction flow directions.

Works on both a regular or irregular grid. Also calculates flow proportions.

KRB Jan 2017
"""

import numpy as np
from landlab.core.utils import as_id_array
from landlab import BAD_INDEX_VALUE
UNDEFINED_INDEX = BAD_INDEX_VALUE

from landlab.grid.raster_gradients import (
         _calc_subtriangle_unit_normals_at_node,
         _calc_subtriangle_slopes_at_node,
         _calc_subtriangle_aspect_at_node)


def flow_directions_dinf(grid, 
                         elevs='topographic_elevation', 
                         baselevel_nodes=None):

    """
    Find Dinfinity-flow-direction flow directions on a grid.

    Finds and returns flow directions and proportions for a given elevation 
    grid by the D infinity method (Tarbotten, 1997). Each node is assigned two 
    flow directions, toward the two neighboring nodes that are on the steepest
    subtriangle. Partitioning of flow is done based on the aspect of the 
    subtriangle. 

    Parameters
    ----------
    grid : ModelGrid
        A grid of type Voroni.
    elevs : field name at node or array of length node
        The surface to direct flow across.       
    baselevel_nodes : array_like, optional
        IDs of open boundary (baselevel) nodes.

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

    >>> import numpy as np
    >>> from landlab.components.flow_routing import flow_directions_dinf

    """

    unit_normals = _calc_subtriangle_unit_normals_at_node(grid, 
                                                          elevs=elevs)
    slopes = np.transpose(_calc_subtriangle_slopes_at_node(grid, 
                                                           elevs=elevs, 
                                                           subtriangle_unit_normals=unit_normals))
    aspects = np.transpose(_calc_subtriangle_aspect_at_node(grid, 
                                                            elevs=elevs, 
                                                            subtriangle_unit_normals=unit_normals))
    unit_normals = np.swapaxes(unit_normals, 0,1)
    
    # next, collect information about each triangle. 
    
    # create list of triangle neighbors at node. Give math-orientation array first always. 
    # has shape, (nnodes, 8 triangles, 2 neighbors)
    n_at_node = grid.neighbors_at_node
    dn_at_node = grid._diagonal_neighbors_at_node
    triangle_neighbors_at_node = np.stack([np.vstack((n_at_node[:,0], dn_at_node[:,0])),
                                           np.vstack((dn_at_node[:,0], n_at_node[:,1])),
                                           np.vstack((n_at_node[:,1], dn_at_node[:,1])),
                                           np.vstack((dn_at_node[:,1], n_at_node[:,2])),
                                           np.vstack((n_at_node[:,2], dn_at_node[:,2])),
                                           np.vstack((dn_at_node[:,2], n_at_node[:,3])),
                                           np.vstack((n_at_node[:,3], dn_at_node[:,3])),
                                           np.vstack((dn_at_node[:,3], n_at_node[:,0]))],
                                          axis=-1)
    triangle_neighbors_at_node = triangle_neighbors_at_node.swapaxes(0,1)
    
    # next create, triangle links at node
    l_at_node = grid.links_at_node
    dl_at_node = grid._diagonal_links_at_node
    triangle_links_at_node = np.stack([np.vstack((l_at_node[:,0], dl_at_node[:,0])),
                                       np.vstack((dl_at_node[:,0], l_at_node[:,1])),
                                       np.vstack((l_at_node[:,1], dl_at_node[:,1])),
                                       np.vstack((dl_at_node[:,1], l_at_node[:,2])),
                                       np.vstack((l_at_node[:,2], dl_at_node[:,2])),
                                       np.vstack((dl_at_node[:,2], l_at_node[:,3])),
                                       np.vstack((l_at_node[:,3], dl_at_node[:,3])),
                                       np.vstack((dl_at_node[:,3], l_at_node[:,0]))],
                                      axis=-1)
    triangle_links_at_node = triangle_links_at_node.swapaxes(0,1)
    
    # next create link directions and active link directions at node
    # link directions
    ld_at_node = grid._link_dirs_at_node
    dld_at_node = grid._diag__link_dirs_at_node
    triangle_link_dirs_at_node = np.stack([np.vstack((ld_at_node[:,0], dld_at_node[:,0])),
                                           np.vstack((dld_at_node[:,0], ld_at_node[:,1])),
                                           np.vstack((ld_at_node[:,1], dld_at_node[:,1])),
                                           np.vstack((dld_at_node[:,1], ld_at_node[:,2])),
                                           np.vstack((ld_at_node[:,2], dld_at_node[:,2])),
                                           np.vstack((dld_at_node[:,2], ld_at_node[:,3])),
                                           np.vstack((ld_at_node[:,3], dld_at_node[:,3])),
                                           np.vstack((dld_at_node[:,3], ld_at_node[:,0]))],
                                          axis=-1)
    triangle_link_dirs_at_node = triangle_link_dirs_at_node.swapaxes(0,1)
    
    # active link directions. 
    ald_at_node = grid.active_link_dirs_at_node
    adld_at_node = grid._diag__active_link_dirs_at_node
    
    triangle_active_link_dirs_at_node = np.stack([np.vstack((ald_at_node[:,0], adld_at_node[:,0])),
                                                  np.vstack((adld_at_node[:,0], ald_at_node[:,1])),
                                                  np.vstack((ald_at_node[:,1], adld_at_node[:,1])),
                                                  np.vstack((adld_at_node[:,1], ald_at_node[:,2])),
                                                  np.vstack((ald_at_node[:,2], adld_at_node[:,2])),
                                                  np.vstack((adld_at_node[:,2], ald_at_node[:,3])),
                                                  np.vstack((ald_at_node[:,3], adld_at_node[:,3])),
                                                  np.vstack((adld_at_node[:,3], ald_at_node[:,0]))],
                                                 axis=-1)
    triangle_active_link_dirs_at_node = triangle_active_link_dirs_at_node.swapaxes(0,1)
    
    # need to create a list of diagonal links since it doesn't exist. 
    diag_links = np.sort(np.unique(grid._diag_links_at_node))
    diag_links = diag_links[diag_links>0]
    
    # calculate graidents across diagonals and orthogonals
    diag_grads = grid._calculate_gradients_at_d8_links(elevs)
    ortho_grads = grid.calc_grad_at_link(elevs)
          
    # finally compile link slopes               
    link_slope = np.hstack((np.arctan(ortho_grads),
                            np.arctan(diag_grads)))
    
    # Calculate the number of nodes.
    num_nodes = len(elevs)
    
    # Create a node array
    node_id = np.arange(num_nodes)
    
    # Set the number of receivers.
    num_receivers = 2
    
    # Initialize receiver and proportion arrays 
    receivers = UNDEFINED_INDEX * np.ones((num_nodes, num_receivers), dtype=int)
    proportions = np.zeros((num_nodes, num_receivers), dtype=float)
    receiver_links = UNDEFINED_INDEX * np.ones((num_nodes, num_receivers), dtype=int)
    receiver_slopes = np.zeros((num_nodes, num_receivers), dtype=float)
    # Construct the array of slope to triangles at node. This also will adjust
    # for the slope convention based on the direction of the links. 
    # this is a (nnodes, 2, 8) array
    slopes_to_triangles_at_node = link_slope[triangle_links_at_node]*triangle_link_dirs_at_node
    
    # determine which triangles have flow going out of the node:
    # these are triangles that have at least one positive slopes 
    # going out of them.
    flow_out_of_node = (np.sum(slopes_to_triangles_at_node>0, axis=1) >= 1)*1                                     
              
    # mask out all possible receivers by those that are
    # closed nodes. 
    #closed_nodes = (triangle_active_link_dirs_at_node==0)
    potential_receiver = triangle_neighbors_at_node
    #potential_receiver[closed_nodes] = -1
    
    # find which triangles have one or more closed nodes attached to them. 
    #flows_to_closed = np.any(potential_receiver == -1, axis=1)                  
                                  
    # find which triangle is steepest
    flow_slopes = slopes * flow_out_of_node
    flow_slopes[np.isnan(flow_slopes)] = 0
    #flow_slopes[flows_to_closed] = 0
               
    # need to choose the triangle with the steepest slope, 
    # and (if there are ties) steepest link.
    tri_numbers = np.arange(8)
    steepest_triangle = np.empty((num_nodes))
    steepest_slope_out_of_triangle = np.max(slopes_to_triangles_at_node, axis=1) * flow_out_of_node
    sum_of_slopes_out_of_triangle = np.sum(slopes_to_triangles_at_node, axis=1) * flow_out_of_node                                       
    #steepest_slope_out_of_triangle[flows_to_closed] = 0
    steepest_triangle = -1 * np.ones((num_nodes), dtype=int) 
                                
    for i in range(num_nodes):         
        steepest = flow_slopes[i,:] == np.max(flow_slopes[i,:])
        if np.sum(steepest)==1:
            steepest_triangle[i] = tri_numbers[steepest]
        else:
            # consider the links.
            # first make sure that at least one of them is going downhill
            if np.all(steepest_slope_out_of_triangle[i,:]<=0):
                steepest_triangle[i] = UNDEFINED_INDEX
            # otherwise, choose the steepest of both. 
            else:
                steepest_links = steepest_slope_out_of_triangle[i,:] == np.max(steepest_slope_out_of_triangle[i,:])
                steepest = steepest_links*steepest
                if np.sum(steepest)==1:
                    steepest_triangle[i] = tri_numbers[steepest]
                else:
                    steepest_sum = sum_of_slopes_out_of_triangle[i,:] == np.max(sum_of_slopes_out_of_triangle[i,:])
       
                    steepest_triangle[i] = np.argmax(steepest_sum*steepest)
            
    # initialize an array to hold the aspect of the steepest triangle
    steepest_aspect = np.empty((num_nodes))
                   
    ### flow slopes must be greater than zero. 
       
    # loop through (can't figure a better way) and find the steepest receievers
    # receiver links, slopes, and aspect. 
    for i in range(num_nodes):
        receivers[i,:] = potential_receiver[i,:,steepest_triangle[i]]
        steepest_aspect[i] = aspects[i,steepest_triangle[i]]
        receiver_links[i, :] = triangle_links_at_node[i,:,steepest_triangle[i]]
        receiver_slopes[i, :] = slopes_to_triangles_at_node[i,:,steepest_triangle[i]]
        
    # for those nodew that flow to themselves, set their receiver id
    # array correctly, and their proportions. 
    flow_to_self = (np.any(receivers == UNDEFINED_INDEX, axis = 1)) + (
                    flow_out_of_node.sum(axis=1) == 0) + (
                    steepest_triangle == UNDEFINED_INDEX)
    
    receivers[flow_to_self, 0] = node_id[flow_to_self]
    receivers[flow_to_self, 1] = -1
    proportions[flow_to_self, 0] = 1
               
    # next determine the flow partitioning. To do this we need the orientation of 
    # the clockwise edge of that triangle if it were flowing INTO the node.  
    orientation_of_cw_edges = np.array([45.,  0., 315., 
                                        270., 225., 180., 135, 90.])
    orientation_of_cw_edge_steepest = np.empty_like(receivers)
    orientation_of_cw_edge_steepest = orientation_of_cw_edges[steepest_triangle]
    
    proportions[flow_to_self==False,0] = (steepest_aspect[flow_to_self==False] - orientation_of_cw_edge_steepest[flow_to_self==False])/45.
    proportions[flow_to_self==False,1] = 1. - proportions[flow_to_self==False,0]
                          
    # by this algorithm, where proportions are greater than on or less than one this 
    # is an indication that all the flow should be along one of the two links. 
    first_col_leq1 = proportions[:,0]<0
    proportions[first_col_leq1,0] = 0                            
    proportions[first_col_leq1,1] = 1
    
    second_col_leq1 = proportions[:,1]<0
    proportions[second_col_leq1,0] = 1                            
    proportions[second_col_leq1,1] = 0           
                    
           
    # mask the receiver_links by where flow doesn't occur to return
    receiver_links[flow_to_self] = UNDEFINED_INDEX
    
    # identify the steepest link so that the steepest receiver, link, and slope
    # can be returned. 
    slope_sort = np.argsort(np.argsort(receiver_slopes, 
                                       axis=1), 
                            axis=1) == (num_receivers-1)
    steepest_slope = receiver_slopes[slope_sort]
    steepest_slope[flow_to_self] = 0.
                  
    ## identify the steepest link and steepest receiever. 
    steepest_link = receiver_links[slope_sort]
    steepest_link[flow_to_self] = UNDEFINED_INDEX
                 
    steepest_receiver = receivers[slope_sort]           
    steepest_receiver[flow_to_self] = UNDEFINED_INDEX
    
               
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
                