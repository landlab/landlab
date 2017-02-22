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
from landlab.utils.decorators import use_field_name_or_array
from landlab import VoronoiDelaunayGrid  # for type tests

from landlab.grid.raster_gradients import (
         _calc_subtriangle_unit_normals_at_node,
         _calc_subtriangle_slopes_at_node,
         _calc_subtriangle_aspect_at_node)

@use_field_name_or_array('node')
def _return_surface(grid, surface):

    """
    Private function to return the surface to direct flow over.

    This function exists to take advantange of the 'use_field_name_or_array
    decorator which permits providing the surface as a field name or array.
    """
    return(surface)

def flow_directions_dinf(grid, 
                         elevs='topographic__elevation', 
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
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_director.flow_direction_dinf import(
    ...                                                  flow_directions_dinf)
    >>> grid = RasterModelGrid((3,3), spacing=(1, 1))
    >>> _ = grid.add_field('topographic__elevation', 
    ...                     grid.node_x*grid.node_y, 
    ...                     at = 'node')
    >>> (receivers, proportions, 
    ... steepest_slope, steepest_receiver, 
    ... sink, receiver_links, steepest_link) = flow_directions_dinf(grid)
    >>> receivers
    array([[ 0, -1],
           [ 1, -1],
           [ 2, -1],
           [ 3, -1],
           [ 5,  8],
           [ 1,  2],
           [ 6, -1],
           [ 6,  3],
           [ 7,  4]])
    >>> proportions

    """
    # grid type testing
    if isinstance(grid, VoronoiDelaunayGrid):
        raise NotImplementedError('Dinfinity is currently implemented for'
                              ' Raster grids only')
    # first set some basic values
    # Calculate the number of nodes.
    num_nodes = len(elevs)
    
    # Create a node array
    node_id = np.arange(num_nodes)
    closed_node = grid.status_at_node == CLOSED_BOUNDARY
    
    # Set the number of receivers and facets.
    num_receivers = 2
    num_facets = 8
    
    # create an array of the triangle numbers
    tri_numbers = np.arange(num_facets)
    
    # create a arrays
    ac = np.array([0., 1., 1., 2., 2., 3., 3., 4.])
    af = np.array([1., -1., 1., -1., 1., -1., 1., -1.])
    
    # Initialize receiver and proportion arrays 
    receivers = UNDEFINED_INDEX * np.ones((num_nodes, num_receivers), dtype=int)
    proportions = np.zeros((num_nodes, num_receivers), dtype=float)
    receiver_links = UNDEFINED_INDEX * np.ones((num_nodes, num_receivers), dtype=int)
    slopes_to_receivers = np.zeros((num_nodes, num_receivers), dtype=float)
    
    # create list of triangle neighbors at node. Give math-orientation array first always. 
    # has shape, (nnodes, 8 triangles, 2 neighbors)
    n_at_node = grid.neighbors_at_node
    dn_at_node = grid._diagonal_neighbors_at_node
    if isinstance(grid, VoronoiDelaunayGrid):
        pass    # construct TNAN array here by circshifing. 
    else: 
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
    
    
    # Construct the array of slope to triangles at node. This also will adjust
    # for the slope convention based on the direction of the links. 
    # this is a (nnodes, 2, 8) array
    slopes_to_triangles_at_node = link_slope[triangle_links_at_node]*triangle_link_dirs_at_node
        
    # construct some arrays that deal with the distances between points on the
    # grid. 
                                            
    # construct d1 and d2, we know these because we know where the orthogonal 
    # links are 
    diag_length = ((grid.dx)**2+(grid.dy)**2)**0.5
                  
    # for irregular grids, d1 and d2 will need to be matricies              
    d1 = np.array([grid.dx, grid.dx , grid.dy, grid.dy, grid.dx, grid.dx, grid.dy, grid.dy])
    d2 = np.array([grid.dx , grid.dy, grid.dy, grid.dx, grid.dx, grid.dy, grid.dy, grid.dx])
    
    thresh = np.arctan(d2/d1)
    
    # begin the algorithm in earnest
    
    # construct e0, e1, e2 for all triangles at all nodes. 
    # will be (nnodes, nfacets=8 for raster or nfacets = max number of patches 
    # for irregular grids. 
    
    # e0 is origin point of the facet
    e0 = elevs[node_id]
    
    # e1 is the point on the side closest to east
    e1 = elevs[triangle_neighbors_at_node[:,0,:]]
    # e2 is the point on the side closest to west. 
    e2 = elevs[triangle_neighbors_at_node[:,1,:]]
    
    # mask out where nodes do not exits (e.g. triangle_neighbors_at_node == -1)
    e2[triangle_neighbors_at_node[:,1,:]==-1] = np.nan
    e1[triangle_neighbors_at_node[:,0,:]==-1] = np.nan
      
    # loop through and calculate s1 and s2
    # this will only loop nfacets times. 
    s1 = np.empty_like(e1)
    s2 = np.empty_like(e2)
    
    for i in range(num_facets):
        s1[:,i] = (e0 - e1[:,i])/d1[i]
        s2[:,i] = (e0 - e2[:,i])/d2[i]
    
    # calculate r and s, the direction and magnitude
    r = np.arctan2(s2, s1)
    s = ((s1**2)+(s2**2))**0.5
        
    # adjust r if it doesn't sit in the realm of (0, arctan(d2,d1))
    too_small = r<0
    radj = r.copy()
    radj[too_small] = 0
    s[too_small] = s1[too_small]
    
    # to consider two big, we need to look by trangle. 
    for i in range(num_facets):
        too_big = r[:,i]>thresh[i]
        radj[too_big, i] = thresh[i]
        s[too_big, i] = (e0[too_big] - e2[too_big,i])/diag_length
    
    # calculate the geospatial version of r based on radj
    rg = np.empty_like(r)
    for i in range(num_facets):
        rg[:, i] = (af[i]*radj[:, i]) + (ac[i]*np.pi/2.)
    
    # set slopes that are nan to zero
    s[np.isnan(s)] = 0     
     
     # sort slopes based on 
    steepest_sort = np.argsort(s)
    
    # determine the steepest triangle
    steepest_triangle = tri_numbers[steepest_sort[:,-1]]
    
    # initialize arrays for the steepest rg and steepest s
    steepest_rg = np.empty_like(node_id, dtype = float)
    steepest_s = np.empty_like(node_id, dtype = float)
    
    for n in node_id:
        steepest_rg[n] = rg[n, steepest_sort[n,-1]]
        steepest_s[n] = s[n, steepest_sort[n,-1]]
        receivers[n,:] = triangle_neighbors_at_node[n,:,steepest_sort[n,-1]]
        receiver_links[n,:] = triangle_links_at_node[n,:,steepest_sort[n,-1]]
        slopes_to_receivers[n,:] = slopes_to_triangles_at_node[n,:,steepest_sort[n,-1]]
    # construct the baseline for proportions
    rg_baseline = np.arange(0, 2, 0.25)*np.pi
    
    # calculate alpha1 and alpha 2
    alpha2 = steepest_rg-rg_baseline[steepest_triangle]
    alpha1 = thresh[steepest_triangle] - alpha2
           
    # calculate proportions from alpha               
    proportions[:,0] = (alpha2)/(alpha1+alpha2)
    proportions[:,1] = (alpha1)/(alpha1+alpha2)
    
    
    drains_to_self = np.isnan(proportions[:,0])
    drains_to_self[steepest_s<=0]=True
    receivers[drains_to_self,0]=node_id[drains_to_self]
    receivers[drains_to_self,1]=-1
    
    proportions[drains_to_self,0]=1.
    proportions[drains_to_self,1]=0.
    
    # mask the receiver_links by where flow doesn't occur to return
    receiver_links[drains_to_self, :] = UNDEFINED_INDEX
    
    # identify the steepest link so that the steepest receiver, link, and slope
    # can be returned. 
    slope_sort = np.argsort(np.argsort(slopes_to_receivers, 
                                       axis=1), 
                            axis=1) == (num_receivers-1)
    steepest_slope = slopes_to_receivers[slope_sort]
    steepest_slope[drains_to_self] = 0.
                  
    ## identify the steepest link and steepest receiever. 
    steepest_link = receiver_links[slope_sort]
    steepest_link[drains_to_self] = UNDEFINED_INDEX
                 
    steepest_receiver = receivers[slope_sort]           
    steepest_receiver[drains_to_self] = UNDEFINED_INDEX
    
               
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
                