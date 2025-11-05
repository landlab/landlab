import numpy as np
import pandas as pd



def _dist_func(x0, x1, y0, y1):
    return np.hypot(x0 - x1, y0 - y1)



def _remove_small_tribs(rmg_nodes_to_nmg_links_mapper, nmg_link_to_rmg_coincident_nodes_mapper):
    """remove channel rmg nodes mapped to link that represent first order channels
    not represented by the network model grid"""
    for link in np.unique(nmg_link_to_rmg_coincident_nodes_mapper['linkID'].values):
        # first get the coincident rmg node that represents the inlet to the link
        mask1 = nmg_link_to_rmg_coincident_nodes_mapper['linkID'] == link 
        # coincident nodes listed in nmg_link_to_rmg_coincident_nodes_mapper are 
        # ordered from outlet to inlet, so last node is coincident with inlet
        # min_area_node = nmg_link_to_rmg_coincident_nodes_mapper['coincident_node'][mask1].iloc[-1] 
        # node with shortest downstream distance from inlet is the inlet node
        min_dist = nmg_link_to_rmg_coincident_nodes_mapper['dist'][mask1].min()
        mask2 = nmg_link_to_rmg_coincident_nodes_mapper['dist'] == min_dist
        min_area_node = nmg_link_to_rmg_coincident_nodes_mapper['coincident_node'][mask1][mask2].iloc[0]
        # now get the contributing area of the rmg channel node mapped to the link 
        # inlet. 
        mask3 = rmg_nodes_to_nmg_links_mapper['coincident_node'] == min_area_node 
        min_area = rmg_nodes_to_nmg_links_mapper['node_drainage_area'][mask3].min() # drainage area of channel inlet node
        # finally, find any nodes assigned to link that have drainage area 
        # less than the inlet node and remove
        mask4 = (rmg_nodes_to_nmg_links_mapper['linkID'] == link) & (rmg_nodes_to_nmg_links_mapper['node_drainage_area']<min_area)
        rmg_nodes_to_nmg_links_mapper = rmg_nodes_to_nmg_links_mapper.drop(rmg_nodes_to_nmg_links_mapper.index[mask4].values)     
    return rmg_nodes_to_nmg_links_mapper



def map_rmg_nodes_to_nmg_links(grid, nmg_link_to_rmg_coincident_nodes_mapper, rmg_nodes, remove_small_tribs = True):
    """Map the nodes representing the channel location in a DEM to the closest
    network model grid location. Network model grid location is described in 
    terms of link id and distance down link, measured from the inlet (tail) node
    of the link.
    
    Parameters
    ----------
    grid : raster model grid
        needs to have node field "drainage_area"
    nmg_link_to_rmg_coincident_nodes_mapper : pandas dataframe
        each row of the dataframe lists the link ID, the coincident node ID, the 
        x and y coordinates and the downstream distance of the coincident node 
        and the drainage area of the link
    rmg_nodes : np.array
        an array of node ids to be mapped to the nmg links
    remove_small_tribs : bool
        If True, first order channels that split from a much higher order channel
        are not matched to a Link. If False, the ndoes representin small, first 
        order channels will be mapped to the same link as the much larger channel 
        they drain into.

    Returns
    -------
    rmg_nodes_to_nmg_links_mapper : pandas dataframe
        each row of the dataframe lists the node ID, the link ID the node has been 
        mapped too, the closest nmg-link-coincident node ID, the drainage area
        of the link and the drainage area of the node

    """

    def dist_between_nmg_and_rmg_nodes(row, xc, yc):
        """distance between channel node and link node"""
        return _dist_func(xc,row['x'],yc,row['y'])
    
    rmg_node_link_id = []
    node_ = []
    dist_ = []
    link_ = []
    for n in rmg_nodes: # for each rmg node
        xc = grid.node_x[n]; yc = grid.node_y[n]
        # compute the distance to all link coincident rmg nodes
        dist = nmg_link_to_rmg_coincident_nodes_mapper.apply(lambda row: dist_between_nmg_and_rmg_nodes(row, xc, yc),axis=1)
        # pick closest coincident node and corresponding link
        # if more than one (which can happen because the confluence between two
        # links overlay the same node), pick link with largest contributing area
        mask = dist == dist.min()
        dist_min_links = nmg_link_to_rmg_coincident_nodes_mapper[['linkID','dist','coincident_node','drainage_area']][mask]
        link = dist_min_links[dist_min_links['drainage_area']==dist_min_links['drainage_area'].max()].head(1)
        link['node_drainage_area'] = grid.at_node['drainage_area'][n] # add node drainage area to attributes
        link_.append(link)
        
    rmg_nodes_to_nmg_links_mapper  = pd.concat(link_)
    rmg_nodes_to_nmg_links_mapper['node'] = rmg_nodes
    # organize column order in mapper
    rmg_nodes_to_nmg_links_mapper =rmg_nodes_to_nmg_links_mapper[['node','linkID','dist','coincident_node','drainage_area','node_drainage_area']].reset_index(drop=True) 

    if remove_small_tribs:# check for small tributary nodes assigned to link and remove them 
        rmg_nodes_to_nmg_links_mapper = _remove_small_tribs(rmg_nodes_to_nmg_links_mapper, 
                                                            nmg_link_to_rmg_coincident_nodes_mapper)
  
    return rmg_nodes_to_nmg_links_mapper 