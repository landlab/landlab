# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 17:50:09 2021

@author: sahrendt
"""

import numpy as np

from ..components import ChannelProfiler, FlowAccumulator, DepressionFinderAndRouter
from ..graph import NetworkGraph
from .raster import RasterModelGrid
from .network import NetworkModelGrid

def create_network_from_raster(
        rmg, min_channel_thresh=10000, outlet_nodes=None, method='variable',
        n_widths=20, a=9.68, b=0.32, d_node_spacing=None, fields=None):
    
    """Create a NetworkModelGrid from a RasterModelGrid. Default behavior
    spaces nodes a certain number of local channel widths apart across the
    network. If method='constant' is specified, the d_node_spacing value
    is used to space nodes a constant distance apart across the network.

    Parameters
    ----------
    rmg : RasterModelGrid object
        A raster grid used to create a network grid
    min_channel_thresh : float, optional
        Value to use for the minimum drainage area associated with a
        plotted channel segment from the ChannelProfiler. Default values 10000.
    outlet_nodes : Single int value in iterable form, optional
        Iterable containing the node ID of nodes to start the channel profiles 
        from in ChannelProfiler. If not provided, the default is the node ID on 
        the model grid boundary with the largest terminal drainage area.
    method : string, 'variable' or 'constant'
        Specifies node-spacing method. 'variable' will dynamically update
        node spacing depending on local channel width. 'constant' will use the
        specified 'node_spacing' value to space nodes evenly across the
        network.
    n_widths : float, optional
        Multiplier to define node spacing as a function of local channel width.
    a : float, optional
        Parameter to be used in the Frasson et al. 2019 (GRL) power
        relationship between drainage area and channel width W=a*A^b. Default
        is value is 9.68 from Frassen et. al 2019
    b : float, optional
        Parameter to be used in the Frasson et al. 2019 (GRL) power
        relationship between drainage area and channel width W=a*A^b. Default
        is value is 0.32 from Frassen et. al 2019        
    d_node_spacing : float, optional
        Distance value for a constant node spacing along channel segments.
        Must be provided if method is 'constant'.
    fields : iterable, optional
        .at_node fields to map from RasterModelGrid to NetworkModelGrid.
        Formatted as strings inside an iterable object
        
    Returns
    -------
    NetworkModelGrid object with .at_node['rmg_node_value'] attribute
    listing the RasterModelGrid node ids at each NetworkModelGrid node.

    """

    if 'drainage_area' not in rmg.at_node:
        
        # run flow accumulator for ChannelProfiler
        fa = FlowAccumulator(rmg, 
                             'topographic__elevation',
                             flow_director='D8',
                             depression_finder='DepressionFinderAndRouter')
        fa.run_one_step()
        
    #delinate channel
    profiler = ChannelProfiler(
        rmg,
        number_of_watersheds=1,
        minimum_channel_threshold=min_channel_thresh,
        outlet_nodes=outlet_nodes,
        main_channel_only=False,
    )
    profiler.run_one_step()

    #obtain watershed key (should only be one)
    wtrshd_key = [k for k in profiler.data_structure.keys()][0]
    # obtain keys for channel segments, keys are raster nodes formated as 
    # tuple for (start, end) for channel segment start and end locations
    channel_segment_keys = profiler.data_structure[wtrshd_key].keys()
    
    # IDENTIFY CHANNEL SEGMENT CONNECTIVITY -----------------------------------
    # obtain node ids for start and end of every channel segments
    seg_starts =[seg[0] for seg in profiler.data_structure[wtrshd_key].keys()] 
    seg_ends = [seg[1] for seg in profiler.data_structure[wtrshd_key].keys()]
    # identify channel connectivity and how to properly link nodes
    # at different channel junctions
    # code does this by identifying the key of the channel seg just downstream 
    # and connects first node of upstream channel seg to downstream channel seg
    for i, seg_key in enumerate(channel_segment_keys):
        #create empty list to store how segment is connected
        connectivity = []
        connectivity_key = None
        #extract a single segment from the profiler data structure
        seg_i = profiler.data_structure[wtrshd_key][seg_key]
        # ask wehther the start of a segment is the end of another segment
        # (i.e. is this going to be connected downstream?)
        if seg_key[0] in seg_ends:
            connectivity.append('connected downstream')
            #find first segment downstream that should connect to last node of
            # segment upstream
            connect_to_channel_idx = np.argmax(seg_key[0]==seg_ends)
            connectivity_key = (seg_starts[connect_to_channel_idx],
                                seg_ends[connect_to_channel_idx])
        # ask whether the end of segment is in start of another segment
        # (i.e. is this going to be connected upstream?) if it isn't going to 
        # be connected upstream, we will use this to prompt the code to take
        # the last node in the profiler datastructure as the end node. (We 
        # don't necessarily need this, if we are okay with some trimming in
        # channel headwaters)
        if seg_key[-1] in seg_starts:
            connectivity.append('connected upstream')
        # note: we do not collect a connectivity_key here, because we
        # connect segments by connecting downstream, not upstream (we don't 
        # need to do it twice)
        
        # lets add this to our segment structure so we can know when and how to
        # connect our segments as we build the grid
        seg_i['connectivity'] = connectivity
        seg_i['connectivity_key'] = connectivity_key
    
    node_xy = [] #empty list to store paired x,y locations of nmg nodes
    rmg_nodes = [] #empty list to store raster model grid node corresponding to each network model grid node    
    links = [] # empty list to store link connections between nodes
    
    # FUNCTION TO ADD LINKS----------------------------------------------------
    # this is called several times in loop below, hoping it
    # makes testing easier to test once as a function
    def add_link(rmg, all_nodes_xy, all_links, head_node_rmg_id,
                 tail_node_rmg_id):
        
        """Add link connections to existing list of NetworkModelGrid nodes
        based upon an upstream and downstream RasterModelGrid node id. Also
        checks whether (x, y) values for upstream and downstream nodes exist
        in list of node locations and adds them if necessary.
        
        Parameters
        ----------
        rmg : RasterModelGrid object
            The RasterModelGrid to which NetworkModelGrid nodes and links will
            be added.
        all_nodes_xy : list of tuples
            List where tuple values for node x and y locations formatted as
            [(x1,y1), (x2,y2)...] already exist or will be stored.
        all_links : list of tuples
            List where tuple values for NetworkModelGrid node ids of upstream
            and downstream nodes for each link already exists or will be
            stored. Formatted as [(id2, id1),(id3, id2)...] where id# '#'
            corresponds to the index of the node entry in all_nodes_xy.
        head_node_rmg_id : int
            Value of the RasterModelGrid node id that corresponds to the
            desired head node of a NetworkModelGrid link.
        tail_node_rmg_id : int
            Value of the RasterModelGrid node id that corresponds to the
            desired head node of a NetworkModelGrid link.

        Returns
        -------
        None.

        """
        # define head node xy value by calling id from raster model grid
        head_node_xy = (rmg.x_of_node[head_node_rmg_id],
                        rmg.y_of_node[head_node_rmg_id])
        # define a tail node xy value by calling id from raster model grid
        tail_node_xy = (rmg.x_of_node[tail_node_rmg_id],
                        rmg.y_of_node[tail_node_rmg_id])
        # if these nodes don't already exist in the array of node xy vals from
        # another channel segment, add them
        if head_node_xy not in all_nodes_xy:
            all_nodes_xy.append(head_node_xy)
        if tail_node_xy not in all_nodes_xy:
            all_nodes_xy.append(tail_node_xy)
        # get the index of the head and tail node from our node_xy list
        # this is important to do in case they were added from a previous
        # channel segment. we need to ensure the order of network nodes is 
        # correct
        head_node__nmg_id = all_nodes_xy.index(head_node_xy)
        tail_node__nmg_id = all_nodes_xy.index(tail_node_xy)
        # append the head and tail network node ids to the link array
        # this if statement should be sufficient since we only connect links
        # one-way (i.e. downstream) between segments, but if there are still
        # issues with duplicate links something to add would be an additional
        # check for the opposite link combo: 
        # if (tail_node__nmg_id, head_node__nmg_id) not in all_links
        if (head_node__nmg_id, tail_node__nmg_id) not in all_links:
            all_links.append((head_node__nmg_id, tail_node__nmg_id))   
        
        
    
    
    # CREATE NETWORK MODEL GRID NODES & LINKS----------------------------------
    # loop over all channel segments and add network model nodes that correspond
    # to a certain cell in the raster model grid
    for i, seg_key in enumerate(channel_segment_keys):
        
        # access data of channel segments
        seg_i = profiler.data_structure[wtrshd_key][seg_key]
        
        # create list to hold rmg node ids where nmg nodes are located
        nmg_nodes = []
        
        # identify rmg value of first node in segment
        # first nodes of channel segments will be included in network 
        # if they don't already exist
        idx_node = 0
        rmg_node = seg_i['ids'][idx_node]
        
        # if seg_i is connected dowstream, add link connecting first node of
        # segment to downstream node
        if seg_i['connectivity_key'] is not None:
            channel_key = seg_i['connectivity_key']
            connecting_seg = profiler.data_structure[wtrshd_key][channel_key]
            
            # check to make sure there are nmg nodes on downstream segment
            # (there might not be any if the downstream segment is shorter than
            # what the user specifies for the desired network grid node spacing)
            if len(connecting_seg['ids_nmg']) > 0:
                connect_node = connecting_seg['ids_nmg'][-1]
            
            # if there are no nmg nodes on the downstream segment
            # it must be too short for calculated node spacing
            # if this is the case, connect upstream segment to first node in 
            # dowsntream connecting seg
            else:
                connect_node = connecting_seg['ids'][0]
            #add a link for this connection if necessary
            add_link(rmg, node_xy, links,
                     head_node_rmg_id=rmg_node,
                     tail_node_rmg_id=connect_node)
        
        # iterate over segment adding new nodes as long as there are upstream nodes
        # that can be placed on network model grid based upon node spacing
        upstrm_node=True
        while upstrm_node is True:

            # if we haven't already stored the rmg id value for this new node
            # add it to our master list of rmg nodes and sub-list of nmg nodes
            if rmg_node not in rmg_nodes:
                rmg_nodes.append(rmg_node)
                nmg_nodes.append(rmg_node)
            
            # Assign node spacing as n_channel_widths or a constant value from
            # input params
            if method is 'variable':
                #calculate drainage area contributing to this node
                da_node = rmg.at_node['drainage_area'][rmg_node]
                #relate drainage area to river width (convert area to km, width in m)
                # from Frasson et al. 2019 GRL
                w_channel = (a*da_node/(1000**2))**b
                #calculate upstream node spacing, n_widths_defines stable node spacing
                node_spacing = n_widths*w_channel
            if method is 'constant':
                node_spacing = d_node_spacing
            
            # if stable node spacing is greater than raster grid resolution
            if node_spacing > rmg.dx:
                #optimal along-channel node location based upon node spacing
                opt_loc = seg_i['distances'][idx_node] + node_spacing
                # define tolerance to not add extra node if opt loc is within half  
                # a node spacing away from end of segment 
                buffer_tol = 0.5*node_spacing
                
                # if we can fit another node on the channel segment
                if opt_loc < (seg_i['distances'][-1] - buffer_tol):
                    # find id of node closest to this location
                    idx_next_node = np.abs(seg_i['distances'] - opt_loc).argmin()
                    # update rmg node with whatever this next node should be
                    rmg_next_node = seg_i['ids'][idx_next_node]
                    
                    # add link from this upstream node to the current node
                    # if necessary
                    add_link(rmg, node_xy, links,
                             head_node_rmg_id=rmg_next_node,
                             tail_node_rmg_id=rmg_node)
                    
                    # update idx_node and rmg node for next loop
                    rmg_node = rmg_next_node
                    idx_node = idx_next_node
                
                # if no more nodes can be placed on this segment, 
                # move to next segment
                else:
                    upstrm_node = False
                    # add last node in segment to list of node xys
                    last_node_xy = (rmg.x_of_node[rmg_node],
                                    rmg.y_of_node[rmg_node])
                    if last_node_xy not in node_xy:
                        node_xy.append(last_node_xy)
    
            # if no more nodes have stable locations on this segment
            # move to next segment
            else:
                upstrm_node = False
                
                # check if segment is connected upstream:
                if 'connected upstream' in seg_i['connectivity']:
                # if we are seeing links on main stem channels that are smaller
                # then raster model grid resolution, flag this as an error
                    raise ValueError(
                        'main stem link lengths are smaller than grid res.'\
                        'try increasing n_widths or changing a and b params')
                
                # add last node in segment to list of node xys
                last_node_xy = (rmg.x_of_node[rmg_node],
                                rmg.y_of_node[rmg_node])
                if last_node_xy not in node_xy:
                    node_xy.append(last_node_xy)
        
        # store location of network nodes as raster model ids in channel profiler
        # datastructure. this will be used for joining channel segments later
        seg_i['ids_nmg'] = np.array(nmg_nodes)
    
    # CREATE NETWORK MODEL GRID OBJECT-----------------------------------------
    x_of_node, y_of_node = zip(*node_xy)
    
    # Maintain sorting by creating an unsorted network graph and sorting.
    # This process is important to ensure that the fields are assigned to the
    # correct links.
    graph_net = NetworkGraph((y_of_node, x_of_node), links=links, sort=False)
    sorted_nodes, sorted_links, sorted_patches = graph_net.sort()
    
    # use the sorting information to make a new network model grid.
    nmg = NetworkModelGrid(
        (np.asarray(y_of_node)[sorted_nodes], np.asarray(x_of_node)[sorted_nodes]),
        np.vstack((graph_net.node_at_link_head, graph_net.node_at_link_tail)).T
    )
    
    #add RMG node locations and extra fields to network model grid from 
    # raster model grid
    nmg.at_node['rmg_node_value'] = np.array(rmg_nodes)[sorted_nodes]
    if fields is None:
        fields = []
    for field in fields:
        nmg.at_node[field] = rmg.at_node[field][nmg.at_node['rmg_node_value']]
        
    return(nmg)
    
