"""
lateral erosion component

needs:




"""

import numpy as np





def FortyfiveNode(donor, i, receiver, link_list, neighbors, diag_neigh):
    debug=0
    print_debug=0
    if (debug):
        print "flow from ", donor, " to ", receiver, " is 45 degrees"
    radcurv_angle=0.67
    if(print_debug):
        print "node is crossing diagonal"
    #OLD WAY: diagonal list goes [SE, SW, NW, NE]. Node list are ordered as [E,S,W,N]
    #LL 2017: diagonal list goes [NE, NW, SW, SE]. Node list are ordered as [E,N,W,S]				
    #if water flows SE-N OR if flow NE-S or E-NW or E-SW, erode west node
    if (donor==diag_neigh[0] and receiver==neighbors[3] or 
    donor==diag_neigh[3] and receiver==neighbors[1]
    or donor==neighbors[0] and receiver==diag_neigh[2] or
    donor==neighbors[0] and receiver==diag_neigh[1]):
        if(print_debug):
            print "flow SE-N or NE-S or E-NW or E-SW, erode west node"
        lat_node=neighbors[2]
        if(print_debug):
            print "lat_node", lat_node
    #if flow is from SW-N or NW-S or W-NE or W-SE, erode east node
    elif (donor==diag_neigh[1] and receiver==neighbors[3] or 
    donor==diag_neigh[2] and receiver==neighbors[1] or
    donor==neighbors[2] and receiver==diag_neigh[3] or
    donor==neighbors[2] and receiver==diag_neigh[0]):
        if(print_debug):
            print "flow from SW-N or NW-S or W-NE or W-SE, erode east node"
        lat_node=neighbors[0]
        if(print_debug):
            print "lat_node", lat_node
    #if flow is from SE-W or SW-E or S-NE or S-NW, erode north node
    elif (donor==diag_neigh[3] and receiver==neighbors[2] or 
    donor==diag_neigh[2] and receiver==neighbors[0] or 
    donor==neighbors[3] and receiver==diag_neigh[0] or
    donor==neighbors[3] and receiver==diag_neigh[1]):
        if(print_debug):
            print "flow from SE-W or SW-E or S-NE or S-NW, erode north node"
            
        lat_node=neighbors[1]
        if(print_debug):
            print "lat_node", lat_node
    #if flow is from NE-W OR NW-E or N-SE or N-SW, erode south node
    elif (donor==diag_neigh[0] and receiver==neighbors[2] or 
    donor==diag_neigh[1] and receiver==neighbors[0] or
    donor==neighbors[1] and receiver==diag_neigh[3] or
    donor==neighbors[1] and receiver==diag_neigh[2]):
        if(print_debug):
            print "flow from NE-W or NW-E or N-SE or N-SW, erode south node"
        lat_node=neighbors[3]
        if(print_debug):
            print "lat_node", lat_node

    return lat_node, radcurv_angle
