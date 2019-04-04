"""
lateral erosion component

needs:
fraction of bed exposed: fe=(1-qsqt)
percent lateral erosion of total erosion : Plat=qsqt
lateral erosion: Elat=Plat*E  (E is total erosion, dzdt)
vertical erosion: Evert=(1-Elat)*E



"""

import numpy as np


def NinetyNode(donor, i, receiver, link_list, neighbors, diag_neigh):
    debug=0
    print_debug=0
    if (debug):
        print 'donor', donor
        print 'i', i
        print 'receiver', receiver

    #if flow is 90 degrees
    if(donor in diag_neigh and receiver in diag_neigh):
        if (debug):
            print "flow is 90 degrees on diagonals from ", donor, " to ", receiver
        radcurv_angle=1.37
        #if flow is NE-SE or NW-SW, erode south node
        if (donor==diag_neigh[0] and receiver==diag_neigh[3] or 
            donor==diag_neigh[1] and receiver==diag_neigh[2]):
            lat_node=neighbors[3]
        #if flow is SW-NW or SE-NE, erode north node
        elif (donor==diag_neigh[2] and receiver==diag_neigh[1] or 
              donor==diag_neigh[3] and receiver==diag_neigh[0]):
            lat_node=neighbors[1]
        #if flow is SW-SE or NW-NE, erode east node
        elif (donor==diag_neigh[2] and receiver==diag_neigh[3] or 
              donor==diag_neigh[1] and receiver==diag_neigh[0]):
            lat_node=neighbors[0]
         #if flow is SE-SW or NE-NW, erode west node
        elif (donor==diag_neigh[3] and receiver==diag_neigh[2] or 
               donor==diag_neigh[0] and receiver==diag_neigh[1]):
             lat_node=neighbors[2]
         #print "lat_node", lat_node
    elif(donor not in diag_neigh and receiver not in diag_neigh):
        if (debug):
            print "flow is 90 degrees (not on diagonal) from ", donor, " to ", receiver
        radcurv_angle=1.37
        
        #if flow is from east, erode west node
        if (donor==neighbors[0]):
            lat_node=neighbors[2]
        #if flow is from north, erode south node
        elif (donor==neighbors[1]):
            lat_node=neighbors[3]
        #if flow is from west, erode east node
        elif (donor==neighbors[2]):
            lat_node=neighbors[0]
        #if flow is from south, erode north node
        elif (donor==neighbors[3]):
            lat_node=neighbors[1]
    if(debug):
        print 'lat_node', lat_node
    return lat_node, radcurv_angle
