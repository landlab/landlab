"""
lateral erosion component

needs:
fraction of bed exposed: fe=(1-qsqt)
percent lateral erosion of total erosion : Plat=qsqt
lateral erosion: Elat=Plat*E  (E is total erosion, dzdt)
vertical erosion: Evert=(1-Elat)*E



"""

import numpy as np
from landlab.components.lateral_erosion.angle_finder import angle_finder
from landlab.components.lateral_erosion.straight_node_finder import StraightNode
from landlab.components.lateral_erosion.fortyfive_node import FortyfiveNode
from landlab.components.lateral_erosion.ninety_node import NinetyNode


def Node_Finder2(grid, i, flowdirs, drain_area):
    debug=0
    print_debug=0


    #receiver node of flow is flowdirs[i]
    receiver=flowdirs[i]

    #find indicies of where flowdirs=i to find donor nodes.
    #will donor nodes always equal the index of flowdir list?
    inflow=np.where(flowdirs==i)
    #if there are more than 1 donors, find the one with largest drainage area
    if len(inflow[0])>1:
        drin=drain_area[inflow]
        drmax=max(drin)
        maxinfl=inflow[0][np.where(drin==drmax)]
        if(debug):
            print( "old inflow", inflow[0])
#            print "max(drin)", max(drin)
#            print "maxinfl", maxinfl
        #inind=np.where(drin==drmax)
        #if donor nodes have same drainage area, choose one randomly
        if len(maxinfl)>1:
            ran_num=np.random.randint(0,len(maxinfl))
            maxinfln=maxinfl[ran_num]
            donor=[maxinfln]
            if(debug):
                print ("random donor", donor)

        else:
            donor=maxinfl
            if(debug):
                print ("donor with larger drainage area", donor)
        #else donor is the only inflow
    else:
        donor=inflow[0]

    if(print_debug):

        print ("donor", donor)
        print ("i", i)
        print ("receiver", receiver)
    #now we have chosen donor cell, next figure out if inflow/outflow lines are
    #straight, 45, or 90 degree angle. and figure out which node to erode

    link_list=grid.links_at_node[i]
    #this gives list of active neighbors for specified node
    # the order of this list is: [E,N,W,S]
    neighbors=grid.active_neighbors_at_node[i]    
    #this gives list of all diagonal neighbors for specified node
    # the order of this list is: [NE,NW,SW,SE]
    diag_neigh=grid.diagonal_adjacent_nodes_at_node[i]
    linkdirs=grid.active_link_dirs_at_node[i]
#    print 'linkdirs', linkdirs

    angle_diff=angle_finder(grid, donor, i, receiver)

    if(debug):
        #print "conlink1", conn_link1
        #print "conlink2", conn_link2
        print ("link_list", link_list)
        print ("neighbors", neighbors)
        print ("diagneighbors", diag_neigh)
        print ("angle_diff", angle_diff)
        print (" ")

    if donor == flowdirs[i]:
        #this is a sink. no lateral ero
        if(debug):
            print ("this is a sink")
        radcurv_angle=0.
        lat_node=0
    if angle_diff==0.0:
        [lat_node, radcurv_angle]=StraightNode(donor, i, receiver, neighbors, diag_neigh)
    if (angle_diff==45.0 or angle_diff==135.0):
        [lat_node, radcurv_angle]=FortyfiveNode(donor, i, receiver, link_list, neighbors, diag_neigh)
    if angle_diff==90.0:
        [lat_node, radcurv_angle]=NinetyNode(donor, i, receiver, link_list, neighbors, diag_neigh)
    if(debug):
        print("lat_node", lat_node)
        print("end of node finder")
    if lat_node > 2e9:
        #print "old latnode", lat_node
        lat_node=0
        radcurv_angle=0.0
        #print "new lat", lat_node
        #print delta
    #
    #
    ##else:
    #    #print "node ", i, " has no upstream donor"
    #    #radcurv_angle=0.
    #
    #
    #lat_node=0
    #radcurv_angle=0.0
#    print(delta)
    #below added 9/18/2014
    dx=grid.dx
    radcurv_angle=radcurv_angle/dx    #May24, 2017: this is actually INVERSE radius of curvature. It works out in the main lateral ero
    return lat_node, radcurv_angle
