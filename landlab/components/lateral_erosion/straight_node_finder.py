"""




"""

import numpy as np




def StraightNode(donor, i, receiver, neighbors, diag_neigh):
    debug=0
    print_debug=0

#####FLOW LINK IS STRAIGHT, NORTH TO SOUTH######
    if ((donor==neighbors[1] or donor==neighbors[3])):
    
        #if (debug):
    #	print "flow is stright, N-S from ", donor, " to ", flowdirs[i]
        radcurv_angle=0.23
        if(print_debug):
            print("erode node to east or west")
       #neighbors are ordered E,N,W, S
        #if the west cell is boundary (neighbors=-1), erode from east node
        if neighbors[2]==-1:
            lat_node=neighbors[0]
            if(print_debug):
                print ("eroding east node, id = ", lat_node)
        elif neighbors[0]==-1:
            lat_node=neighbors[2]
            if(print_debug):
                print( "eroding west node, id = ", lat_node)

        else:
            #if could go either way, choose randomly. 0 goes East, 1 goes west
            ran_num=np.random.randint(0,2)
            if ran_num==0:
                lat_node=neighbors[0]
                if(print_debug):
                    print ("eroding east node (random)", lat_node)
            if ran_num==1:
                lat_node=neighbors[2]
                if(print_debug):
                    print( "eroding west node (random)", lat_node)

    #####FLOW LINK IS STRAIGHT, EAST-WEST#####	
    elif (donor==neighbors[0] or donor==neighbors[2]):
        if (debug):
            print ("flow is stright, E-W")
        radcurv_angle=0.23
        if(print_debug):
            print ("erode node to north or south")
    #  Node list are ordered as [E,N,W,S]
    #if the north cell is boundary (neighbors=-1), erode from south node
        if neighbors[1]==-1:
            lat_node=neighbors[3]
            if(print_debug):  
                print ("eroding south node, id = ", lat_node)
        elif neighbors[3]==-1:
            lat_node=neighbors[1]
            if(print_debug):
                print ("eroding north node, id = ", lat_node)
        else:
    #if could go either way, choose randomly. 0 goes south, 1 goes north
            ran_num=np.random.randint(0,2)
            if ran_num==0:
                lat_node=neighbors[1]
                if(print_debug):
                    print ("eroding north node (random), id = ", lat_node)
            if ran_num==1:
                lat_node=neighbors[3]
                if(print_debug):
                    print ("eroding south node (random), id = ", lat_node)

    #if flow is straight across diagonal, choose node to erode at random
    elif(donor in diag_neigh and receiver in diag_neigh):
        radcurv_angle=0.23
        if (debug):
            print ("flow is straight across diagonal")
        if receiver==diag_neigh[0]:
            if(print_debug):
                print ("erode east or north")
            poss_diag_nodes=neighbors[0:1+1]
            if(print_debug):
                print ("poss_diag_nodes", poss_diag_nodes)
        elif receiver==diag_neigh[1]:
            if(print_debug):
                print ("erode north or west")
            poss_diag_nodes=neighbors[1:2+1]
            if(print_debug):
                print ("poss_diag_nodes", poss_diag_nodes)
        elif receiver==diag_neigh[2]:
            if(print_debug):
                print ("erode west or south")
            poss_diag_nodes=neighbors[2:3+1]
            if(print_debug):
                print( "poss_diag_nodes", poss_diag_nodes)
        elif receiver==diag_neigh[3]:
            if(print_debug):
                print ("erode south or east")
            poss_diag_nodes=[neighbors[3], neighbors[0]]
            if(print_debug):
                print( "poss_diag_nodes", poss_diag_nodes)
        ran_num=np.random.randint(0,2)
        if ran_num==0:
            lat_node=poss_diag_nodes[0]
            if(print_debug):
                print ("eroding first poss diag node (random), id = ", lat_node)
        if ran_num==1:
            lat_node=poss_diag_nodes[1]
            if(print_debug):
                print ("eroding second poss diag node (random), id = ", lat_node)
    return lat_node, radcurv_angle
