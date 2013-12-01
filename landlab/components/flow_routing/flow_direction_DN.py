#! /usr/env/python

"""
flow_direction_DN.py: calculates single-direction flow directions on a regular
or irregular grid.

GT Nov 2013
"""


import numpy


def flow_directions(elev, fromnode, tonode, link_slope, baselevel_nodes=None):
    """
    Finds and returns flow directions for a given elevation grid. Each node is
    assigned a single direction, toward one of its N neighboring nodes (or
    itself, if none of its neighbors are lower).
    
    Inputs:
        elev = array or list of elevations at nodes
        active_links = array or list of links (ordered pairs of nodes)
        fromnode = array or list containing ID of the "from" node for each link
        tonode = array or list containing ID of the "to" node for each link
        link_slope = slope of each link, defined POSITIVE DOWNHILL (i.e., a
                     negative value means the link runs uphill from the fromnode
                     to the tonode)
    
    Returns:
        receiver = Numpy array containing, for each node, the ID of the node
                   that receives its flow. Defaults to the node itself if no 
                   other receiver is assigned.
        steepest_slope = the slope value (positive downhill) in the direction of
                         flow
    
    The example below assigns elevations to the 10-node example network in
    Braun and Willett (2012), so that their original flow pattern should be
    re-created.
    
    Example:
        >>> z = numpy.array([2.4, 1.0, 2.2, 3.0, 0.0, 1.1, 2.0, 2.3, 3.1, 3.2])
        >>> fn = numpy.array([1,4,4,0,1,2,5,1,5,6,7,7,8,6,3,3,2,0])
        >>> tn = numpy.array([4,5,7,1,2,5,6,5,7,7,8,9,9,8,8,6,3,3])
        >>> s = z[fn] - z[tn]  # slope with unit link length, positive downhill
        >>> r, ss, snk = flow_directions(z, fn, tn, s)
        >>> r
        array([1, 4, 1, 6, 4, 4, 5, 4, 6, 7])
        >>> ss
        array([ 1.4,  1. ,  1.2,  1. ,  0. ,  1.1,  0.9,  2.3,  1.1,  0.9])
        >>> snk
        array([4])
    """
    
    # Setup
    num_nodes = len(elev)
    steepest_slope = numpy.zeros(num_nodes)
    receiver = numpy.arange(num_nodes)
    
    # For each link, find the higher of the two nodes. The higher is the
    # potential donor, and the lower is the potential receiver. If the slope
    # from donor to receiver is steeper than the steepest one found so far for
    # the donor, then assign the receiver to the donor and record the new slope.
    # (Note the minus sign when looking at slope from "t" to "f").
    for i in range(len(fromnode)):
        f = fromnode[i]
        t = tonode[i]
        #print 'link from',f,'to',t,'with slope',link_slope[i]
        if elev[f]>elev[t] and link_slope[i]>steepest_slope[f]:
            receiver[f] = t
            steepest_slope[f] = link_slope[i]
            #print ' flows from',f,'to',t
        elif elev[t]>elev[f] and -link_slope[i]>steepest_slope[t]:
            receiver[t] = f
            steepest_slope[t] = -link_slope[i]
            #print ' flows from',t,'to',f
        #else:
            #print ' is flat'
            
    node_id = numpy.arange(num_nodes)

    # Optionally, handle baselevel nodes: they are their own receivers
    if baselevel_nodes is not None:
        receiver[baselevel_nodes] = node_id[baselevel_nodes]
    
    # The sink nodes are those that are their own receivers (this will normally
    # include boundary nodes as well as interior ones; "pits" would be sink
    # nodes that are also interior nodes).
    (sink, ) = numpy.where(node_id==receiver)
    
    return receiver, steepest_slope, sink
    
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
