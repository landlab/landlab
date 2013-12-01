#!/usr/env/python

"""
flow_accum_bw.py:
    
Implementation of Braun & Willett (2012) algorithm for calculating drainage
area and (optionally) water discharge. Assumes each node has only one downstream
receiver. If water discharge is calculated, the result assumes steady flow 
(that is, hydrologic equilibrium).

The main public function is 

    a, q, s = flow_accumulation(r, b)

which takes an array of receiver-node IDs, r (the nodes that "receive" the flow 
from a each node; this array would be returned by the flow_routing component's
calc_flowdirs() method), and an array of baselevel nodes, b. It returns Numpy 
arrays with the drainage area (a) and discharge (q) at each node, along with an 
array (s) that contains the IDs of the nodes in downstream-to-upstream order.

If you simply want the ordered list by itself, use:
    
    s = make_ordered_node_array(r, b)

Created: GT Nov 2013
"""

import numpy


class _DrainageStack():
    """
    The _DrainageStack() class implements Braun & Willett's add_to_stack
    function (as a method) and also keeps track of the counter (j) and the
    stack (s). It is used by the make_ordered_node_array() function.
    """
    def __init__(self, delta, D):
        """
        Initializes the index counter j to zero, creates the stack array s,
        and stores references to delta and D.
        """
        self.j = 0
        self.s = numpy.zeros(len(D), dtype=int)
        self.delta = delta
        self.D = D

    def add_to_stack(self, l):
        """
        Adds node l to the stack and increments the current index (j).
        
        Example:
            >>> delta = numpy.array([ 0,  0,  2,  2,  2,  6,  7,  9, 10, 10, 10])
            >>> D = numpy.array([0, 2, 1, 4, 5, 7, 6, 3, 8, 9])
            >>> ds = _DrainageStack(delta, D)
            >>> ds.add_to_stack(4)
            >>> ds.s
            array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
        """
        self.s[self.j] = l
        self.j += 1
        for n in range(self.delta[l], self.delta[l+1]):
            m = self.D[n]
            if m != l:
                self.add_to_stack(m)
                
                
def _make_number_of_donors_array(r):
    """
    Creates and returns an array containing the number of donors for each node.
    
    Inputs: r = Numpy array with ID of receiver for each node.
    Returns: Numpy array containing number of donors for each node.
    
    The example below is from Braun and Willett (2012); nd corresponds to their
    d_i in Table 1.
    
    Example:
        >>> r = numpy.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
        >>> nd = _make_number_of_donors_array(r)
        >>> nd
        array([0, 2, 0, 0, 4, 1, 2, 1, 0, 0])
    """
    np = len(r)
    nd = numpy.zeros(np, dtype=int)
    for i in range(np):
        nd[r[i]] += 1
    return nd
    
    
def _make_delta_array(nd):
    """
    Creates and returns the "delta" array, which is a list containing, for each 
    node, the array index where that node's donor list begins.
    
    Inputs: nd = array containing number of donors for each node
    Returns: delta array
    
    The example below is from Braun and Willett (2012), and represents 
    \delta_i in their Table 1. Here, the numbers are all one less than in their
    table because here we number indices from 0 rather than 1.
    
    Example:
        >>> nd = numpy.array([0, 2, 0, 0, 4, 1, 2, 1, 0, 0])
        >>> delta = _make_delta_array(nd)
        >>> delta
        array([ 0,  0,  2,  2,  2,  6,  7,  9, 10, 10, 10])
    """
    np = len(nd)
    delta = numpy.zeros(np+1, dtype=int)
    delta[np] = np   # not np+1 as in B&W because here we number from 0
    for i in range(np-1, -1, -1):
        delta[i] = delta[i+1] - nd[i]
    return delta
    
    
def _make_array_of_donors(r, delta):
    """
    Creates and returns an array containing the IDs of donors for each node.
    Essentially, the array is a series of lists (not in the Python list object
    sense) of IDs for each node. See Braun & Willett (2012) for details.
    
    The example below is from Braun & Willett (2012), and produces D_i in their
    Table 1 (except that here the ID numbers are one less, because we number
    indices from zero).
    
    Example:
        >>> r = numpy.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
        >>> delta = numpy.array([ 0,  0,  2,  2,  2,  6,  7,  9, 10, 10, 10])
        >>> D = _make_array_of_donors(r, delta)
        >>> D
        array([0, 2, 1, 4, 5, 7, 6, 3, 8, 9])
    """    
    np = len(r)
    w = numpy.zeros(np, dtype=int)
    D = numpy.zeros(np, dtype=int)
    for i in range(0, np):
        ri = r[i]
        D[delta[ri]+w[ri]] = i
        w[ri] += 1
    return D


def make_ordered_node_array(receiver_nodes, baselevel_nodes):
    """
    Creates and returns an array of node IDs that is arranged in order from
    downstream to upstream. 
    
    The lack of a leading underscore is meant to signal that this operation
    could be useful outside of this module!
    
    Example:
        >>> r = numpy.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
        >>> b = numpy.array([4])
        >>> s = make_ordered_node_array(r, b)
        >>> s
        array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
    """
    nd = _make_number_of_donors_array(receiver_nodes)
    delta = _make_delta_array(nd)
    D = _make_array_of_donors(receiver_nodes, delta)
    dstack = _DrainageStack(delta, D)
    for k in baselevel_nodes:
        dstack.add_to_stack(k)
    return dstack.s
    
    
def find_drainage_area_and_discharge(s, r, node_cell_area=1.0, runoff=1.0):
    """
    Calculates and returns the drainage area and water discharge at each node.
    
    Inputs: s = ordered (downstream to upstream) array of node IDs
            r = array of receiver IDs for each node
            node_cell_area = scalar or numpy array of cell surface areas for
                             each node. If it's an array, must have same length
                             as s (that is, the number of nodes).
            runoff = scalar or numpy array of local runoff rate at each cell
                     (in water depth per time). If it's an array, must have same
                     length as s (that is, the number of nodes).
                     
    Returns: drainage area and discharge as Numpy arrays
    
    Notes:
        - If node_cell_area not given, the output drainage area is equivalent
          to the number of nodes/cells draining through each point, including
          the local node itself.
        - Give node_cell_area as a scalar when using a regular raster grid.
        - If runoff is not given, the discharge returned will be the same as
          drainage area (i.e., drainage area times unit runoff rate).
        - If using an unstructured Landlab grid, make sure that the input
          argument for node_cell_area is the cell area at each NODE rather than
          just at each CELL. This means you need to include entries for the
          perimeter nodes too. They can be zeros.
          
    Example:
        >>> r = numpy.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
        >>> s = numpy.array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
        >>> a, q = find_drainage_area_and_discharge(s, r)
        >>> a
        array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
        >>> q
        array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
        
    """
    
    # Number of points
    np = len(s)
    
    # Initialize the drainage_area and discharge arrays. Drainage area starts
    # out as the area of the cell in question, then (unless the cell has no
    # donors) grows from there. Discharge starts out as the cell's local runoff
    # rate times the cell's surface area.
    drainage_area = numpy.zeros(np) + node_cell_area
    discharge = numpy.zeros(np) + node_cell_area*runoff
    
    # Iterate backward through the list, which means we work from upstream to
    # downstream.
    for i in range(np-1, -1, -1):
        donor = s[i]
        recvr = r[donor]
        if donor != recvr:
            drainage_area[recvr] += drainage_area[donor]
            discharge[recvr] += discharge[donor]
            
    return drainage_area, discharge
    

def flow_accumulation(receiver_nodes, baselevel_nodes, node_cell_area=1.0,
                      runoff_rate=1.0):
    """
    Calculates and returns the drainage area and (steady) discharge at each
    node, along with a downstream-to-upstream ordered list (array) of node IDs.
    
    Example:
        >>> r = numpy.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8])-1
        >>> b = numpy.array([4])
        >>> a, q, s = flow_accumulation(r, b)
        >>> a
        array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
        >>> q
        array([  1.,   3.,   1.,   1.,  10.,   4.,   3.,   2.,   1.,   1.])
        >>> s
        array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
    """
    
    s = make_ordered_node_array(receiver_nodes, baselevel_nodes)
    
    a, q = find_drainage_area_and_discharge(s, receiver_nodes, node_cell_area,
                                            runoff_rate)
    
    return a, q, s
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()
