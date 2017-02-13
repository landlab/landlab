#! /usr/env/python

"""
flow_direction_two_two.py: calculates two-direction flow directions.

Works on both a regular or irregular grid.

KRB Jan 2017
"""
from six.moves import range
import numpy as np
from landlab import RasterModelGrid, BAD_INDEX_VALUE, CLOSED_BOUNDARY

def flow_directions_dinfinity(elev, 
                        neighbors_at_node,
                        links_at_node,
                        active_link_dir_at_node,
                        tail_node, 
                        head_node, 
                        link_slope, 
                        baselevel_nodes=None):

    """
    

    """
    
    dum=1
    
    return dum

def flow_directions_dtrig(elev, 
                        neighbors_at_node,
                        links_at_node,
                        active_link_dir_at_node,
                        tail_node, 
                        head_node, 
                        link_slope, 
                        baselevel_nodes=None):
    """
    stuff.
    
    more stuff.
    """
    
    dum=1

    return dum


if __name__ == '__main__':
    import doctest
    doctest.testmod()