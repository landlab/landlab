# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 08:10:22 2015

@author: gtucker
"""

from landlab import HexModelGrid
import numpy as np



def make_links_at_node_array(grid):
    """Make array with links at each node"""
    
    # Create arrays for link-at-node information
    grid.gt_links_at_node = -np.ones((6, grid.number_of_nodes), dtype=np.int32)
    grid.gt_link_dirs_at_node = np.zeros((6, grid.number_of_nodes), dtype=np.int8)
    grid.gt_active_link_dirs_at_node = np.zeros((6, grid.number_of_nodes), dtype=np.int8)
    grid.gt_num_links_at_node = np.zeros(grid.number_of_nodes, dtype=np.uint8)  # assume <256 links at any node
    grid.gt_num_active_links_at_node = np.zeros(grid.number_of_nodes, dtype=np.uint8)  # assume <256 links at any node
    
    # Sweep over all links
    for lk in xrange(grid.number_of_links):
        
        # Find the ID of the tail node
        t = grid.node_at_link_tail[lk]
        
        # Its row in the 2D array is equal to the number of links we've found
        # so far for this node
        index = grid.gt_num_links_at_node[t]
        
        # Add this link to the list for this node, and increment the number
        # found so far
        grid.gt_links_at_node[index][t] = lk
        grid.gt_num_links_at_node[t] += 1
        

def testing_flux_divergence_with_hex():
    """Test flux divergence function(s).

    Notes
    -----
    Test grid looks like this:

        (7)--1-.(8)-17-.(9)
        . .     . .     . .
       0   2  18  11   9  16
      /     \ /     \ /     \ 
    (3)--4-.(4)-12-.(5)-10-.(6)
      .     . .     . .     .
       3   6  15   8   7  13
        \ /     \ /     \ /
        (0)--5-.(1)-14-.(2)

    Node numbers in parentheses; others are link numbers; period indicates
    link head.
    """
    hmg = HexModelGrid(3, 3, reorient_links=True)
    
    f = hmg.add_zeros('link', 'test_flux')
    f[:] = np.arange(hmg.number_of_links)
    
    print('Nodes:')
    for n in range(hmg.number_of_nodes):
        print(str(n)+' '+str(hmg.node_x[n])+' '+str(hmg.node_y[n]))
    
    print('Links:')
    for lk in range(hmg.number_of_links):
        print(str(lk)+' '+str(hmg.node_at_link_tail[lk])+' '+
              str(hmg.node_at_link_head[lk])+' '+str(f[lk]))
              
    make_links_at_node_array(hmg)
    

if __name__=='__main__':
    testing_flux_divergence_with_hex()
    