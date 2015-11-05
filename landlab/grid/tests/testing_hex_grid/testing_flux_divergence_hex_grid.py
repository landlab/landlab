# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 08:10:22 2015

@author: gtucker
"""

from landlab import HexModelGrid
import numpy as np

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
    

if __name__=='__main__':
    testing_flux_divergence_with_hex()
    