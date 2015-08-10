#! /usr/bin/env python

import numpy as np

''' 
Functions which can map data contained on links  onto their neighboring nodes.

Each link has a "from" and "to" node. The "from" nodes are located at the link 
tail, while the "to" nodes are located at link heads.

Below, the numbering scheme for links in RasterModelGrid is illustrated with an
example of a four-row by five column grid (4x5). In this example, each * (or X)
is a node, the lines represent links, and the ^ and > symbols indicate the
direction and "head" of each link. Link heads in the RasterModelGrid always
point in the cardinal directions North (N) or East (E).

    *--27-->*--28-->*--29-->*--30-->*
    ^       ^       ^       ^       ^
   10      11      12      13      14
    |       |       |       |       |
    *--23-->*--24-->*--25-->*--26-->*
    ^       ^       ^       ^       ^
    5       6       7       8       9
    |       |       |       |       |
    *--19-->*--20-->X--21-->*--22-->*
    ^       ^       ^       ^       ^
    0       1       2       3       4
    |       |       |       |       |
    *--15-->*--16-->*--17-->*--18-->*
  

For example, node 'X' has four link-neighbors. From south and going clockwise,
these neighbors are [2, 20, 7, 21]. Both link 2 and link 20 have node 'X' as their
'to' node, while links 7 and 21 have node 'X' as their from node. 
'''
def map_values_from_link_head_node_to_link(mg, var_name):
    '''
    map_values_from_link_head_node_to_link iterates across the grid and
    identifies the node at the "head", or the "to" node for each link. For
    each link, the value of 'var_name' at the "to" node is mapped to the corresponding
    link. 
    
    In a RasterModelGrid, each one node has two adjacent "link heads". This means
    each node value is mapped to two corresponding links. 
    '''
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    values_at_links[:] = values_at_nodes[mg.node_index_at_link_tail]


def map_values_from_link_tail_node_to_link(mg, var_name):
    '''
    map_values_from_link_tail_node_to_link iterates across the grid and
    identifies the node at the "tail", or the "from" node for each link. For
    each link, the value of 'var_name' at the "from" node is mapped to the 
    corresponding link. 
    
    In a RasterModelGrid, each one node has two adjacent "link tails". This means
    each node value is mapped to two corresponding links. 
    '''
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    values_at_links[:] = values_at_nodes[mg.node_index_at_link_head]


def map_link_end_node_min_value_to_link(mg, var_name):
    '''
    map_link_end_node_min_value_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.  
    This function evaluates the value of 'var_name' at both the "to" and "from" node.
    The minimum value of the two node values is then mapped to the link.
    '''
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    np.minimum(values_at_nodes[mg.node_index_at_link_head],
               values_at_nodes[mg.node_index_at_link_tail],
               out=values_at_links)


def map_link_end_node_max_value_to_link(mg, var_name):
    '''
    map_link_end_node_max_value_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.  
    This function evaluates the value of 'var_name' at both the "to" and "from" node.
    The maximum value of the two node values is then mapped to the link.
    '''
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    np.maximum(values_at_nodes[mg.node_index_at_link_head],
               values_at_nodes[mg.node_index_at_link_tail],
               out=values_at_links)


def map_values_from_link_end_nodes_to_link(mg, var_name):
    '''
    map_values_from_link_end_nodes_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.  
    This function takes the sum of the two values of 'var_name' at both the "to"
    and "from" node. The average value of the two node values of 'var_name'
    is then mapped to the link.
    '''
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('link', var_name)
    values_at_links = mg.at_link[var_name]
    values_at_links[:] = 0.5 * (values_at_nodes[mg.node_index_at_link_head] +
                                values_at_nodes[mg.node_index_at_link_tail])

def map_values_from_cell_node_to_cell(mg, var_name):
    '''
    map_values_from_cell_node_to_cell iterates across the grid and
    identifies the all node values of 'var_name'.  
    
    This function takes node values of 'var_name' and mapes that value to the
    corresponding cell area for each node.
    '''
    values_at_nodes = mg.at_node[var_name]
    mg.add_empty('cell', var_name)
    values_at_cells = mg.at_cell[var_name]
    values_at_cells[:] = values_at_nodes[mg.node_index_at_cells]
