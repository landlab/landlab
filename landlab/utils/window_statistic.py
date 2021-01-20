# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:35:49 2020

@author: laure
"""

import numpy as np

from landlab import FieldError

def calculate_window_statistic(grid,field,func,search_radius=100,
                               calc_on_closed_nodes=True,**kwargs):
    """Calculate a statistic using a function within a search window.
    This works on any ModelGrid type.

    Currently this utility is not optimized to run quickly on large grids
    (e.g. 250x250 or larger), as the distance to every node on the grid is
    calculated at each iteration of the for loop.
    A future version will act only on a subset of nodes.

    Parameters
    ----------
    grid : ModelGrid
        Landlab ModelGrid object.
    field : Field attributed to ModelGrid nodes
        User-input field on which to calculate the statistic of interest.
        Must exist in grid.
    func : function
        User-input function that acts on the field.
    search_radius : float, optional (defaults to 100 m)
        Radius of window within which the statistic is calculated.
    calc_on_closed_nodes : boolean, optional (defaults to True)
        Toggle calculation over all nodes including closed nodes (True) or all
        nodes except closed nodes (False)
    **kwargs : optional
        Keyword arguments that must be passed to func in addition to the field.
    
    Returns
    -------
    output : array
        Output array containing the calculated values of the statistic.
        Same length as input field. 
    
    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((5, 5), xy_spacing=10.0)
    >>> grid.set_closed_boundaries_at_grid_edges(False,True,False,True)
    >>> z = grid.add_zeros("topographic__elevation", at="node")
    >>> z += np.arange(len(z))
    
    # Calculate relief using np.ptp function.
    >>> relief = calculate_window_statistic(grid,'topographic__elevation',
                                            np.ptp,search_radius=15)
    >>> grid.at_node['topographic__elevation']
    array([0.,   1.,   2.,   3.,   4.,
           5.,   6.,   7.,   8.,   9.,
           10.,  11.,  12.,  13.,  14.,
           15.,  16.,  17.,  18.,  19.,
           20.,  21.,  22.,  23.,  24.])
    >>> relief
    array([6.,   7.,   7.,   7.,   6.,
           11.,  12.,  12.,  12.,  11.,
           11.,  12.,  12.,  12.,  11.,
           11.,  12.,  12.,  12.,  11.,
           6.,   7.,   7.,   7.,   6.])
    
    # Calculate relief using np.ptp function excluding closed nodes.
    >>> relief = calculate_window_statistic(grid,'topographic__elevation',
                                        np.ptp,search_radius=15,
                                        calc_on_closed_nodes=False)
    >>> grid.at_node['topographic__elevation']
    array([0.,   1.,   2.,   3.,   4.,
           5.,   6.,   7.,   8.,   9.,
           10.,  11.,  12.,  13.,  14.,
           15.,  16.,  17.,  18.,  19.,
           20.,  21.,  22.,  23.,  24.])
    >>> relief
    array([nan,  nan,  nan,  nan,  nan,
           6.,   7.,   7.,   7.,   6.,
           11.,  12.,  12.,  12.,  11.,
           6.,   7.,   7.,   7.,   6.,
           nan,  nan,  nan,  nan,  nan])
    
    # Calculate 90th percentile using np.percentile function and **kwargs.
    >>> perc_90 = calculate_window_statistic(grid,'topographic__elevation',
                                             np.percentile,search_radius=15,
                                             calc_on_closed_nodes=False,
                                             q=90)
    >>> grid.at_node['topographic__elevation']
    array([0.,   1.,   2.,   3.,   4.,
           5.,   6.,   7.,   8.,   9.,
           10.,  11.,  12.,  13.,  14.,
           15.,  16.,  17.,  18.,  19.,
           20.,  21.,  22.,  23.,  24.])
    >>> perc_90
    array([nan,   nan,   nan,   nan,   nan,
           10.7,  11.5,  12.5,  13.5,  13.7,
           15.5,  16.2,  17.2,  18.2,  18.5,
           15.7,  16.5,  17.5,  18.5,  18.7,
           nan,   nan,   nan,   nan,   nan])
    """

    if field not in grid.at_node:
        raise FieldError(
            f"A {field} field is required at the nodes of the input grid."
        )
    
    # Create output array
    output = np.zeros(grid.number_of_nodes)
   
    # Create arrays of x and y coords for input to "distance to point' calc
    x_coord = grid.x_of_node
    y_coord = grid.y_of_node
    
    nodes_in_loop = grid.nodes.flatten()
    nodes_to_include = np.ones(grid.number_of_nodes,dtype=bool)
    
    if calc_on_closed_nodes is False:
        closed_nodes = grid.status_at_node == grid.BC_NODE_IS_CLOSED
        nodes_in_loop = nodes_in_loop[~closed_nodes]
        nodes_to_include[closed_nodes] = False
        output[closed_nodes] = np.NaN
        
    # Calculate "dist to point" then local relief at nodes within radius
    for node in nodes_in_loop:
        node_dist_to_point = grid.calc_distances_of_nodes_to_point(
            (x_coord[node], y_coord[node])
        )
        nodes_in_window = np.all([node_dist_to_point <= search_radius,
                                  nodes_to_include],0)
        values_in_window = grid.at_node[field][nodes_in_window]
        output[node] = func(values_in_window, **kwargs)
        
    return output
