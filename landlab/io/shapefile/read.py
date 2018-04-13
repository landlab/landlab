#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 09:10:56 2018

@author: barnhark
"""
import os
import numpy as np

import shapefile as ps

from landlab import NetworkModelGrid


def read_shapefile(file):
    """Read shapefile and create a NetworkModelGrid.
    
    Text here about assumptions   


    Parameters
    ----------
    file : str
        File path to a valid shapefile
        
    Examples 
    --------
    >>> # an example will be here. 
    """
    
    if os.path.exists(file) == False:
        raise ValueError(('landlab.io.shapefile was passed a filepath that does '
                          'not exist.'))
    
    sf = ps.Reader(file)
    
    if sf.shapeType != 3:
        raise ValueError(('landlab.io.shapefile read requires a polyline '
                          'type shapefile. The provided shapefile does '
                          'not meet these requirements.'))
    
    # get record information, the firste element is # ('DeletionFlag', 'C', 1, 0) 
    # which we will ignore. 
    records = sf.fields[1:]
    
    # initialize data structures for node (x,y) tuples, 
    # link (head_node_id, tail_node_id) tuples, and a dictionary of at-link fields. 
    # besides the at-link fields on the shapefile, we'll also store an array of 
    # x and y of the full polyline segment'. 
    
    node_xy = []
    links = []
    fields = {rec[0]:[] for rec in records}
    fields['x_of_polyline'] = []
    fields['y_of_polyline'] = []
        
    record_order = [rec[0] for rec in records]
        
    # itterate through shapes and records
    shapeRecs = sf.shapeRecords()
    for sr in shapeRecs:
        
        # if not a multi-part polyline:
        if len(sr.shape.parts) == 1:
            
            # get all the points on the polyline and deconstruct into x and y
            points = sr.shape.points
            x, y = zip(*points)
            
            # construct the (x,y) tuples of the head and tail nodes of each 
            # polyline. Note here, that head and tail just refer to starting and 
            # ending, they will be re-oriented if necessary by landlab. 
            
            head_node_xy = (x[0], y[0])
            tail_node_xy = (x[-1], y[-1])
            
            # we should expect that the head node and tail node of later links will
            # already be part of the model grid. So we check, and add the nodes, 
            # if they don't already exist. 
            
            if head_node_xy not in node_xy:
                node_xy.append(head_node_xy)
            
            if tail_node_xy not in node_xy:
                node_xy.append(tail_node_xy)
        
            # get the index of the head and tail node index. 
            head_node__node_id = node_xy.index(head_node_xy)
            tail_node__node_id = node_xy.index(tail_node_xy)
            
            # append the head and tail node ids to the link array
            links.append((head_node__node_id, tail_node__node_id))
            
            for i in range(len(sr.record)):
                field_name = record_order[i]
                fields[field_name].append(sr.record[i])
            
            fields['x_of_polyline'].append(x)
            fields['y_of_polyline'].append(y)
            
        else:
            raise ValueError(('landlab.io.shapefile currently does not support ',
                              'reading multipart polyline shapefiles.'))
        
    
    ## Create a Network Model Grid
    x_of_node, y_of_node = zip(*node_xy)    
    grid = NetworkModelGrid((y_of_node, x_of_node), links)  
    for field_name in fields:
        grid.at_link[field_name] = fields[field_name]
        
    return grid

#%%

file = '/Users/barnhark/Downloads/Methow4Allison/Methow/Methow_Network.shp'
grid = read_shapefile(file)

#%%

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# for plotting and testing
x_of_polylines = grid['link']['x_of_polyline']
y_of_polylines = grid['link']['y_of_polyline']

segments = []

for i in range(len(x_of_polylines)):
    x = np.array(x_of_polylines[i])
    y = np.array(y_of_polylines[i])
    segment = np.array((x, y)).T
    segments.append(segment)
        

from landlab.plot import graph
fig, ax = plt.subplots(figsize=(8,8), dpi=300)


graph.plot_links(grid, color='c', linestyle='solid', with_id=False,
               as_arrow=False, linewidth=1)
graph.plot_nodes(grid, color='r', with_id=False, markersize=1)

line_segments = LineCollection(segments, color='b', linewidth=0.5)
ax.add_collection(line_segments)

plt.savefig('test.png') 
