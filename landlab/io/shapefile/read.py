#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 09:10:56 2018

@author: barnhark
"""

import numpy as np

import shapefile as ps


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


file = '/Users/barnhark/Downloads/Methow4Allison/Methow/Methow_Network.shp'

#def read_shapefile(file):
#    """Read shapefile and create a NetworkModelGrid."""

sf = ps.Reader(file)

if sf.shapeType != 3:
    raise ValueError(('landlab.io.shapefile read requires a polyline '
                      'type shapefile. The provided shapefile does '
                      'not meet these requirements.'))

numRecords = sf.numRecords
fields = sf.fields

nodes = []
links = []

fields = {'node': [],
          'link': []}


    
shapes = sf.shapes()

segments = []
for sr in sf.iterShapeRecords():
    
    #%%
    if len(shape.parts) == 1:
        points = shape.points
        x, y = zip(*points)
        x = np.array(x)
        y = np.array(y)
        segment = np.array((x, y)).T
        segments.append(segment)
    else:
        raise ValueError(('landlab.io.shapefile currently does not support ',
                          'reading multipart polyline shapefiles.'))
    
#    sf.fields
    
#%%
fig, ax = plt.subplots()
line_segments = LineCollection(segments, linewidth=1)
ax.add_collection(line_segments)
plt.xlim(sf.bbox[0], sf.bbox[2])
plt.ylim(sf.bbox[1], sf.bbox[3])

plt.show()