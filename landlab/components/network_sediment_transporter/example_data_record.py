# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:20:19 2018

@author: pfeif
"""

import numpy as np

from landlab import RasterModelGrid
from landlab.data_record import DataRecord
from landlab.plot import graph

grid = RasterModelGrid((3, 3))

graph.plot_graph(grid, at="node,link")


# Example from DataRecord doc
dr1 = DataRecord(
    grid,
    time=[0.],
    data_vars={"mean_elevation": (["time"], np.array([100]))},
    attrs={"time_units": "y"},
)


# Vaguely how we'll be using it


my_items = {
    "grid_element": np.array(("node", "link"), dtype=str),
    "element_id": np.array([1, 3]),
}

dr = DataRecord(grid, items=my_items)
