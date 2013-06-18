#! /usr/env/python
"""
Runs RasterModelGrid's unit_test() method
"""

from landlab import model_grid

reload(model_grid)
    
mg = model_grid.RasterModelGrid()

mg.unit_test()
