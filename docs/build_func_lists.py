#! /usr/env/python
"""
This script auto-constructs strings for Landlab methods, based around
Sphinx and the new LLCATS type declaration system.

It saves as text files strings of this format, as seen at the tops of the grid
files, and used by sphinx to build our docs. These can be pasted manually into
each file.

'''
Information about the grid as a whole
+++++++++++++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.raster.RasterModelGrid.axis_name
    ~landlab.grid.raster.RasterModelGrid.axis_units
    ~landlab.grid.raster.RasterModelGrid.cell_grid_shape
    ...

'''
"""
from landlab.core.utils import get_categories_from_grid_methods
from copy import copy
import re

grid_types = ('ModelGrid', 'RasterModelGrid', 'VoronoiDelaunayGrid',
              'HexModelGrid', 'RadialModelGrid')
str_sequence = ('Raster', 'Irregular Voronoi-cell', 'Hexagonal', 'Radial')
paths = ('raster', 'voronoi', 'hex', 'radial')

autosummary = '\n\n.. autosummary::\n    :toctree: generated/\n\n'


all_methods_for_cat_allgrid = {}
all_cats_for_method_allgrid = {}
fails_allgrid = {}  # not needed here, but good for manual checking
for grid_type in grid_types:
    catdict, griddict, fails = get_categories_from_grid_methods(grid_type)
    all_methods_for_cat_allgrid[grid_type] = copy(catdict)
    all_cats_for_method_allgrid[grid_type] = copy(griddict)
    fails_allgrid[grid_type] = copy(fails)

grid_info__whole_grid_0 = (
'''
Information about the grid as a whole
+++++++++++++++++++++++++++++++++++++

''')
for grid, print_name, path in zip(grid_types[1:], str_sequence, paths):
    grid_info__whole_grid_0 += (
        '\n' + '.. _GINF_' + grid + ':\n\n' +
        print_name + '\n' + '-'*len(print_name) + autosummary)
    allmethsGINF = all_methods_for_cat_allgrid[grid]['GINF']
    allmethsGINF.sort()
    for meth in allmethsGINF:
        if 'DEPR' not in all_cats_for_method_allgrid[grid][meth]:
            grid_info__whole_grid_0 += (
                '    ~landlab.grid.' + path + '.' + grid + '.' + meth + '\n'
            )
    grid_info__whole_grid_0 += '\n\n'

f = open('./auto_grid_info__whole_grid_0.rst', "wb")
f.write(grid_info__whole_grid_0)
f.close()


