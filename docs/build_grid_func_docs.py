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
import numpy as np

# grid_types = ('ModelGrid', 'RasterModelGrid', 'VoronoiDelaunayGrid',
#               'HexModelGrid', 'RadialModelGrid')
str_sequence = ('Base class', 'Raster', 'Irregular Voronoi-cell', 'Hexagonal',
                'Radial')
paths = ('base', 'raster', 'voronoi', 'hex', 'radial')

autosummary = '\n\n.. autosummary::\n    :toctree: generated/\n\n'

LLCATS = ('GINF', 'NINF', 'LINF', 'CINF', 'PINF', 'FINF', 'CNINF', 'GRAD',
          'MAP', 'BC', 'SUBSET', 'SURF')
grid_name_to_class = ('base': 'ModelGrid', )
#                       'raster': 'RasterModelGrid',
#                       'voronoi': 


def create_dicts_of_cats():
    '''
    Create the dicts that record grid methods by grid and LLCAT.

    Returns
    -------
    all_methods_for_cat_allgrid : dict of dicts of lists
        First key is grid type, second key is LLCAT. Values are lists of method
        name strings.
    all_cats_for_method_allgrid : dict of dicts of lists
        First key is grid type, second key is method name strings. Values are
        lists of LLCATS assigned to each method.
    fails_allgrid : dict of lists
        Key is grid type, value is list of methods with no LLCATS.
    '''
    all_methods_for_cat_allgrid = {}
    all_cats_for_method_allgrid = {}
    fails_allgrid = {}
    for grid_type in grid_types:
        catdict, griddict, fails = get_categories_from_grid_methods(grid_type)
        all_methods_for_cat_allgrid[grid_type] = copy(catdict)
        all_cats_for_method_allgrid[grid_type] = copy(griddict)
        fails_allgrid[grid_type] = copy(fails['MISSING'])
    return (all_methods_for_cat_allgrid, all_cats_for_method_allgrid,
            fails_allgrid)

(all_methods_for_cat_allgrid, all_cats_for_method_allgrid,
 fails_allgrid) = create_dicts_of_cats()

for LLCAT, grid_to_modify in zip(LLCATS, grid_name_to_class.keys()):
    f = open('./text_for_' + grid_to_modify + '.py.txt', "rb")
    text = f.read()
    f.close()
    for print_name, path in zip(str_sequence, paths):
        grid = grid_name_to_class[grid_to_modify]
        text += (
            '\n' + '.. _' + LLCAT + '_' + grid + ':\n\n' +
            print_name + '\n' + '-'*len(print_name) + autosummary)
        try:
            allmeths = all_methods_for_cat_allgrid[grid][LLCAT]
            allmeths.sort()
            for meth in allmeths:
                if LLCAT in ('GINF', 'NINF', 'LINF', 'CINF', 'PINF', 'FINF',
                             'CNINF'):
                    exclude = ('DEPR', 'MAP', 'GRAD', 'SURF')
                else:
                    exclude = ('DEPR', )
                if np.in1d(exclude,
                           all_cats_for_method_allgrid[grid][meth]).sum() == 0:
                    text += (
                        '    ~landlab.grid.' + path + '.' + grid + '.' + meth +
                        '\n'
                    )
            text += '\n\n'
        except KeyError:
            print('Failed at ' + grid + ' looking for ' + LLCAT)

    f = open('./' + outfile_name, "wb")
    f.write(text)
    f.close()

'(None are available for this grid type)'
