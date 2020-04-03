#! /usr/env/python
"""
This script builds several of Landlab's key grid documentation files, based
around Sphinx and the new LLCATS type declaration system. The files affected
are:

landlab.grid.base.rst
landlab.grid.raster.rst
landlab.grid.voronoi.rst
landlab.grid.radial.rst
landlab.grid.hex.rst

It takes the files named text_for_XXXX.py.txt, then uses these text blocks
coupled with lists of grid methods according to the LLCATS system to build the
above files, automatically.

This script is designed to be run as part of the commit process for LL.

Any changes made directly to the above files will be lost whenever this script
is run.
"""
import re
from copy import copy

import numpy as np

from landlab.core.utils import get_categories_from_grid_methods

grid_types = ('ModelGrid', 'RasterModelGrid', 'VoronoiDelaunayGrid',
              'HexModelGrid', 'RadialModelGrid', 'NetworkModelGrid')
str_sequence = ('Base class', 'Raster', 'Irregular Voronoi-cell', 'Hexagonal',
                'Radial', 'Network')
paths = ('base', 'raster', 'voronoi', 'hex', 'radial', 'network')

autosummary = '.. currentmodule:: landlab \n\n.. autosummary::\n\n'

LLCATS = ('GINF', 'NINF', 'LINF', 'CINF', 'PINF', 'FINF', 'CNINF', 'GRAD',
          'MAP', 'BC', 'SUBSET', 'SURF')
grid_name_to_class = {'base': 'ModelGrid',
                      'hex': 'HexModelGrid',
                      'radial': 'RadialModelGrid',
                      'raster': 'RasterModelGrid',
                      'voronoi': 'VoronoiDelaunayGrid',
                      'network': 'NetworkModelGrid'}


def create_dicts_of_cats():
    """Create the dicts that record grid methods by grid and LLCAT.

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
    """
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

for grid_to_modify in grid_name_to_class.keys():
    f = open('source/reference/grid/text_for_' + grid_to_modify + '.py.txt', "rt")
    text = f.read()
    f.close()
    for LLCAT in LLCATS:
        text_to_add = ''
        grid = grid_name_to_class[grid_to_modify]
        try:
            allmeths = all_methods_for_cat_allgrid[grid][LLCAT]
            allmeths.sort()
            text_to_add += autosummary
            for meth in allmeths:
                if LLCAT in ('GINF', 'NINF', 'LINF', 'CINF', 'PINF', 'FINF',
                             'CNINF'):
                    exclude = ('DEPR', 'MAP', 'GRAD', 'SURF')
                else:
                    exclude = ('DEPR', )
                if np.in1d(exclude,
                           all_cats_for_method_allgrid[grid][meth]).sum() == 0:
                    text_to_add += (
                        '    ~landlab.grid.' + grid_to_modify + '.' +
                        grid + '.' + meth + '\n'
                    )
        except KeyError:
            # print('For ' + grid + ' found no ' + LLCAT)
            text_to_add += '(None are available for this grid type)\n'

        text = text.replace('LLCATKEY: ' + LLCAT, text_to_add)

    f = open('source/reference/grid/' + grid_to_modify + '.rst', "wt")
    f.write(text)
    f.close()
