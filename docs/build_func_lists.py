#! /usr/env/python
"""This script auto-constructs strings for Landlab methods, based around Sphinx
and the new LLCATS type declaration system.

It saves as text files strings of this format, as seen at the tops of the grid
files, and used by sphinx to build our docs. These can be pasted manually into
each file.

'''
Information about the grid as a whole
+++++++++++++++++++++++++++++++++++++

.. autosummary::

    ~landlab.grid.raster.RasterModelGrid.axis_name
    ~landlab.grid.raster.RasterModelGrid.axis_units
    ~landlab.grid.raster.RasterModelGrid.cell_grid_shape
    ...

'''
"""
import re
from copy import copy

import numpy as np

from landlab.core.utils import get_categories_from_grid_methods

grid_types = ('ModelGrid', 'RasterModelGrid', 'VoronoiDelaunayGrid',
              'HexModelGrid', 'RadialModelGrid')
str_sequence = ('Base class', 'Raster', 'Irregular Voronoi-cell', 'Hexagonal',
                'Radial')
paths = ('base', 'raster', 'voronoi', 'hex', 'radial')

autosummary = '\n\n.. currentmodule:: landlab \n\n.. autosummary::\n\n'


all_methods_for_cat_allgrid = {}
all_cats_for_method_allgrid = {}
fails_allgrid = {}  # not needed here, but good for manual checking
for grid_type in grid_types:
    catdict, griddict, fails = get_categories_from_grid_methods(grid_type)
    all_methods_for_cat_allgrid[grid_type] = copy(catdict)
    all_cats_for_method_allgrid[grid_type] = copy(griddict)
    fails_allgrid[grid_type] = copy(fails)


LLCATS = ('GINF', 'NINF', 'LINF', 'CINF', 'PINF', 'FINF', 'CNINF', 'GRAD',
          'MAP', 'BC', 'SUBSET', 'SURF')
outfile_names = ('whole.rst', 'nodes.rst',
                 'links.rst', 'cells.rst',
                 'patches.rst', 'faces.rst',
                 'corners.rst', 'grads.rst',
                 'mappers.rst', 'BCs.rst',
                 'subsets.rst', 'surf.rst')

text_heads = (
'''
=====================================
Information about the grid as a whole
=====================================

''',
'''
=======================
Information about nodes
=======================

''',
'''
=======================
Information about links
=======================

''',
'''
=======================
Information about cells
=======================

''',
'''
=========================
Information about patches
=========================

''',
'''
=======================
Information about faces
=======================

''',
'''
=========================
Information about corners
=========================

''',
'''
==============================================
Gradients, fluxes, and divergences on the grid
==============================================

''',
'''
=======
Mappers
=======

''',
'''
==========================
Boundary condition control
==========================

''',
'''
========================
Identifying node subsets
========================

''',
'''
================
Surface analysis
================

''')

for LLCAT, outfile_name, text in zip(LLCATS, outfile_names, text_heads):
    for grid, print_name, path in zip(grid_types, str_sequence, paths):
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

    f = open('source/reference/grid/auto/' + outfile_name, "w")
    f.write(text)
    f.close()
