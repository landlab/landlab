#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def imshow_grid(grid, values, var_name=None, var_units=None,
                grid_units=(None, None), symmetric_cbar=False,
                limits = None, cmap='jet'):
    data = values.view()
    data.shape = grid.shape

    y = np.arange(data.shape[0] + 1) - grid.dx * .5
    x = np.arange(data.shape[1] + 1) - grid.dx * .5

    if symmetric_cbar:
        (var_min, var_max) = (data.min(), data.max())
        limit = max(abs(var_min), abs(var_max))
        limits = (-limit, limit)
    else:
        limits = (None, None)

    plt.pcolormesh(x, y, data, vmin=limits[0], vmax=limits[1], cmap=cmap)

    plt.gca().set_aspect(1.)
    plt.autoscale(tight=True)

    plt.colorbar()

    plt.xlabel('X (%s)' % grid_units[1])
    plt.ylabel('Y (%s)' % grid_units[0])

    if var_name is not None:
        plt.title('%s (%s)' % (var_name, var_units))

    #plt.show()


def imshow_field(field, name, **kwds):
    imshow_grid(field, field.field_values('node', name), var_name=name,
                var_units=field.field_units('node', name), **kwds)

###
# Added by Sai Nudurupati 29Oct2013 
# This function is exactly the same as imshow_grid but this function plots
# arrays spread over cells rather than nodes

def imshow_active_cells(grid, values, var_name=None, var_units=None,
                grid_units=(None, None), symmetric_cbar=False,
                cmap='jet'):
    data = values.view()
    data.shape = (grid.shape[0]-2, grid.shape[1]-2)

    y = np.arange(data.shape[0]) - grid.dx * .5
    x = np.arange(data.shape[1]) - grid.dx * .5

    if symmetric_cbar:
        (var_min, var_max) = (data.min(), data.max())
        limit = max(abs(var_min), abs(var_max))
        limits = (-limit, limit)
    else:
        limits = (None, None)

    plt.pcolormesh(x, y, data, vmin=limits[0], vmax=limits[1], cmap=cmap)

    plt.gca().set_aspect(1.)
    plt.autoscale(tight=True)

    plt.colorbar()

    plt.xlabel('X (%s)' % grid_units[1])
    plt.ylabel('Y (%s)' % grid_units[0])

    if var_name is not None:
        plt.title('%s (%s)' % (var_name, var_units))

    plt.show()

###