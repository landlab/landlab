#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


def imshow_grid(grid, values, var_name=None, var_units=None,
                grid_units=(None, None), symmetric_cbar=False):
    data = values.view()
    data.shape = grid.shape

    y = np.linspace(0, grid.get_grid_ydimension(), data.shape[0])
    x = np.linspace(0, grid.get_grid_xdimension(), data.shape[1])

    if symmetric_cbar:
        (var_min, var_max) = (data.min(), data.max())
        limit = max(abs(var_min), abs(var_max))
        limits = (-limit, limit)
    else:
        limits = (None, None)

    plt.pcolormesh(y, x, data, vmin=limits[0], vmax=limits[1])
    plt.colorbar()

    plt.xlabel('X (%s)' % grid_units[1])
    plt.ylabel('Y (%s)' % grid_units[0])

    if var_name is not None:
        plt.title('%s (%s)' % (var_name, var_units))

    plt.show()


def imshow_field(field, name, **kwds):
    imshow_grid(field, field[name], var_name=name, var_units=field.units[name],
                **kwds)
