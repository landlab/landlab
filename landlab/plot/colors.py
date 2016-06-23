# -*- coding: utf-8 -*-
"""
colors.py

Created on Mon Jan 18 13:28:17 2016

@author: gtucker
"""

from matplotlib.colors import LinearSegmentedColormap


def water_colormap():
    """Return matplotlib colormap with 'water' theme."""
    cdict = { 'red'   : ((0.0, 0.0, 169.0/255.0),
                         (1.0,  38.0/255.0, 1.0)),
              'green' : ((0.0, 0.0, 222.0/255.0),
                         (1.0, 39.0/255.0, 1.0)),
              'blue'  : ((0.0, 0.0, 242.0/255.0),
                         (1.0, 23.0/255.0, 1.0)) }
    return LinearSegmentedColormap('landlab_water', cdict)


def earth_colormap():
    """Return matplotlib colormap with 'earth' theme."""
    cdict = { 'red'   : ((0.0, 0.0, 252.0/255.0),
                         (1.0,  33.0/255.0, 1.0)),
              'green' : ((0.0, 0.0, 237.0/255.0),
                         (1.0, 38.0/255.0, 1.0)),
              'blue'  : ((0.0, 0.0, 179.0/255.0),
                         (1.0, 24.0/255.0, 1.0)) }
    return LinearSegmentedColormap('landlab_earth', cdict)
    

def colormap(name):
    """Return named Landlab colormap as a matplotlib colormap.

    Parameters
    ----------
    name : str
        Name of colormap

    Currently available maps are:
        'water': black to light blue
        'earth': dark olive to light sand color
    """
    colormap_fns = { 'water' : water_colormap(),
                     'earth' : earth_colormap() }
    try:
        return colormap_fns[name]
    except:
        print('Warning: colormap "' + name + '" does not exist')
        return None