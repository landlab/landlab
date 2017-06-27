#! /usr/bin/env python
"""
Functions to interact with figures that plot Landlab grid data.
"""

from landlab import RasterModelGrid
from pprint import pprint


def query_grid_on_button_press(event, grid):
    """Print node information using an imshow plot.

    This function is triggered when a mouse button is pressed on the matplotlib
    figure connected by event. Coordinates of grid and its node attributes are
    printed only when the event location is within the axes of the figure.
    The node whose attributes are printed is the node at the center of the
    patch that envelopes the event coordinates.

    This function works only with raster model grids.

    Parameters
    ----------
    event : matplotlib event
        Event associated with a figure canvas using mpl_connect().
    grid : RasterModelGrid
        Grid containing the attributes to print.
    """
    if type(grid) is not RasterModelGrid:
        raise TypeError('Only raster grids can be queried.')
   
    if all([event.xdata, event.ydata]):        
        x_pressed = int(round(event.xdata / grid.dx))
        y_pressed = int(round(event.ydata / grid.dy))
        
        attributes = {
                'grid location': {
                        'x_coord': event.xdata,
                        'y_coord': event.ydata
                        },
                'node': {
                        'ID': grid.grid_coords_to_node_id(y_pressed,
                                                          x_pressed),
                        'column': x_pressed,
                        'row': y_pressed
                        }
                }
        
        print('\nGrid query results:\n')
        pprint(attributes, width=1)
