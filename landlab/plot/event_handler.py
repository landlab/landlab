#! /usr/bin/env python
"""
Functions to interact with figures that plot Landlab grid data.
"""

from pprint import pprint

from landlab import RasterModelGrid


def query_grid_on_button_press(event, grid):
    """Print and returns node information using an imshow plot.

    This function is triggered when a mouse button is pressed on the matplotlib
    figure connected by event. Coordinates of grid and its node attributes are
    queried only when the event location is within the axes of the figure.
    The node whose attributes are queried is the node at the center of the
    cell that contains the event coordinates.

    This function only works with raster model grids.

    Parameters
    ----------
    event : matplotlib event
        Event associated with a figure canvas using mpl_connect().
    grid : RasterModelGrid
        Grid containing the attributes to print.

    Returns
    -------
    query_results :
        Dictionary containing grid query results.
    """
    if type(grid) is not RasterModelGrid:
        raise TypeError("Only raster grids can be queried.")

    if all([event.xdata, event.ydata]):
        x_pressed = int(round(event.xdata / grid.dx))
        y_pressed = int(round(event.ydata / grid.dy))

        query_results = {
            "grid location": {"x_coord": event.xdata, "y_coord": event.ydata},
            "node": {
                "ID": grid.grid_coords_to_node_id(y_pressed, x_pressed),
                "column": x_pressed,
                "row": y_pressed,
            },
        }

        print("\nGrid query results:\n")
        pprint(query_results, width=1)

    return query_results
