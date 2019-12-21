#! /usr/bin/env python
"""
Methods to plot data defined on Landlab grids.

Plotting functions
++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.plot.imshow.imshow_grid
    ~landlab.plot.imshow.imshow_grid_at_cell
    ~landlab.plot.imshow.imshow_grid_at_node
"""


import inspect

import numpy as np

from landlab.field.scalar_data_fields import FieldError
from landlab.grid import CLOSED_BOUNDARY
from landlab.grid.raster import RasterModelGrid
from landlab.grid.voronoi import VoronoiDelaunayGrid
from landlab.plot.event_handler import query_grid_on_button_press

try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings

    warnings.warn("matplotlib not found", ImportWarning)


def imshow_grid_at_node(grid, values, **kwds):
    """imshow_grid_at_node(grid, values, plot_name=None, var_name=None,
                           var_units=None, grid_units=None,
                           symmetric_cbar=False, cmap='pink',
                           limits=(values.min(), values.max()),
                           vmin=values.min(), vmax=values.max(),
                           allow_colorbar=True,
                           norm=[linear], shrink=1.,
                           color_for_closed='black',
                           color_for_background=None,
                           show_elements=False, output=None)

    Prepare a map view of data over all nodes in the grid.

    Data is plotted as cells shaded with the value at the node at its center.
    Outer edges of perimeter cells are extrapolated. Closed elements are
    colored uniformly (default black, overridden with kwd 'color_for_closed');
    other open boundary nodes get their actual values.

    *values* can be a field name, a regular array, or a masked array. If a
    masked array is provided, masked entries will be treated as if they were
    Landlab CLOSED_BOUNDARYs. Used together with the color_at_closed=None
    keyword (i.e., "transparent"), this can allow for construction of overlay
    layers in a figure (e.g., only defining values in a river network, and
    overlaying it on another landscape).

    Use matplotlib functions like xlim, ylim to modify your plot after calling
    :func:`imshow_grid`, as desired.

    Node coordinates are printed when a mouse button is pressed on a cell in
    the plot.

    This function happily works with both regular and irregular grids.

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like, masked_array, or str
        Node values, or a field name as a string from which to draw the data.
    plot_name : str, optional
        String to put as the plot title.
    var_name : str, optional
        Variable name, to use as a colorbar label.
    var_units : str, optional
        Units for the variable being plotted, for the colorbar.
    grid_units : tuple of str, optional
        Units for y, and x dimensions. If None, component will look to the
        gri property `axis_units` for this information. If no units are
        specified there, no entry is made.
    symmetric_cbar : bool
        Make the colormap symetric about 0.
    cmap : str
        Name of a colormap
    limits : tuple of float
        Minimum and maximum of the colorbar.
    vmin, vmax: floats
        Alternatives to limits.
    allow_colorbar : bool
        If True, include the colorbar.
    colorbar_label : str or None
        The string with which to label the colorbar.
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed nodes (default 'black'). If None, closed
        (or masked) nodes will be transparent.
    color_for_background : color str or other color declaration, or None
        Color to use for closed elements (default None). If None, the
        background will be transparent, and appear white.
    show_elements : bool
        If True, and grid is a Voronoi, the faces will be plotted in black
        along with just the colour of the cell, defining the cell outlines
        (defaults False).
    output : None, string, or bool
        If None (or False), the image is sent to the imaging buffer to await
        an explicit call to show() or savefig() from outside this function.
        If a string, the string should be the path to a save location, and the
        filename (with file extension). The function will then call
        plt.savefig([string]) itself. If True, the function will call
        plt.show() itself once plotting is complete.
    """
    if isinstance(values, str):
        values_at_node = grid.at_node[values]
    else:
        values_at_node = values.reshape((-1,))

    if values_at_node.size != grid.number_of_nodes:
        raise ValueError("number of values does not match number of nodes")

    values_at_node = np.ma.masked_where(
        grid.status_at_node == CLOSED_BOUNDARY, values_at_node
    )

    try:
        shape = grid.shape
    except AttributeError:
        shape = (-1,)

    _imshow_grid_values(grid, values_at_node.reshape(shape), **kwds)

    if isinstance(values, str):
        plt.title(values)

    plt.gcf().canvas.mpl_connect(
        "button_press_event", lambda event: query_grid_on_button_press(event, grid)
    )


def imshow_grid_at_cell(grid, values, **kwds):
    """imshow_grid_at_cell(grid, values, plot_name=None, var_name=None,
                           var_units=None, grid_units=None,
                           symmetric_cbar=False, cmap='pink',
                           limits=(values.min(), values.max()),
                           vmin=values.min(), vmax=values.max(),
                           allow_colorbar=True, colorbar_label=None,
                           norm=[linear], shrink=1.,
                           color_for_closed='black',
                           color_for_background=None,
                           show_elements=False, output=None)

    Map view of grid data over all grid cells.

    Prepares a map view of data over all cells in the grid.
    Method can take any of the same ``**kwds`` as :func:`imshow_grid_at_node`.

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like, masked_array, or str
        Values at the cells on the grid. Alternatively, can be a field name
        (string) from which to draw the data from the grid.
    plot_name : str, optional
        String to put as the plot title.
    var_name : str, optional
        Variable name, to use as a colorbar label.
    var_units : str, optional
        Units for the variable being plotted, for the colorbar.
    grid_units : tuple of str, optional
        Units for y, and x dimensions. If None, component will look to the
        gri property `axis_units` for this information. If no units are
        specified there, no entry is made.
    symmetric_cbar : bool
        Make the colormap symetric about 0.
    cmap : str
        Name of a colormap
    limits : tuple of float
        Minimum and maximum of the colorbar.
    vmin, vmax: floats
        Alternatives to limits.
    allow_colorbar : bool
        If True, include the colorbar.
    colorbar_label : str or None
        The string with which to label the colorbar.
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed elements (default 'black'). If None, closed
        (or masked) elements will be transparent.
    color_for_background : color str or other color declaration, or None
        Color to use for closed elements (default None). If None, the
        background will be transparent, and appear white.
    show_elements : bool
        If True, and grid is a Voronoi, the faces will be plotted in black
        along with just the colour of the cell, defining the cell outlines
        (defaults False).
    output : None, string, or bool
        If None (or False), the image is sent to the imaging buffer to await
        an explicit call to show() or savefig() from outside this function.
        If a string, the string should be the path to a save location, and the
        filename (with file extension). The function will then call
        plt.savefig([string]) itself. If True, the function will call
        plt.show() itself once plotting is complete.

    Raises
    ------
    ValueError
        If input grid is not uniform rectilinear.
    """
    if isinstance(values, str):
        try:
            values_at_cell = grid.at_cell[values]
        except FieldError:
            values_at_cell = grid.at_node[values]
    else:
        values_at_cell = values

    if values_at_cell.size == grid.number_of_nodes:
        values_at_cell = values_at_cell[grid.node_at_cell]

    if values_at_cell.size != grid.number_of_cells:
        raise ValueError(
            "number of values must match number of cells or " "number of nodes"
        )

    values_at_cell = np.ma.asarray(values_at_cell)
    values_at_cell.mask = True
    values_at_cell.mask[grid.core_cells] = False

    myimage = _imshow_grid_values(
        grid, values_at_cell.reshape(grid.cell_grid_shape), **kwds
    )

    if isinstance(values, str):
        plt.title(values)

    return myimage


def _imshow_grid_values(
    grid,
    values,
    plot_name=None,
    var_name=None,
    var_units=None,
    grid_units=(None, None),
    symmetric_cbar=False,
    cmap="pink",
    limits=None,
    colorbar_label=None,
    allow_colorbar=True,
    vmin=None,
    vmax=None,
    norm=None,
    shrink=1.0,
    color_for_closed="black",
    color_for_background=None,
    show_elements=False,
    output=None,
):

    gridtypes = inspect.getmro(grid.__class__)

    cmap = plt.get_cmap(cmap)
    if color_for_closed is not None:
        cmap.set_bad(color=color_for_closed)
    else:
        cmap.set_bad(alpha=0.0)

    if isinstance(grid, RasterModelGrid):
        if values.ndim != 2:
            raise ValueError("values must have ndim == 2")

        y = np.arange(values.shape[0] + 1) * grid.dy - grid.dy * 0.5
        x = np.arange(values.shape[1] + 1) * grid.dx - grid.dx * 0.5

        kwds = dict(cmap=cmap)
        (kwds["vmin"], kwds["vmax"]) = (values.min(), values.max())
        if (limits is None) and ((vmin is None) and (vmax is None)):
            if symmetric_cbar:
                (var_min, var_max) = (values.min(), values.max())
                limit = max(abs(var_min), abs(var_max))
                (kwds["vmin"], kwds["vmax"]) = (-limit, limit)
        elif limits is not None:
            (kwds["vmin"], kwds["vmax"]) = (limits[0], limits[1])
        else:
            if vmin is not None:
                kwds["vmin"] = vmin
            if vmax is not None:
                kwds["vmax"] = vmax

        if np.isclose(grid.dx, grid.dy):
            if values.size == grid.number_of_nodes:
                myimage = plt.imshow(
                    values.reshape(grid.shape),
                    origin="lower",
                    extent=(x[0], x[-1], y[0], y[-1]),
                    **kwds
                )
            else:  # this is a cell grid, and has been reshaped already...
                myimage = plt.imshow(
                    values, origin="lower", extent=(x[0], x[-1], y[0], y[-1]), **kwds
                )
        myimage = plt.pcolormesh(x, y, values, **kwds)
        myimage.set_rasterized(True)
        plt.gca().set_aspect(1.0)
        plt.autoscale(tight=True)

        if allow_colorbar:
            cb = plt.colorbar(norm=norm, shrink=shrink)
            if colorbar_label:
                cb.set_label(colorbar_label)
    elif VoronoiDelaunayGrid in gridtypes:
        # This is still very much ad-hoc, and needs prettifying.
        # We should save the modifications needed to plot color all the way
        # to the diagram edge *into* the grid, for faster plotting.
        # (see http://stackoverflow.com/questions/20515554/...
        # colorize-voronoi-diagram)
        # (This technique is not implemented yet)
        from scipy.spatial import voronoi_plot_2d
        import matplotlib.colors as colors
        import matplotlib.cm as cmx

        cm = plt.get_cmap(cmap)

        if (limits is None) and ((vmin is None) and (vmax is None)):
            # only want to work with NOT CLOSED nodes
            open_nodes = grid.status_at_node != 4
            if symmetric_cbar:
                (var_min, var_max) = (
                    values.flat[open_nodes].min(),
                    values.flat[open_nodes].max(),
                )
                limit = max(abs(var_min), abs(var_max))
                (vmin, vmax) = (-limit, limit)
            else:
                (vmin, vmax) = (
                    values.flat[open_nodes].min(),
                    values.flat[open_nodes].max(),
                )
        elif limits is not None:
            (vmin, vmax) = (limits[0], limits[1])
        else:
            open_nodes = grid.status_at_node != 4
            if vmin is None:
                vmin = values.flat[open_nodes].min()
            if vmax is None:
                vmax = values.flat[open_nodes].max()

        cNorm = colors.Normalize(vmin, vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        colorVal = scalarMap.to_rgba(values)

        if show_elements:
            myimage = voronoi_plot_2d(grid.vor, show_vertices=False, show_points=False)

        # show_points to be supported in scipy0.18, but harmless for now
        mycolors = (i for i in colorVal)
        for order in grid.vor.point_region:
            region = grid.vor.regions[order]
            colortouse = next(mycolors)
            if -1 not in region:
                polygon = [grid.vor.vertices[i] for i in region]
                plt.fill(*zip(*polygon), color=colortouse)

        plt.gca().set_aspect(1.0)
        # plt.autoscale(tight=True)
        # Tempting though it is to move the boundary outboard of the outermost
        # nodes (e.g., to the outermost corners), this is a bad idea, as the
        # outermost cells tend to have highly elongated shapes which make the
        # plot look stupid
        plt.xlim((np.min(grid.node_x), np.max(grid.node_x)))
        plt.ylim((np.min(grid.node_y), np.max(grid.node_y)))

        scalarMap.set_array(values)
        if allow_colorbar:
            cb = plt.colorbar(scalarMap, shrink=shrink)

    if grid_units[1] is None and grid_units[0] is None:
        grid_units = grid.axis_units
        if grid_units[1] == "-" and grid_units[0] == "-":
            plt.xlabel("X")
            plt.ylabel("Y")
        else:
            plt.xlabel("X (%s)" % grid_units[1])
            plt.ylabel("Y (%s)" % grid_units[0])
    else:
        plt.xlabel("X (%s)" % grid_units[1])
        plt.ylabel("Y (%s)" % grid_units[0])

    if plot_name is not None:
        plt.title("%s" % (plot_name))

    if var_name is not None or var_units is not None:
        if var_name is not None:
            assert type(var_name) is str
            if var_units is not None:
                assert type(var_units) is str
                colorbar_label = var_name + " (" + var_units + ")"
            else:
                colorbar_label = var_name
        else:
            assert type(var_units) is str
            colorbar_label = "(" + var_units + ")"
        assert type(colorbar_label) is str
        assert allow_colorbar
        cb.set_label(colorbar_label)

    if color_for_background is not None:
        plt.gca().set_facecolor(color_for_background)

    if output is not None:
        if type(output) is str:
            plt.savefig(output)
            plt.clf()
        elif output:
            plt.show()


def imshow_grid(grid, values, **kwds):
    """imshow_grid(grid, values, plot_name=None, var_name=None,
                   var_units=None, grid_units=None,
                   symmetric_cbar=False, cmap='pink',
                   limits=(values.min(), values.max()),
                   vmin=values.min(), vmax=values.max(),
                   allow_colorbar=True, colorbar_label=None,
                   norm=[linear], shrink=1.,
                   color_for_closed='black',
                   color_for_background=None,
                   show_elements=False)

    Prepare a map view of data over all nodes or cells in the grid.

    Data is plotted as colored cells. If at='node', the surrounding cell is
    shaded with the value at the node at its center. If at='cell', the cell
    is shaded with its own value. Outer edges of perimeter cells are
    extrapolated. Closed elements are colored uniformly (default black,
    overridden with kwd 'color_for_closed'); other open boundary nodes get
    their actual values.

    *values* can be a field name, a regular array, or a masked array. If a
    masked array is provided, masked entries will be treated as if they were
    Landlab CLOSED_BOUNDARYs. Used together with the color_for_closed=None
    keyword (i.e., "transparent"), this can allow for construction of overlay
    layers in a figure (e.g., only defining values in a river network, and
    overlaying it on another landscape).

    Use matplotlib functions like xlim, ylim to modify your plot after calling
    :func:`imshow_grid`, as desired.

    This function happily works with both regular and irregular grids.

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like, masked_array, or str
        Node or cell values, or a field name as a string from which to draw
        the data.
    at : str, {'node', 'cell'}
        Tells plotter where values are defined.
    plot_name : str, optional
        String to put as the plot title.
    var_name : str, optional
        Variable name, to use as a colorbar label.
    var_units : str, optional
        Units for the variable being plotted, for the colorbar.
    grid_units : tuple of str, optional
        Units for y, and x dimensions. If None, component will look to the
        gri property `axis_units` for this information. If no units are
        specified there, no entry is made.
    symmetric_cbar : bool
        Make the colormap symetric about 0.
    cmap : str
        Name of a colormap
    limits : tuple of float
        Minimum and maximum of the colorbar.
    vmin, vmax: floats
        Alternatives to limits.
    allow_colorbar : bool
        If True, include the colorbar.
    colorbar_label : str or None
        The string with which to label the colorbar.
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed elements (default 'black'). If None, closed
        (or masked) elements will be transparent.
    color_for_background : color str or other color declaration, or None
        Color to use for closed elements (default None). If None, the
        background will be transparent, and appear white.
    show_elements : bool
        If True, and grid is a Voronoi, the faces will be plotted in black
        along with just the colour of the cell, defining the cell outlines
        (defaults False).
    output : None, string, or bool
        If None (or False), the image is sent to the imaging buffer to await
        an explicit call to show() or savefig() from outside this function.
        If a string, the string should be the path to a save location, and the
        filename (with file extension). The function will then call
        plt.savefig([string]) itself. If True, the function will call
        plt.show() itself once plotting is complete.
    """
    show = kwds.pop("show", False)
    values_at = kwds.pop("values_at", "node")
    values_at = kwds.pop("at", values_at)

    if isinstance(values, str):
        values = grid.field_values(values_at, values)

    if isinstance(values, str):
        values = grid.field_values(values_at, values)

    if values_at == "node":
        imshow_grid_at_node(grid, values, **kwds)
    elif values_at == "cell":
        imshow_grid_at_cell(grid, values, **kwds)
    else:
        raise TypeError("value location %s not understood" % values_at)

    # retained for backwards compatibility:
    if show:
        plt.show()
