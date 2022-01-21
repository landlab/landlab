#! /usr/bin/env python
"""Methods to plot data defined on Landlab grids.

Plotting functions
++++++++++++++++++

.. autosummary::

    ~landlab.plot.imshow.imshow_grid
    ~landlab.plot.imshow.imshowhs_grid
    ~landlab.plot.imshow.imshow_grid_at_cell
    ~landlab.plot.imshow.imshow_grid_at_node
    ~landlab.plot.imshow.imshowhs_grid_at_node
"""
import numpy as np
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from landlab.grid.raster import RasterModelGrid
from landlab.plot.event_handler import query_grid_on_button_press

from ..field import FieldError

try:
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.patches import Polygon
except ImportError:
    import warnings

    warnings.warn("matplotlib not found", ImportWarning)


def imshow_grid_at_node(grid, values, **kwds):
    """imshow_grid_at_node(grid, values, plot_name=None, var_name=None,
    var_units=None, grid_units=None, symmetric_cbar=False, cmap='pink',
    limits=(values.min(), values.max()), vmin=values.min(), vmax=values.max(),
    allow_colorbar=True, norm=[linear], shrink=1., color_for_closed='black',
    color_for_background=None, show_elements=False, output=None)

    Prepare a map view of data over all nodes in the grid.

    Data is plotted as cells shaded with the value at the node at its center.
    Outer edges of perimeter cells are extrapolated. Closed elements are
    colored uniformly (default black, overridden with kwd 'color_for_closed');
    other open boundary nodes get their actual values.

    *values* can be a field name, a regular array, or a masked array. If a
    masked array is provided, masked entries will be treated as if they were
    Landlab BC_NODE_IS_CLOSED. Used together with the color_at_closed=None
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
        grid property `axis_units` for this information. If no units are
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
        grid.status_at_node == grid.BC_NODE_IS_CLOSED, values_at_node
    )

    if isinstance(grid, RasterModelGrid):
        shape = grid.shape
    else:
        shape = (-1,)

    _imshow_grid_values(grid, values_at_node.reshape(shape), **kwds)

    if isinstance(values, str):
        plt.title(values)

    plt.gcf().canvas.mpl_connect(
        "button_press_event", lambda event: query_grid_on_button_press(event, grid)
    )


def imshowhs_grid_at_node(grid, values, **kwds):
    """imshowhs_grid_at_node(grid, values, **kwds)

    Prepare a map view of data over all nodes in the grid using a hillshade
    topography map in the background.

    Data is plotted as cells shaded with the value at the node at its center.
    Outer edges of perimeter cells are extrapolated. Closed elements are
    colored uniformly (default black, overridden with kwd 'color_for_closed');
    other open boundary nodes get their actual values.

    *values* can be a field name, a regular array, or a masked array. If a
    masked array is provided, masked entries will be treated as if they were
    Landlab BC_NODE_IS_CLOSED. Used together with the color_at_closed=None
    keyword (i.e., "transparent"), this can allow for construction of overlay
    layers in a figure (e.g., only defining values in a river network, and
    overlaying it on another landscape).

    Use matplotlib functions like xlim, ylim to modify your plot after calling
    :func:`imshowhs_grid`, as desired.

    Node coordinates are printed when a mouse button is pressed on a cell in
    the plot.

    For now, this function only works with regular grids.

    Developed by: Benjamin Campforts

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
    var_name_two : str, optional
        Variable name of second layer, to use as a colorbar label.
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
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    thicks_km : bool, optional
        Display thicks in km instead of m
        Default: False
    allow_colorbar : bool
        If True, include the colorbar.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed nodes (default 'black'). If None, closed
        (or masked) nodes will be transparent.
    color_for_background : color str or other color declaration, or None
        Color to use for closed elements (default None). If None, the
        background will be transparent, and appear white.
    output : None, string, or bool
        If None (or False), the image is sent to the imaging buffer to await
        an explicit call to show() or savefig() from outside this function.
        If a string, the string should be the path to a save location, and the
        filename (with file extension). The function will then call
        plt.savefig([string]) itself. If True, the function will call
        plt.show() itself once plotting is complete.
    fontweight_xlabel : str, optional
        weight of x label. The default is 'bold'.
    fontweight_ylabel : str, optional
        weight of y label. The default is 'bold'.
    plot_type : str, optional
        There are four options:
        * 'DEM': Display a digital elevation map underlain by a shaded relief,
        based on the same DEM ('topographic__elevation')
        * 'Hillshade': Display the shaded relief, of the provided DEM ('topographic__elevation')
        * 'Drape1': Display any kind of provided layer on top of a shaded
        relief provided in the 'topographic__elevation' field
        * 'Drape2': Display two layers on top of a shaded relief provided in
        the 'topographic__elevation' field
        The default is "DEM".
    drape1 : array_like, masked_array
        Node values to plot on top of a hillshade map. The default is None.
    drape2 : array_like, masked_array
        Node values to plot on top of drape1 and a hillshade map. The default is None.
    cmap2 : str
        Name of a colormap for drape 2. The default is None.
    vertical_exa : float, optional
        vertical exageration of hillshade map. The default is None.
    azdeg : float, optional
        azimuth of light source. The default is 315.
    altdeg : float, optional
        elevation of light source. The default is 65.
    thres_drape1 : float, optional
        threshold below which drape1 is made transparant. The default is None.
    alpha : float (0-1), optional
        transparency of DEM/Drape1 . The default is None.
    thres_drape2 : float, optional
        threshold below which drape2 is made transparant. The default is None.
    alpha2 : float (0-1), optional
        transparency of Drape2 . The default is None.
    add_double_colorbar : bool, optional
        add a double colorbar when two drapes are plotted. The default is False.
    plt_contour : bool, optional
        Add contour lines to elevation plot . The default is False.
    contour_nb : int, optional
        number of contour lines. The default is 50.
    default_fontsize : float, optional
        Default font size of plot labels. The default is 10.
    cbar_height : percentage, optional
        height of colorbar in percentage of figure. The default is "5%".
    cbar_width : percentage, optional
        width of colorbar in percentage of figure. The default is "30%".
    cbar_or : str, optional
        orientation of colorbar. The default is "horizontal".
    cbar_loc : str, optional
        location of colorbar. The default is "lower right".
    bbox_to_anchor : vector, optional
        bbox to anchor. The default is (0, 0, 1, 1).
    cbar_ticks_position : str, optional
        location of colorbar ticks (below or on top of the colorbar). The default is "top".
    cbar_ticks_position2 : str, optional
        location of colorbar ticks for colorbar of Drape2 (below or on top of the colorbar). The default is "bottom".
    colorbar_label_y : float, optional
        location of colorbar label with respect to the colorbar in y direction. The default is -40.
    colorbar_label_x : float , optional
        location of colorbar label with respect to the colorbar in x direction. The default is 0.5.
    cbar_tick_size : float, optional
        colorbar tick size. The default is 10.
    cbar_label_color : str, optional
        colorbar tick color. The default is 'black'.
    cbar_label_fontweight : str, optional
        colorbar font weight. The default is 'bold'.
    add_label_bbox : bool, optional
        Add a bbox surrounding the colorbar label. The default is False.

    Returns
    -------
    ax : figure ax
        return ax if output == True.
    """
    if isinstance(values, str):
        values_at_node = grid.at_node[values]
    else:
        values_at_node = values.reshape((-1,))

    if values_at_node.size != grid.number_of_nodes:
        raise ValueError("number of values does not match number of nodes")

    values_at_node = np.ma.masked_where(
        grid.status_at_node == grid.BC_NODE_IS_CLOSED, values_at_node
    )

    if isinstance(grid, RasterModelGrid):
        shape = grid.shape
    else:
        raise NotImplementedError(
            "For now, only RasterModelGrids are supported in the imshowhs functions"
        )

    ax = _imshowhs_grid_values(grid, values_at_node.reshape(shape), **kwds)

    if isinstance(values, str):
        plt.title(values)

    plt.gcf().canvas.mpl_connect(
        "button_press_event", lambda event: query_grid_on_button_press(event, grid)
    )
    # plt.show()
    return ax


def imshow_grid_at_cell(grid, values, **kwds):
    """imshow_grid_at_cell(grid, values, plot_name=None, var_name=None,
    var_units=None, grid_units=None, symmetric_cbar=False, cmap='pink',
    limits=(values.min(), values.max()), vmin=values.min(), vmax=values.max(),
    allow_colorbar=True, colorbar_label=None, norm=[linear], shrink=1.,
    color_for_closed='black', color_for_background=None, show_elements=False,
    output=None)

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
    kwds.setdefault("color_for_closed", None)

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

    values_at_node = np.ma.masked_array(grid.empty(at="node"))
    values_at_node.mask = True
    values_at_node[grid.node_at_cell] = values_at_cell
    values_at_node.mask[grid.node_at_cell] = False

    if isinstance(grid, RasterModelGrid):
        values_at_node = values_at_node.reshape(grid.shape)

    myimage = _imshow_grid_values(grid, values_at_node, **kwds)

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
    cmap = plt.get_cmap(cmap)

    if color_for_closed is not None:
        cmap.set_bad(color=color_for_closed)
    else:
        cmap.set_bad(alpha=0.0)

    if isinstance(grid, RasterModelGrid):
        if values.ndim != 2:
            raise ValueError("values must have ndim == 2")

        y = (
            np.arange(values.shape[0] + 1) * grid.dy
            - grid.dy * 0.5
            + grid.xy_of_lower_left[1]
        )
        x = (
            np.arange(values.shape[1] + 1) * grid.dx
            - grid.dx * 0.5
            + grid.xy_of_lower_left[0]
        )

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

        myimage = plt.pcolormesh(x, y, values, **kwds)
        myimage.set_rasterized(True)
        plt.gca().set_aspect(1.0)
        plt.autoscale(tight=True)

        if allow_colorbar:
            cb = plt.colorbar(norm=norm, shrink=shrink)
            if colorbar_label:
                cb.set_label(colorbar_label)
    else:
        import matplotlib.cm as cmx
        import matplotlib.colors as colors

        if limits is not None:
            (vmin, vmax) = (limits[0], limits[1])
        else:
            if vmin is None:
                vmin = values.min()
            if vmax is None:
                vmax = values.max()
            if symmetric_cbar:
                vmin, vmax = -max(abs(vmin), abs(vmax)), max(abs(vmin), abs(vmax))

        cNorm = colors.Normalize(vmin, vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
        colorVal = scalarMap.to_rgba(values)[grid.node_at_cell]

        patches = []

        for corners in grid.corners_at_cell:
            valid_corners = corners[corners != grid.BAD_INDEX]
            closed_loop_corners = np.concatenate([valid_corners, [valid_corners[0]]])

            x = grid.x_of_corner[closed_loop_corners]
            y = grid.y_of_corner[closed_loop_corners]
            xy = np.vstack((x, y)).T
            patches.append(Polygon(xy, closed=True, fill=True))

        patchcollection = PatchCollection(
            patches, facecolor=colorVal, edgecolor=colorVal
        )

        ax = plt.gca()
        ax.add_collection(patchcollection)

        if show_elements:
            x = grid.x_of_corner[grid.corners_at_face]
            y = grid.y_of_corner[grid.corners_at_face]

            segs = np.dstack((x, y))
            line_segments = LineCollection(segs)
            line_segments.set_color("black")
            ax.add_collection(line_segments)

        ax.set_aspect(1.0)
        ax.set_rasterized(True)

        plt.xlim((np.min(grid.x_of_node), np.max(grid.x_of_node)))
        plt.ylim((np.min(grid.y_of_node), np.max(grid.y_of_node)))

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


def _imshowhs_grid_values(
    grid,
    values,
    plot_name=None,
    var_name=None,
    var_name_two=None,
    var_units=None,
    fontweight_xlabel="bold",
    fontweight_ylabel="bold",
    grid_units=(None, None),
    symmetric_cbar=False,
    cmap="pink",
    limits=None,
    allow_colorbar=True,
    vmin=None,
    vmax=None,
    norm=None,
    thicks_km=False,
    shrink=1.0,
    color_for_closed=None,
    color_for_background=None,
    output=None,
    plot_type="DEM",
    drape1=None,
    drape2=None,
    cmap2=None,
    vertical_exa=None,
    azdeg=315,
    altdeg=65,
    thres_drape1=None,
    alpha=None,
    thres_drape2=None,
    alpha2=None,
    add_double_colorbar=False,
    plt_contour=False,
    contour_nb=50,
    default_fontsize=10,
    cbar_height="5%",
    cbar_width="30%",
    cbar_or="horizontal",
    cbar_loc="lower right",
    bbox_to_anchor=(0, 0, 1, 1),
    cbar_ticks_position="top",
    cbar_ticks_position2="bottom",
    colorbar_label_y=-40,
    colorbar_label_x=0.5,
    cbar_tick_size=10,
    cbar_label_color="black",
    cbar_label_fontweight="bold",
    add_label_bbox=False,
):
    plot_type_options = ["DEM", "Hillshade", "Drape1", "Drape2"]
    if plot_type not in plot_type_options:
        raise ValueError(
            "plot_type should be one of the following: "
            + ", ".join(map(str, plot_type_options))
        )
    if plot_type == "Drape1" and drape1 is None:
        raise ValueError(
            "if plot_type is Drape1, 'drape1' input argument cannot be None. \
                         Provide at least one array with the size of the number of grid nodes as drape1='field_to_be_plotted'"
        )
    if plot_type == "Drape2" and (drape1 is None or drape2 is None):
        raise ValueError(
            "if plot_type is Drape2, 'drape1' and 'drape2' input arguments cannot be None. \
                         Provide an array for both with the size of the number of grid nodes as drape1='field1_to_be_plotted' and drape2='field2_to_be_plotted' "
        )

    # Poperties of bounding box of colorbar label, if used:
    if add_label_bbox:
        bbox_prop = dict(
            boxstyle="round", pad=0.1, facecolor="white", alpha=0.7, edgecolor="white"
        )
    else:
        bbox_prop = None

    cmap = plt.get_cmap(cmap)

    if color_for_closed is not None:
        cmap.set_bad(color=color_for_closed)
    else:
        cmap.set_bad(alpha=0.0)

    if isinstance(grid, RasterModelGrid):
        # somethingToPlot is a flag indicating if any pixels should be plotted.
        somethingToPlot = True

        if values.ndim != 2:
            raise ValueError("values must have ndim == 2")

        y = (
            np.arange(values.shape[0] + 1) * grid.dy
            - grid.dy * 0.5
            + grid.xy_of_lower_left[1]
        )
        x = (
            np.arange(values.shape[1] + 1) * grid.dx
            - grid.dx * 0.5
            + grid.xy_of_lower_left[0]
        )

        ls = LightSource(azdeg=azdeg, altdeg=altdeg)
        if cmap is None:
            cmap = plt.cm.terrain

        dx = x[1] - x[0]
        dy = y[1] - y[0]

        if vertical_exa is not None:
            ve = vertical_exa
        else:
            ve = 3
        extent = np.array([x[0] - dx, x[-1] + dx, y[-1] + dy, y[0] - dy])
        if thicks_km:
            extent /= 1e3

        ax1 = plt.gca()
        if alpha is None:
            alpha = 1
        if alpha2 is None:
            alpha2 = 1
        blend_modes = ["hsv", "overlay", "soft"]
        if plot_type == "DEM":

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

            val = values.data
            rgb = ls.shade(
                val,
                cmap=cmap,
                blend_mode=blend_modes[0],
                vert_exag=ve,
                dx=dx,
                dy=dy,
                fraction=0.4,
            )
            ima = ax1.imshow(rgb, extent=extent, **kwds)

        elif plot_type == "Hillshade":
            ima = plt.imshow(
                ls.hillshade(values, vert_exag=ve, dx=dx, dy=dy),
                cmap="gray",
                extent=extent,
            )
            allow_colorbar = False

        elif plot_type == "Drape1" or plot_type == "Drape2":
            # Process values from first drape
            if isinstance(drape1, str):
                values_at_node_drape1 = grid.at_node[drape1]
            else:
                values_at_node_drape1 = drape1.reshape((-1,))

            if values_at_node_drape1.size != grid.number_of_nodes:
                raise ValueError("number of values does not match number of nodes")

            values_at_node_drape1 = np.ma.masked_where(
                grid.status_at_node == grid.BC_NODE_IS_CLOSED, values_at_node_drape1
            )

            # Add mask if thres_drape1 is given
            if thres_drape1 is not None:
                # check if any value exceeds threshold
                if not np.any(values_at_node_drape1 > thres_drape1):
                    somethingToPlot = False

                values_at_node_drape1 = np.ma.masked_where(
                    values_at_node_drape1 < thres_drape1, values_at_node_drape1
                )

            if isinstance(grid, RasterModelGrid):
                shape = grid.shape
            else:
                shape = (-1,)
            val1 = values_at_node_drape1.reshape(shape)

            kwds = dict(cmap=cmap)
            (kwds["vmin"], kwds["vmax"]) = (val1.min(), val1.max())
            if (limits is None) and ((vmin is None) and (vmax is None)):
                if symmetric_cbar:
                    (var_min, var_max) = (val1.min(), val1.max())
                    limit = max(abs(var_min), abs(var_max))
                    (kwds["vmin"], kwds["vmax"]) = (-limit, limit)
            elif limits is not None:
                (kwds["vmin"], kwds["vmax"]) = (limits[0], limits[1])
            else:
                if vmin is not None:
                    kwds["vmin"] = vmin
                if vmax is not None:
                    kwds["vmax"] = vmax
            plt.imshow(
                ls.hillshade(values, vert_exag=ve, dx=dx, dy=dy),
                cmap="gray",
                extent=extent,
            )
            ima = ax1.imshow(val1, extent=extent, alpha=alpha, **kwds)
            if plt_contour:
                plt.contour(
                    x[0:-1] * 1e-3,
                    y[0:-1] * 1e-3,
                    val1,
                    contour_nb,
                    colors="black",
                    linewidths=0.2,
                )
        if somethingToPlot:
            # To cartezian  coordinates   (not if other layers has to be plotted on top!)
            if plot_type != "Drape2":
                ax1.invert_yaxis()
            plt.xticks(fontsize=default_fontsize)
            plt.yticks(fontsize=default_fontsize)

            # if Drape2, default behavior is to add colorbar of first layer if add_double_colorbar == False
            if allow_colorbar and (
                plot_type == "DEM"
                or plot_type == "Drape1"
                or (plot_type == "Drape2" and not add_double_colorbar)
            ):

                cb_or = cbar_or
                cb_ticks_position = cbar_ticks_position

                axins1 = inset_axes(
                    ax1,
                    width=cbar_width,  # width = 50% of parent_bbox width
                    height=cbar_height,  # height : 5%
                    loc=cbar_loc,
                    bbox_transform=ax1.transAxes,
                    borderpad=0,
                    bbox_to_anchor=bbox_to_anchor,
                )

                maxV = kwds["vmax"]
                minV = kwds["vmin"]
                cb_length = maxV - minV
                if maxV <= 10:
                    cb = plt.colorbar(
                        ima,
                        ax=ax1,
                        cax=axins1,
                        orientation=cb_or,
                        ticks=[
                            np.round(minV + 0.2 * cb_length, 1),
                            np.round(minV + 0.8 * cb_length, 1),
                        ],
                    )
                elif maxV <= 100:
                    cb = plt.colorbar(
                        ima,
                        ax=ax1,
                        cax=axins1,
                        orientation=cb_or,
                        ticks=[
                            np.round(minV + 0.2 * cb_length, 0),
                            np.round(minV + 0.8 * cb_length, 0),
                        ],
                    )
                else:
                    cb = plt.colorbar(
                        ima,
                        ax=ax1,
                        cax=axins1,
                        orientation=cb_or,
                        ticks=[
                            np.round(0.1 * (minV + 0.2 * cb_length)) * 10,
                            np.round(0.1 * (minV + 0.8 * cb_length)) * 10,
                        ],
                    )
                axins1.xaxis.set_ticks_position(cb_ticks_position)
                cb.ax.tick_params(
                    labelsize=cbar_tick_size,
                    color=cbar_label_color,
                    labelcolor=cbar_label_color,
                )

                # if colorbar_label:
                #     cb.set_label(colorbar_label, rotation=270)
                #     # ax1.xaxis.set_label_coords(0,2.5)

            if plot_type == "Drape2":

                # Process values from first drape

                if isinstance(drape2, str):
                    values_at_node_drape2 = grid.at_node[drape2]
                else:
                    values_at_node_drape2 = drape2.reshape((-1,))

                if values_at_node_drape2.size != grid.number_of_nodes:
                    raise ValueError("number of values does not match number of nodes")

                values_at_node_drape2 = np.ma.masked_where(
                    grid.status_at_node == grid.BC_NODE_IS_CLOSED, values_at_node_drape2
                )

                # Add mask if thres_drape1 is given
                if thres_drape2 is not None:
                    values_at_node_drape2 = np.ma.masked_where(
                        values_at_node_drape2 < thres_drape2, values_at_node_drape2
                    )

                if isinstance(grid, RasterModelGrid):
                    shape = grid.shape
                else:
                    shape = (-1,)
                val2 = values_at_node_drape2.reshape(shape)

                if cmap2 is None:
                    cmap2 = plt.cm.terrain
                kwds = dict(cmap=cmap2)
                (kwds["vmin"], kwds["vmax"]) = (val2.min(), val2.max())
                if (limits is None) and ((vmin is None) and (vmax is None)):
                    if symmetric_cbar:
                        (var_min, var_max) = (val2.min(), val2.max())
                        limit = max(abs(var_min), abs(var_max))
                        (kwds["vmin"], kwds["vmax"]) = (-limit, limit)
                elif limits is not None:
                    (kwds["vmin"], kwds["vmax"]) = (limits[0], limits[1])
                else:
                    if vmin is not None:
                        kwds["vmin"] = vmin
                    if vmax is not None:
                        kwds["vmax"] = vmax

                ima2 = ax1.imshow(val2, extent=extent, alpha=alpha2, **kwds)
                ax1.invert_yaxis()

                # Add colorbars
                if add_double_colorbar:
                    axins1 = inset_axes(
                        ax1,
                        width=cbar_width,  # width = 50% of parent_bbox width
                        height=cbar_height,  # height : 5%
                        loc=cbar_loc,
                        bbox_to_anchor=(-0.005, 0.25, 1, 1),
                        bbox_transform=ax1.transAxes,
                        borderpad=0,
                    )

                    cb_or = cbar_or
                    cb_ticks_position = cbar_ticks_position
                    maxV = np.max(val1)
                    minV = np.min(val1)
                    cb_length = maxV - minV
                    if maxV <= 10:
                        cb = plt.colorbar(
                            ima,
                            ax=ax1,
                            cax=axins1,
                            orientation=cb_or,
                            ticks=[
                                np.round(minV + 0.2 * cb_length, 1),
                                np.round(minV + 0.8 * cb_length, 1),
                            ],
                        )
                    elif maxV <= 100:
                        cb = plt.colorbar(
                            ima,
                            ax=ax1,
                            cax=axins1,
                            orientation=cb_or,
                            ticks=[
                                np.round(minV + 0.2 * cb_length, 0),
                                np.round(minV + 0.8 * cb_length, 0),
                            ],
                        )
                    else:
                        cb = plt.colorbar(
                            ima,
                            ax=ax1,
                            cax=axins1,
                            orientation=cb_or,
                            ticks=[
                                np.round(0.1 * (minV + 0.2 * cb_length)) * 10,
                                np.round(0.1 * (minV + 0.8 * cb_length)) * 10,
                            ],
                        )
                    cb.ax.tick_params(
                        labelsize=cbar_tick_size,
                        color=cbar_label_color,
                        labelcolor=cbar_label_color,
                    )
                    axins1.xaxis.set_ticks_position(cb_ticks_position)

                    axins1.set_xlabel(
                        var_name,
                        usetex=True,
                        fontsize=default_fontsize,
                        rotation=0,
                        color=cbar_label_color,
                        fontweight=cbar_label_fontweight,
                        bbox=bbox_prop,
                    )
                    axins1.xaxis.set_label_coords(0.5, 3.5)

                    axins2 = inset_axes(
                        ax1,
                        width=cbar_width,  # width = 50% of parent_bbox width
                        height=cbar_height,  # height : 5%
                        loc=cbar_loc,
                        bbox_to_anchor=(-0.005, 0.15, 1, 1),
                        bbox_transform=ax1.transAxes,
                        borderpad=0,
                    )
                    cb_or = cbar_or
                    cb_ticks_position = cbar_ticks_position2
                    maxV = np.max(val2)
                    minV = np.min(val2)
                    cb_length = maxV - minV
                    if maxV <= 10:
                        cb = plt.colorbar(
                            ima2,
                            ax=ax1,
                            cax=axins2,
                            orientation=cb_or,
                            ticks=[
                                np.round(minV + 0.2 * cb_length, 1),
                                np.round(minV + 0.8 * cb_length, 1),
                            ],
                        )
                    elif maxV <= 100:
                        cb = plt.colorbar(
                            ima2,
                            ax=ax1,
                            cax=axins2,
                            orientation=cb_or,
                            ticks=[
                                np.round(minV + 0.2 * cb_length, 0),
                                np.round(minV + 0.8 * cb_length, 0),
                            ],
                        )
                    else:
                        cb = plt.colorbar(
                            ima2,
                            ax=ax1,
                            cax=axins2,
                            orientation=cb_or,
                            ticks=[
                                np.round(0.1 * (minV + 0.2 * cb_length)) * 10,
                                np.round(0.1 * (minV + 0.8 * cb_length)) * 10,
                            ],
                        )
                    cb.ax.tick_params(
                        labelsize=cbar_tick_size,
                        color=cbar_label_color,
                        labelcolor=cbar_label_color,
                    )

                    axins2.xaxis.set_ticks_position(cb_ticks_position)
                    axins2.set_xlabel(
                        var_name_two,
                        usetex=True,
                        fontsize=default_fontsize,
                        rotation=0,
                        color=cbar_label_color,
                        fontweight=cbar_label_fontweight,
                        bbox=bbox_prop,
                    )
                    axins2.xaxis.set_label_coords(0.5, -1.75)

    if grid_units[1] is None and grid_units[0] is None:
        grid_units = grid.axis_units
        if grid_units[1] == "-" and grid_units[0] == "-":
            ax1.set_xlabel(
                "Easting", fontweight=fontweight_xlabel, fontsize=default_fontsize
            )
            ax1.set_ylabel(
                "Northing", fontweight=fontweight_ylabel, fontsize=default_fontsize
            )
        else:
            ax1.set_xlabel(
                "Easting, %s" % grid_units[1],
                fontweight=fontweight_xlabel,
                fontsize=default_fontsize,
            )
            ax1.set_ylabel(
                "Northing, %s" % grid_units[1],
                fontweight=fontweight_ylabel,
                fontsize=default_fontsize,
            )
    else:
        ax1.set_xlabel(
            "Easting, %s" % grid_units[1],
            fontweight=fontweight_xlabel,
            fontsize=default_fontsize,
        )
        ax1.set_ylabel(
            "Northing, %s" % grid_units[1],
            fontweight=fontweight_ylabel,
            fontsize=default_fontsize,
        )

    if plot_name is not None:
        plt.title("%s" % (plot_name))

    if (
        somethingToPlot
        and (var_name is not None or var_units is not None)
        and plot_type != "Drape2"
    ):
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
        if allow_colorbar:
            cb.set_label(
                colorbar_label,
                fontsize=default_fontsize,
                labelpad=colorbar_label_y,
                color=cbar_label_color,
                x=colorbar_label_x,
                fontweight=cbar_label_fontweight,
                bbox=bbox_prop,
            )

    if color_for_background is not None:
        plt.gca().set_facecolor(color_for_background)

    if output is not None:
        if type(output) is str:
            plt.savefig(output)
            plt.clf()
        elif output:
            plt.show()

    return ax1


def imshow_grid(grid, values, **kwds):
    """imshow_grid(grid, values, plot_name=None, var_name=None, var_units=None,
    grid_units=None, symmetric_cbar=False, cmap='pink', limits=(values.min(),
    values.max()), vmin=values.min(), vmax=values.max(), allow_colorbar=True,
    colorbar_label=None, norm=[linear], shrink=1., color_for_closed='black',
    show_elements=False, color_for_background=None)

    Prepare a map view of data over all nodes or cells in the grid.

    Data is plotted as colored cells. If at='node', the surrounding cell is
    shaded with the value at the node at its center. If at='cell', the cell
    is shaded with its own value. Outer edges of perimeter cells are
    extrapolated. Closed elements are colored uniformly (default black,
    overridden with kwd 'color_for_closed'); other open boundary nodes get
    their actual values.

    *values* can be a field name, a regular array, or a masked array. If a
    masked array is provided, masked entries will be treated as if they were
    Landlab BC_NODE_IS_CLOSED. Used together with the color_for_closed=None
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


def imshowhs_grid(grid, values, **kwds):
    """imshowhs_grid(grid, values, **kwds)

    Prepare a map view of data over all nodes in the grid using a hillshade
    topography map in the background.

    Data is plotted as cells shaded with the value at the node at its center.
    Outer edges of perimeter cells are extrapolated. Closed elements are
    colored uniformly (default black, overridden with kwd 'color_for_closed');
    other open boundary nodes get their actual values.

    *values* can be a field name, a regular array, or a masked array. If a
    masked array is provided, masked entries will be treated as if they were
    Landlab BC_NODE_IS_CLOSED. Used together with the color_at_closed=None
    keyword (i.e., "transparent"), this can allow for construction of overlay
    layers in a figure (e.g., only defining values in a river network, and
    overlaying it on another landscape).

    Use matplotlib functions like xlim, ylim to modify your plot after calling
    :func:`imshowhs_grid`, as desired.

    Node coordinates are printed when a mouse button is pressed on a cell in
    the plot.

    For now, this function only works with regular grids.

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
    var_name_two : str, optional
        Variable name of second layer, to use as a colorbar label.
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
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    thicks_km : bool, optional
        Display thicks in km instead of m
    allow_colorbar : bool
        If True, include the colorbar.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed nodes (default 'black'). If None, closed
        (or masked) nodes will be transparent.
    color_for_background : color str or other color declaration, or None
        Color to use for closed elements (default None). If None, the
        background will be transparent, and appear white.
    output : None, string, or bool
        If None (or False), the image is sent to the imaging buffer to await
        an explicit call to show() or savefig() from outside this function.
        If a string, the string should be the path to a save location, and the
        filename (with file extension). The function will then call
        plt.savefig([string]) itself. If True, the function will call
        plt.show() itself once plotting is complete.
    fontweight_xlabel : str, optional
        weight of x label. The default is 'bold'.
    fontweight_ylabel : str, optional
        weight of y label. The default is 'bold'.
    plot_type : str, optional
        The type of plot that will be plotted.
        There are four options:
        * 'DEM': Display a digital elevation map underlain by a shaded relief,
        based on the same DEM ('topographic__elevation')
        * 'Hillshade': Display the shaded relief, of the provided DEM ('topographic__elevation')
        * 'Drape1': Display any kind of provided layer on top of a shaded
        relief provided in the 'topographic__elevation' field
        * 'Drape2': Display two layers on top of a shaded relief provided in
        the 'topographic__elevation' field
        The default is "DEM".
    drape1 : array_like, masked_array
        Node values to plot on top of a hillshade map. The default is None.
    drape2 : array_like, masked_array
        Node values to plot on top of drape1 and a hillshade map. The default is None.
    cmap2 : str
        Name of a colormap for drape 2. The default is None.
    vertical_exa : float, optional
        vertical exageration of hillshade map. The default is None.
    azdeg : float, optional
        azimuth of light source. The default is 315.
    altdeg : float, optional
        elevation of light source. The default is 65.
    thres_drape1 : float, optional
        threshold below which drape1 is made transparant. The default is None.
    alpha : float (0-1), optional
        transparency of DEM/Drape1 . The default is None.
    thres_drape2 : float, optional
        threshold below which drape2 is made transparant. The default is None.
    alpha2 : float (0-1), optional
        transparency of Drape2 . The default is None.
    add_double_colorbar : bool, optional
        add a double colorbar when two drapes are plotted. The default is False.
    plt_contour : bool, optional
        Add contour lines to elevation plot . The default is False.
    contour_nb : int, optional
        number of contour lines. The default is 50.
    default_fontsize : float, optional
        Default font size of plot labels. The default is 10.
    cbar_height : percentage, optional
        height of colorbar in percentage of figure. The default is "5%".
    cbar_width : percentage, optional
        width of colorbar in percentage of figure. The default is "30%".
    cbar_or : str, optional
        orientation of colorbar. The default is "horizontal".
    cbar_loc : str, optional
        location of colorbar. The default is "lower right".
    bbox_to_anchor : vector, optional
        bbox to anchor. The default is (0, 0, 1, 1).
    cbar_ticks_position : str, optional
        location of colorbar ticks (below or on top of the colorbar). The default is "top".
    cbar_ticks_position2 : str, optional
        location of colorbar ticks for colorbar of Drape2 (below or on top of the colorbar). The default is "bottom".
    colorbar_label_y : float, optional
        location of colorbar label with respect to the colorbar in y direction. The default is -40.
    colorbar_label_x : float , optional
        location of colorbar label with respect to the colorbar in x direction. The default is 0.5.
    cbar_tick_size : float, optional
        colorbar tick size. The default is 10.
    cbar_label_color : str, optional
        colorbar tick color. The default is 'black'.
    cbar_label_fontweight : str, optional
        colorbar font weight. The default is 'bold'.
    add_label_bbox : bool, optional
        Add a bbox surrounding the colorbar label. The default is False.

    Returns
    -------
    ax : figure ax
        return ax if output == True.
    """
    values_at = kwds.pop("values_at", "node")
    values_at = kwds.pop("at", values_at)

    if isinstance(values, str):
        values = grid.field_values(values_at, values)

    if values_at == "node":
        ax = imshowhs_grid_at_node(grid, values, **kwds)
    elif values_at == "cell":
        raise NotImplementedError(
            "For now, only values at nodes can be displayed using the in the imshowhs functions"
        )
    else:
        raise TypeError("value location %s not understood" % values_at)

    return ax
