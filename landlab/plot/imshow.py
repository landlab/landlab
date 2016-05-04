#! /usr/bin/env python

import numpy as np
import inspect
from landlab.field.scalar_data_fields import FieldError
try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)
from landlab.grid import CLOSED_BOUNDARY
from landlab.grid.raster import RasterModelGrid
from landlab.grid.voronoi import VoronoiDelaunayGrid
from landlab.utils.decorators import deprecated


def imshow_grid_at_node(grid, values, **kwds):
    """Prepare a map view of data over all nodes in the grid.

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

    This function happily works with both regular and irregular grids.

    Construction ::

        imshow_grid_at_node(grid, values, var_name=[field_name],
                            var_units=[field_units], grid_units=None,
                            symmetric_cbar=False, cmap='pink',
                            limits=(values.min(), values.max()),
                            vmin=values.min(), vmax=values.max(),
                            allow_colorbar=True, norm=[linear], shrink=1.,
                            color_for_closed='black', show_elements=False)

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like, masked_array, or str
        Node values, or a field name as a string from which to draw the data.
    var_name : str, optional
        Name of the variable to put in plot title.
    var_units : str, optional
        Units for the variable being plotted.
    grid_units : tuple of str, optional
        Units for y, and x dimensions.
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
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed nodes (default 'black'). If None, closed
        (or masked) nodes will be transparent.
    show_elements : bool
        If True, and grid is a Voronoi, extra grid elements (nodes, faces,
        corners) will be plotted along with just the colour of the cell
        (defaults False).
    """
    if isinstance(values, str):
        values_at_node = grid.at_node[values]
    else:
        values_at_node = values

    if values_at_node.size != grid.number_of_nodes:
        raise ValueError('number of values does not match number of nodes')

    values_at_node = np.ma.masked_where(
        grid.status_at_node == CLOSED_BOUNDARY, values_at_node)

    try:
        shape = grid.shape
    except AttributeError:
        shape = (-1, )

    _imshow_grid_values(grid, values_at_node.reshape(shape), **kwds)

    if isinstance(values, str):
        plt.title(values)


@deprecated(use='imshow_grid_at_node', version='0.5')
def imshow_node_grid(grid, values, **kwds):
    imshow_grid_at_node(grid, values, **kwds)


def imshow_grid_at_cell(grid, values, **kwds):
    """Map view of grid data over all grid cells.

    Prepares a map view of data over all cells in the grid.
    Method can take any of the same ``**kwds`` as :func:`imshow_grid_at_node`.

    Construction ::

        imshow_grid_at_cell(grid, values, var_name=[field_name],
                            var_units=[field_units], grid_units=None,
                            symmetric_cbar=False, cmap='pink',
                            limits=(values.min(), values.max()),
                            vmin=values.min(), vmax=values.max(),
                            allow_colorbar=True, norm=[linear], shrink=1.,
                            color_for_closed='black', show_elements=False)

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like, masked_array, or str
        Values at the cells on the grid. Alternatively, can be a field name
        (string) from which to draw the data from the grid.
    var_name : str, optional
        Name of the variable to put in plot title.
    var_units : str, optional
        Units for the variable being plotted.
    grid_units : tuple of str, optional
        Units for y, and x dimensions.
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
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed elements (default 'black'). If None, closed
        (or masked) elements will be transparent.
    show_elements : bool
        If True, and grid is a Voronoi, extra grid elements (nodes, faces,
        corners) will be plotted along with just the colour of the cell
        (defaults False).

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
        raise ValueError('number of values must match number of cells or '
                         'number of nodes')

    values_at_cell = np.ma.asarray(values_at_cell)
    values_at_cell.mask = True
    values_at_cell.mask[grid.core_cells] = False

    myimage = _imshow_grid_values(grid,
                                  values_at_cell.reshape(grid.cell_grid_shape),
                                  **kwds)

    if isinstance(values, str):
        plt.title(values)

    return myimage


@deprecated(use='imshow_grid_at_cell', version='0.5')
def imshow_cell_grid(grid, values, **kwds):
    imshow_grid_at_cell(grid, values, **kwds)


def _imshow_grid_values(grid, values, var_name=None, var_units=None,
                        grid_units=(None, None), symmetric_cbar=False,
                        cmap='pink', limits=None, allow_colorbar=True,
                        vmin=None, vmax=None,
                        norm=None, shrink=1., color_for_closed='black',
                        show_elements=False):

    gridtypes = inspect.getmro(grid.__class__)

    cmap = plt.get_cmap(cmap)
    if color_for_closed is not None:
        cmap.set_bad(color=color_for_closed)
    else:
        cmap.set_bad(alpha=0.)

    if isinstance(grid, RasterModelGrid):
        if values.ndim != 2:
            raise ValueError('values must have ndim == 2')

        y = np.arange(values.shape[0] + 1) * grid.dy - grid.dy * .5
        x = np.arange(values.shape[1] + 1) * grid.dx - grid.dx * .5

        kwds = dict(cmap=cmap)
        (kwds['vmin'], kwds['vmax']) = (values.min(), values.max())
        if (limits is None) and ((vmin is None) and (vmax is None)):
            if symmetric_cbar:
                (var_min, var_max) = (values.min(), values.max())
                limit = max(abs(var_min), abs(var_max))
                (kwds['vmin'], kwds['vmax']) = (- limit, limit)
        elif limits is not None:
            (kwds['vmin'], kwds['vmax']) = (limits[0], limits[1])
        else:
            if vmin is not None:
                kwds['vmin'] = vmin
            if vmax is not None:
                kwds['vmax'] = vmax

        myimage = plt.pcolormesh(x, y, values, **kwds)

        plt.gca().set_aspect(1.)
        plt.autoscale(tight=True)

        if allow_colorbar:
            plt.colorbar(norm=norm, shrink=shrink)

        plt.xlabel('X (%s)' % grid_units[1])
        plt.ylabel('Y (%s)' % grid_units[0])

        if var_name is not None:
            plt.title('%s (%s)' % (var_name, var_units))

        # plt.show()

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
        if limits is None:
            # only want to work with NOT CLOSED nodes
            open_nodes = grid.status_at_node != 4
            if symmetric_cbar:
                (var_min, var_max) = (values.flat[
                    open_nodes].min(), values.flat[open_nodes].max())
                limit = max(abs(var_min), abs(var_max))
                (vmin, vmax) = (- limit, limit)
            else:
                (vmin, vmax) = (values.flat[
                    open_nodes].min(), values.flat[open_nodes].max())
        else:
            (vmin, vmax) = (limits[0], limits[1])
        cNorm = colors.Normalize(vmin, vmax)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        colorVal = scalarMap.to_rgba(values)

        if show_elements:
            myimage = voronoi_plot_2d(grid.vor)
        mycolors = (i for i in colorVal)
        for order in grid.vor.point_region:
            region = grid.vor.regions[order]
            colortouse = next(mycolors)
            if -1 not in region:
                polygon = [grid.vor.vertices[i] for i in region]
                plt.fill(*zip(*polygon), color=colortouse)

        plt.gca().set_aspect(1.)
        # plt.autoscale(tight=True)
        plt.xlim((np.min(grid.node_x), np.max(grid.node_x)))
        plt.ylim((np.min(grid.node_y), np.max(grid.node_y)))

        scalarMap.set_array(values)
        plt.colorbar(scalarMap)

        plt.xlabel('X (%s)' % grid_units[1])
        plt.ylabel('Y (%s)' % grid_units[0])

        if var_name is not None:
            plt.title('%s (%s)' % (var_name, var_units))


def imshow_grid(grid, values, **kwds):
    """Prepare a map view of data over all nodes or cells in the grid.

    Data is plotted as colored cells. If at='node', the surrounding cell is
    shaded with the value at the node at its center. If at='cell', the cell
    is shaded with its own value. Outer edges of perimeter cells are
    extrapolated. Closed elements are colored uniformly (default black,
    overridden with kwd 'color_for_closed'); other open boundary nodes get
    their actual values.

    *values* can be a field name, a regular array, or a masked array. If a
    masked array is provided, masked entries will be treated as if they were
    Landlab CLOSED_BOUNDARYs. Used together with the color_at_closed=None
    keyword (i.e., "transparent"), this can allow for construction of overlay
    layers in a figure (e.g., only defining values in a river network, and
    overlaying it on another landscape).

    Use matplotlib functions like xlim, ylim to modify your plot after calling
    :func:`imshow_grid`, as desired.

    This function happily works with both regular and irregular grids.

    Construction ::

        imshow_grid(grid, values, at='node', var_name=[field_name],
                    var_units=[field_units], grid_units=None,
                    symmetric_cbar=False, cmap='pink',
                    limits=(values.min(), values.max()), vmin=values.min(),
                    vmax=values.max(), allow_colorbar=True, norm=[linear],
                    shrink=1., color_for_closed='black', show_elements=False,
                    show=False)

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
    var_name : str, optional
        Name of the variable to put in plot title.
    var_units : str, optional
        Units for the variable being plotted.
    grid_units : tuple of str, optional
        Units for y, and x dimensions.
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
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    shrink : float
        Fraction by which to shrink the colorbar.
    color_for_closed : str or None
        Color to use for closed elements (default 'black'). If None, closed
        (or masked) elements will be transparent.
    show_elements : bool
        If True, and grid is a Voronoi, extra grid elements (nodes, faces,
        corners) will be plotted along with just the colour of the cell
        (defaults False).
    show : bool
        If True, plot the figure immediately after calling this method.
    """
    show = kwds.pop('show', False)
    values_at = kwds.pop('values_at', 'node')
    values_at = kwds.pop('at', values_at)

    if isinstance(values, str):
        values = grid.field_values(values_at, values)

    if isinstance(values, str):
        values = grid.field_values(values_at, values)

    if values_at == 'node':
        imshow_node_grid(grid, values, **kwds)
    elif values_at == 'cell':
        imshow_cell_grid(grid, values, **kwds)
    else:
        raise TypeError('value location %s not understood' % values_at)

    if show:
        plt.show()
