#!/usr/bin/env python3

"""
Methods to 3d plot data defined on Landlab grids.
"""

import numpy as np
import inspect
from landlab.field.scalar_data_fields import FieldError
try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)

from mpl_toolkits.mplot3d import Axes3D

from landlab.grid import CLOSED_BOUNDARY
from landlab.grid.raster import RasterModelGrid
#from landlab.grid.voronoi import VoronoiDelaunayGrid
from landlab.plot.event_handler import query_grid_on_button_press

def _plot_3d_surface(grid, values, plot_name=None, var_name=None,
                        var_units=None, grid_units=(None, None),
                        cmap=None, norm=None, vmin=None,
                        vmax=None, alpha=0.5,
                        colorbar_label = None, allow_colorbar=True, shrink=1.,
                        rstride=1, cstride=1,
                        init_view=(None, -130),output=None):


#    gridtypes = inspect.getmro(grid.__class__) # to be used for Voronoi

    cmap = plt.get_cmap(cmap)


    if isinstance(grid, RasterModelGrid):
        if values.ndim != 2:
            raise ValueError('values must have ndim == 2')


        core_mask = grid.node_is_core(grid.nodes)
        number_of_core_lines = np.max(core_mask.sum(axis=0))
        number_of_core_rows = np.max(core_mask.sum(axis=1))
        x_plot3D = grid.node_x[grid.core_nodes].reshape(number_of_core_lines, number_of_core_rows)
        y_plot3D = grid.node_y[grid.core_nodes].reshape(number_of_core_lines, number_of_core_rows)
        z_plot3D = grid.at_node['topographic__elevation'][grid.core_nodes].reshape(number_of_core_lines, number_of_core_rows)

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        

        kwds = dict(cmap=cmap)
        (kwds['vmin'], kwds['vmax']) = (values.min(), values.max())
        if norm is not None:
            kwds['norm'] = norm
        if vmin is not None:
            kwds['vmin'] = vmin
        if vmax is not None:
            kwds['vmax'] = vmax
        if cstride is not None:
            kwds['cstride'] = cstride
        if rstride is not None:
            kwds['rstride'] = rstride
        if alpha is not None:
            kwds['alpha'] = alpha

        myimage = ax.plot_surface(x_plot3D, y_plot3D, z_plot3D, **kwds)


        if allow_colorbar:
            cbar = fig.colorbar(myimage, shrink=shrink, aspect=5)
            if colorbar_label:
                cbar.set_label(colorbar_label)


    else:
        raise TypeError('Landlab Plot3d not implemented for grid types other'\
                        'than RasterModelGrid')

    # Figure finish: labels, title, etc.
    if grid_units[1] is None and grid_units[0] is None:
        grid_units = grid.axis_units
        if grid_units[1] == '-' and grid_units[0] == '-':
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
        else:
            ax.set_xlabel('X (%s)' % grid_units[1])
            ax.set_ylabel('Y (%s)' % grid_units[0])
    else:
        ax.set_xlabel('X (%s)' % grid_units[1])
        ax.set_ylabel('Y (%s)' % grid_units[0])

    if plot_name is not None:
        plt.title('%s' % (plot_name))

    if var_name is not None or var_units is not None:
        if var_name is not None:
            assert type(var_name) is str
            if var_units is not None:
                assert type(var_units) is str
                colorbar_label = var_name + ' (' + var_units + ')'
            else:
                colorbar_label = var_name
        else:
            assert type(var_units) is str
            colorbar_label = '(' + var_units + ')'
        assert type(colorbar_label) is str
        assert allow_colorbar
        cbar.set_label(colorbar_label)

    if init_view[1] is None and init_view[0] is None:
        ax.view_init(elev=None, azim=-130)
    else:
        ax.view_init(elev=init_view[0], azim=init_view[1])

    if output is not None:
        if type(output) is str:
            plt.savefig(output)
            plt.clf()
        elif output:
            plt.show()
    
    
    
    
    
    
    
    
    

def plot_3d_surface_at_node(grid, values, **kwds):
    """plot3D_surface_at_node(grid, values, plot_name=None, var_name=None,
                        var_units=None, grid_units=(None, None),
                        cmap=None, norm=None, vmin=None,
                        vmax=None, alpha=0.5,
                        colorbar_label = None, allow_colorbar=True, shrink=1.,
                        rstride=1, cstride=1,
                        init_view=(None, -130),output=None)

    Prepare a 3-dimension view of data over core nodes in the grid.

    Data is plotted as a 3d surface interpolated between values at nodes, which
    make the corners of the surface elements (i.e. the surface's 'cells' are
    Landlab patches). By default the surface elements are colored in shades of
    a solid color but color mapping by supplying cmap argument is also
    supported.
    Boundary nodes are excluded from representation for better clarity: this
    means that the outer edges of the figure are core nodes.

    *values* can be a field name or a regular array.

    Node coordinates are printed when a mouse button is pressed on a cell in
    the plot.

    This function only works with regular grids at the moment.

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like or str
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
    cmap : str
        Name of a colormap
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    vmin, vmax: floats
        Minimum and maximum of the colorbar.
    alpha: float
        Transparency factor [0, 1]. Default is 0.5.
    allow_colorbar : bool
        If True, include the colorbar. Default is True.
    colorbar_label : str or None
        The string with which to label the colorbar.
    shrink : float
        Fraction by which to shrink the colorbar. Default is 1 (no shrink).
    rstride, cstride: integers
        Array row and column strides, respectively (step sizes). Set the
        stride used to sample the input data to generate the graph. If 1k by
        1k arrays are passed in, the default values for the strides will
        result in a 1000x1000 grid being plotted. Defaults to 1.
    init_view: tuple
        Elevation and azimuth of initial plot view. Default is (None, -130).
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
        values_at_node = values

    if values_at_node.size != grid.number_of_nodes:
        raise ValueError('number of values does not match number of nodes')

    values_at_node = np.ma.masked_where(
        grid.status_at_node == CLOSED_BOUNDARY, values_at_node)

    try:
        shape = grid.shape
    except AttributeError:
        shape = (-1, )

    _plot_3d_surface(grid, values_at_node.reshape(shape), **kwds)


    plt.gcf().canvas.mpl_connect('button_press_event',
       lambda event: query_grid_on_button_press(event, grid))

def plot_3d_surface_at_cell(grid, values, **kwds):
    """plot3D_surface_at_node(grid, values, plot_name=None, var_name=None,
                        var_units=None, grid_units=(None, None),
                        cmap='Oranges', norm=None, vmin=None,
                        vmax=None, alpha=0.5,
                        colorbar_label = None, allow_colorbar=True, shrink=1.,
                        rstride=1, cstride=1,
                        init_view=(None, -130),output=None)

    Prepare a 3-dimension view of data over core cells in the grid.

    Data is plotted as a 3d surface interpolated between values at cells, which
    make the corners of the surface elements (i.e. the surface's 'cells' are
    Landlab patches). By default the surface elements are colored in shades of
    a solid color but color mapping by supplying cmap argument is also
    supported.
    Boundary cells are excluded from representation for better clarity: this
    means that the outer edges of the figure are core cell values.

    *values* can be a field name or a regular array.

    Cell coordinates are printed when a mouse button is pressed on a cell in
    the plot.

    This function only works with regular grids at the moment.

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like or str
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
    cmap : str
        Name of a colormap
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    vmin, vmax: floats
        Minimum and maximum of the colorbar.
    alpha: float
        Transparency factor [0, 1]. Default is 0.5.
    allow_colorbar : bool
        If True, include the colorbar. Default is True.
    colorbar_label : str or None
        The string with which to label the colorbar.
    shrink : float
        Fraction by which to shrink the colorbar. Default is 1 (no shrink).
    rstride, cstride: integers
        Array row and column strides, respectively (step sizes). Set the
        stride used to sample the input data to generate the graph. If 1k by
        1k arrays are passed in, the default values for the strides will
        result in a 1000x1000 grid being plotted. Defaults to 1.
    init_view: tuple
        Elevation and azimuth of initial plot view. Default is (None, -130).
    output : None, string, or bool
        If None (or False), the image is sent to the imaging buffer to await
        an explicit call to show() or savefig() from outside this function.
        If a string, the string should be the path to a save location, and the
        filename (with file extension). The function will then call
        plt.savefig([string]) itself. If True, the function will call
        plt.show() itself once plotting is complete.
    """
    if isinstance(values, str):
        values_at_cell = grid.at_cell[values]
    else:
        values_at_cell = values

    if values_at_cell.size != grid.number_of_cells:
        raise ValueError('number of values does not match number of nodes')

    values_at_cell = np.ma.asarray(values_at_cell)
    values_at_cell.mask = True
    values_at_cell.mask[grid.core_cells] = False

    my_image = _plot_3d_surface(grid,
                                values_at_cell.reshape(grid.cell_grid_shape),
                                **kwds)

    return my_image

    plt.gcf().canvas.mpl_connect('button_press_event',
       lambda event: query_grid_on_button_press(event, grid))






def plot_3d_surface(grid, values, **kwds):
    """plot_3d_surface(grid, values, plot_name=None, var_name=None,
                        var_units=None, grid_units=(None, None),
                        cmap='Oranges', norm=None, vmin=None,
                        vmax=None, alpha=0.5,
                        colorbar_label = None, allow_colorbar=True, shrink=1.,
                        rstride=1, cstride=1,
                        init_view=(None, -130),output=None)

    Prepare a 3-dimension view of data over core nodes or cells in the grid.

    Data is plotted as an interpolated 3d surface. If at='node', the surface
    is interpolated between values at core nodes, which then make the corners
    of the surface elements (i.e. the surface's 'cells' are Landlab patches).
    If at='cell', the surface is interpolated between values at core nodes,
    which then make the corners of the surface elements (i.e. the surface's
    'cells' are Landlab patches).

    By default the surface elements are colored in shades of
    a solid color but color mapping by supplying cmap argument is also
    supported.
    Boundary nodes/cells are excluded from representation for better clarity:
    this means that the outer edges of the figure are core nodes (or core cell
    values).

    *values* can be a field name or a regular array.

    Node coordinates are printed when a mouse button is pressed on a cell in
    the plot.

    This function only works with regular grids at the moment.

    Parameters
    ----------
    grid : ModelGrid
        Grid containing the field to plot, or describing the geometry of the
        provided array.
    values : array_like or str
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
    cmap : str
        Name of a colormap
    norm : matplotlib.colors.Normalize
        The normalizing object which scales data, typically into the interval
        [0, 1]. Ignore in most cases.
    vmin, vmax: floats
        Minimum and maximum of the colorbar.
    alpha: float
        Transparency factor [0, 1]. Default is 0.5.
    allow_colorbar : bool
        If True, include the colorbar. Default is True.
    colorbar_label : str or None
        The string with which to label the colorbar.
    shrink : float
        Fraction by which to shrink the colorbar. Default is 1 (no shrink).
    rstride, cstride: integers
        Array row and column strides, respectively (step sizes). Set the
        stride used to sample the input data to generate the graph. If 1k by
        1k arrays are passed in, the default values for the strides will
        result in a 1000x1000 grid being plotted. Defaults to 1.
    init_view: tuple
        Elevation and azimuth of initial plot view. Default is (None, -130).
    output : None, string, or bool
        If None (or False), the image is sent to the imaging buffer to await
        an explicit call to show() or savefig() from outside this function.
        If a string, the string should be the path to a save location, and the
        filename (with file extension). The function will then call
        plt.savefig([string]) itself. If True, the function will call
        plt.show() itself once plotting is complete.
    """
    
    show = kwds.pop('show', False)
    values_at = kwds.pop('values_at', 'node')
    values_at = kwds.pop('at', values_at)
    
    if isinstance(values, str):
        values = grid.field_values(values_at, values)

    if isinstance(values, str):
        values = grid.field_values(values_at, values)

    if values_at == 'node':
        plot_3d_surface_at_node(grid, values, **kwds)
    elif values_at == 'cell':
        plot_3d_surface_at_cell(grid, values, **kwds)
    else:
        raise TypeError('value location %s not understood' % values_at)

    # retained for backwards compatibility:
    if show:
        plt.show()

