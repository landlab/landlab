"""
This code was designed to plot outputs of the NetworkSedimentTransporter
landlab component.

This code plots:
    - the network, with option to color each link according to a link attribute.
    - the parcels, with option to color and size each parcel according to
    parcel attributes.

Authors: Katy Barnhart, Jon Czuba, Allison Pfeiffer
"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

from landlab.plot.network_sediment_transporter.locate_parcel_xy import locate_parcel_xy
from landlab.utils.return_array import return_array_at_link


def plot_network_and_parcels(
    grid,
    parcels,
    parcel_time_index=None,
    map_buffer=0.1,
    parcel_filter=None,
    network_color=None,
    link_attribute=None,
    link_attribute_title=None,
    network_cmap="cividis",
    network_norm=None,
    network_linewidth=None,
    parcel_color=None,
    parcel_color_attribute=None,
    parcel_color_attribute_title=None,
    parcel_color_cmap="plasma",
    parcel_color_norm=None,
    parcel_size=None,
    parcel_size_attribute=None,
    parcel_size_attribute_title=None,
    parcel_size_norm=None,
    parcel_size_min=5,
    parcel_size_max=40,
    parcel_alpha=0.5,
    fig=None,
    **kwargs,
):
    """Plot a river network and parcels on the river network. Intended to
    display the results of the NetworkSedimentTransporter component.

    The river network (an instance of NetworkModelGrid) is plotted either as
    straight links between grid nodes, or (if the network was created using a
    shapefile to set network topology) as sinuous lines representing the actual
    link geometry.

    The parcels (an instance of DataRecord) are represented as dot markers
    along the links, with the marker location set by parcel attribute
    `location_at_link`. The default is to plot the parcel locations at the
    last timestep in DataRecord, though any time index may be specified.

    Use of this plotting tool is described in detail in a landlab tutorial.

    Parameters
    ----------
    grid : NetworkModelGrid
        Instance of NetworkModelGrid.
    parcels : DataRecord
        Instance of Landlab DataRecord, with the same attribute requirements as
        NetworkSedimentTransporter.
    parcel_time_index : int, time index of parcels DataRecord
        Parcel time index to plot. Default is last timestep in parcels
        DataRecord.
    map_buffer : 0.1
        Increase the plot extent by at least this much (default 0.1). Note,
        b/c of axis equal, may be more.
    parcel_filter : boolean array of shape (number_of_parcels, )
        Filter to plot only a selection of the parcels.

    ########################
    ## Part 1: Network. To set the link colors provide either:

    network_color="k"
        Uniform color for network links.

    # or

    link_attribute : array or field name at link
        Value used to set link color. Categorical options not supported. Must
        be continuous.
    link_attribute_title : str
        String to use as the title, if link_attribute is a string, it is
        used as the default.
    network_cmap : "cividis"
        Name of colormap for network.
    network_norm : matplotlib color normalizer.
        https://matplotlib.org/3.1.1/tutorials/colors/colormapnorms.html
        Default is linear between min and max of link_attribute.

    # linewidth will be recognized by either link coloring option.

    network_linewidth : float
        Width of network lines (default 0.5).

    ########################
    ## Part 2: Parcels. To set the parcel color, provide either:

    parcel_color : color str
        Constant color used for parcel markers (default "k").

    # or

    parcel_color_attribute : parcel attribute name.
        Categorical options not supported. Must be continuous.
    parcel_color_attribute_title : str
        String to use as the legend title. If parcel_color_attribute is a
        string, it is used as the default.
    parcel_color_cmap : cmap str
        Name of colormap for variable parcel color (default "plasma").
    parcel_color_norm : matplotlib color normalizer.
        https://matplotlib.org/3.1.1/tutorials/colors/colormapnorms.html
        Default is linear between min and max of parcel_color_attribute.

    # for parcel size use either:

    parcel_size : float
        Marker size, in points.

    # or

    parcel_size_attribute: parcel atribute name.
        Categorical options not supported. Must be continuous.
    parcel_size_attribute_title : str
        string to use as the title, if parcel_size_attribute is a string, it is
        used as the default.
    parcel_size_norm : par
        matplotlib color normalizer.
        https://matplotlib.org/3.1.1/tutorials/colors/colormapnorms.html
        Default is linear between min and max of parcel_size_attribute.,
    parcel_size_min : float
        Specify the smallest size of the dot markers plotted, in
        units of points (default 5). Use with parcel_size_max. They will be
        aligned with the limits of parcel_size_norm.
    parcel_size_max : float
        Specify the largest size of the dot markers plotted, in
        units of points (default 40). Use with parcel_size_min. They will be
        aligned with the limits of parcel_size_norm.

    # with constant or attribute, can set parcel transparency

    parcel_alpha : float, between 0 and 1
        Specify parcel marker transparency (default 0.5).

    ################### Miscellaneous

    fig :  figure object
        Default is to create a new figure object.

    **kwargs Anything else to pass to figure creation.

    Returns
    -------
    fig :
        Figure object.

    """
    # part 0 checking and default setting.

    # only network color/linewidth provided OR link attribute.
    if (link_attribute is not None) and (network_color is not None):
        raise ValueError(
            "Only one of link_attribute and network_color can be provided."
        )

    if link_attribute is None:
        network_color = network_color or "c"
        network_linewidth = network_linewidth or 0.5
        legend_link = False
    else:
        legend_link = True
        if link_attribute_title is None:
            if isinstance(link_attribute, str):
                link_attribute_title = link_attribute
            else:
                link_attribute_title = ""

    # only parcel color OR parcel_color_attribute.
    if (parcel_color_attribute is not None) and (parcel_color is not None):
        raise ValueError(
            "Only one of parcel_color_attribute and parcel_color can be provided."
        )

    if parcel_color_attribute is None:
        parcel_color = parcel_color or "k"
        legend_parcel_color = False
    else:
        legend_parcel_color = True
        if parcel_color_attribute_title is None:
            parcel_color_attribute_title = parcel_color_attribute

    # only parcel size or parcel_size_attribute
    if (parcel_size_attribute is not None) and (parcel_size is not None):
        raise ValueError(
            "Only one of parcel_size_attribute and parcel_size can be provided."
        )

    if parcel_size_attribute is None:
        parcel_size = parcel_size or 1.0
        legend_parcel_size = False
    else:
        legend_parcel_size = True

        if parcel_size_attribute_title is None:
            parcel_size_attribute_title = parcel_size_attribute

    # parcel time:
    # cant use standard value or default because a value of 0 is valid.
    if parcel_time_index is None:
        parcel_time_index = -1

    # Figure out whether the legend will have one, two, or three
    # parts (linewidth, parcel size, parcel color)
    n_legends = legend_link + legend_parcel_size + legend_parcel_color

    # set up figure, label and legend gridspecs.
    if fig is None:
        fig = plt.figure(**kwargs)

    spec = gridspec.GridSpec(
        ncols=1,
        nrows=3,
        left=0,
        right=1,
        top=1,
        bottom=0,
        figure=fig,
        height_ratios=[1, 0.1, 0.2],
    )
    ax = fig.add_subplot(spec[0, 0])
    if n_legends > 0:
        label_spec = spec[1, 0].subgridspec(
            ncols=2 * n_legends - 1,
            nrows=1,
            wspace=0.1,
            hspace=1,
            width_ratios=[1] + (n_legends - 1) * [0.2, 1],
        )
        legend_spec = spec[2, 0].subgridspec(
            ncols=2 * n_legends - 1,
            nrows=1,
            wspace=0.1,
            hspace=1,
            width_ratios=[1] + (n_legends - 1) * [0.2, 1],
        )
        legend_idx = 0

    # SET UP LINESEGMENTS FOR NETWORK. If polylines exist use, otherwise use
    # endpoints. Also get the ranges so a buffer can be placed around the
    # network.

    if "x_of_polyline" in grid.at_link:
        xy_of_polylines = _get_xy_of_polylines(
            grid.at_link["x_of_polyline"], grid.at_link["y_of_polyline"]
        )
    else:
        xy_of_polylines = grid.xy_of_node[grid.nodes_at_link]

    xlim, ylim = _calc_xy_limits(xy_of_polylines, buffer_frac=map_buffer)

    # Add Linesegments and Configure.

    # if there is a link attribute.
    if link_attribute is not None:
        line_segments = LineCollection(
            xy_of_polylines,
            cmap=network_cmap,
            norm=network_norm,
            linewidth=network_linewidth,
            zorder=1,
        )
        line_segments.set_array(return_array_at_link(grid, link_attribute))
        ax.add_collection(line_segments)

        # create label
        lax = fig.add_subplot(label_spec[0, legend_idx])
        lax.text(
            0.5,
            0.0,
            "Line Color",
            transform=lax.transAxes,
            color="k",
            ha="center",
            va="top",
            size=plt.rcParams["font.size"] + 2,
        )
        lax.axis("off")

        # add legend.
        lax = fig.add_subplot(legend_spec[0, legend_idx])
        legend_idx += 2
        fig.colorbar(
            line_segments, cax=lax, orientation="horizontal", label=link_attribute_title
        )

    # if link values are constant.
    else:
        line_segments = LineCollection(
            xy_of_polylines, colors=network_color, linewidth=network_linewidth, zorder=1
        )
        ax.add_collection(line_segments)

    # Part 2: Add Parcels.
    X = np.empty(len(parcels.dataset.element_id))
    Y = np.empty(len(parcels.dataset.element_id))

    # Locate parcel XY for each parcel at a particular time
    # some aspects of this may be possible to speed up, but at minimum
    # locate_parcel_xy must be called for each link (since calculating location)
    # along link requires interpoloation.
    # if this occurs we must also ensure the parcel order is maintained b/c of
    # color and shape formatting.

    for parcel_idx in range(parcels.dataset.item_id.size):
        XY = locate_parcel_xy(grid, parcels, parcel_time_index, parcel_idx)
        X[parcel_idx] = XY[0]
        Y[parcel_idx] = XY[1]

    # plot X,Y point on delineated network and color/size point according to a
    # certain attribute of the parcel or the link in which the parcel resides

    # if a parcel color attribute is provided.
    if parcel_color_attribute is not None:
        # if this is true, then instead of none we need to get the right
        # values from the parcels and scale/normalize them correctly. At present
        # plan to support only continuous values. Can be extended to strs as
        # categorical.
        if parcel_color_attribute in parcels.dataset:
            if "time" in parcels.dataset[parcel_color_attribute].sizes:
                parcel_color = parcels.dataset[parcel_color_attribute].values[
                    :, parcel_time_index
                ]
            else:
                parcel_color = parcels.dataset[parcel_color_attribute].values
        else:
            raise ValueError(
                f"Parcel color attribute {parcel_color_attribute} not present in "
                "parcels."
            )

        if parcel_filter is not None:
            parcel_color = parcel_color[parcel_filter]

    if parcel_size_attribute is not None:
        # if this is true, then instead of none we need to get the right
        # values from the parcels and scale/normalize them correctly. At present
        # plan to support only continuous values. Can be extended to strs as
        # categorical.
        if parcel_size_attribute in parcels.dataset:
            if "time" in parcels.dataset[parcel_size_attribute].sizes:
                parcel_size_values = parcels.dataset[parcel_size_attribute].values[
                    :, parcel_time_index
                ]
            else:
                parcel_size_values = parcels.dataset[parcel_size_attribute].values

            if parcel_size_norm is None:
                parcel_size_norm = Normalize(
                    vmin=parcel_size_values.min(), vmax=parcel_size_values.max()
                )

            parcel_size = parcel_size_min + (
                parcel_size_max - parcel_size_min
            ) * parcel_size_norm(parcel_size_values)
        else:
            raise ValueError(
                f"Parcel size attribute {parcel_size_attribute} not present in "
                "parcels."
            )

        if parcel_filter is not None:
            parcel_size = parcel_size[parcel_filter]

    # add scatter, filter x and y if necessary.
    if parcel_filter is not None:
        X = X[parcel_filter]
        Y = Y[parcel_filter]

    scatter = ax.scatter(
        X,
        Y,
        s=parcel_size,
        c=parcel_color,
        alpha=parcel_alpha,
        cmap=parcel_color_cmap,
        norm=parcel_color_norm,
        zorder=2,
    )

    # create legends.
    if legend_parcel_color:
        lax = fig.add_subplot(label_spec[0, legend_idx])
        lax.text(
            0.5,
            0.0,
            "Parcel Color",
            transform=lax.transAxes,
            color="k",
            ha="center",
            va="top",
            size=plt.rcParams["font.size"] + 2,
        )
        lax.axis("off")

        lax = fig.add_subplot(legend_spec[0, legend_idx])
        legend_idx += 2
        fig.colorbar(
            scatter,
            cax=lax,
            orientation="horizontal",
            label=parcel_color_attribute_title,
        )

    if legend_parcel_size:
        lax = fig.add_subplot(label_spec[0, legend_idx])
        lax.text(
            0.5,
            0.0,
            "Parcel Size",
            transform=lax.transAxes,
            color="k",
            ha="center",
            va="top",
            size=plt.rcParams["font.size"] + 2,
        )
        lax.axis("off")

        lax = fig.add_subplot(legend_spec[0, legend_idx])
        handles, _ = scatter.legend_elements(prop="sizes", alpha=0.6)

        if len(handles) - 1 == 0:
            han = handles
            lab = [parcel_size_values.min()]
        else:
            han = handles[:: len(handles) - 1]
            lab = [parcel_size_values.min(), parcel_size_values.max()]

        lax.legend(
            han,
            lab,
            title=parcel_size_attribute_title,
            loc="center",
            frameon=False,
        )

        plt.axis("off")

    # Set the plot limits
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # make axes equal
    ax.axis("equal")

    return fig


def _get_xy_of_polylines(x_of_polylines, y_of_polylines):
    """Zip together x and y coordinate arrays.

    Parameters
    ----------
    x_of_polylines : ndarray
        x coordinates of a series of polyline segments.
    y_of_polylines : ndarray
        y coordinates of a series of polyline segments.

    Returns
    -------
    ndarray
        An ndarray of zipped polyline coordinates.

    Examples
    --------
    >>> x = [[0, 1, 2], [3, 4], [4, 3, 2, 1]]
    >>> y = [[5, 7, 6], [9, 8], [4, 5, 6, 7]]
    >>> xy_of_polylines = _get_xy_of_polylines(x, y)
    >>> len(xy_of_polylines)
    3
    >>> xy_of_polylines[0]
    array([[0, 5],
           [1, 7],
           [2, 6]])
    >>> xy_of_polylines[1]
    array([[3, 9],
           [4, 8]])
    >>> xy_of_polylines[2]
    array([[4, 4],
           [3, 5],
           [2, 6],
           [1, 7]])
    """
    return [np.stack(xy, axis=1) for xy in zip(x_of_polylines, y_of_polylines)]


def _calc_xy_limits(xy_of_segment, buffer_frac=0.0):
    """Calculate xy limits with an optional buffer.

    Parameters
    ----------
    xy_of_segment : iterable of ndarray
        xy coordinates of each segment.
    buffer_frac : float, optional
        Size of buffer as a fraction of the size of the bounding
        box of all the segments. A value of zero mean the limits
        will be 'tight'.

    Returns
    -------
    x_limits, y_limits
        x and y limits.

    Examples
    --------
    >>> xy_of_segments = ([[0, 1], [1, 2]], [[4, 5], [2, 3], [6, 6]], [[2, 9]])
    >>> _calc_xy_limits(xy_of_segments)
    ((0.0, 6.0), (1.0, 9.0))
    >>> _calc_xy_limits(xy_of_segments, buffer_frac=0.5)
    ((-3.0, 9.0), (-3.0, 13.0))
    """
    segments = np.concatenate(xy_of_segment)

    x_limits = _calc_limits(segments[:, 0], buffer_frac=buffer_frac)
    y_limits = _calc_limits(segments[:, 1], buffer_frac=buffer_frac)

    return x_limits, y_limits


def _calc_limits(values, buffer_frac=0.0):
    """Calculate min and max limits with a buffer on each end.

    Parameters
    ----------
    values : iterable
        Values to find limits of.
    buffer_frac : float, optional
        Size of buffer, as a fraction of peak-to-peak, to add to the
        upper and lower limit.

    Returns
    -------
    limits : tuple
        Lower and upper limits.

    Examples
    --------
    >>> _calc_limits([1, 3, 5, 4, 2])
    (1.0, 5.0)
    >>> _calc_limits([1, 3, 5, 4, 2], buffer_frac=0.25)
    (0.0, 6.0)
    """
    values = np.asarray(values)
    min_value, max_value = values.min(), values.max()
    buffer_width = buffer_frac * (max_value - min_value)
    return min_value - buffer_width, max_value + buffer_width
