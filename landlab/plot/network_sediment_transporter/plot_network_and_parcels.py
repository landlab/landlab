# -*- coding: utf-8 -*-
"""
This code plots the network and colors each link according to a link attribute.

This code plots the network and parcels and colors each
parcel according to a link or parcel attribute.

Authors: Jon Czuba, Allison Pfeiffer, Katy Barnhart
"""
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from landlab.plot.network_sediment_transporter.locate_parcel_xy import locate_parcel_xy
from landlab.utils.return_array import return_array_at_link


def plot_network_and_parcels(
    grid,
    parcels,
    parcel_time_index=None,
    network_color=None,
    network_linewidth=None,
    link_attribute=None,
    network_cmap="blues",
    parcel_color=None,
    parcel_size=None,
    parcel_color_attribute=None,
    parcel_cmap=None,
    parcel_size_attribute=None,
    parcel_alpha=0.5,
    fig=None,
    **kwargs
):
    """Descriptive text.

    More descriptive text.


    Parameters
    ----------
    grid
    parcels
    parcel_time_index = None. Default is last timestep.

    # to set the link colors provide either:
    network_color="k"
    network_linewidth=0.1

    # or
    link_attribute, array or field name at link.
    network_cmap = "blues"

    # to set the parcel size, shape, either provide constant values or attributes.

    #for color:

    parcel_color = "g"

    # or

    parcel_color_attribute : parcel attribute name.
    parcel_color_cmap = "viridis"

    # for size:

    parcel_size = 1

    # or

    parcel_size_attribute: parcel atribute name.


    parcel_alpha

    # figure information.

    fig, figure object to use.

    **kwargs Anything else to pass to figure creation.

    Returns
    -------
    fig


    """
    # part 0 checking and default setting.

    # only network color/linewidth provided OR link attribute.
    if (link_attribute is not None) and (
        (network_color is not None) or (network_linewidth is not None)
    ):
        assert ValueError

    if link_attribute is None:
        network_color = network_color or "b"
        network_linewidth = network_linewidth or 0.5
        legend_link = False
    else:
        legend_link = True

    # only parcel color OR parcel_color_attribute.
    if (parcel_color_attribute is not None) and (parcel_color is not None):
        assert ValueError

    if parcel_color_attribute is None:
        parcel_color = parcel_color or "c"
        legend_parcel_color = False
    else:
        legend_parcel_color = True

    # only parcel size or parcel_size_attribute
    if (parcel_size_attribute is not None) and (parcel_size is not None):
        assert ValueError

    if parcel_size_attribute is None:
        parcel_size = parcel_size or 1.0
        legend_parcel_size = False
    else:
        legend_parcel_size = True

    # parcel time:
    parcel_time_index = parcel_time_index or -1
    # also figure out whether the legend will have one, two, or three
    # parts (linewidth, parcel size, parcel color)
    n_legends = legend_link + legend_parcel_size + legend_parcel_color

    # SET UP FIGURE
    if fig is None:
        fig = plt.figure(**kwargs)

    spec = gridspec.GridSpec(
        ncols=1, nrows=1 + n_legends, left=0, right=1, top=1, bottom=0, figure=fig
    )

    ax = fig.add_subplot(spec[0, 0])
    legend_idx = 1
    # SET UP LINESEGMENTS FOR NETWORK. If polylines exist use, otherwise use
    # endpoints.
    segments = []
    if "x_of_polyline" in grid.at_link:
        x_of_polylines = grid["link"]["x_of_polyline"]
        y_of_polylines = grid["link"]["y_of_polyline"]
        for x, y in zip(x_of_polylines, y_of_polylines):
            segment = np.array((np.array(x), np.array(y))).T
            segments.append(segment)
    else:
        for i in range(grid.size("links")):
            nal = grid.nodes_at_link[i]
            segment = np.array(
                (np.array(grid.x_of_node[nal]), np.array(grid.y_of_node[nal]))
            ).T
            segments.append(segment)

    # Add Linesegments and Configure.
    if link_attribute:
        line_segments = LineCollection(segments)
        line_segments.set_array(return_array_at_link(grid, link_attribute))
        ax.add_collection(line_segments)
        lax = fig.add_subplot(spec[legend_idx, 0])
        legend_idx += 1
        fig.colorbar(line_segments, cax=lax)
    else:
        line_segments = LineCollection(
            segments, colors=network_color, linewidth=network_linewidth
        )
        ax.add_collection(line_segments)

    # Part 2: Add Parcels.
    X = np.ones(len(parcels.dataset.element_id)) - 1
    Y = np.ones(len(parcels.dataset.element_id)) - 1

    # Locate parcel XY for each parcel at a particular time
    for parcel_number in range(parcels.dataset.item_id.size):
        XY = locate_parcel_xy(grid, parcels, parcel_time_index, parcel_number)
        X[parcel_number] = XY[0]
        Y[parcel_number] = XY[1]

    # plot X,Y point on delineated network and color/size point according to a
    # certain attribute of the parcel or the link in which the parcel resides

    # Plot parcels
    if parcel_color_attribute is not None:
        parcel_color = None
        # if this is true, then instead of none we need to get the right
        # values from the parcels and scale/normalize them correctly.

    if parcel_size_attribute is not None:
        parcel_size = None
        # if this is true, then instead of none we need to get the right
        # values from the parcels and scale/normalize them correctly.

    ax.scatter(X, Y, s=parcel_size, c=parcel_color, alpha=parcel_alpha)

    # add parcel size and parcel color legends.

    # We need to set the plot limits, they will not autoscale
    ax.set_xlim(grid.x_of_node.min(), grid.x_of_node.max())
    ax.set_ylim(grid.y_of_node.min(), grid.y_of_node.max())

    return fig
