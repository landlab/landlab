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
from matplotlib.colors import Normalize

from landlab.plot.network_sediment_transporter.locate_parcel_xy import locate_parcel_xy
from landlab.utils.return_array import return_array_at_link


def plot_network_and_parcels(
    grid,
    parcels,
    parcel_time_index=None,
    network_color=None,
    network_linewidth=None,
    link_attribute=None,
    link_attribute_title=None,
    network_cmap="cividis",
    network_norm=None,
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
    link_attribute_title
    network_cmap = "cividis"
    network_norm

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
        if link_attribute_title is None and isinstance(link_attribute, str):
            link_attribute_title = link_attribute

    # only parcel color OR parcel_color_attribute.
    if (parcel_color_attribute is not None) and (parcel_color is not None):
        assert ValueError

    if parcel_color_attribute is None:
        parcel_color = parcel_color or "c"
        legend_parcel_color = False
    else:
        legend_parcel_color = True
        if parcel_color_attribute_title is None and isinstance(
            parcel_color_attribute, str
        ):
            parcel_color_attribute_title = parcel_color_attribute

    # only parcel size or parcel_size_attribute
    if (parcel_size_attribute is not None) and (parcel_size is not None):
        assert ValueError

    if parcel_size_attribute is None:
        parcel_size = parcel_size or 1.0
        legend_parcel_size = False
    else:
        legend_parcel_size = True

        if parcel_size_attribute_title is None and isinstance(
            parcel_size_attribute, str
        ):
            parcel_size_attribute_title = parcel_size_attribute

    # parcel time:
    # cant use standard or because a value of 0 is valid.
    if parcel_time_index is None:
        parcel_time_index = -1

    # also figure out whether the legend will have one, two, or three
    # parts (linewidth, parcel size, parcel color)
    n_legends = legend_link + legend_parcel_size + legend_parcel_color

    # SET UP FIGURE
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
    ax = fig.add_subplot(spec[0, 0])
    legend_idx = 0
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
        for i in range(grid.size("link")):
            nal = grid.nodes_at_link[i]
            segment = np.array(
                (np.array(grid.x_of_node[nal]), np.array(grid.y_of_node[nal]))
            ).T
            segments.append(segment)

    # Add Linesegments and Configure.
    if link_attribute:
        line_segments = LineCollection(
            segments, cmap=network_cmap, norm=network_norm, zorder=1
        )
        line_segments.set_array(return_array_at_link(grid, link_attribute))
        ax.add_collection(line_segments)

        # create legend and label.
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

        lax = fig.add_subplot(legend_spec[0, legend_idx])
        legend_idx += 2
        fig.colorbar(
            line_segments, cax=lax, orientation="horizontal", label=link_attribute_title
        )

        # Todo: add ability to specify color limits.

    else:
        line_segments = LineCollection(
            segments, colors=network_color, linewidth=network_linewidth, zorder=1
        )
        ax.add_collection(line_segments)

    # Part 2: Add Parcels.
    X = np.ones(len(parcels.dataset.element_id)) - 1
    Y = np.ones(len(parcels.dataset.element_id)) - 1

    # Locate parcel XY for each parcel at a particular time
    for parcel_idx in range(parcels.dataset.item_id.size):
        XY = locate_parcel_xy(grid, parcels, parcel_time_index, parcel_idx)
        X[parcel_idx] = XY[0]
        Y[parcel_idx] = XY[1]

    # plot X,Y point on delineated network and color/size point according to a
    # certain attribute of the parcel or the link in which the parcel resides

    # Plot parcels
    if parcel_color_attribute is not None:
        if parcel_color_attribute in parcels.dataset:
            if "time" in parcels.dataset[parcel_color_attribute].dims:
                parcel_color = parcels.dataset[parcel_color_attribute].values[
                    :, parcel_time_index
                ]
            else:
                parcel_color = parcels.dataset[parcel_color_attribute].values
        else:
            raise ValueError

        # if this is true, then instead of none we need to get the right
        # values from the parcels and scale/normalize them correctly. At present
        # plan to support only continuous values. Can be extended to strs as
        # categorical.

    if parcel_size_attribute is not None:
        if parcel_size_attribute in parcels.dataset:
            if "time" in parcels.dataset[parcel_size_attribute].dims:
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
            raise ValueError

        # if this is true, then instead of none we need to get the right
        # values from the parcels and scale/normalize them correctly. At present
        # plan to support only continuous values. Can be extended to strs as
        # categorical.

    # Todo: add ability to specify color limits.

    # Todo: add ability to specify size limits (and size edges).

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
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)

        new_labels = [parcel_size_values.min(), parcel_size_values.max()]
        lax.legend(
            handles[:: len(handles) - 1],
            new_labels,
            title=parcel_size_attribute_title,
            loc="center",
            frameon=False,
        )

        plt.axis("off")
    # We need to set the plot limits, they will not autoscale
    ax.set_xlim(grid.x_of_node.min(), grid.x_of_node.max())
    ax.set_ylim(grid.y_of_node.min(), grid.y_of_node.max())

    return fig
