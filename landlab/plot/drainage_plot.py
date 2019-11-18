"""Plot drainage network.

"""
import matplotlib.pylab as plt
import numpy as np

# KRB, FEB 2017.
import six

from landlab import (
    CLOSED_BOUNDARY,
    CORE_NODE,
    FIXED_GRADIENT_BOUNDARY,
    FIXED_VALUE_BOUNDARY,
)
from landlab.plot.imshow import imshow_grid


def drainage_plot(
    mg,
    surface="topographic__elevation",
    receivers=None,
    proportions=None,
    surf_cmap="gray",
    quiver_cmap="viridis",
    title="Drainage Plot",
):

    if isinstance(surface, six.string_types):
        colorbar_label = surface
    else:
        colorbar_label = "topographic_elevation"
    imshow_grid(mg, surface, cmap=surf_cmap, colorbar_label=colorbar_label)

    if receivers is None:
        receivers = mg.at_node["flow__receiver_node"]
        if proportions is None:
            if "flow__receiver_proportions" in mg.at_node:
                proportions = mg.at_node["flow__receiver_proportions"]
    else:
        receivers = np.asarray(receivers)

    if receivers.ndim == 1:
        receivers = np.expand_dims(receivers, axis=1)

    nreceievers = receivers.shape[-1]

    propColor = plt.get_cmap(quiver_cmap)
    for j in range(nreceievers):
        rec = receivers[:, j]
        is_bad = rec == -1

        xdist = -0.8 * (mg.x_of_node - mg.x_of_node[rec])
        ydist = -0.8 * (mg.y_of_node - mg.y_of_node[rec])

        if proportions is None:
            proportions = np.ones_like(receivers, dtype=float)

        is_bad[proportions[:, j] == 0.0] = True

        xdist[is_bad] = np.nan
        ydist[is_bad] = np.nan

        prop = proportions[:, j] * 256.0
        lu = np.floor(prop)
        colors = propColor(lu.astype(int))

        shape = (mg.number_of_nodes, 1)

        plt.quiver(
            mg.x_of_node.reshape(shape),
            mg.y_of_node.reshape(shape),
            xdist.reshape(shape),
            ydist.reshape(shape),
            color=colors,
            angles="xy",
            scale_units="xy",
            scale=1,
            zorder=3,
        )

    # Plot differen types of nodes:
    o, = plt.plot(
        mg.x_of_node[mg.status_at_node == CORE_NODE],
        mg.y_of_node[mg.status_at_node == CORE_NODE],
        "b.",
        label="Core Nodes",
        zorder=4,
    )
    fg, = plt.plot(
        mg.x_of_node[mg.status_at_node == FIXED_VALUE_BOUNDARY],
        mg.y_of_node[mg.status_at_node == FIXED_VALUE_BOUNDARY],
        "c.",
        label="Fixed Gradient Nodes",
        zorder=5,
    )
    fv, = plt.plot(
        mg.x_of_node[mg.status_at_node == FIXED_GRADIENT_BOUNDARY],
        mg.y_of_node[mg.status_at_node == FIXED_GRADIENT_BOUNDARY],
        "g.",
        label="Fixed Value Nodes",
        zorder=6,
    )
    c, = plt.plot(
        mg.x_of_node[mg.status_at_node == CLOSED_BOUNDARY],
        mg.y_of_node[mg.status_at_node == CLOSED_BOUNDARY],
        "r.",
        label="Closed Nodes",
        zorder=7,
    )

    node_id = np.arange(mg.number_of_nodes)
    flow_to_self = receivers[:, 0] == node_id

    fts, = plt.plot(
        mg.x_of_node[flow_to_self],
        mg.y_of_node[flow_to_self],
        "kx",
        markersize=6,
        label="Flows To Self",
        zorder=8,
    )

    ax = plt.gca()

    ax.legend(
        labels=[
            "Core Nodes",
            "Fixed Gradient Nodes",
            "Fixed Value Nodes",
            "Closed Nodes",
            "Flows To Self",
        ],
        handles=[o, fg, fv, c, fts],
        numpoints=1,
        loc="center left",
        bbox_to_anchor=(1.7, 0.5),
    )
    sm = plt.cm.ScalarMappable(cmap=propColor, norm=plt.Normalize(vmin=0, vmax=1))
    sm._A = []
    cx = plt.colorbar(sm)
    cx.set_label("Proportion of Flow")
    plt.title(title)
