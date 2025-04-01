import matplotlib.pyplot as plt
import numpy as np


def plot_nodes(graph, color="r", with_id=True, markersize=4):
    for node in range(len(graph.x_of_node)):
        x, y = graph.x_of_node[node], graph.y_of_node[node]
        plt.plot(
            graph.x_of_node[node],
            graph.y_of_node[node],
            "o",
            color=color,
            markersize=markersize,
        )
        if with_id:
            plt.text(x, y, node, color=color, size=16)


def plot_links(
    graph, color="b", linestyle="solid", with_id=True, as_arrow=True, linewidth=None
):
    if as_arrow:
        head_width = 0.1
    else:
        head_width = 0.0
    for link, nodes in enumerate(graph.nodes_at_link):
        x, y = graph.x_of_node[nodes[0]], graph.y_of_node[nodes[0]]
        dx, dy = graph.x_of_node[nodes[1]] - x, graph.y_of_node[nodes[1]] - y
        plt.arrow(
            x,
            y,
            dx,
            dy,
            head_width=head_width,
            linewidth=linewidth,
            length_includes_head=True,
            color=color,
            linestyle=linestyle,
        )
        if with_id:
            plt.text(x + dx * 0.5, y + dy * 0.5, link, size=16, color=color)


def plot_patches(graph, color="g", with_id=False):
    from matplotlib.patches import Polygon

    for patch, nodes in enumerate(graph.nodes_at_patch):
        nodes = nodes[nodes >= 0]
        x, y = np.mean(graph.x_of_node[nodes]), np.mean(graph.y_of_node[nodes])
        plt.gca().add_patch(
            Polygon(graph.xy_of_node[nodes], ec=color, fc=None, alpha=0.5)
        )
        if with_id:
            plt.text(
                x,
                y,
                patch,
                color=color,
                size=16,
                horizontalalignment="center",
                verticalalignment="center",
            )


def plot_graph(graph, at="node,link,patch", with_id=True, axes=None):
    """Plot elements of a graph.

    Parameters
    ----------
    graph : graph-like
        A landlab graph-like object.
    at : str or iterable of str
        Comma-separated list of elements to plot.
    with_id : str, iterable of str or bool
        Indicate which elements should be plotted with their corresponding id.
        Either a comma-separated list of grid elements or ``True`` to include
        ids for all elements of ``False`` for no elements.
    axes : , optional
        Add the plot to an existing matplotlib ``Axes``, otherwise, create a new one.

    Returns
    -------
    ``Axes``
        The ``Axes`` containing the plot.
    """
    EVERYWHERE = {"node", "link", "patch", "corner", "face", "cell"}

    if isinstance(with_id, bool):
        with_id = EVERYWHERE if with_id else set()
    else:
        with_id = _parse_locations_as_set(with_id)
    locs = _parse_locations_as_set(at)

    ax = plt.axes() if axes is None else axes

    ax.set_xlim([min(graph.x_of_node) - 0.5, max(graph.x_of_node) + 0.5])
    ax.set_ylim([min(graph.y_of_node) - 0.5, max(graph.y_of_node) + 0.5])

    if "node" in locs:
        plot_nodes(graph, with_id="node" in with_id, markersize=4)
    if "link" in locs:
        plot_links(graph, with_id="link" in with_id, linewidth=None, as_arrow=True)
    if "patch" in locs:
        plot_patches(graph, with_id="patch" in with_id)

    if "corner" in locs:
        plot_nodes(graph.dual, color="c", with_id="corner" in with_id)
    if "face" in locs:
        plot_links(graph.dual, linestyle="dotted", color="k", with_id="face" in with_id)
    if "cell" in locs and graph.number_of_cells > 0:
        plot_patches(graph.dual, color="m", with_id="cell" in with_id)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect(1.0)

    return ax


def _parse_locations_as_set(locations):
    """Parse grid element locations as a set.

    Parameters
    ----------
    locations : str or iterable of str
        Grid locations.

    Returns
    -------
    set
        Grid locations as strings.

    Raises
    ------
    ValueError
        If any of the locations are invalid.
    """
    EVERYWHERE = {"node", "link", "patch", "corner", "face", "cell"}

    if isinstance(locations, str):
        as_set = set(locations.split(","))
    else:
        as_set = set(locations)

    as_set = {item.strip() for item in as_set}

    unknown = sorted(as_set - EVERYWHERE)
    if unknown:
        unknown = [repr(item) for item in unknown]
        raise ValueError(
            f"unknown location{'s' if len(unknown) > 1 else ''} ({', '.join(unknown)})"
        )

    return as_set
