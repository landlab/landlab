from collections.abc import Iterable
from collections.abc import Mapping
from collections.abc import Sequence
from functools import partial
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.typing import ColorType


def plot_nodes(
    graph,
    *,
    color: ColorType = "r",
    with_id: bool = True,
    markersize: int = 4,
    text: dict[str, Any] | None = None,
):
    text_kwds = merge_text_kwds(text, defaults={"size": 16, "color": color})

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
            plt.text(x, y, node, **text_kwds)


def plot_links(
    graph,
    *,
    color: ColorType = "b",
    linestyle="solid",
    with_id=True,
    as_arrow=True,
    linewidth=None,
    text: dict[str, Any] | None = None,
):
    if as_arrow:
        head_width = 0.1
    else:
        head_width = 0.0
    text_kwds = merge_text_kwds(text, defaults={"size": 16, "color": color})

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
            plt.text(x + dx * 0.5, y + dy * 0.5, link, **text_kwds)


def plot_patches(
    graph,
    *,
    color: ColorType = "g",
    with_id: bool = False,
    text: dict[str, Any] | None = None,
):
    from matplotlib.patches import Polygon

    text_kwds = merge_text_kwds(
        text,
        defaults={
            "size": 16,
            "color": color,
            "horizontalalignment": "center",
            "verticalalignment": "center",
        },
    )

    for patch, nodes in enumerate(graph.nodes_at_patch):
        nodes = nodes[nodes >= 0]
        x, y = np.mean(graph.x_of_node[nodes]), np.mean(graph.y_of_node[nodes])
        plt.gca().add_patch(
            Polygon(graph.xy_of_node[nodes], ec=color, fc=None, alpha=0.5)
        )
        if with_id:
            plt.text(x, y, patch, **text_kwds)


DEFAULT_TEXT_PROPS = (("fontsize", 16),)


def plot_graph(
    graph,
    *,
    at="node,link,patch",
    with_id=True,
    axes=None,
    text_props: Mapping[str, Any] | None = None,
):
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

    text_props = merge_text_kwds(text_props, defaults=DEFAULT_TEXT_PROPS)

    if "node" in locs:
        plot_nodes(
            graph,
            with_id="node" in with_id,
            text=text_props,
        )
    if "link" in locs:
        plot_links(
            graph,
            with_id="link" in with_id,
            as_arrow=True,
            text=text_props,
        )
    if "patch" in locs:
        plot_patches(
            graph,
            with_id="patch" in with_id,
            text=text_props,
        )

    if "corner" in locs:
        plot_nodes(
            graph.dual,
            with_id="corner" in with_id,
            text=text_props,
        )
    if "face" in locs:
        plot_links(
            graph.dual,
            with_id="face" in with_id,
            text=text_props,
        )
    if "cell" in locs and getattr(graph, "number_of_cells", 0) > 0:
        plot_patches(graph.dual, with_id="cell" in with_id, text=text_props)

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


def _merge_kwds(
    base: Mapping[str, Any] | None = None,
    *,
    defaults: Mapping[str, Any] | Iterable[tuple[str, Any]] | None = None,
    aliases: Sequence[tuple[str, str]] | None = None,
) -> dict[str, Any]:
    return {
        **_norm_dict_with_aliases(defaults or {}, aliases=aliases),
        **_norm_dict_with_aliases(base or {}, aliases=aliases),
    }


def _norm_dict_with_aliases(
    base: Mapping[str, Any] | Iterable[tuple[str, Any]],
    *,
    aliases: Sequence[tuple[str, str]] | None = None,
) -> dict[str, Any]:
    base_map = dict(base or {})
    if not aliases:
        return base_map

    alias_map = {alias: primary for primary, alias in aliases}
    for alias, primary in alias_map.items():
        if alias in base_map and primary in base_map:
            raise ValueError(
                f"conflicting keys. Both {primary!r} and its alias, {alias!r}, provided"
            )

    return {alias_map.get(k, k): v for k, v in base_map.items()}


merge_text_kwds = partial(
    _merge_kwds,
    aliases=[
        ("color", "c"),
        ("fontsize", "size"),
        ("fontstyle", "style"),
    ],
)
