import matplotlib.pyplot as plt
import numpy as np


def plot_nodes(graph, color="r", with_id=True, markersize=10):
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


def plot_patches(graph, color="g"):
    for patch, nodes in enumerate(graph.nodes_at_patch):
        x, y = np.mean(graph.x_of_node[nodes]), np.mean(graph.y_of_node[nodes])
        plt.text(x, y, patch, color=color, size=16)


def plot_graph(graph, at="node,link,patch", with_id=True):
    locs = [loc.strip() for loc in at.split(",")]
    for loc in locs:
        if loc not in ("node", "link", "patch", "corner", "face", "cell"):
            raise ValueError('{at}: "at" element not understood'.format(at=loc))

    plt.plot(graph.x_of_node, graph.y_of_node, ".", color="r")
    plt.xlim([min(graph.x_of_node) - 0.5, max(graph.x_of_node) + 0.5])
    plt.ylim([min(graph.y_of_node) - 0.5, max(graph.y_of_node) + 0.5])

    if "node" in locs:
        plot_nodes(graph, with_id=with_id, markersize=10)
    if "link" in locs:
        plot_links(graph, with_id=with_id, linewidth=None, as_arrow=False)
    if "patch" in locs:
        plot_patches(graph)

    if "corner" in locs:
        plot_nodes(graph.dual, color="c")
    if "face" in locs:
        plot_links(graph.dual, linestyle="dotted", color="k")
    if "cell" in locs:
        plot_patches(graph.dual, color="m")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect(1.0)

    plt.show()
