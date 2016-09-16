import numpy as np
import matplotlib.pyplot as plt


def plot_nodes(graph, color='r'):
    for node in range(len(graph.x_of_node)):
        x, y = graph.x_of_node[node], graph.y_of_node[node]
        plt.plot(graph.x_of_node[node], graph.y_of_node[node], 'o',
                 color=color)
        plt.text(x, y, node, color=color, size=16)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect(1.)


def plot_links(graph, color='b', linestyle='solid'):
    for link, nodes in enumerate(graph.nodes_at_link):
        x, y = graph.x_of_node[nodes[0]], graph.y_of_node[nodes[0]]
        dx, dy = graph.x_of_node[nodes[1]] - x, graph.y_of_node[nodes[1]] - y
        plt.arrow(x, y, dx, dy, head_width=.1, length_includes_head=True,
                  color=color, linestyle=linestyle)
        plt.text(x + dx * .5, y + dy * .5, link, size=16, color=color)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect(1.)


def plot_patches(graph, color='g'):
    for patch, nodes in enumerate(graph.nodes_at_patch):
        x, y = np.mean(graph.x_of_node[nodes]), np.mean(graph.y_of_node[nodes])
        plt.text(x, y, patch, color=color, size=16)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect(1.)


def plot_graph(graph, at='node,link,patch'):
    locs = [loc.strip() for loc in at.split(',')]
    for loc in locs:
        if loc not in ('node', 'link', 'patch', 'corner', 'face', 'cell'):
            raise ValueError(
                '{at}: "at" element not understood'.format(at=loc))

    plt.plot(graph.x_of_node, graph.y_of_node, '.', color='r')

    if 'node' in locs:
        plot_nodes(graph)
    if 'link' in locs:
        plot_links(graph)
    if 'patch' in locs:
        plot_patches(graph)

    if 'corner' in locs:
        plot_nodes(graph.dual, color='c')
    if 'face' in locs:
        plot_links(graph.dual, linestyle='dotted', color='k')
    if 'cell' in locs:
        plot_patches(graph.dual, color='m')

    plt.show()
