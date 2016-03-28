import numpy as np


def get_nodes_at_patch(graph):
    """Set up data structure that describes node-patch connectivity.

    Parameters
    ----------
    links_at_patch: ndarray
        Links that define each patch.
    nodes_at_link: ndarray
        Nodes that define each link.

    Returns
    -------
    ndarray
        Nodes that define each patch.
    """
    nodes_at_patch = np.full(graph.links_at_patch.shape, -1, dtype=int)

    for patch, links in enumerate(graph.links_at_patch):
        unique_nodes = np.unique(graph.nodes_at_link[links])
        nodes_at_patch[patch, :len(unique_nodes)] = unique_nodes

    return nodes_at_patch
