import numpy as np

from .ext.at_patch import get_nodes_at_patch as _get_nodes_at_patch


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

    _get_nodes_at_patch(graph.links_at_patch, graph.nodes_at_link, nodes_at_patch)

    return nodes_at_patch
