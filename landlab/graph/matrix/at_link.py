import numpy as np


def nodes(links, sort=True, xy_of_node=None):
    """Get nodes at link.

    Parameters
    ----------
    links : array_like
        Node identifiers for each link tail and head.
    sort : bool
        
    """
    nodes_at_link = np.asarray(links, dtype=np.int)

    if sort:
        from ..ext.remap_element import reorient_links
        from ..sort import sort_links

        if xy_of_node is None and not isinstance(xy_of_node, np.ndarray):
            raise ValueError('xy_of_node must be ndarray')

        reorient_links(nodes_at_link, xy_of_node)
        sort_links(nodes_at_link, (xy_of_node[:, 1], xy_of_node[:, 0]))

    return nodes_at_link
