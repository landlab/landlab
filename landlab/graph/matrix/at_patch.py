import numpy as np

from ...utils.jaggedarray import unravel


def links_at_patch(patches, sort=True, nodes_at_link=None):
    """Construct as links_at_patch array for a graph.

    Parameters
    ----------
    patches : tuple of ndarray of int
        Links that define patches as `(links, offset_to_patch)`.
    sort : bool, optional
        Sort the links.
    nodes_at_link : ndarray of int, shape `(n_links, 2)`
        Nodes at link tail and head.

    Examples
    --------

    """
    from ..sort.ext.remap_element import reorder_links_at_patch
    from ..quantity.ext.of_link import calc_midpoint_of_link

    links_at_patch, offset_to_patch = patches

    sort = False
    # print 'SORTING TURNED OFF'
    if sort:
        xy_of_link = np.empty((len(nodes_at_link), 2), dtype=float)
        calc_midpoint_of_link(nodes_at_link, x_of_node, y_of_node, xy_of_link)
        reorder_links_at_patch(links_at_patch, offset_to_patch, xy_of_link)

        sort_patches(links_at_patch, offset_to_patch, xy_of_link)

    return unravel(links_at_patch, offset_to_patch, pad=-1)
