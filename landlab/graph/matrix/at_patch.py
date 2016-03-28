import numpy as np

from ...utils.jaggedarray import unravel


def new_links_at_patch(patches, nodes_at_link=None, nodes_at_patch=None):
    links_at_patch, offset_to_patch = patches

    connect_links_at_patch(links_at_patch, offset_to_patch, nodes_at_link,
                           nodes_at_patch)


def links_at_patch(patches, sort=True, nodes_at_link=None):
    """Set up patch data structures."""
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


def old_links(patches):
    """Set up data structure that describes link-patch connectivity.

    Parameters
    ----------
    patches: iterable of iterables
        List of links for each patch.

    Returns
    -------
    ndarray
        Links for each patch.
    """
    from .ext.at_patch import fill_links_at_patch

    max_links_per_patch = np.diff(patches[1]).max()
    n_patches = len(patches[1]) - 1
    links_at_patch = np.full((n_patches, max_links_per_patch), -1, dtype=int)

    fill_links_at_patch(patches[0], patches[1], links_at_patch)

    return links_at_patch
