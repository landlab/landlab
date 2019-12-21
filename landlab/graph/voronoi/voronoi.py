import numpy as np
import xarray as xr
from scipy.spatial import Delaunay

from ...core.utils import as_id_array
from ..graph import Graph
from ..ugrid import (
    MESH_ATTRS,
    update_links_at_patch,
    update_node_coords,
    update_nodes_at_link,
)


def remove_bad_patches(max_node_spacing, nodes_at_patch, neighbors_at_patch):
    from .ext.delaunay import remove_tris

    max_node_dist = np.ptp(nodes_at_patch, axis=1)
    bad_patches = as_id_array(np.where(max_node_dist > max_node_spacing)[0])

    if len(bad_patches) > 0:
        remove_tris(nodes_at_patch, neighbors_at_patch, bad_patches)
        nodes_at_patch = nodes_at_patch[: -len(bad_patches), :]
        neighbors_at_patch = neighbors_at_patch[: -len(bad_patches), :]

    return nodes_at_patch, neighbors_at_patch


def setup_links_and_patches(node_y_and_x, max_node_spacing=None):
    from .ext.delaunay import _setup_links_at_patch

    delaunay = Delaunay(list(zip(node_y_and_x[1], node_y_and_x[0])))

    nodes_at_patch = np.asarray(delaunay.simplices, dtype=int)
    neighbors_at_patch = np.asarray(delaunay.neighbors, dtype=int)

    if max_node_spacing is not None:
        nodes_at_patch, neighbors_at_patch = remove_bad_patches(
            max_node_spacing, nodes_at_patch, neighbors_at_patch
        )

    n_patches = len(nodes_at_patch)
    n_shared_links = np.count_nonzero(neighbors_at_patch > -1)
    n_links = 3 * n_patches - n_shared_links // 2

    links_at_patch = np.empty((n_patches, 3), dtype=int)
    nodes_at_link = np.empty((n_links, 2), dtype=int)

    _setup_links_at_patch(
        nodes_at_patch, neighbors_at_patch, nodes_at_link, links_at_patch
    )

    return nodes_at_link, links_at_patch


def ugrid_from_voronoi(node_y_and_x, max_node_spacing=None):
    ugrid = xr.Dataset({"mesh": xr.DataArray(data=1, attrs=MESH_ATTRS)})

    nodes_at_link, links_at_patch = setup_links_and_patches(
        node_y_and_x, max_node_spacing=max_node_spacing
    )

    update_node_coords(ugrid, node_y_and_x)
    update_nodes_at_link(ugrid, nodes_at_link)
    update_links_at_patch(ugrid, links_at_patch)

    return ugrid


class VoronoiGraph(Graph):

    """Graph of a voronoi grid.

    Examples
    --------
    >>> from landlab.graph import VoronoiGraph
    """

    def __init__(self, node_y_and_x, **kwds):
        """Create a voronoi grid.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of every node. First *y*, then *x*.

        Examples
        --------
        >>> from landlab.graph import VoronoiGraph
        >>> node_x = [0, 1, 2,
        ...           1, 2, 3]
        >>> node_y = [0, 0, 0,
        ...           2, 2, 2]
        >>> graph = VoronoiGraph((node_y, node_x))
        >>> graph.x_of_node
        array([ 0.,  1.,  2.,  1.,  2.,  3.])
        >>> graph.y_of_node
        array([ 0.,  0.,  0.,  2.,  2.,  2.])
        >>> graph.nodes_at_link # doctest: +NORMALIZE_WHITESPACE
        array([[0, 1], [1, 2],
               [0, 3], [1, 3], [1, 4], [2, 4], [2, 5],
               [3, 4], [4, 5]])
        >>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  2, -1, -1], [ 1,  4,  3,  0], [ 6,  5,  1, -1],
               [ 7,  2,  3, -1], [ 8,  7,  4,  5], [ 8,  6, -1, -1]])
        >>> graph.links_at_patch # doctest: +NORMALIZE_WHITESPACE
        array([[3, 2, 0], [5, 4, 1], [4, 7, 3], [6, 8, 5]])
        >>> graph.nodes_at_patch # doctest: +NORMALIZE_WHITESPACE
        array([[3, 0, 1], [4, 1, 2], [4, 3, 1], [5, 4, 2]])
        """
        max_node_spacing = kwds.pop("max_node_spacing", None)
        mesh = ugrid_from_voronoi(node_y_and_x, max_node_spacing=max_node_spacing)
        Graph.__init__(self, mesh, **kwds)
