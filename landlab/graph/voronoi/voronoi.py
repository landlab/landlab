import numpy as np
import xarray as xr
from scipy.spatial import Delaunay, Voronoi

from ...core.utils import as_id_array
from ..graph import Graph
from ..ugrid import (
    MESH_ATTRS,
    update_links_at_patch,
    update_node_coords,
    update_nodes_at_link,
)
from .voronoi_helpers import VoronoiConverter
from .voronoi_to_graph import VoronoiDelaunayToGraph


def remove_bad_patches(
    max_node_spacing,
    nodes_at_patch,
    neighbors_at_patch,
    boundary_nodes=None,
):
    from .ext.delaunay import remove_tris

    bad_patches = []
    if boundary_nodes is not None and len(boundary_nodes) > 0:
        boundary_nodes = set(boundary_nodes)
        for patch, nodes in enumerate(nodes_at_patch):
            if set(nodes).issubset(boundary_nodes):
                bad_patches.append(patch)
    bad_patches = as_id_array(bad_patches)

    if len(bad_patches) > 0:
        remove_tris(nodes_at_patch, neighbors_at_patch, bad_patches)
        nodes_at_patch = nodes_at_patch[: -len(bad_patches), :]
        neighbors_at_patch = neighbors_at_patch[: -len(bad_patches), :]

    return nodes_at_patch, neighbors_at_patch


def setup_links_and_patches(node_y_and_x, max_node_spacing=None, boundary_nodes=None):
    from .ext.delaunay import _setup_links_at_patch

    delaunay = Delaunay(list(zip(node_y_and_x[1], node_y_and_x[0])))

    nodes_at_patch = np.asarray(delaunay.simplices, dtype=int)
    neighbors_at_patch = np.asarray(delaunay.neighbors, dtype=int)

    if boundary_nodes is not None:
        nodes_at_patch, neighbors_at_patch = remove_bad_patches(
            max_node_spacing, nodes_at_patch, neighbors_at_patch, boundary_nodes=boundary_nodes
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


def ugrid_from_voronoi(node_y_and_x, max_node_spacing=None, boundary_nodes=None):
    ugrid = xr.Dataset({"mesh": xr.DataArray(data=1, attrs=MESH_ATTRS)})

    nodes_at_link, links_at_patch = setup_links_and_patches(
        node_y_and_x, max_node_spacing=max_node_spacing, boundary_nodes=boundary_nodes,
    )

    update_node_coords(ugrid, node_y_and_x)
    update_nodes_at_link(ugrid, nodes_at_link)
    update_links_at_patch(ugrid, links_at_patch)

    return ugrid


class DelaunayGraph(Graph):

    """Graph of a voronoi grid.

    Examples
    --------
    >>> from landlab.graph import DelaunayGraph
    """
    def __init__(self, node_y_and_x, max_node_spacing=None, sort=False, perimeter_links=None):
        """Create a voronoi grid.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of every node. First *y*, then *x*.

        Examples
        --------
        >>> from landlab.graph import DelaunayGraph
        >>> node_x = [0.0, 1.0, 2.0,
        ...           0.9, 1.9, 2.9]
        >>> node_y = [0, 0, 0,
        ...           2, 2, 2]
        >>> graph = DelaunayGraph((node_y, node_x), sort=True)
        >>> graph.x_of_node
        array([ 0. ,  1. ,  2. ,  0.9,  1.9,  2.9])
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
        array([[3, 2, 0], [5, 4, 1], [7, 3, 4], [8, 5, 6]])

        >>> graph.nodes_at_patch # doctest: +NORMALIZE_WHITESPACE
        array([[3, 0, 1], [4, 1, 2], [4, 3, 1], [5, 4, 2]])
        """
        mesh = VoronoiDelaunayToGraph(
            np.vstack((node_y_and_x[1], node_y_and_x[0])).T,
            perimeter_links=perimeter_links,
            # perimeter_links=kwds.get("perimeter_links", None),
        )

        # nodes_at_link, links_at_patch = DelaunayGraph.get_links_and_patches(
        #     node_y_and_x, max_node_spacing=max_node_spacing
        # )
        # mesh.links_at_patch

        Graph.__init__(
            self,
            node_y_and_x,
            links=mesh.nodes_at_link,
            patches=mesh.links_at_patch,
            sort=sort,
        )

    @staticmethod
    def get_links_and_patches(node_y_and_x, max_node_spacing=None):
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
