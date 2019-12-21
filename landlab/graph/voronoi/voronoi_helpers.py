"""Helper functions that work with Voronoi data structures.

This module provides functions to help work with the Voronoi data
structures and to map them to more sensible data structures.

The terminology used within this module can sometimes be confusing.

region :
    Polygon definied by Voronoi. These polygons can be unbound (infinite in
    area).
vertex :
    The points that define the Voronoi regions.
ridge :
    Lines that join adjacent vertices. Ridges are the edges of Voronoi
    regions.
point :
    Each region has an associated point, which lies inside the region. Each
    region has only one point, and each point has only on region.
patch :
    Similar to a region but is always bound.
link :
    Similar to a ridge but is always part of at least one patch.
node :
    Similar to a vertex but bounds only links.
corner :
    Same as a point.
"""
import numpy as np
from six.moves import range

from ...core.utils import as_id_array
from ..sort.sort import remap


def flatten_vertices_at_region(regions):
    """
    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import flatten_vertices_at_region
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> (vertices, count) = flatten_vertices_at_region(voronoi.regions)
    >>> vertices
    array([ 1,  0, -1,  4,  2, -1,  3,  3, -1,  4,  0, -1,  3,  7,  5, -1,  6,
            5,  2, -1,  6, -1,  7,  1,  0,  4,  2,  5,  7,  1, -1,  6])
    >>> count
    array([3, 0, 4, 2, 4, 4, 3, 2, 6, 4])
    """
    vertices_at_region = np.array(np.concatenate(regions), dtype=int)
    vertices_per_region = np.array([len(region) for region in regions], dtype=int)

    return vertices_at_region, vertices_per_region


class VoronoiConverter(object):
    def __init__(self, voronoi, min_patch_size=3):
        self._voronoi = voronoi
        self._min_patch_size = min_patch_size

    @property
    def n_regions(self):
        return len(self._voronoi.regions)

    @property
    def regions(self):
        return self._voronoi.regions

    @property
    def max_patch_size(self):
        try:
            return self._max_patch_size
        except AttributeError:
            self._max_patch_size = max([len(r) for r in self.regions])
            return self._max_patch_size

    def get_finite_regions(self):
        """Get regions of finite area.

        Examples
        --------
        >>> from scipy.spatial import Voronoi
        >>> from landlab.graph.voronoi.voronoi_helpers import VoronoiConverter

        >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
        ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
        ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
        >>> voronoi = Voronoi(points)

        >>> converter = VoronoiConverter(voronoi)
        >>> converter.get_finite_regions()
        array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
        """
        try:
            return self._is_finite_region
        except AttributeError:
            from .ext.voronoi import _is_finite_region

            n_regions = self.n_regions
            vertices_at_region, vertices_per_region = flatten_vertices_at_region(
                self.regions
            )

            is_finite_region = np.empty(n_regions, dtype=int)
            _is_finite_region(
                vertices_at_region,
                vertices_per_region,
                is_finite_region,
                self._min_patch_size,
            )

            self._is_finite_region = is_finite_region
            return self._is_finite_region

    def is_patch(self, region):
        """Check if a voronoi region in bounded.

        Parameters
        ----------
        region : array_like of int
            Identifiers for the vertices of a region.

        Examples
        --------
        >>> from landlab.graph.voronoi.voronoi_helpers import VoronoiConverter

        >>> converter = VoronoiConverter(None)
        >>> converter.is_patch([1, 2, 3])
        True
        >>> converter.is_patch([1, 2, 3, -1])
        False
        >>> converter.is_patch([])
        False

        Specify the minimum number sides for a valid patch.

        >>> converter = VoronoiConverter(None, min_patch_size=4)
        >>> converter.is_patch([1, 2, 3])
        False
        """
        return len(region) >= self._min_patch_size and -1 not in region

    def is_link(self, ridge):
        """Check if a voronoi ridge is a valid link.

        A ridge is a valid link if at least one of its neighbor regions is
        bound.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.
        ridge : int
            Identifier of a ridge in the voronoi graph.

        Returns
        -------
        bool
            True if the ridge is a link.
        """
        for point in self._voronoi.ridge_points[ridge]:
            region = self._voronoi.point_region[point]
            if self.is_patch(self._voronoi.regions[region]):
                return True
        return False

    def get_patch_at_region(self):
        """Map voronoi regions to patches.

        Regions that do not have corresponding patches, are -1.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_regions, 1)`
            Array of patch identifiers for every voronoi region.
        """
        patch_at_region = np.empty(len(self._voronoi.regions), dtype=int)

        patch = 0
        for n, region in enumerate(self._voronoi.regions):
            if self.is_patch(region):
                patch_at_region[n] = patch
                patch += 1
            else:
                patch_at_region[n] = -1

        return patch_at_region

    def get_link_at_ridge(self):
        """Map voronoi ridges to links.

        Ridges that do not have corresponding links, are -1.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_ridges, 1)`
            Array of link identifiers for every voronoi ridge.
        """
        link_at_ridge = np.empty(len(self._voronoi.ridge_vertices), dtype=int)

        link = 0
        for n, ridge in enumerate(self._voronoi.ridge_vertices):
            if -1 not in ridge and self.is_link(n):
                link_at_ridge[n] = link
                link += 1
            else:
                link_at_ridge[n] = -1

        return link_at_ridge

    def get_patches_at_link(self):
        """Get patches on either side of each link.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_links, 2)`
            Array of patch identifiers for every link.
        """
        points_at_ridge = self._voronoi.ridge_points
        region_at_point = self._voronoi.point_region

        regions_at_ridge = as_id_array(region_at_point[points_at_ridge])

        patch_at_region = self.get_patch_at_region()
        link_at_ridge = self.get_link_at_ridge()

        patch_at_ridge = remap(regions_at_ridge, patch_at_region, inplace=True)

        return patch_at_ridge[link_at_ridge >= 0]

    def get_node_at_vertex(self):
        """Map voronoi vertices to nodes.

        Vertices that do not have corresponding nodes, are -1.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_vertices, 1)`
            Array of node identifiers for every voronoi vertex.
        """
        node_at_vertex = np.full(len(self._voronoi.vertices), -1, dtype=int)

        node = 0
        for region in self._voronoi.regions:
            if self.is_patch(region):
                for vertex in region:
                    if node_at_vertex[vertex] == -1:
                        node_at_vertex[vertex] = node
                        node += 1

        return node_at_vertex

    def get_nodes_at_link(self):
        """Get end nodes for links of a voronoi graph.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_links, 2)`
            Array of node identifiers for every link.
        """
        node_at_vertex = self.get_node_at_vertex()
        link_at_ridge = self.get_link_at_ridge()
        vertices_at_ridge = np.array(self._voronoi.ridge_vertices)

        nodes_at_link = node_at_vertex[vertices_at_ridge[link_at_ridge >= 0].flat]
        nodes_at_link.shape = (-1, 2)

        return nodes_at_link

    def get_nodes(self):
        """Get nodes of a voronoi graph.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of float, shape `(n_nodes, 2)`
            Array of coordinates `(x, y)` for every node.
        """
        node_at_vertex = self.get_node_at_vertex()
        nodes = np.argsort(node_at_vertex[node_at_vertex >= 0])

        return self._voronoi.vertices[node_at_vertex >= 0][nodes]

    def get_ridges_at_region(self):
        """Get ridges that bound each voronoi region.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_regions, max_ridges_per_region)`
            Ridge identifiers that define each region (padded with -1).
        """
        n_regions = self.n_regions
        max_links_per_patch = self.max_patch_size

        ridges_at_region = np.full((n_regions, max_links_per_patch), -1, dtype=int)

        ridge_at_vertices = {}
        for ridge, vertices in enumerate(self._voronoi.ridge_vertices):
            ridge_at_vertices[tuple(vertices)] = ridge

        for region, vertices in enumerate(self.regions):
            if self.is_patch(vertices):
                for n in range(len(vertices) - 1):
                    pair = [vertices[n], vertices[n + 1]]
                    pair.sort()

                    ridge = ridge_at_vertices[tuple(pair)]

                    ridges_at_region[region, n] = ridge
                pair = [vertices[len(vertices) - 1], vertices[0]]
                pair.sort()

                ridge = ridge_at_vertices[tuple(pair)]

                ridges_at_region[region, len(vertices) - 1] = ridge

        return ridges_at_region

    def get_links_at_patch(self):
        """Get links that bound each patch from a voronoi graph.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_patches, max_links_per_region)`
            Link identifiers that define each patch (padded with -1).
        """
        from ..sort.ext.remap_element import remap_graph_element_ignore

        ridges_at_region = self.get_ridges_at_region()
        link_at_ridge = self.get_link_at_ridge()

        patches = self.get_patch_at_region() >= 0
        links_at_patch = ridges_at_region[patches]

        remap_graph_element_ignore(links_at_patch.reshape((-1,)), link_at_ridge, -1)
        return links_at_patch

    def get_corner_at_patch(self):
        """Get corner that is contained within each patch from a voronoi graph.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_patches, )`
            Corner identifiers contained in patches.
        """
        finite_regions = self.get_finite_regions()
        n_patches = finite_regions.sum()

        region_at_patch = np.argsort(self.get_patch_at_region())[-n_patches:]
        point_at_region = np.argsort(self._voronoi.point_region)

        return point_at_region[region_at_patch - 1]

    def get_corners_at_link(self):
        """Get corners on either side of links.

        Parameters
        ----------
        voronoi : Voronoi
            A voronoi graph.

        Returns
        -------
        ndarray of int, shape `(n_links, 2)`
            Corner identifiers on either side of each link.
        """
        points_at_ridge = self._voronoi.ridge_points
        link_at_ridge = self.get_link_at_ridge()
        return points_at_ridge[link_at_ridge >= 0]
