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
from ...utils.decorators import cache_result_in_object
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


def link_at_ridge(ridge_points, ridge_vertices, boundary_points):
    link_at_ridge = np.full(len(ridge_points), -1, dtype=int)
    link = 0
    for ridge, points in enumerate(ridge_points):
        if tuple(points) in boundary_points or -1 in ridge_vertices[ridge]:
            pass
        else:
            link_at_ridge[ridge] = link
            link += 1

    ridge_at_link = np.full(link, -1, dtype=int)
    for ridge, link in enumerate(link_at_ridge):
        if link != -1:
            ridge_at_link[link] = ridge
    return link_at_ridge, ridge_at_link


def node_at_vertex(ridge_points, ridge_vertices, boundary_points):
    link_at_ridge_, _ = link_at_ridge(ridge_points, ridge_vertices, boundary_points)
    good_ridges = np.where(link_at_ridge_ != -1)[0]

    # print link_at_ridge_

    # node_at_vertex = np.full(len(good_ridges), -1, dtype=int)
    node_at_vertex = np.full(np.max(ridge_vertices) + 1, -1, dtype=int)
    # print ridge_vertices, np.max(ridge_vertices)
    node = 0
    for ridge in good_ridges:
        for vertex in ridge_vertices[ridge]:
            if node_at_vertex[vertex] == -1:
                node_at_vertex[vertex] = node
                node += 1

    vertex_at_node = np.full(node, -1, dtype=int)
    for vertex, node in enumerate(node_at_vertex):
        if node != -1:
            vertex_at_node[node] = vertex

    return node_at_vertex, vertex_at_node


def patch_at_region(point_region, boundary_points):
    # patch_at_region = np.full(len(point_region), -1,
    # patch_at_region = np.full(len(point_region) - len(boundary_points), -1,
    patch_at_region = np.full(np.max(point_region) + 1, -1,
                              dtype=int)
    patch = 0
    for point, region in enumerate(point_region):
        if point in boundary_points:
            pass
        else:
            patch_at_region[region] = patch
            patch += 1

    region_at_patch = np.full(patch, -1, dtype=int)
    for region, patch in enumerate(patch_at_region):
        if patch != -1:
            region_at_patch[patch] = region

    return patch_at_region, region_at_patch


def corner_at_region(point_region, boundary_points):
    corner_at_region = np.full(len(point_region), -1, dtype=int)
    corner = 0
    for point in range(len(boundary_points)):
        if point not in boundary_points:
            patch_at_region[point] = corner
            corner += 1
    return corner_at_region


class VoronoiConverter(object):

    """Convert scipy.spatial.Voronoi data structures to landlab.

    points == xy_of_node
    vertices == xy_of_corner
    regions == corners_at_cell
    ridge_vertices == corners_at_face
    ridge_points == nodes_at_face
    point_region == node_at_cell

        # clean up:
        # * ridge_points
        # * ridge_vertices
        # * point_region
        # * regions (region_vertices)
        # * vertices

    corners = converter.get_nodes()
    corners = (corners[:, 1], corners[:, 0])
    faces = converter.get_nodes_at_link()
    cells = converter.get_links_at_patch()
    if len(cells) > 0:
        cells = [cell for cell in JaggedArray(cells) if -1 not in cell]

    node_at_cell = converter.get_corner_at_patch()
    nodes_at_face = converter.get_corners_at_link()
    """

    def __init__(self, voronoi, min_patch_size=3, boundary_nodes=None):
        if boundary_nodes is None:
            boundary_nodes = []

        self._voronoi = voronoi
        self._min_patch_size = min_patch_size
        self._boundary_points = set()
        for pair in boundary_nodes:
            self._boundary_points.add(tuple(pair))
            self._boundary_points.add(tuple(pair[::-1]))

            # NOTE: Fix this
            self._boundary_points.add(pair[0])
            self._boundary_points.add(pair[1])

        self.link_at_ridge, self.ridge_at_link = link_at_ridge(
            self._voronoi.ridge_points,
            self._voronoi.ridge_vertices,
            self._boundary_points)
        self.node_at_vertex, self.vertex_at_node = node_at_vertex(
            self._voronoi.ridge_points, self._voronoi.ridge_vertices,
            self._boundary_points)
        self.patch_at_region, self.region_at_patch = patch_at_region(
            self._voronoi.point_region,
            self._boundary_points)

    @property
    @cache_result_in_object()
    def vertices_at_ridge(self):
        return np.array(self._voronoi.ridge_vertices)

    @property
    @cache_result_in_object()
    def points_at_ridge(self):
        return np.array(self._voronoi.ridge_points)

    @property
    @cache_result_in_object()
    def xy_at_vertex(self):
        return np.array(self._voronoi.vertices)

    @property
    @cache_result_in_object()
    def xy_at_node(self):
        return self.xy_at_vertex[self.vertex_at_node]

    @property
    @cache_result_in_object()
    def nodes_at_link(self):
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
        self.ridge_at_link
        # print(self.ridge_at_link)
        # print(self.vertices_at_ridge)
        self.vertices_at_ridge[self.ridge_at_link]
        self.node_at_vertex[self.vertices_at_ridge[self.ridge_at_link]]
        return self.node_at_vertex[self.vertices_at_ridge[self.ridge_at_link]]

    # def nodes_at_patch(self):
    #     return self.node_at_vertex[self._voronoi.regions[self.region_at_patch]]

    @property
    @cache_result_in_object()
    def nodes_at_patch(self): 
        nodes_at_patch = []
        for region in self.region_at_patch:
            nodes_at_patch.append(self.node_at_vertex[self.regions[region]])
        return nodes_at_patch

    def corner_at_patch(self):
        return self.patch_at_region[self._voronoi.point_region[self.point_at_corner]]

    def corners_at_link(self):
        return self.point_at_corner[self._voronoi.ridge_points[self.ridge_at_link]]

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
        >>> import numpy as np
        >>> from scipy.spatial import Voronoi
        >>> from landlab.graph.voronoi.voronoi_helpers import VoronoiConverter

        >>> node_x = (0., 1., 2., 3., .1, 1.1, 2.1, 3.1, .2, 1.2, 2.2, 3.2)
        >>> node_y = (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.)
        >>> voronoi = Voronoi(list(zip(node_x, node_y)))

        >>> converter = VoronoiConverter(voronoi)
        >>> converter.is_patch([1, 2, 3])
        True
        >>> converter.is_patch([1, 2, 3, -1])
        False
        >>> converter.is_patch([])
        False

        Specify the minimum number sides for a valid patch.

        >>> converter = VoronoiConverter(voronoi, min_patch_size=4)
        >>> converter.is_patch([1, 2, 3])
        True
        """
        # print len(region) >= self._min_patch_size, region
        return len(region) > 0 and -1 not in region
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
        return tuple(self._voronoi.ridge_points[ridge]) not in self._boundary_points

        # return True
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
        # print 'get_patch_at_region'
        patch = 0
        for n, region in enumerate(self._voronoi.regions):
            if self.is_patch(region):
                patch_at_region[n] = patch
                patch += 1
            else:
                patch_at_region[n] = -1
        # print 'patch_at_region', patch_at_region
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
            # if -1 not in ridge and self.is_link(n):
            if -1 not in ridge:
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
        for vertices in self._voronoi.ridge_vertices:
            for vertex in vertices:
                if vertex != -1 and node_at_vertex[vertex] == -1:
                    node_at_vertex[vertex] = node
                    node += 1
        return node_at_vertex

        node = 0
        for region in self._voronoi.regions:
            # if self.is_patch(region):
            if True:
                for vertex in region:
                    if node_at_vertex[vertex] == -1:
                        node_at_vertex[vertex] = node
                        node += 1

        return node_at_vertex

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
        return self.xy_at_vertex[self.vertex_at_node]

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
            if len(vertices) > 0 and self.is_patch(vertices):
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

    @property
    @cache_result_in_object()
    def link_at_nodes(self):
        link_at_nodes = {}
        for link, nodes in enumerate(self.nodes_at_link):
            nodes = list(nodes)
            nodes.sort()
            link_at_nodes[tuple(nodes)] = link
        return link_at_nodes

    @property
    @cache_result_in_object()
    def links_at_patch(self):
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
        link_at_nodes = self.link_at_nodes
        nodes_at_patch = self.nodes_at_patch

        # nodes_at_patch = []
        # for region in self.region_at_patch:
        #     nodes_at_patch.append(self.node_at_vertex[self._voronoi.regions[region]])

        links_at_patch = []
        for patch in nodes_at_patch:
            links = []
            for t, h in zip(patch[:-1], patch[1:]):
                pair = [t, h]
                pair.sort()
                links.append(link_at_nodes[tuple(pair)])
            pair = [patch[0], patch[-1]]
            pair.sort()
            links.append(link_at_nodes[tuple(pair)])
            links_at_patch.append(links)

        return links_at_patch

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

        if n_patches > 0:
            region_at_patch = np.argsort(self.get_patch_at_region())[- n_patches:]
            point_at_region = np.argsort(self._voronoi.point_region)

            return point_at_region[region_at_patch - 1]
        else:
            return np.array([], dtype=int)

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
        return self.points_at_ridge[self.ridge_at_link]
        # return self._voronoi.ridge_points[self.ridge_at_link]


        points_at_ridge = self._voronoi.ridge_points
        link_at_ridge = self.get_link_at_ridge()
        return points_at_ridge[link_at_ridge >= 0]
