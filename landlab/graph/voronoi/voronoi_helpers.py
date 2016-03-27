import numpy as np


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
    vertices_per_region = np.array(
        [len(region) for region in regions], dtype=int)

    return vertices_at_region, vertices_per_region


def get_finite_regions(voronoi):
    """Get regions of finite area.

    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import get_finite_regions
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> get_finite_regions(voronoi)
    array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
    """
    from .ext.voronoi import _is_finite_region

    n_regions = len(voronoi.regions)
    vertices_at_region, vertices_per_region = flatten_vertices_at_region(
        voronoi.regions)

    is_finite_region = np.empty(n_regions, dtype=int)
    _is_finite_region(vertices_at_region, vertices_per_region,
                      is_finite_region)

    return is_finite_region


def get_neighbor_regions(voronoi):
    """Get regions on either side of each ridge.

    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import get_neighbor_regions
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> get_neighbor_regions(voronoi)
    array([[-1, -1],
           [-1,  8],
           [-1, -1],
           [-1, -1],
           [-1, -1],
           [-1,  8],
           [-1, -1],
           [-1, -1],
           [-1,  8],
           [-1, -1],
           [-1, -1],
           [-1,  8],
           [-1, -1],
           [-1,  8],
           [-1, -1],
           [ 8, -1]])
    """
    from .ext.voronoi import _get_neighbor_regions

    n_ridges = len(voronoi.ridge_vertices)
    is_finite_region = get_finite_regions(voronoi)

    regions_at_ridge = np.empty((n_ridges, 2), dtype=int)
    _get_neighbor_regions(np.array(voronoi.ridge_points, dtype=int),
                          np.array(voronoi.point_region, dtype=int),
                          is_finite_region, regions_at_ridge)

    return regions_at_ridge


def get_ridges_at_cell(voronoi):
    """Get the ridges that define a cell.

    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import get_ridges_at_cell
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> ridges_at_cell, cell_at_region = get_ridges_at_cell(voronoi)
    >>> ridges_at_cell
    array([[ 1,  5,  8, 11, 13, 15]])
    >>> cell_at_region
    array([-1, -1, -1, -1, -1, -1, -1, -1,  0, -1])
    """
    from .ext.voronoi import _get_cell_at_region

    n_regions = len(voronoi.regions)

    regions_at_ridge = get_neighbor_regions(voronoi)

    max_vertices_per_region = max([len(region) for region in voronoi.regions])
    ridges_at_cell = np.empty((n_regions, max_vertices_per_region), dtype=int)
    cell_at_region = np.empty(n_regions, dtype=int)

    n_cells = _get_cell_at_region(regions_at_ridge, ridges_at_cell,
                                  cell_at_region)

    return ridges_at_cell[:n_cells, :], cell_at_region


def get_faces_at_cell(voronoi):
    """Find the faces that make up each cell.

    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import get_faces_at_cell
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> faces_at_cell, face_at_ridge, cell_at_region = get_faces_at_cell(voronoi)
    >>> faces_at_cell
    array([[0, 1, 2, 3, 4, 5]])
    >>> face_at_ridge
    array([-1,  0, -1, -1, -1,  1, -1, -1,  2, -1, -1,  3, -1,  4, -1,  5])
    >>> cell_at_region
    array([-1, -1, -1, -1, -1, -1, -1, -1,  0, -1])
    """
    from .ext.voronoi import _get_faces_at_cell

    n_regions = len(voronoi.regions)
    n_ridges = len(voronoi.ridge_vertices)
    max_vertices_per_region = max([len(region) for region in voronoi.regions])

    ridges_at_cell, cell_at_region = get_ridges_at_cell(voronoi)
    n_cells = len(ridges_at_cell)

    faces_at_cell = np.empty((n_cells, max_vertices_per_region), dtype=int)
    face_at_ridge = np.empty(n_ridges, dtype=int)
    _get_faces_at_cell(ridges_at_cell, faces_at_cell, face_at_ridge)

    return faces_at_cell, face_at_ridge, cell_at_region


def get_corners_at_face(voronoi, face_at_ridge):
    """Get the two corners at either end of each face.

    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import (get_faces_at_cell,
    ...                                            get_corners_at_face)
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.], [3. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.], [3.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.], [3.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> _, face_at_ridge, _ = get_faces_at_cell(voronoi)
    >>> corners_at_face, corner_at_vertex = get_corners_at_face(voronoi, face_at_ridge)
    >>> corners_at_face
    array([[0, 1],
           [3, 1],
           [0, 6],
           [8, 9],
           [6, 8],
           [3, 9],
           [0, 2],
           [4, 5],
           [2, 5],
           [4, 7],
           [6, 7]])
    >>> corner_at_vertex
    array([ 3,  0,  2,  1,  4, -1,  5,  6,  7,  8, -1,  9])
    """
    from .ext.voronoi import _get_corners_at_face

    n_vertices = len(voronoi.vertices)
    n_faces = max(face_at_ridge) + 1

    corner_at_vertex = np.empty(n_vertices, dtype=int)
    corners_at_face = np.empty((n_faces, 2), dtype=int)

    _get_corners_at_face(face_at_ridge,
                         np.array(voronoi.ridge_vertices, dtype=int),
                         corner_at_vertex, corners_at_face)

    return corners_at_face, corner_at_vertex


def get_xy_of_corner(voronoi, corner_at_vertex):
    """Get the x and y position of each corner.

    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import (get_faces_at_cell,
    ...                                            get_corners_at_face,
    ...                                            get_xy_of_corner)
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> _, face_at_ridge, _ = get_faces_at_cell(voronoi)
    >>> _, corner_at_vertex = get_corners_at_face(voronoi, face_at_ridge)
    >>> get_xy_of_corner(voronoi, corner_at_vertex)
    array([[ 0.6  ,  1.455],
           [ 0.7  ,  1.545],
           [ 1.5  ,  0.455],
           [ 0.6  ,  0.545],
           [ 1.6  ,  0.545],
           [ 1.6  ,  1.455]])
    """
    from .ext.voronoi import _get_xy_at_corners

    n_corners = max(corner_at_vertex) + 1
    xy_of_corner = np.empty((n_corners, 2), dtype=float)
    n_corners = _get_xy_at_corners(voronoi.vertices, corner_at_vertex, xy_of_corner)

    return xy_of_corner


def get_node_at_cell(voronoi, cell_at_region):
    """Get node for each cell.

    Examples
    --------
    >>> from scipy.spatial import Voronoi
    >>> from landlab.graph.voronoi.voronoi_helpers import (get_faces_at_cell,
    ...                                            get_node_at_cell)
    >>> points = [[0. , 0.], [1. , 0.], [2. , 0.],
    ...           [0.1, 1.], [1.1, 1.], [2.1, 1.],
    ...           [0.2, 2.], [1.2, 2.], [2.2, 2.]]
    >>> voronoi = Voronoi(points)
    >>> _, _, cell_at_region = get_faces_at_cell(voronoi)
    >>> get_node_at_cell(voronoi, cell_at_region)
    array([4])
    """
    from .ext.voronoi import _get_node_at_cell

    n_cells = max(cell_at_region) + 1
    node_at_cell = np.empty(n_cells, dtype=int)

    _get_node_at_cell(np.array(voronoi.ridge_points, dtype=int),
                      np.array(voronoi.point_region, dtype=int),
                      cell_at_region, node_at_cell)

    return node_at_cell


def setup_voronoi_connectivity(voronoi):
    faces_at_cell, face_at_ridge, cell_at_region = get_faces_at_cell(voronoi)

    corners_at_face, corner_at_vertex = get_corners_at_face(voronoi,
                                                            face_at_ridge)

    xy_at_corner = get_xy_of_corner(voronoi, corner_at_vertex)

    node_at_cell = get_node_at_cell(voronoi, cell_at_region)

    return faces_at_cell, corners_at_face, xy_at_corner, node_at_cell
