import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free


DTYPE = np.int
ctypedef np.int_t DTYPE_t
ctypedef np.uint8_t uint8


cdef _id_array_contains(long *array, long size, long bad_id):
    cdef long n
    for n in range(size):
        if array[n] == bad_id:
            return 1
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
def id_array_contains(
    np.ndarray[long, ndim=2, mode="c"] corners_at_cell not None,
    np.ndarray[long, ndim=1, mode="c"] n_corners_at_cell not None,
    long bad_val,
    np.ndarray[uint8, ndim=1, mode="c"] out not None,
):
    cdef long n_cells = corners_at_cell.shape[0]
    cdef long cell

    for cell in range(n_cells):
        out[cell] = _id_array_contains(
            &corners_at_cell[cell, 0],
            n_corners_at_cell[cell],
            bad_val,
        )


def _is_finite_region(np.ndarray[DTYPE_t, ndim=1] vertices_at_region,
                      np.ndarray[DTYPE_t, ndim=1] vertices_per_region,
                      np.ndarray[DTYPE_t, ndim=1] is_finite_region,
                      DTYPE_t min_patch_size):
    """Test if each region if finite.

    Parameters
    ----------
    vertices_at_region : ndarray of int
        Indices of the Voronoi vertices forming each Voronoi region.
    vertices_per_region : ndarray of int, shape `(n_regions, )`
        Number of vertices that form each voronoi region.
    is_finite_region : ndarray of int, shape `(n_regions, )`
        Output buffer. 1 if a region is finite, otherwise 0.
    """
    cdef int i
    cdef int region
    cdef int vertex
    cdef int n_vertices
    cdef int offset = 0
    cdef int n_regions = len(vertices_per_region)

    for region in range(n_regions):
        n_vertices = vertices_per_region[region]
        if n_vertices == 0 or n_vertices < min_patch_size:
            is_finite_region[region] = 0
        else:
            is_finite_region[region] = 1
            for i in range(n_vertices):
                vertex = vertices_at_region[offset + i]
                if vertex == -1:
                    is_finite_region[region] = 0
                    break
        offset += n_vertices


def _get_neighbor_regions(np.ndarray[DTYPE_t, ndim=2] ridge_points,
                          np.ndarray[DTYPE_t, ndim=1] point_region,
                          np.ndarray[DTYPE_t, ndim=1] is_finite_region,
                          np.ndarray[DTYPE_t, ndim=2] regions_at_ridge):
    """Get voronoi regions on either side of ridges.

    Parameters
    ----------
    ridge_points : ndarray of int, shape `(n_points, 2)`
        Indices of the points between which each Voronoi ridge lies (as
        provided by the `scipy.spatial.Voronoi` class).
    point_region : ndarray of int, shape `(n_points, )`
        Index of the Voronoi region for each input point (as provided by the
        `scipy.spatial.Voronoi` class).
    is_finite_region : ndarray of int, shape `(n_regions, )`
        1 if a region is finite, otherwise 0.
    regions_at_ridge : ndarray of int, shape `(n_ridges, 2)`
        Indices of the voronoi regions on either side of each ridge.
    """
    cdef int i
    cdef int ridge
    cdef int region
    cdef int n_ridges = len(ridge_points)

    for ridge in range(n_ridges):
        for i in range(2):
            region = point_region[ridge_points[ridge, i]]
            if is_finite_region[region]:
                regions_at_ridge[ridge, i] = region
            else:
                regions_at_ridge[ridge, i] = -1


@cython.boundscheck(False)
def _get_cell_at_region(np.ndarray[DTYPE_t, ndim=2] regions_at_ridge,
                        np.ndarray[DTYPE_t, ndim=2] ridges_at_cell,
                        np.ndarray[DTYPE_t, ndim=1] cell_at_region):
    """Get cell corresponding to each voronoi region.

    Parameters
    ----------
    regions_at_ridge : ndarray of int, shape `(n_regions, 2)`
        Indices of the voronoi regions on either side of each ridge.
    ridges_at_cell : ndarray of int, shape `(n_cells, 2)`
        Output buffer. The voronoi regions on either side of each cell.
    cell_at_region : ndarray of int, shape `(n_regions, )`
        Output buffer. ID of the cell corresponding to each voronoi region.
        Regions without cells are set to -1.
    """
    cdef int i
    cdef int ridge
    cdef int cell
    cdef int region
    cdef int n_cells = 0
    cdef int n_ridges = len(regions_at_ridge)
    cdef int *ridges_per_cell = <int *>malloc(n_ridges * sizeof(int))

    if not ridges_per_cell:
        raise MemoryError('unable to allocate {bytes} bytes'.format(
            bytes=n_ridges * sizeof(int)))

    try:
        for ridge in range(n_ridges):
            ridges_per_cell[ridge] = 0

        for region in range(cell_at_region.shape[0]):
            cell_at_region[region] = -1

        for cell in range(ridges_at_cell.shape[0]):
            for ridge in range(ridges_at_cell.shape[1]):
                ridges_at_cell[cell, ridge] = -1

        for ridge in range(n_ridges):
            for i in range(2):
                region = regions_at_ridge[ridge, i]
                if region >= 0:
                    if cell_at_region[region] == -1:
                        cell_at_region[region] = n_cells
                        n_cells += 1
                    cell = cell_at_region[region]

                    ridges_at_cell[cell, ridges_per_cell[cell]] = ridge
                    ridges_per_cell[cell] += 1

    finally:
        free(ridges_per_cell)

    return n_cells


@cython.boundscheck(False)
def _get_faces_at_cell(np.ndarray[DTYPE_t, ndim=2] ridges_at_cell,
                       np.ndarray[DTYPE_t, ndim=2] faces_at_cell,
                       np.ndarray[DTYPE_t, ndim=1] face_at_ridge):
    """Get faces that define each cell.

    Parameters
    ----------
    ridges_at_cell : ndarray of int, shape `(n_cells, max_ridges)`
        Indices of the voronoi ridges that define each cell. Cells with less
        than `max_ridges` ridges are padded with -1.
    faces_at_cell : ndarray of int, shape `(n_cells, max_ridges)`
        Output buffer. IDs of faces that define each cell. Cells with less
        than `max_ridges` faces are padded with -1.
    face_at_ridge : ndarray of int, shape `(n_ridges, )`
        Output buffer. ID of face corresponding to each voronoi ridge. Ridges
        without a face are set to -1.
    """
    cdef int i
    cdef int face
    cdef int cell
    cdef int ridge
    cdef int n_faces
    cdef int n_ridges = len(face_at_ridge)
    cdef int n_cells = ridges_at_cell.shape[0]
    cdef int max_faces_per_cell = ridges_at_cell.shape[1]

    for i in range(n_ridges):
        face_at_ridge[i] = -1

    face = 0
    n_faces = 0
    for cell in range(n_cells):
        for i in range(max_faces_per_cell):
            ridge = ridges_at_cell[cell, i]
            if ridge >= 0:
                face = face_at_ridge[ridge]

                if face == -1:
                  face_at_ridge[ridge] = n_faces
                  n_faces += 1
                face = face_at_ridge[ridge]

                faces_at_cell[cell, i] = face
                face_at_ridge[ridge] = face
            else:
                faces_at_cell[cell, i] = -1


@cython.boundscheck(False)
def _get_corners_at_face(np.ndarray[DTYPE_t, ndim=1] face_at_ridge,
                         np.ndarray[DTYPE_t, ndim=2] vertices_at_ridge,
                         np.ndarray[DTYPE_t, ndim=1] corner_at_vertex,
                         np.ndarray[DTYPE_t, ndim=2] corners_at_face):
    """Get corners for each face.

    Parameters
    ----------
    face_at_ridge : ndarray of int, shape `(n_ridges, )`
        ID of the face corresponding to each ridge (-1 for ridges without
        a face). A ridge will not have a face if it is appears only in infinite
        regions.
    vertices_at_ridge : ndarray of int, shape `(n_ridges, 2)`
        Indices of the Voronoi vertices forming each Voronoi ridge (as provided
        by the `ridge_vertices` attribute of `scipy.spatial.Voronoi.class`.
    corner_at_vertex : ndarray of int, shape `(n_vertices, )`
        Output buffer. ID of corner corresponding to each voronoi vertex.
    corners_at_face : ndarray of int, shape `(n_vertices, 2)`
        Output buffer. IDs of two corners that define each face.
    """
    cdef int i
    cdef int ridge
    cdef int corner
    cdef int vertex
    cdef int face
    cdef int n_ridges = len(face_at_ridge)
    cdef int n_vertices = len(corner_at_vertex)
    cdef int n_faces = len(corners_at_face)
    cdef int n_corners = 0

    for vertex in range(n_vertices):
      corner_at_vertex[vertex] = -1

    for face in range(n_faces):
      for i in range(2):
        corners_at_face[face, i] = -1

    for ridge in range(n_ridges):
        for i in range(2):
            vertex = vertices_at_ridge[ridge, i]
            if vertex >= 0 and face_at_ridge[ridge] >= 0:
              if corner_at_vertex[vertex] == -1:
                  corner_at_vertex[vertex] = n_corners
                  n_corners += 1

              corners_at_face[face_at_ridge[ridge], i] = corner_at_vertex[vertex]


@cython.boundscheck(False)
def _get_xy_at_corners(np.ndarray[np.double_t, ndim=2] vertices,
                       np.ndarray[DTYPE_t, ndim=1] corner_at_vertex,
                       np.ndarray[np.double_t, ndim=2] xy_at_corner):
    """Get x and y coordinates for each corner.

    Parameters
    ----------
    vertices : ndarray of double, shape `(n_vertices, 2)`
        Coordinates of the Voronoi vertices (as provided by the
        `scipy.spacial.Voronoi` class.
    corner_at_vertex : ndarray of int, shape `(n_vertices, )`
        ID of the corner corresponding to each voronoi vertex.
    xy_at_corner : ndarray of double, shape `(n_vertices, 2)`
        Output buffer. Coordinates of the corners (first column is *x*,
        second column is *y*).
    """
    cdef int vertex
    cdef int corner
    cdef int n_vertices = len(vertices)

    for vertex in range(n_vertices):
        corner = corner_at_vertex[vertex]

        if corner >= 0:
          xy_at_corner[corner, 0] = vertices[vertex, 0]
          xy_at_corner[corner, 1] = vertices[vertex, 1]


@cython.boundscheck(False)
def _get_node_at_cell(np.ndarray[DTYPE_t, ndim=2] ridge_points,
                      np.ndarray[DTYPE_t, ndim=1] point_region,
                      np.ndarray[DTYPE_t, ndim=1] cell_at_region,
                      np.ndarray[DTYPE_t, ndim=1] node_at_cell):
    """Get node-to-cell connectivity.

    Parameters
    ----------
    ridge_points : ndarray of int, shape `(n_ridges, 2)`
        Indices of the points between which each Voronoi ridge lies (as
        provided by the `scipy.spatial.Voronoi` class).
    point_region : ndarray of int, shape `(n_points, )`
        Index of the Voronoi region for each input point (as provided by the
        `scipy.spatial.Voronoi` class).
    cell_at_region : ndarray of int, shape `(n_regions, )`
        ID of the cell corresponding to each voronoi region.
    node_at_cell : ndarray of int, shape `(n_cells, )`
        Output buffer. The node corresponding to each cell.
    """
    cdef int point
    cdef int region
    cdef int cell
    cdef int n_regions = len(point_region)

    for point in range(n_regions):
        region = point_region[point]
        cell = cell_at_region[region]
        if cell >= 0:
            node_at_cell[cell] = point
