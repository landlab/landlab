import numpy as np

from landlab.grid.base import BAD_INDEX_VALUE

from . import nodes


def number_of_cells(shape):
    """Number of cells is a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of cells in the grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.cells import number_of_cells
    >>> number_of_cells((3, 4))
    2
    """
    return np.prod(np.array(shape) - 2)


def shape_of_cells(shape):
    """Shape of cells in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Shape of cells in grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.cells import shape_of_cells
    >>> shape_of_cells((3, 4))
    (1, 2)
    """
    return tuple(np.array(shape) - 2)


def node_id_at_cells(shape):
    """Node ID at each cell.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        ID of node associated with each cell.

    Examples
    --------
    >>> from landlab.grid.structured_quad.cells import node_id_at_cells
    >>> node_id_at_cells((3, 4))
    array([[5, 6]])
    """
    node_ids = nodes.node_ids(shape)
    return node_ids[1:-1, 1:-1].copy().reshape(shape_of_cells(shape))


def cell_id_at_nodes(shape, bad=BAD_INDEX_VALUE):
    """Cell ID at each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        ID of cell associated with each node.

    Examples
    --------
    >>> from landlab.grid.structured_quad.cells import cell_id_at_nodes
    >>> cell_id_at_nodes((4, 5), bad=-1) # doctest: +NORMALIZE_WHITESPACE
    array([[-1, -1, -1, -1, -1],
           [-1,  0,  1,  2, -1],
           [-1,  3,  4,  5, -1],
           [-1, -1, -1, -1, -1]])
    """
    cells = np.empty(shape, dtype=int)
    cells[1:-1, 1:-1] = cell_ids(shape)
    cells[(0, -1), :] = bad
    cells[:, (0, -1)] = bad
    return cells


def cell_ids(shape):
    """IDs of all cells.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        ID of node associated with each cell.

    Examples
    --------
    >>> from landlab.grid.structured_quad.cells import cell_ids
    >>> cell_ids((3, 4))
    array([[0, 1]])
    """
    return np.arange(number_of_cells(shape)).reshape(shape_of_cells(shape))


class StructuredQuadCellGrid(object):
    def __init__(self, shape):
        self._shape = shape_of_cells(shape)
        self._number_of_cells = np.prod(self._shape)
        self._node_at_cell = node_id_at_cells(shape)

    @property
    def number_of_cells(self):
        return self._number_of_cells

    @property
    def shape(self):
        return self._shape

    @property
    def number_of_cell_rows(self):
        return self.shape[0]

    @property
    def number_of_cell_columns(self):
        return self.shape[1]

    @property
    def node_at_cell(self):
        return self._node_at_cell

    def number_of_vertices_at_cell(self, cell):
        return 4
