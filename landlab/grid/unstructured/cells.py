import numpy as np
from six.moves import range

from ...utils.jaggedarray import JaggedArray


class CellGrid(object):
    """
    Parameters
    ----------
    vertices : array-like
        Vertex IDs at each grid cell
    vertices_per_cell : array-like
        Number of vertices per grid cell

    Returns
    -------
    CellGrid :
        A newly-created CellGrid

    Examples
    --------
    Create a grid of two cells where the first cell has four vertices and
    the second has three.

    >>> from landlab.grid.unstructured.cells import CellGrid
    >>> cgrid = CellGrid([0, 1, 3, 2, 1, 4, 3], [4, 3])
    >>> cgrid.number_of_cells
    2
    >>> cgrid.number_of_vertices_at_cell(0)
    4
    >>> cgrid.number_of_vertices_at_cell(1)
    3
    >>> cgrid.vertices_at_cell(0)
    array([0, 1, 3, 2])
    >>> cgrid.vertices_at_cell(1)
    array([1, 4, 3])

    Associate nodes with each cell.

    >>> cgrid = CellGrid([0, 1, 2, 3, 1, 3, 4], [4, 3], node_at_cell=[10, 11])
    >>> cgrid.node_at_cell
    array([10, 11])
    >>> cgrid.cell_at_node[10]
    0
    >>> cgrid.cell_at_node[11]
    1
    """

    def __init__(self, vertices, vertices_per_cell, node_at_cell=None):
        """
        Parameters
        ----------
        vertices : array-like
            Vertex IDs at each grid cell
        vertices_per_cell : array-like
            Number of vertices per grid cell

        Returns
        -------
        CellGrid :
            A newly-created CellGrid

        Examples
        --------
        Create a grid of two cells where the first cell has four vertices and
        the second has three.

        >>> from landlab.grid.unstructured.cells import CellGrid
        >>> cgrid = CellGrid([0, 1, 3, 2, 1, 4, 3], [4, 3])
        >>> cgrid.number_of_cells
        2
        >>> cgrid.number_of_vertices_at_cell(0)
        4
        >>> cgrid.number_of_vertices_at_cell(1)
        3
        >>> cgrid.vertices_at_cell(0)
        array([0, 1, 3, 2])
        >>> cgrid.vertices_at_cell(1)
        array([1, 4, 3])

        Associate nodes with each cell.

        >>> cgrid = CellGrid([0, 1, 2, 3, 1, 3, 4], [4, 3], node_at_cell=[10, 11])
        >>> cgrid.node_at_cell
        array([10, 11])
        >>> cgrid.cell_at_node[10]
        0
        >>> cgrid.cell_at_node[11]
        1
        """
        self._vertices_at_cell = JaggedArray(vertices, vertices_per_cell)
        self._number_of_cells = len(vertices_per_cell)

        if node_at_cell:
            self._node_at_cell = np.array(node_at_cell)
            self._cell_at_node = np.ma.masked_all(
                max(node_at_cell) + 1, dtype=int)
            self._cell_at_node[self._node_at_cell] = range(len(node_at_cell))
            #self._cell_id_map = dict(zip(node_at_cell, range(len(node_at_cell))))

    @property
    def number_of_cells(self):
        return self._number_of_cells

    @property
    def node_at_cell(self):
        return self._node_at_cell

    @property
    def cell_at_node(self):
        return self._cell_at_node

    def number_of_vertices_at_cell(self, cell):
        return self._vertices_at_cell.length_of_row(cell)

    def vertices_at_cell(self, cell):
        return self._vertices_at_cell.row(cell)

    def iter(self):
        """Iterate over the cells of the grid.

        Returns
        -------
        ndarray :
            Nodes entering and leaving each node
        """
        for cell in range(self.number_of_cells):
            yield self.vertices_at_cell(cell)

    @property
    def cell_id(self):
        try:
            return self._cell_ids
        except AttributeError:
            return np.arange(self.number_of_cells)

    def vertices_at_cell_id(self, cell_id):
        try:
            return self.vertices_at_cell(self._cell_id_map[cell_id])
        except AttributeError:
            return self.vertices_at_cell(cell_id)
