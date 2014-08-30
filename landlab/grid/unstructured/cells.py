from ...utils.jaggedarray import JaggedArray


class CellGrid(object):
    def __init__(self, nodes_at_cell, nodes_per_cell):
        """
        Parameters
        ----------
        nodes_at_cell : array-like
            Node IDs at each grid cell
        nodes_per_cell : array-like
            Number of nodes per grid cell

        Returns
        -------
        CellGrid :
            A newly-created CellGrid

        Examples
        --------
        >>> cgrid = CellGrid([0, 1, 2 ,3, 1, 3, 4], [4, 3])
        >>> cgrid.number_of_cells
        2
        >>> cgrid.number_of_nodes_at_cell(0)
        4
        >>> cgrid.number_of_nodes_at_cell(1)
        3
        >>> cgrid.nodes_at_cell(0)
        array([0, 1, 2, 3])
        >>> cgrid.nodes_at_cell(1)
        array([1, 3, 4])
        """
        self._nodes_at_cell = JaggedArray(nodes_at_cell, nodes_per_cell)
        self._number_of_cells = len(nodes_per_cell)

    @property
    def number_of_cells(self):
        return self._number_of_cells

    def number_of_nodes_at_cell(self, cell):
        return self._nodes_at_cell.length_of_row(cell)

    def nodes_at_cell(self, cell):
        return self._nodes_at_cell.row(cell)

    def iter(self):
        """Iterate over the cells of the grid.

        Returns
        -------
        ndarray :
            Nodes entering and leaving each node
        """
        for cell in xrange(self.number_of_cells):
            yield self.nodes_at_cell(cell)
