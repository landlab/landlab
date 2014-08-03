#! /usr/bin/env python
"""
Examples
--------
>>> import numpy as np
>>> (y, x) = np.meshgrid(np.arange(4.), np.arange(5.), indexing='ij')
>>> grid = StructuredQuadModelGrid((y, x))
>>> grid.number_of_nodes
20
>>> grid.number_of_core_nodes
6
>>> grid.number_of_node_rows
4
>>> grid.number_of_node_columns
5
>>> grid.corner_nodes
array([ 0,  4, 15, 19])
>>> grid.number_of_cells
6
"""
import numpy as np

from ..base import FIXED_VALUE_BOUNDARY
from . import links, nodes, cells, faces


class StructuredQuadModelGrid(object):
    def __init__(self, node_coord, shape=None, node_status=None):
        """
        Parameters
        ----------
        node_coord : tuple
            Coordinates of all grid nodes.
        shape : tuple, optional
            Shape of the grid of nodes.
        """
        if len(node_coord) != 2:
            raise ValueError('only 2d grids are supported')

        if shape is None:
            shape = node_coord[0].shape

        self._shape = shape

        self._num_nodes = nodes.number_of_nodes(self.shape)
        self._num_cells = cells.number_of_cells(self.shape)
        self._num_links = links.number_of_links(self.shape)
        self._num_faces = faces.number_of_faces(self.shape)

        self._num_core_nodes = nodes.number_of_core_nodes(self.shape)
        self._num_core_cells = self._num_cells

        self._node_x, self._node_y = (
            np.ravel(node_coord[0]),
            np.ravel(node_coord[1]),
        )

        assert(node_status is None or node_status.size == self._num_nodes)

        if node_status is None:
            self._status = nodes.status_with_perimeter_as_boundary(
                self.shape, status_on_perimeter=FIXED_VALUE_BOUNDARY)
        else:
            self._status = node_status

        self._node_id_at_cells = cells.node_id_at_cells(self.shape)
        self._cell_id_at_nodes = cells.cell_ids(self.shape)

        self._cell_node = cells.node_id_at_cells(self.shape)

        self._in_link_id_at_nodes = links.node_in_link_ids(self.shape)
        self._out_link_id_at_nodes = links.node_out_link_ids(self.shape)

        self._node_id_at_link_start = links.node_id_at_link_start(self.shape)
        self._node_id_at_link_end = links.node_id_at_link_end(self.shape)

        self._active_link_ids = links.active_link_ids(self.shape, self._status)

    @property
    def shape(self):
        """Shape of the grid as rows, columns.
        """
        return self._shape

    @property
    def number_of_cells(self):
        """Number of cells.
        """
        return self._num_cells

    @property
    def number_of_nodes(self):
        """Number of nodes.
        """
        return self._num_nodes

    @property
    def number_of_core_nodes(self):
        """Number of core nodes.
        """
        return self._num_core_nodes

    @property
    def number_of_node_columns(self):
        """Number of node columns.

        Returns the number of columns, including boundaries.
        """
        return self.shape[1]

    @property
    def number_of_node_rows(self):
        """Number of node rows.

        Returns the number of rows, including boundaries.
        """
        return self.shape[0]

    @property
    def corner_nodes(self):
        """Nodes in grid corners.

        Return the IDs to the corner nodes of the grid, sorted by ID.

        Returns
        -------
        (4, ) ndarray
            Array of corner node IDs.

        Examples
        --------
        >>> import numpy as np
        >>> (x, y) = np.meshgrid(np.arange(4.), np.arange(5.), indexing='ij')
        >>> grid = StructuredQuadModelGrid((x, y))
        >>> grid.corner_nodes
        array([ 0,  4, 15, 19])
        """
        return nodes.corners(self.shape)
