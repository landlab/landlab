#! /usr/bin/env python
"""
Examples
--------
>>> import numpy as np
>>> from landlab.grid.structured_quad.structured import StructuredQuadGrid
>>> (y, x) = np.meshgrid(np.arange(4.), np.arange(5.), indexing='ij')
>>> grid = StructuredQuadGrid((y, x))
>>> grid.number_of_nodes
20

>>> grid.number_of_node_rows == 4
True
>>> grid.number_of_node_columns == 5
True
>>> grid.nodes_at_corners_of_grid
array([ 0,  4, 15, 19])
>>> grid.number_of_cells
6
"""
import numpy as np

from landlab.utils.decorators import deprecated

from ..base import FIXED_VALUE_BOUNDARY
from ..unstructured.base import BaseGrid
from . import cells as quad_cells, faces as quad_faces, links as quad_links, nodes
from .cells import StructuredQuadCellGrid
from .links import StructuredQuadLinkGrid


class StructuredQuadGrid(BaseGrid):
    """
    Parameters
    ----------
    node_coord : tuple
        Coordinates of all grid nodes.
    shape : tuple, optional
        Shape of the grid of nodes.
    """

    def __init__(
        self,
        node_coord,
        shape=None,
        axis_name=None,
        axis_units=None,
        links=True,
        cells=True,
        node_status=None,
    ):
        """
        Parameters
        ----------
        node_coord : tuple
            Coordinates of all grid nodes.
        shape : tuple, optional
            Shape of the grid of nodes.
        """
        if len(node_coord) != 2:
            raise ValueError("only 2d grids are supported")

        self._shape = shape or node_coord[0].shape

        if node_status is not None:
            if node_status.size != nodes.number_of_nodes(self.shape):
                raise ValueError("incorrect size for node_status array")

        if node_status is None:
            self._status = nodes.status_with_perimeter_as_boundary(
                self.shape, node_status=FIXED_VALUE_BOUNDARY
            )
        else:
            self._status = node_status

        if links:
            # links = (node_id_at_link_start(self.shape),
            #         node_id_at_link_end(self.shape))
            link_grid = StructuredQuadLinkGrid(self.shape)
        if cells:
            cell_grid = StructuredQuadCellGrid(self.shape)

        # super(StructuredQuadGrid, self).__init__(node_status=node_status)
        BaseGrid.__init__(
            self,
            (node_coord[0].flatten(), node_coord[1].flatten()),
            links=link_grid,
            cells=cell_grid,
        )

        self._num_nodes = nodes.number_of_nodes(self.shape)
        self._num_cells = quad_cells.number_of_cells(self.shape)
        self._num_links = quad_links.number_of_links(self.shape)
        self._num_faces = quad_faces.number_of_faces(self.shape)

        self._num_core_nodes = nodes.number_of_core_nodes(self.shape)
        self._num_core_cells = self._num_cells

        self._node_x, self._node_y = (np.ravel(node_coord[0]), np.ravel(node_coord[1]))

        self._node_id_at_cells = quad_cells.node_id_at_cells(self.shape)
        self._cell_id_at_nodes = quad_cells.cell_ids(self.shape)

        self._cell_node = quad_cells.node_id_at_cells(self.shape)

        self._in_link_id_at_nodes = quad_links.node_in_link_ids(self.shape)
        self._out_link_id_at_nodes = quad_links.node_out_link_ids(self.shape)

        self._node_id_at_link_start = quad_links.node_id_at_link_start(self.shape)
        self._node_id_at_link_end = quad_links.node_id_at_link_end(self.shape)

        self._active_link_ids = quad_links.active_link_ids(self.shape, self._status)

    @property
    def shape(self):
        """Shape of the grid as rows, columns.
        """
        return self._shape

    # @property
    # def number_of_core_nodes(self):
    #    """Number of core nodes.
    #    """
    #    return self._num_core_nodes

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
    @deprecated(use="nodes_at_corners_of_grid", version=1.0)
    def corner_nodes(self):
        return self.nodes_at_corners_of_grid

    @property
    def nodes_at_corners_of_grid(self):
        """Nodes in grid corners.

        Return the IDs to the corner nodes of the grid, sorted by ID.

        Returns
        -------
        (4, ) ndarray
            Array of corner node IDs.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.grid.structured_quad.structured import StructuredQuadGrid
        >>> (x, y) = np.meshgrid(np.arange(4.), np.arange(5.), indexing='ij')
        >>> grid = StructuredQuadGrid((x, y))
        >>> grid.nodes_at_corners_of_grid
        array([ 0,  4, 15, 19])
        """
        return nodes.corners(self.shape)
