#! /usr/bin/env python
"""
Examples
--------
>>> import numpy as np
>>> (y, x) = np.meshgrid(np.arange(4.), np.arange(5.), indexing='ij')
>>> grid = StructuredQuadModelGrid((y, x))
>>> grid.number_of_nodes
20
>>> grid.number_of_interior_nodes
6
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
from numpy import ravel

from ..base import ModelGrid, FIXED_VALUE_BOUNDARY, CLOSED_BOUNDARY
from ...utils import structured_grid as sgrid
from . import links
from . import nodes


class StructuredQuadModelGrid(ModelGrid):
    def __init__(self, node_coord, shape=None):
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
        self._num_active_nodes = self.number_of_nodes

        self._num_cells = sgrid.cell_count(self.shape)
        self._num_active_cells = self.number_of_cells

        self._num_core_nodes = self.number_of_cells
        self._num_core_cells = self.number_of_cells

        self._num_links = links.number_of_links(self.shape)
        self._num_active_links = sgrid.active_link_count(self.shape)

        self._num_faces = sgrid.face_count(self.shape)
        self._num_active_faces = sgrid.active_face_count(self.shape)

        self._node_x, self._node_y = (
            ravel(node_coord[0]),
            ravel(node_coord[1]),
        )

        # Node boundary/active status:
        # Next, we set up an array of "node status" values, which indicate
        # whether a given node is an active, non-boundary node, or some type of
        # boundary. Here we default to having all perimeter nodes be active
        # fixed-value boundaries.
        self.node_status = sgrid.node_status(
            self.shape, boundary_status=FIXED_VALUE_BOUNDARY)

        self.cell_node = sgrid.node_index_at_cells(self.shape)
        self.node_corecell = sgrid.core_cell_index_at_nodes(self.shape)
        self.core_cells = sgrid.core_cell_index(self.shape)
        self.corecell_node = self.cell_node

        (self.link_fromnode,
         self.link_tonode) = sgrid.node_index_at_link_ends(self.shape)

        # set up in-link and out-link matrices and numbers
        self._setup_inlink_and_outlink_matrices()
        #(self.node_inlink_matrix,
        # self.node_numinlink) = sgrid.setup_inlink_matrix(shape)

        #(self.node_outlink_matrix,
        # self.node_numoutlink) = sgrid.setup_outlink_matrix(shape)

        # Flag indicating whether we have created diagonal links.
        self._diagonal_links_created = False

        # Set up the list of active links
        self._reset_list_of_active_links()

        #   set up link faces
        #
        #   Here we assume that we've already created a list of active links
        # in which all 4 boundaries are "open", such that each boundary node
        # (except the 4 corners) is connected to an adjacent interior node. In
        # this case, there will be the same number of faces as active links,
        # and the numbering of faces will be the same as the corresponding
        # active links. We start off creating a list of all None values. Only
        # those links that cross a face will have this None value replaced with
        # a face ID.
        self.link_face = sgrid.face_index_at_links(self.shape,
                                                   actives=self.active_link_ids)

        # List of neighbors for each cell: we will start off with no
        # list. If a caller requests it via get_neighbor_list or
        # create_neighbor_list, we'll create it if necessary.
        self.neighbor_list_created = False

        # List of diagonal neighbors. As with the neighbor list, we'll only
        # create it if requested.
        self.diagonal_list_created = False

        super(StructuredModelGrid, self).__init__()

    def _setup_inlink_and_outlink_matrices(self):
        (self.node_inlink_matrix,
         self.node_numinlink) = sgrid.setup_inlink_matrix(self.shape)

        (self.node_outlink_matrix,
         self.node_numoutlink) = sgrid.setup_outlink_matrix(self.shape)

    def _setup_active_inlink_and_outlink_matrices(self):
        """
        Creates data structures to record the numbers of active inlinks and
        active outlinks for each node. These data structures are equivalent to
        the "regular" inlink and outlink matrices, except that it uses the IDs
        of active links (only).
        """

        node_status = self.node_status != CLOSED_BOUNDARY

        (self.node_active_inlink_matrix, self.node_numactiveinlink) = (
            sgrid.setup_active_inlink_matrix(self.shape, node_status=node_status))

        (self.node_active_outlink_matrix, self.node_numactiveoutlink) = (
            sgrid.setup_active_outlink_matrix(self.shape, node_status=node_status))

    def _reset_list_of_active_links(self):
        """
        Assuming the active link list has already been created elsewhere, this
        helper method checks link statuses (active/inactive) for internal
        consistency after the BC status of some nodes has been changed.
        """
        super(StructuredModelGrid, self)._reset_list_of_active_links()
        #if self._diagonal_links_created:
        #    self._reset_list_of_active_diagonal_links()

    @property
    def shape(self):
        """Shape of the grid as rows, columns.
        """
        return self._shape

    @property
    def number_of_interior_nodes(self):
        """Number of interior nodes.

        Returns the number of interior nodes on the grid, i.e., non-perimeter
        nodes. Compare self.number_of_core_nodes.
        """
        return sgrid.interior_node_count(self.shape)

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
        >>> grid = StructuredModelGrid((x, y))
        >>> grid.corner_nodes
        array([ 0,  4, 15, 19])
        """
        return sgrid.corners(self.shape)
