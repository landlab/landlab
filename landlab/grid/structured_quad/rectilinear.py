#! /usr/bin/env python

import numpy as np

from .structured import StructuredQuadModelGrid


class RectilinearModelGrid(StructuredQuadModelGrid):
    """
    Parameters
    ----------
    coord : tuple
        Coordinates of node rows and node columns.

    Examples
    --------
    >>> import numpy as np
    >>> (y, x) = np.arange(4.), np.arange(5.)
    >>> grid = RectilinearModelGrid((y, x))
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
    >>> grid.node_row_coord
    array([ 0.,  1.,  2.,  3.])
    >>> grid.node_col_coord
    array([ 0.,  1.,  2.,  3.,  4.])
    """
    def __init__(self, coord):
        """
        Parameters
        ----------
        coord : tuple
            Coordinates of node rows and node columns.
        """
        if len(coord) != 2:
            raise ValueError('only 2d grids are supported')

        shape = (
            len(coord[0]),
            len(coord[1]),
        )
        node_coord = np.meshgrid(*coord, indexing='ij')

        super(RectilinearModelGrid, self).__init__(node_coord, shape)

        self._coord = (coord[0], coord[1])

    @property
    def coord(self):
        """Node row and column coordinates.
        """
        return self._coord

    @property
    def node_row_coord(self):
        return self._coord[0]

    @property
    def node_col_coord(self):
        return self._coord[1]


class UniformRectilinearModelGrid(RectilinearModelGrid):
    """
    Examples
    --------
    >>> grid = UniformRectilinearModelGrid((4, 5), spacing=(2, 3), origin=(-1, 1))
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
    >>> grid.node_row_coord
    array([-1.,  1.,  3.,  5.])
    >>> grid.node_col_coord
    array([  1.,   4.,   7.,  10.,  13.])
    """
    def __init__(self, shape, spacing=(1., 1.), origin=(0., 0.)):
        """
        Parameters
        ----------
        shape : tuple
            Shape of the grid in nodes.
        spacing : tuple, optional
            Spacing between rows and columns.
        origin : tuple, optional
            Coordinates of grid origin.
        """
        if len(shape) != 2:
            raise ValueError('only 2d grids are supported')

        coords = (
            np.arange(origin[0], origin[0] + shape[0] * spacing[0], spacing[0], dtype=np.float64),
            np.arange(origin[1], origin[1] + shape[1] * spacing[1], spacing[1], dtype=np.float64),
        )

        super(UniformRectilinearModelGrid, self).__init__(coords)

        self._spacing = tuple(spacing)

    @property
    def spacing(self):
        """Spacing of rows and columns of grid nodes.
        """
        return self._spacing

    @property
    def dy(self):
        """Spacing between rows of grid nodes.
        """
        return self._spacing[0]

    @property
    def dx(self):
        """Spacing between columns of grid nodes.
        """
        return self._spacing[1]
