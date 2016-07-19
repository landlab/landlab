#! /usr/bin/env python

import numpy as np

from .structured import StructuredQuadGrid


class RectilinearGrid(StructuredQuadGrid):
    """
    Parameters
    ----------
    coord : tuple
        Coordinates of node rows and node columns.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.structured_quad.rectilinear import RectilinearGrid
    >>> (y, x) = np.arange(4.), np.arange(5.)
    >>> grid = RectilinearGrid((y, x))
    >>> grid.number_of_nodes
    20
    >>> grid.number_of_node_rows
    4
    >>> grid.number_of_node_columns
    5
    >>> grid.nodes_at_corners_of_grid
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

        super(RectilinearGrid, self).__init__(node_coord, shape, cells=True)

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


class UniformRectilinearGrid(RectilinearGrid):
    """
    Parameters
    ----------
    shape : tuple
        Shape of the grid in nodes.
    spacing : tuple, optional
        Spacing between rows and columns.
    origin : tuple, optional
        Coordinates of grid origin.

    Examples
    --------
    >>> from landlab.grid.structured_quad.rectilinear import UniformRectilinearGrid
    >>> grid = UniformRectilinearGrid((4, 5), spacing=(2, 3), origin=(-1, 1))
    >>> grid.number_of_nodes
    20

    #>>> grid.number_of_core_nodes
    6
    >>> grid.number_of_node_rows
    4
    >>> grid.number_of_node_columns
    5
    >>> grid.nodes_at_corners_of_grid
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
            np.arange(origin[0], origin[0] + shape[0] * spacing[0], spacing[0],
                      dtype=np.float64),
            np.arange(origin[1], origin[1] + shape[1] * spacing[1], spacing[1],
                      dtype=np.float64),
        )

        super(UniformRectilinearGrid, self).__init__(coords)

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


class RasterGrid(UniformRectilinearGrid):
    """
    Parameters
    ----------
    shape : tuple
        Shape of the grid in nodes.
    spacing : float, optional
        Spacing between rows and columns.
    origin : tuple, optional
        Coordinates of grid origin.
        
    Examples
    --------
    >>> from landlab.grid.structured_quad.rectilinear import RasterGrid
    >>> grid = RasterGrid((4, 5), spacing=2, origin=(-1, 1))
    >>> grid.number_of_nodes
    20

    #>>> grid.number_of_core_nodes
    6
    >>> grid.number_of_node_rows
    4
    >>> grid.number_of_node_columns
    5
    >>> grid.nodes_at_corners_of_grid
    array([ 0,  4, 15, 19])
    >>> grid.number_of_cells
    6
    >>> grid.node_row_coord
    array([-1.,  1.,  3.,  5.])
    >>> grid.node_col_coord
    array([ 1.,  3.,  5.,  7.,  9.])
    """

    def __init__(self, shape, spacing=1., origin=0.):
        """
        Parameters
        ----------
        shape : tuple
            Shape of the grid in nodes.
        spacing : float, optional
            Spacing between rows and columns.
        origin : tuple, optional
            Coordinates of grid origin.
        """
        super(RasterGrid, self).__init__(shape, origin=origin,
                                         spacing=(spacing, ) * len(shape))
