import numpy as np


class NodeGrid(object):
    """__init__((coord0, coord1))
    Create a grid of nodes.

    Parameters
    ----------
    coord0, coord1 : sequence of array-like
        Coordinates of grid nodes

    Returns
    -------
    NodeGrid :
        A newly-created NodeGrid

    Examples
    --------
    >>> from landlab.grid.unstructured.nodes import NodeGrid
    >>> ngrid = NodeGrid(([0, 0, 1, 1], [0, 1, 0, 1]))
    >>> ngrid.ndim == 2
    True
    >>> ngrid.number_of_nodes == 4
    True
    >>> ngrid.x
    array([ 0.,  1.,  0.,  1.])
    >>> ngrid.y
    array([ 0.,  0.,  1.,  1.])

    Create a 1D grid.

    >>> ngrid = NodeGrid(((0, 1, 3), ))
    >>> ngrid.ndim == 1
    True
    >>> ngrid.number_of_nodes
    3
    >>> ngrid.x
    array([ 0.,  1.,  3.])
    >>> ngrid.y
    Traceback (most recent call last):
    AttributeError: Grid has no y-coordinate
    """

    def __init__(self, nodes):
        """__init__((coord0, coord1))
        Create a grid of nodes.

        Parameters
        ----------
        coord0, coord1 : sequence of array-like
            Coordinates of grid nodes

        Returns
        -------
        NodeGrid :
            A newly-created NodeGrid

        Examples
        --------
        >>> from landlab.grid.unstructured.nodes import NodeGrid
        >>> ngrid = NodeGrid(([0, 0, 1, 1], [0, 1, 0, 1]))
        >>> ngrid.ndim == 2
        True
        >>> ngrid.number_of_nodes == 4
        True
        >>> ngrid.x
        array([ 0.,  1.,  0.,  1.])
        >>> ngrid.y
        array([ 0.,  0.,  1.,  1.])

        Create a 1D grid.

        >>> ngrid = NodeGrid(((0, 1, 3), ))
        >>> ngrid.ndim == 1
        True
        >>> ngrid.number_of_nodes
        3
        >>> ngrid.x
        array([ 0.,  1.,  3.])
        >>> ngrid.y
        Traceback (most recent call last):
        AttributeError: Grid has no y-coordinate
        """
        self._coords = np.vstack(
            [np.array(coord, dtype=float) for coord in nodes]
        )
        self._coords.flags['WRITEABLE'] = False

        self._number_of_nodes = len(nodes[0])

    @property
    def ndim(self):
        return self._coords.shape[0]

    @property
    def number_of_nodes(self):
        """
        Examples
        --------
        >>> from landlab.grid.unstructured.nodes import NodeGrid
        >>> ngrid = NodeGrid(([0, 0, 1], [0, 1, 0]))
        >>> ngrid.number_of_nodes
        3
        """
        return self._number_of_nodes

    @property
    def x(self):
        """Node coordinates of the "fast" dimension.

        Examples
        --------
        >>> from landlab.grid.unstructured.nodes import NodeGrid
        >>> ngrid = NodeGrid(([0, 0, 1], [0, 1, 0]))
        >>> ngrid.x
        array([ 0.,  1.,  0.])
        """
        return self._coords[-1]

    @property
    def y(self):
        """Node y-coordinates.

        Examples
        --------
        >>> from landlab.grid.unstructured.nodes import NodeGrid
        >>> ngrid = NodeGrid(([0, 0, 1], [0, 1, 0]))
        >>> ngrid.y
        array([ 0.,  0.,  1.])
        """
        try:
            return self._coords[-2]
        except IndexError:
            raise AttributeError('Grid has no y-coordinate')

    @property
    def z(self):
        """Node z-coordinates.

        Examples
        --------
        >>> from landlab.grid.unstructured.nodes import NodeGrid
        >>> ngrid = NodeGrid(([0, 0, 1], [0, 1, 0]))
        >>> ngrid.y
        array([ 0.,  0.,  1.])
        """
        try:
            return self._coords[-2]
        except IndexError:
            raise AttributeError('Grid has no z-coordinate')

    @property
    def coord(self):
        """
        Examples
        --------
        >>> from landlab.grid.unstructured.nodes import NodeGrid
        >>> ngrid = NodeGrid(([0, 0, 1], [0, 1, 0]))
        >>> ngrid.coord[0]
        array([ 0.,  0.,  1.])
        >>> ngrid.coord[1]
        array([ 0.,  1.,  0.])
        """
        return self._coords

    @property
    def point(self):
        """
        Examples
        --------
        >>> from landlab.grid.unstructured.nodes import NodeGrid
        >>> ngrid = NodeGrid(([0, 0, 1], [0, 1, 0]))
        >>> ngrid.point
        array([[ 0.,  0.],
               [ 0.,  1.],
               [ 1.,  0.]])
        >>> ngrid.point[1]
        array([ 0.,  1.])
        """
        return self._coords.T
