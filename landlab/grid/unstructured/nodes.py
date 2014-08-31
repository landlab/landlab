import numpy as np


class NodeGrid(object):
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
        >>> ngrid = NodeGrid(([0, 0, 1, 1], [0, 1, 0, 1]))
        >>> ngrid.number_of_nodes
        4
        >>> ngrid.x
        array([ 0.,  1.,  0.,  1.])
        >>> ngrid.y
        array([ 0.,  0.,  1.,  1.])
        """
        if len(nodes) > 2:
            raise ValueError('only 2D grids supported')

        self._coords = np.vstack((np.array(nodes[0], dtype=float),
                                  np.array(nodes[1], dtype=float)))

        self._number_of_nodes = len(nodes[0])

    @property
    def number_of_nodes(self):
        return self._number_of_nodes

    @property
    def x(self):
        return self._coords[-1]

    @property
    def y(self):
        return self._coords[-2]


