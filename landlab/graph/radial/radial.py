import numpy as np

from ..voronoi.voronoi import VoronoiGraph


def number_of_nodes(shape):
    return np.sum(np.arange(1, shape[0] + 1)) * shape[1] + 1


def create_xy_of_node(shape, spacing=1., xy_of_center=(0., 0.)):
    n_shells, n_points = shape
    n_nodes = number_of_nodes(shape)

    x = np.empty((n_nodes,), dtype=float)
    y = np.empty((n_nodes,), dtype=float)

    x[0] = y[0] = 0.
    offset = 1
    for shell in range(1, n_shells + 1):
        rho = spacing * shell
        d_theta = np.pi * 2 / (shell * shape[1])
        theta = np.arange(shell * shape[1]) * d_theta

        y[offset : offset + len(theta)] = rho * np.sin(theta)
        x[offset : offset + len(theta)] = rho * np.cos(theta)

        offset += len(theta)

    x = np.round(x, decimals=6)
    y = np.round(y, decimals=6)

    x += xy_of_center[0]
    y += xy_of_center[1]

    return (x, y)


class RadialGraph(VoronoiGraph):

    """Graph of a series of points on concentric circles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import RadialGraph
    >>> graph = RadialGraph((1, 4))
    >>> graph.number_of_nodes == 5
    True
    >>> graph.y_of_node
    array([-1.,  0.,  0.,  0.,  1.])
    >>> graph.x_of_node
    array([ 0., -1.,  0.,  1.,  0.])
    """

    def __init__(self, shape, spacing=1., xy_of_center=(0., 0.)):
        """Create a structured grid of triangles arranged radially.

        Parameters
        ----------
        shape : tuple of int
            Shape of the graph as number of rings and number of points
            in the first ring.
        spacing : float, optional
            Spacing between rings.
        xy_of_center : tuple of float, optional
            Coordinates of the node at the center of the grid.
        """
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError("spacing must be a float")

        xy_of_center = tuple(np.broadcast_to(xy_of_center, 2))

        x_of_node, y_of_node = create_xy_of_node(
            shape, spacing=spacing, xy_of_center=xy_of_center
        )

        super(RadialGraph, self).__init__(
            (y_of_node, x_of_node), xy_sort=True, rot_sort=True
        )

        self._shell_spacing = spacing
        self._shape = tuple(shape)
        self._xy_of_center = xy_of_center

    @property
    def xy_of_center(self):
        return self._xy_of_center

    @property
    def number_of_shells(self):
        """Number of node shells in grid.

        Returns
        -------
        int
            The number of node shells in the radial grid (not counting the
            center node).

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((1, 4))
        >>> graph.number_of_shells
        1

        LLCATS: GINF
        """
        return self._shape[0]

    @property
    def spacing_of_shells(self):
        """Fixed distance between shells.

        Returns
        -------
        ndarray of float
            The distance from the center node of each node.

        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((2, 6), spacing=2.)
        >>> graph.spacing_of_shells
        2.0

        LLCATS: GINF MEAS
        """
        return self._shell_spacing

    @property
    def radius_at_node(self):
        """Distance for center node to each node.

        Returns
        -------
        ndarray of float
            The distance from the center node of each node.

        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((2, 6))
        >>> np.round(graph.radius_at_node, 3)
        array([ 2.,  2.,  2.,  2.,  2.,  1.,  1.,  2.,  1.,  0.,  1.,  2.,  1.,
                1.,  2.,  2.,  2.,  2.,  2.])

        LLCATS: NINF MEAS
        """
        return np.sqrt(
            np.square(self.x_of_node - self._xy_of_center[0])
            + np.square(self.y_of_node - self._xy_of_center[1])
        )

    @property
    def number_of_nodes_in_shell(self):
        """Number of nodes in each shell.

        Returns
        -------
        ndarray of int
            Number of nodes in each shell, excluding the center node.

        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((4, 6))
        >>> graph.number_of_nodes_in_shell
        array([ 6, 12, 18, 24])

        LLCATS: NINF MEAS
        """
        return np.arange(1, self._shape[0] + 1) * self._shape[1]
