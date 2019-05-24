import numpy as np

from ..voronoi.voronoi import VoronoiGraph


def number_of_nodes(shape):
    return np.sum(np.arange(1, shape[0] + 1)) * shape[1] + 1


def create_xy_of_node(shape, spacing=1.0, origin=(0.0, 0.0)):
    n_shells, n_points = shape
    n_nodes = number_of_nodes(shape)

    x = np.empty((n_nodes,), dtype=float)
    y = np.empty((n_nodes,), dtype=float)

    x[0] = y[0] = 0.0
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

    y += origin[0]
    x += origin[1]

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

    def __init__(self, shape, spacing=1.0, origin=(0.0, 0.0)):
        """Create a structured grid of triangles arranged radially.

        Parameters
        ----------
        shape : tuple of int
            Shape of the graph as number of rings and number of points
            in the first ring.
        spacing : float, optional
            Spacing between rings.
        origin : tuple of float, optional
            Coordinates of the center of the grid.
        """
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError("spacing must be a float")

        x_of_node, y_of_node = create_xy_of_node(shape, spacing=spacing, origin=origin)

        super(RadialGraph, self).__init__(
            (y_of_node, x_of_node), xy_sort=True, rot_sort=True
        )
