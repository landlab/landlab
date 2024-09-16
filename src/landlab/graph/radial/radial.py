from functools import cached_property

import numpy as np

from ...core.utils import as_id_array
from ...utils.decorators import read_only_array
from ..voronoi.voronoi import DelaunayGraph


class RadialGraphLayout:
    @staticmethod
    def number_of_nodes(shape):
        return np.sum(np.arange(1, shape[0] + 1)) * shape[1] + 1

    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_center=(0.0, 0.0)):
        """Create the node layout for a radial grid.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph.radial.radial import RadialGraphLayout
        >>> x, y = RadialGraphLayout.xy_of_node((1, 6))
        >>> x
        array([ 0. ,  1. ,  0.5, -0.5, -1. , -0.5,  0.5])
        >>> np.round(y / np.sin(np.pi / 3.0))
        array([ 0.,  0.,  1.,  1.,  0., -1., -1.])
        """
        n_rings, n_points = shape
        n_nodes = RadialGraphLayout.number_of_nodes(shape)

        x = np.empty((n_nodes,), dtype=float)
        y = np.empty((n_nodes,), dtype=float)

        x[0] = y[0] = 0.0
        offset = 1
        for ring in range(1, n_rings + 1):
            rho = spacing * ring
            d_theta = np.pi * 2 / (ring * shape[1])
            theta = np.arange(ring * shape[1]) * d_theta

            y[offset : offset + len(theta)] = rho * np.sin(theta)
            x[offset : offset + len(theta)] = rho * np.cos(theta)

            offset += len(theta)

        x = np.round(x, decimals=6)
        y = np.round(y, decimals=6)

        x += xy_of_center[0]
        y += xy_of_center[1]

        return (x, y)


class RadialGraphExtras:
    @property
    def shape(self):
        return self._shape

    @property
    def spacing(self):
        return self._spacing

    @property
    def origin(self):
        return self._origin

    @property
    def number_of_rings(self):
        return self.shape[0]

    @property
    def spacing_of_rings(self):
        return self.spacing

    @cached_property
    @read_only_array
    def radius_of_ring(self):
        return np.arange(0, self.number_of_rings, dtype=float) * self.spacing_of_rings

    @cached_property
    @read_only_array
    def angle_spacing_of_ring(self):
        return 2.0 * np.pi / self.nodes_per_ring

    @cached_property
    @read_only_array
    def nodes_per_ring(self):
        nodes_per_ring = np.empty(self.number_of_rings, dtype=int)
        nodes_per_ring[0] = 1
        nodes_per_ring[1:] = np.round(2.0 * np.pi * np.arange(1, self.number_of_rings))
        return nodes_per_ring

    @cached_property
    @read_only_array
    def ring_at_node(self):
        return np.repeat(np.arange(self.number_of_rings), self.nodes_per_ring)

    @cached_property
    @read_only_array
    def radius_at_node(self):
        return self.radius_of_ring[self.ring_at_node]

    @cached_property
    @read_only_array
    def angle_at_node(self):
        angle_at_node = np.empty(self.nodes_per_ring.sum(), dtype=float)
        angle_at_node[0] = 0.0
        offset = 1
        for n_nodes in self.nodes_per_ring[1:]:
            angles, step = np.linspace(
                0.0, 2 * np.pi, n_nodes, endpoint=False, retstep=True, dtype=float
            )
            angle_at_node[offset : offset + n_nodes] = np.add(
                angles, 0.5 * step, out=angles
            )
            offset += n_nodes
        return angle_at_node


class RadialGraph(RadialGraphExtras, DelaunayGraph):
    """Graph of a series of points on concentric circles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import RadialGraph
    >>> graph = RadialGraph((1, 4), sort=True)
    >>> graph.number_of_nodes
    5
    >>> graph.y_of_node
    array([-1.,  0.,  0.,  0.,  1.])
    >>> graph.x_of_node
    array([ 0., -1.,  0.,  1.,  0.])
    """

    def __init__(self, shape, spacing=1.0, xy_of_center=(0.0, 0.0), sort=False):
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
        except TypeError as exc:
            raise TypeError("spacing must be a float") from exc

        xy_of_center = tuple(np.broadcast_to(xy_of_center, 2))

        x_of_node, y_of_node = RadialGraphLayout.xy_of_node(
            shape, spacing=spacing, xy_of_center=xy_of_center
        )

        self._ring_spacing = spacing
        self._shape = tuple(shape)
        self._xy_of_center = xy_of_center

        DelaunayGraph.__init__(self, (y_of_node, x_of_node))

        if sort:
            self.sort()

    @property
    def xy_of_center(self):
        return self._xy_of_center

    @property
    def number_of_rings(self):
        """Number of node rings in grid.

        Returns
        -------
        int
            The number of node rings in the radial grid (not counting the
            center node).

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((1, 4))
        >>> graph.number_of_rings
        1

        :meta landlab: info-grid
        """
        return self._shape[0]

    @property
    def spacing_of_rings(self):
        """Fixed distance between rings.

        Returns
        -------
        ndarray of float
            The distance from the center node of each node.

        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((2, 6), spacing=2.0)
        >>> graph.spacing_of_rings
        2.0

        :meta landlab: info-grid, quantity
        """
        return self._ring_spacing

    @cached_property
    def radius_at_node(self):
        """Distance for center node to each node.

        Returns
        -------
        ndarray of float
            The distance from the center node of each node.

        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((2, 6), sort=True)
        >>> np.round(graph.radius_at_node, 3)
        array([2., 2., 2., 2., 2., 1., 1., 2., 1., 0., 1., 2., 1.,
               1., 2., 2., 2., 2., 2.])

        :meta landlab: info-node, quantity
        """
        return np.sqrt(
            np.square(self.x_of_node - self._xy_of_center[0])
            + np.square(self.y_of_node - self._xy_of_center[1])
        )

    @cached_property
    def number_of_nodes_in_ring(self):
        """Number of nodes in each ring.

        Returns
        -------
        ndarray of int
            Number of nodes in each ring, excluding the center node.

        >>> from landlab.graph import RadialGraph
        >>> graph = RadialGraph((4, 6))
        >>> graph.number_of_nodes_in_ring
        array([ 6, 12, 24, 48])

        :meta landlab: info-node, quantity
        """
        return as_id_array(self._shape[1] * 2 ** np.arange(self.number_of_rings))
