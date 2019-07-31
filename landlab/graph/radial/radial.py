import numpy as np

from ...core.utils import argsort_points_by_x_then_y
from ..voronoi.voronoi import DelaunayGraph
from ...utils.decorators import store_result_in_grid, read_only_array


def number_of_nodes(shape):
    return np.sum(np.arange(1, shape[0] + 1)) * shape[1] + 1


def create_xy_of_node(shape, spacing=1.0, xy_of_center=(0.0, 0.0)):
    n_shells, n_points = shape
    n_nodes = number_of_nodes(shape)

    x = np.empty((n_nodes,), dtype=float)
    y = np.empty((n_nodes,), dtype=float)
    # nodes_per_shell = np.round(2. * np.pi * np.arange(1, n_shells + 1)).astype(int)
    # n_nodes = np.sum(nodes_per_shell) + 1

    x[0] = y[0] = 0.0
    offset = 1
    # for shell in range(0, n_shells):
    for shell in range(1, n_shells + 1):
        # rho = spacing * (shell + 1)
        rho = spacing * shell
        d_theta = np.pi * 2 / (shell * shape[1])
        theta = np.arange(shell * shape[1]) * d_theta

        y[offset : offset + len(theta)] = rho * np.sin(theta)
        x[offset : offset + len(theta)] = rho * np.cos(theta)
        # d_theta = 2. * np.pi / nodes_per_shell[shell]
        # theta = np.arange(nodes_per_shell[shell]) * d_theta

        offset += len(theta)

    x = np.round(x, decimals=6)
    y = np.round(y, decimals=6)

    x += xy_of_center[0]
    y += xy_of_center[1]

    return (x, y)


class RadialNodeLayout(object):

    def __init__(self, n_shells, spacing=1., origin=0.):
        self._n_shells = n_shells

        self._n_shells = int(n_shells)
        self._spacing_of_shells = float(spacing)
        self._origin = tuple(np.broadcast_to(origin, (2, )).astype(float))

        # self._shell_at_node = np.repeat(np.arange(self.number_of_shells),
        #                                 self.nodes_per_shell)

        # y_of_node = graph.radius_at_node * np.sin(graph.angle_at_node) - origin[0]
        # x_of_node = graph.radius_at_node * np.cos(graph.angle_at_node) - origin[1]

        # sorted_nodes = argsort_points_by_x_then_y((x_of_node, y_of_node))

    @property
    def origin(self):
        return self._origin

    @property
    def y_of_node(self):
        return self.radius_at_node * np.sin(self.angle_at_node) - self.origin[0]

    @property
    def x_of_node(self):
        return self.radius_at_node * np.cos(self.angle_at_node) - self.origin[1]

    @property
    def xy_of_node(self):
        return self.x_of_node, self.y_of_node

    @property
    def number_of_shells(self):
        return self._n_shells

    @property
    def spacing_of_shells(self):
        return self._spacing_of_shells

    @property
    @read_only_array
    def radius_of_shell(self):
        return np.arange(0, self.number_of_shells, dtype=float) * self.spacing_of_shells

    @property
    @read_only_array
    def shell_at_node(self):
        return np.repeat(np.arange(self.number_of_shells), self.nodes_per_shell)

    @property
    @read_only_array
    def radius_at_node(self):
        return self.radius_of_shell[self.shell_at_node]

    @property
    @read_only_array
    def angle_at_node(self):
        angle_at_node = np.empty(self.nodes_per_shell.sum(), dtype=float)
        angle_at_node[0] = 0.
        offset = 1
        for n_nodes in self.nodes_per_shell[1:]:
            angles, step = np.linspace(0., 2 * np.pi, n_nodes, endpoint=False,
                                       retstep=True, dtype=float)
            angle_at_node[offset:offset + n_nodes] = np.add(angles, .5 * step, out=angles)
            offset += n_nodes
        return angle_at_node

    @property
    @read_only_array
    def nodes_per_shell(self):
        # nodes_per_shell = np.arange(self.number_of_shells, dtype=int) * self.shape[1]
        # nodes_per_shell[0] = 1
        # return nodes_per_shell
        nodes_per_shell = np.empty(self.number_of_shells, dtype=int)
        nodes_per_shell[0] = 1
        nodes_per_shell[1:] = np.round(2. * np.pi * np.arange(1, self.number_of_shells))
        return nodes_per_shell


class RadialGraphExtras(object):

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
    def number_of_shells(self):
        return self.shape[0]

    @property
    def spacing_of_shells(self):
        return self.spacing

    @property
    # @store_result_in_grid()
    @read_only_array
    def radius_of_shell(self):
        return np.arange(0, self.number_of_shells, dtype=float) * self.spacing_of_shells

    @property
    # @store_result_in_grid()
    @read_only_array
    def angle_spacing_of_shell(self):
        return 2. * np.pi / self.nodes_per_shell

    @property
    # @store_result_in_grid()
    @read_only_array
    def nodes_per_shell(self):
        # nodes_per_shell = np.arange(self.number_of_shells, dtype=int) * self.shape[1]
        # nodes_per_shell[0] = 1
        # return nodes_per_shell
        nodes_per_shell = np.empty(self.number_of_shells, dtype=int)
        nodes_per_shell[0] = 1
        nodes_per_shell[1:] = np.round(2. * np.pi * np.arange(1, self.number_of_shells))
        return nodes_per_shell

    @property
    # @store_result_in_grid()
    @read_only_array
    def shell_at_node(self):
        return np.repeat(np.arange(self.number_of_shells), self.nodes_per_shell)

    @property
    # @store_result_in_grid()
    @read_only_array
    def radius_at_node(self):
        return self.radius_of_shell[self.shell_at_node]

    @property
    # @store_result_in_grid()
    @read_only_array
    def angle_at_node(self):
        angle_at_node = np.empty(self.nodes_per_shell.sum(), dtype=float)
        angle_at_node[0] = 0.
        offset = 1
        for n_nodes in self.nodes_per_shell[1:]:
            angles, step = np.linspace(0., 2 * np.pi, n_nodes, endpoint=False,
                                       retstep=True, dtype=float)
            angle_at_node[offset:offset + n_nodes] = np.add(angles, .5 * step, out=angles)
            offset += n_nodes
        return angle_at_node

    def empty_cache(self):
        for attr in ('_angle_at_node', '_radius_at_node', '_shell_at_node',
                     '_nodes_per_shell', '_angle_spacing_of_shell',
                     '_radius_of_shell', ):
            try:
                del self.__dict__[attr]
            except KeyError:
                pass


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
        # x_of_node, y_of_node = create_xy_of_node(shape, spacing=spacing,
        #                                          origin=origin)
        # graph = RadialGraphMaker(shape, spacing=spacing, origin=origin)

        # self._shape = graph._shape
        # self._spacing = graph.spacing
        # self._origin = graph.origin
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError("spacing must be a float")

        xy_of_center = tuple(np.broadcast_to(xy_of_center, 2))

        x_of_node, y_of_node = create_xy_of_node(
            shape, spacing=spacing, xy_of_center=xy_of_center
        )

        # super(RadialGraph, self).__init__(
        #     (y_of_node, x_of_node), xy_sort=True, rot_sort=True
        # )

        self._shell_spacing = spacing
        self._shape = tuple(shape)
        self._xy_of_center = xy_of_center

        # nodes = RadialNodeLayout(self.shape[0], spacing=spacing, origin=xy_of_center)

        # VoronoiGraph.__init__(
        #     self, (y_of_node, x_of_node), xy_sort=True, rot_sort=True
        # )

        DelaunayGraph.__init__(self, (y_of_node, x_of_node))

        if sort:
            self.sort()

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
        >>> graph = RadialGraph((2, 6), sort=True)
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
