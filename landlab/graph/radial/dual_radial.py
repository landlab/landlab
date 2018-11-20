import numpy as np

from ...core.utils import argsort_points_by_x_then_y
from ..voronoi import DualVoronoiGraph
from .radial import create_xy_of_node, RadialGraphExtras, RadialNodeLayout, RadialGraph


# class DualRadialGraph(RadialGraphExtras, DualVoronoiGraph):
class DualRadialGraph(RadialGraph, DualVoronoiGraph):

    """Graph of a series of points on concentric circles.

    Examples
    --------
    >>> from landlab.graph import DualRadialGraph
    >>> graph = DualRadialGraph((1, 4))
    >>> graph.number_of_corners
    4
    >>> graph.y_of_corner
    array([-0.5, -0.5,  0.5,  0.5])
    >>> graph.x_of_corner
    array([-0.5,  0.5, -0.5,  0.5])
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
            Coordinates of the center of the grid.
        """
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError("spacing must be a float")

        # x_of_node, y_of_node = create_xy_of_node(shape, spacing=spacing, xy_of_center=xy_of_center)

        # super(DualRadialGraph, self).__init__(
        #     (y_of_node, x_of_node), xy_sort=True, rot_sort=True
        # )

        # self._shape = tuple(np.asarray(shape).astype(int))
        # self._origin = tuple(np.broadcast_to(origin, (2, )).astype(float))
        # self._spacing = float(spacing)
        self._shape = int(shape[0]), int(shape[1])
        self._spacing = float(spacing)
        self._xy_of_center = tuple(np.broadcast_to(xy_of_center, (2, )).astype(float))

        # graph = RadialGraphExtras()
        # graph._shape = self._shape
        # graph._spacing = self._spacing
        # graph._origin = self._origin

        # y_of_node = graph.radius_at_node * np.sin(graph.angle_at_node) - graph.origin[0]
        # x_of_node = graph.radius_at_node * np.cos(graph.angle_at_node) - graph.origin[1]

        # sorted_nodes = argsort_points_by_x_then_y((x_of_node, y_of_node))

        ###   nodes = RadialNodeLayout(shape[0], spacing=spacing, origin=origin)

        RadialGraph.__init__(self, shape, spacing=spacing, xy_of_center=xy_of_center)
        # super(DualRadialGraph, self).__init__(
        # DualVoronoiGraph.__init__(self, (y_of_node[sorted_nodes], x_of_node[sorted_nodes]), xy_sort=True,
        # DualVoronoiGraph.__init__(self, (nodes.y_of_node, nodes.x_of_node),
        DualVoronoiGraph.__init__(self, (self.y_of_node, self.x_of_node),
                                  xy_sort=True, rot_sort=True)

    @property
    def shape(self):
        return self._shape

    @property
    def spacing(self):
        return self._spacing

    @property
    def origin(self):
        return self._xy_of_center

    @property
    def xy_of_center(self):
        return self._xy_of_center
