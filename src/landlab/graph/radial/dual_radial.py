import numpy as np

from ..dual import DualGraph
from ..voronoi.dual_voronoi import DualVoronoiGraph
from .radial import RadialGraph
from .radial import RadialGraphLayout


class DualRadialGraph(DualGraph, RadialGraph):
    """Graph of a series of points on concentric circles.

    Examples
    --------
    >>> from landlab.graph import DualRadialGraph
    >>> graph = DualRadialGraph((1, 4), sort=True)
    >>> graph.number_of_corners
    4
    >>> graph.y_of_corner
    array([-0.5, -0.5,  0.5,  0.5])
    >>> graph.x_of_corner
    array([-0.5,  0.5, -0.5,  0.5])
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
            Coordinates of the center of the grid.
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

        DualVoronoiGraph.__init__(self, (y_of_node, x_of_node), sort=False)

        if sort:
            self.sort()

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
