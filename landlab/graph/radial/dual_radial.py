from ..voronoi import DualVoronoiGraph
from .radial import create_xy_of_node


class DualRadialGraph(DualVoronoiGraph):

    """Graph of a series of points on concentric circles.

    Examples
    --------
    >>> from landlab.graph import DualRadialGraph
    >>> graph = DualRadialGraph((1, 4))
    >>> graph.number_of_corners == 4
    True
    >>> graph.y_of_corner
    array([-0.5, -0.5,  0.5,  0.5])
    >>> graph.x_of_corner
    array([-0.5,  0.5, -0.5,  0.5])
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

        super(DualRadialGraph, self).__init__(
            (y_of_node, x_of_node), xy_sort=True, rot_sort=True
        )
