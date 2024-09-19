import numpy as np

from ..dual import DualGraph
from ..voronoi.dual_voronoi import DualVoronoiGraph
from .hex import HorizontalHexTriGraph
from .hex import HorizontalRectTriGraph
from .hex import TriGraph
from .hex import VerticalHexTriGraph
from .hex import VerticalRectTriGraph


class DualHexGraph(DualGraph, TriGraph):
    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import DualHexGraph

    >>> graph = DualHexGraph((3, 2), node_layout="hex")
    >>> graph.number_of_nodes
    7
    >>> graph.number_of_corners
    6

    >>> np.round(graph.y_of_node * 2.0 / np.sqrt(3))
    array([0., 0., 1., 1., 1., 2., 2.])
    >>> graph.x_of_node
    array([0.5, 1.5, 0. , 1. , 2. , 0.5, 1.5])
    """

    def __init__(
        self,
        shape,
        spacing=1.0,
        xy_of_lower_left=(0.0, 0.0),
        orientation="horizontal",
        node_layout="rect",
        sort=False,
    ):
        """Create a structured grid of triangles.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of the hex grid. The first value
            is the number of nodes in the first column and the second the
            number of nodes in the first column.
        spacing : float, optional
            Length of links.
        xy_of_lower_left : tuple of float, optional
            Coordinates of lower-left corner of the grid.
        orientation: {'horizontal', 'vertical'}
            Specify if triangles should be laid out in rows or columns.
        node_layout: {'rect', 'hex'}
            Specify the overall layout of the nodes. Use *rect* for
            the layout to approximate a rectangle and *hex* for
            a hexagon.
        """
        if node_layout not in ("rect", "hex"):
            raise ValueError("node_layout not understood")

        if orientation not in ("horizontal", "vertical"):
            raise ValueError("orientation not understood")

        layouts = {
            "horizontal_hex": HorizontalHexTriGraph,
            "vertical_hex": VerticalHexTriGraph,
            "horizontal_rect": HorizontalRectTriGraph,
            "vertical_rect": VerticalRectTriGraph,
        }
        layout = layouts["_".join([orientation, node_layout])]

        try:
            spacing = float(spacing)
        except TypeError as exc:
            raise TypeError("spacing must be a float") from exc

        self._shape = tuple(shape)
        self._spacing = spacing
        self._orientation = orientation
        self._node_layout = node_layout

        x_of_node, y_of_node = layout.xy_of_node(
            shape, spacing=spacing, xy_of_lower_left=xy_of_lower_left
        )
        self._perimeter_nodes = layout.perimeter_nodes(shape)

        perimeter_links = np.empty((len(self._perimeter_nodes), 2), dtype=int)
        perimeter_links[:, 0] = self._perimeter_nodes
        perimeter_links[:-1, 1] = self._perimeter_nodes[1:]
        perimeter_links[-1, 1] = self._perimeter_nodes[0]

        DualVoronoiGraph.__init__(
            self, (y_of_node, x_of_node), perimeter_links=perimeter_links, sort=False
        )

        if sort:
            self.sort()
