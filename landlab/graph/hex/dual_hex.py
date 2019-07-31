from ..dual import DualGraph
from ..voronoi.dual_voronoi import create_dual_graph, DualVoronoiGraph
from .hex import TriGraph, setup_xy_of_node
from .perimeternodes import perimeter_links


class DualHexGraph(DualGraph, TriGraph):

    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import DualHexGraph

    >>> graph = DualHexGraph((3, 2), node_layout='hex')
    >>> graph.number_of_nodes
    7
    >>> graph.number_of_corners
    6

    >>> np.round(graph.y_of_node * 2. / np.sqrt(3))
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  0.,  1.,  1.,  1.,  2.,  2.])
    >>> graph.x_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5,  1.5,  0. ,  1. ,  2. ,  0.5,  1.5])
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
        if 0:
            TriGraph.__init__(
                self,
                shape,
                spacing=spacing,
                xy_of_lower_left=xy_of_lower_left,
                orientation=orientation,
                node_layout=node_layout,
                sort=False,
            )

            if node_layout == "hex":
                max_node_spacing = None
            else:
                max_node_spacing = shape[1] + 1

            dual_graph, node_at_cell, nodes_at_face = create_dual_graph(
                (self.y_of_node, self.x_of_node),
                min_cell_size=6,
                max_node_spacing=max_node_spacing,
            )

            self.merge(dual_graph, node_at_cell=node_at_cell, nodes_at_face=nodes_at_face)

            if sort:
                self.sort()

            return

        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError("spacing must be a float")

        self._shape = tuple(shape)
        self._spacing = spacing

        if node_layout not in ("rect", "hex"):
            raise ValueError("node_layout not understood")
        else:
            self._node_layout = node_layout

        if orientation not in ("horizontal", "vertical"):
            raise ValueError("orientation not understood")
        else:
            self._orientation = orientation

        x_of_node, y_of_node = setup_xy_of_node(
            shape,
            spacing=spacing,
            xy_of_lower_left=xy_of_lower_left,
            orientation=orientation,
            node_layout=node_layout,
        )

        _perimeter_links = perimeter_links(
            self.shape, orientation=self.orientation, node_layout=self.node_layout
        )
        self._perimeter_nodes = _perimeter_links[:, 0].copy()

        DualVoronoiGraph.__init__(
            self,
            (y_of_node, x_of_node),
            perimeter_links=_perimeter_links,
            sort=False,
        )

        if sort:
            self.sort()
