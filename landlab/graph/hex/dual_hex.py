from ..dual import DualGraph
from ..voronoi.dual_voronoi import create_dual_graph
from .hex import TriGraph


class DualHexGraph(DualGraph, TriGraph):
    def __init__(
        self,
        shape,
        spacing=1.0,
        origin=(0.0, 0.0),
        orientation="horizontal",
        node_layout="rect",
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
        origin : tuple of float, optional
            Coordinates of lower-left corner of the grid.
        orientation: {'horizontal', 'vertical'}
            Specify if triangles should be laid out in rows or columns.
        node_layout: {'rect', 'hex', 'rect1'}
            Specify the overall layout of the nodes. Use *rect* for
            the layout to approximate a rectangle and *hex* for
            a hexagon.
        """
        TriGraph.__init__(
            self,
            shape,
            spacing=spacing,
            origin=origin,
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

        self.sort()
