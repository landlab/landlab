import numpy as np

from ..voronoi.voronoi import VoronoiGraph


def number_of_nodes(shape, node_layout="rect"):
    """Get the number of nodes in a hex graph.

    Parameters
    ----------
    shape : tuple of int
        Number of rows and columns of the hex grid. The first value
        is the number of nodes in the first column and the second the
        number of nodes in the first column.

    Examples
    --------
    >>> from landlab.graph.hex.hex import number_of_nodes
    >>> number_of_nodes((3, 2))
    6
    >>> number_of_nodes((2, 4))
    8
    >>> number_of_nodes((4, 2))
    8

    >>> number_of_nodes((3, 2), node_layout='hex')
    7
    >>> number_of_nodes((2, 4), node_layout='hex')
    9
    >>> number_of_nodes((4, 2), node_layout='hex')
    12

    >>> number_of_nodes((3, 2), node_layout='rect1')
    7
    >>> number_of_nodes((2, 4), node_layout='rect1')
    9
    >>> number_of_nodes((4, 2), node_layout='rect1')
    10
    """
    if node_layout not in ("rect", "hex", "rect1"):
        raise ValueError("node_layout not understood")

    n_rows, n_cols = shape

    if node_layout == "rect":
        return n_rows * n_cols
    elif node_layout == "hex":
        return n_rows * n_cols + (n_rows // 2) ** 2
    elif node_layout == "rect1":
        return (2 * n_cols + 1) * (n_rows // 2) + n_cols * (n_rows % 2)


def setup_perimeter_nodes(shape, orientation="horizontal", node_layout="rect"):
    from .ext.hex import fill_perimeter_nodes, fill_hex_perimeter_nodes

    n_perimeter_nodes = 2 * shape[0] + 2 * (shape[1] - 2)
    if node_layout in ("hex", "rect1"):
        n_perimeter_nodes += (shape[0] + 1) % 2
    perimeter_nodes = np.empty(n_perimeter_nodes, dtype=int)

    if node_layout == "hex":
        fill_hex_perimeter_nodes(shape, perimeter_nodes)
    else:
        fill_perimeter_nodes(shape, perimeter_nodes)

    return perimeter_nodes


def setup_xy_of_node(
    shape, spacing=1.0, origin=(0.0, 0.0), orientation="horizontal", node_layout="rect"
):
    """Create arrays of coordinates of a node on a hex grid.

    Parameters
    ----------
    shape : tuple of int
        Number of rows and columns of the hex grid. The first value
        is the number of nodes in the first column and the second the
        number of nodes in the first column.
    spacing : float, optional
        Length of links.
    origin : tuple of float, optional
        (y, x) coordinates of lower-left corner of the grid.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph.hex.hex import setup_xy_of_node

    >>> x, y = setup_xy_of_node((3, 2))
    >>> x
    array([ 0. ,  1. ,  0.5,  1.5,  0. ,  1. ])
    >>> y / (np.sqrt(3) / 2.)
    array([ 0.,  0.,  1.,  1.,  2.,  2.])
    >>> x, y = setup_xy_of_node((2, 2))
    >>> x
    array([ 0. ,  1. ,  0.5,  1.5])
    >>> y / (np.sqrt(3) / 2.)
    array([ 0.,  0.,  1.,  1.])

    >>> x, y = setup_xy_of_node((2, 2), spacing=2, origin=(1, 2))
    >>> x
    array([ 2.,  4.,  3.,  5.])
    >>> (y - 1) / (np.sqrt(3) / 2.)
    array([ 0.,  0.,  2.,  2.])
    """
    from .ext.hex import get_xy_of_node, fill_xy_of_node, fill_hex_xy_of_node

    if orientation == "vertical":
        return setup_xy_of_node(
            (shape[1], shape[0]),
            spacing=spacing,
            origin=(origin[1], origin[0]),
            node_layout=node_layout,
            orientation="horizontal",
        )[::-1]

    n_nodes = number_of_nodes(shape, node_layout=node_layout)

    x_of_node = np.empty((n_nodes,), dtype=float)
    y_of_node = np.empty((n_nodes,), dtype=float)

    if node_layout == "rect":
        fill_xy_of_node(shape, x_of_node, y_of_node)
    elif node_layout == "hex":
        fill_hex_xy_of_node(shape, x_of_node, y_of_node)
    elif node_layout == "rect1":
        get_xy_of_node(shape, x_of_node, y_of_node)

    x_of_node *= spacing
    x_of_node += origin[1]
    y_of_node *= spacing * np.sin(np.pi / 3.0)
    y_of_node += origin[0]

    return (x_of_node, y_of_node)


class HexGraph(VoronoiGraph):

    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import HexGraph

    >>> graph = HexGraph((3, 2))
    >>> graph.number_of_nodes == 6
    True
    >>> np.round(graph.y_of_node * 2. / np.sqrt(3))
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  0.,  1.,  1.,  2.,  2.])
    >>> graph.x_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0. ,  1. ,  0.5,  1.5,  0. ,  1. ])
    """

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
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError("spacing must be a float")

        self._shape = tuple(shape)

        if node_layout not in ("rect", "hex", "rect1"):
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
            origin=origin,
            orientation=orientation,
            node_layout=node_layout,
        )
        if node_layout == "hex":
            max_node_spacing = shape[1] + shape[0] / 2 + 2
            max_node_spacing = None
        elif node_layout == "rect":
            max_node_spacing = shape[1] + 1
        elif node_layout == "rect1":
            max_node_spacing = shape[1] + 1

        VoronoiGraph.__init__(
            self,
            (y_of_node, x_of_node),
            xy_sort=True,
            rot_sort=True,
            max_node_spacing=max_node_spacing,
        )

    @property
    def shape(self):
        return self._shape

    @property
    def orientation(self):
        return self._orientation

    @property
    def node_layout(self):
        return self._node_layout

    @property
    # @store_result_in_grid()
    def perimeter_nodes(self):
        return setup_perimeter_nodes(self.shape, self.orientation, self.node_layout)
