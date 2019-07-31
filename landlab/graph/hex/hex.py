"""
Examples
--------

::

        * - *
       / \ / \
      * - * - *
     / \ / \ / \
    * - * - * - *
     \ / \ / \ /
      * - * - *
       \ / \ /
        * - *

>>> from landlab.graph import TriGraph
>>> graph = TriGraph((5, 2), node_layout="hex", sort=True)
>>> graph.number_of_nodes
14
>>> graph.x_of_node
array([ 1. ,  2. ,
        0.5,  1.5,  2.5,
        0. ,  1. ,  2. ,  3. ,
        0.5,  1.5, 2.5,
        1. ,  2. ])
>>> graph.number_of_links
29
>>> graph.number_of_patches
16

::

    * - * - * - *
     \ / \ / \ / \
      * - * - * - *
     / \ / \ / \ /
    * - * - * - *

>>> from landlab.graph import TriGraph
>>> graph = TriGraph((3, 4), orientation="horizontal", node_layout="rect", sort=True)
>>> graph.number_of_nodes
12
>>> graph.x_of_node.reshape((3, 4))
array([[ 0. ,  1. ,  2. ,  3. ],
       [ 0.5,  1.5,  2.5,  3.5],
       [ 0. ,  1. ,  2. ,  3. ]])
>>> graph.number_of_links
23
>>> graph.number_of_patches
12
"""


import numpy as np

from ...utils.decorators import cache_result_in_object, make_return_array_immutable
from ..voronoi.voronoi import DelaunayGraph
from .ext.hex import (
    fill_xy_of_node_hex_horizontal,
    fill_xy_of_node_hex_vertical,
    fill_xy_of_node_rect_horizontal,
    fill_xy_of_node_rect_vertical,
)
from .perimeternodes import perimeter_links, perimeter_nodes


def number_of_nodes(shape, node_layout="rect", orientation="horizontal"):
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
    if orientation == "vertical":
        return number_of_nodes(
            shape[::-1], node_layout=node_layout, orientation="horizontal"
        )
    if node_layout not in ("rect", "hex", "rect1"):
        raise ValueError("node_layout not understood")

    n_rows, n_cols = shape

    if node_layout == "rect":
        return n_rows * n_cols
    elif node_layout == "hex":
        return n_rows * n_cols + (n_rows // 2) ** 2
    elif node_layout == "rect1":
        return (2 * n_cols + 1) * (n_rows // 2) + n_cols * (n_rows % 2)


def setup_xy_of_node(
    shape,
    spacing=1.0,
    xy_of_lower_left=(0.0, 0.0),
    orientation="horizontal",
    node_layout="rect",
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
    xy_of_lower_left : tuple of float, optional
        (x, y) coordinates of lower-left corner of the grid.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph.hex.hex import setup_xy_of_node

    >>> x, y = setup_xy_of_node((3, 2))
    >>> x
    array([ 0. ,  1. ,  0.5,  1.5,  0. ,  1. ])
    >>> np.round(y / (np.sqrt(3) / 2.))
    array([ 0.,  0.,  1.,  1.,  2.,  2.])
    >>> x, y = setup_xy_of_node((2, 2))
    >>> x
    array([ 0. ,  1. ,  0.5,  1.5])
    >>> np.round(y / (np.sqrt(3) / 2.))
    array([ 0.,  0.,  1.,  1.])

    >>> x, y = setup_xy_of_node((2, 2), spacing=2, xy_of_lower_left=(2, 1))
    >>> x
    array([ 2.,  4.,  3.,  5.])
    >>> np.round((y - 1) / (np.sqrt(3) / 2.))
    array([ 0.,  0.,  2.,  2.])
    """
    fill_xy_of_node = {
        "rect": {
            "horizontal": fill_xy_of_node_rect_horizontal,
            "vertical": fill_xy_of_node_rect_vertical,
        },
        "hex": {
            "horizontal": fill_xy_of_node_hex_horizontal,
            "vertical": fill_xy_of_node_hex_vertical,
        },
    }

    n_nodes = number_of_nodes(shape, orientation=orientation, node_layout=node_layout)

    x_of_node = np.empty((n_nodes,), dtype=float)
    y_of_node = np.empty((n_nodes,), dtype=float)

    fill_xy_of_node[node_layout][orientation](shape, x_of_node, y_of_node)

    if orientation == "horizontal":
        x_of_node *= spacing
        y_of_node *= spacing * np.round(np.sin(np.pi / 3.0), 5)
    else:
        x_of_node *= spacing * np.sin(np.pi / 3.0)
        y_of_node *= spacing
    x_of_node += xy_of_lower_left[0]
    y_of_node += xy_of_lower_left[1]

    return (x_of_node, y_of_node)


class HexGraphExtras:
    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def nodes_at_right_edge(self):
        """Get nodes along the right edge.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import TriGraph
        >>> graph = TriGraph((3, 4), node_layout='rect')
        >>> graph.nodes_at_right_edge
        array([ 3,  7, 11])
        """
        return np.arange(
            self.shape[1] - 1, self.shape[0] * self.shape[1], self.shape[1], dtype=int
        )

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def nodes_at_top_edge(self):
        """Get nodes along the top edge.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import TriGraph
        >>> graph = TriGraph((3, 4), node_layout='rect')
        >>> graph.nodes_at_top_edge
        array([ 8,  9, 10, 11])
        """
        return np.arange(
            self.number_of_nodes - self.shape[1], self.number_of_nodes, dtype=int
        )

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def nodes_at_left_edge(self):
        """Get nodes along the left edge.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import TriGraph
        >>> graph = TriGraph((3, 4), node_layout='rect')
        >>> graph.nodes_at_left_edge
        array([0, 4, 8])
        """
        return np.arange(0, self.shape[0] * self.shape[1], self.shape[1], dtype=int)

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def nodes_at_bottom_edge(self):
        """Get nodes along the bottom edge.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.graph import TriGraph
        >>> graph = TriGraph((3, 4), node_layout='rect')
        >>> graph.nodes_at_bottom_edge
        array([0, 1, 2, 3])
        """
        return np.arange(self.shape[1], dtype=int)

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def perimeter_nodes(self):
        return setup_perimeter_nodes(self.shape, self.orientation, self.node_layout)

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def length_of_link(self):
        return np.full(self.number_of_links, self.spacing, dtype=float)


class TriGraph(HexGraphExtras, DelaunayGraph):

    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import TriGraph

    >>> graph = TriGraph((3, 2))
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
        xy_of_lower_left=(0.0, 0.0),
        orientation="horizontal",
        node_layout="rect",
        sort=True,
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

        DelaunayGraph.__init__(
            self,
            (y_of_node, x_of_node),
            perimeter_links=_perimeter_links,
            sort=False,
        )

        if sort:
            self.sort()

    @property
    def shape(self):
        return self._shape

    @property
    def spacing(self):
        return self._spacing

    @property
    def orientation(self):
        return self._orientation

    @property
    def node_layout(self):
        return self._node_layout

    @property
    def perimeter_nodes(self):
        return self._perimeter_nodes
