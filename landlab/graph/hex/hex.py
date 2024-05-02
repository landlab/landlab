r"""
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
array([1. , 2. ,
       0.5, 1.5, 2.5,
       0. , 1. , 2. , 3. ,
       0.5, 1.5, 2.5,
       1. , 2. ])
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
array([[0. , 1. , 2. , 3. ],
       [0.5, 1.5, 2.5, 3.5],
       [0. , 1. , 2. , 3. ]])
>>> graph.number_of_links
23
>>> graph.number_of_patches
12
"""

from functools import cached_property

import numpy as np

from ...core.utils import as_id_array
from ...grid.linkorientation import LinkOrientation
from ...utils.decorators import cache_result_in_object
from ...utils.decorators import make_return_array_immutable
from ..graph import Graph
from ..voronoi.voronoi import DelaunayGraph


class HorizontalRectTriGraphCython:
    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        from .ext.hex import fill_xy_of_node_rect_horizontal

        x_of_node = np.empty(shape[0] * shape[1], dtype=float)
        y_of_node = np.empty(shape[0] * shape[1], dtype=float)
        fill_xy_of_node_rect_horizontal(shape, x_of_node, y_of_node)

        x_of_node[:] *= spacing
        y_of_node[:] *= spacing * np.sin(np.pi / 3.0)
        x_of_node[:] += xy_of_lower_left[0]
        y_of_node[:] += xy_of_lower_left[1]

        return x_of_node, y_of_node


class VerticalRectTriGraphCython:
    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        from .ext.hex import fill_xy_of_node_rect_vertical

        n_rows, n_cols = shape
        x_spacing = np.sin(np.pi / 3.0) * spacing
        y_spacing = spacing

        x_of_node = np.empty(n_rows * n_cols, dtype=float)
        y_of_node = np.empty(n_rows * n_cols, dtype=float)

        fill_xy_of_node_rect_vertical(shape, x_of_node, y_of_node)

        x_of_node *= x_spacing
        y_of_node *= y_spacing
        x_of_node += xy_of_lower_left[0]
        y_of_node += xy_of_lower_left[1]

        return x_of_node, y_of_node

        x_of_node[:, (n_cols + 1) // 2 :] = (
            np.arange(n_cols // 2) * x_spacing * 2.0 + x_spacing + xy_of_lower_left[0]
        )
        x_of_node[:, : (n_cols + 1) // 2] = (
            np.arange((n_cols + 1) // 2) * x_spacing * 2.0 + xy_of_lower_left[0]
        )

        y_of_node[:, : (n_cols + 1) // 2] = (
            np.arange(n_rows) * y_spacing + xy_of_lower_left[1]
        ).reshape((n_rows, 1))
        y_of_node[:, (n_cols + 1) // 2 :] = (
            np.arange(n_rows) * y_spacing + xy_of_lower_left[1] + y_spacing * 0.5
        ).reshape((n_rows, 1))

        return x_of_node.reshape(-1), y_of_node.reshape(-1)


class HorizontalHexTriGraphCython:
    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        from .ext.hex import fill_xy_of_node_hex_horizontal

        n_rows, n_cols = shape
        n_nodes = n_rows * n_cols + (n_rows // 2) ** 2

        x_of_node = np.empty(n_nodes, dtype=float)
        y_of_node = np.empty(n_nodes, dtype=float)

        fill_xy_of_node_hex_horizontal(shape, x_of_node, y_of_node)

        x_of_node[:] *= spacing
        y_of_node[:] *= spacing * np.sin(np.pi / 3.0)
        x_of_node[:] += xy_of_lower_left[0]
        y_of_node[:] += xy_of_lower_left[1]

        return x_of_node, y_of_node


class VerticalHexTriGraphCython:
    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        from .ext.hex import fill_xy_of_node_hex_vertical

        n_rows, n_cols = shape
        n_nodes = n_cols * n_rows + (n_cols // 2) ** 2

        x_of_node = np.empty(n_nodes, dtype=float)
        y_of_node = np.empty(n_nodes, dtype=float)

        fill_xy_of_node_hex_vertical(shape, x_of_node, y_of_node)

        x_of_node[:] *= spacing * np.sin(np.pi / 3.0)
        y_of_node[:] *= spacing
        x_of_node[:] += xy_of_lower_left[0]
        y_of_node[:] += xy_of_lower_left[1]

        return x_of_node, y_of_node


class HorizontalRectTriGraph:
    @staticmethod
    def number_of_nodes(shape):
        n_rows, n_cols = shape
        return n_rows * n_cols

    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        n_rows, n_cols = shape
        x_spacing, y_spacing = spacing, spacing * np.sin(np.pi / 3.0)

        x_of_node, y_of_node = np.meshgrid(
            np.arange(n_cols) * x_spacing + xy_of_lower_left[0],
            np.arange(n_rows) * y_spacing + xy_of_lower_left[1],
        )
        x_of_node[1::2] += spacing * 0.5

        return x_of_node.reshape(-1), y_of_node.reshape(-1)

    @staticmethod
    def corner_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import HorizontalRectTriGraph
        >>> HorizontalRectTriGraph.corner_nodes((3, 4))
        (11, 8, 0, 3)
        >>> HorizontalRectTriGraph.corner_nodes((3, 2))
        (5, 4, 0, 1)
        >>> HorizontalRectTriGraph.corner_nodes((7, 1))
        (6, 6, 0, 0)
        >>> HorizontalRectTriGraph.corner_nodes((1, 3))
        (2, 0, 0, 2)
        """
        n_rows, n_cols = shape
        return (n_rows * n_cols - 1, n_cols * (n_rows - 1), 0, n_cols - 1)

    @staticmethod
    def number_of_perimeter_nodes(shape):
        if 1 in shape:
            return np.prod(shape)
        return 2 * shape[0] + 2 * (shape[1] - 2)

    @staticmethod
    def perimeter_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import HorizontalRectTriGraph
        >>> HorizontalRectTriGraph.perimeter_nodes((3, 2))
        array([1, 3, 5, 4, 2, 0])
        """
        return np.concatenate(HorizontalRectTriGraph.nodes_at_edge(shape))

    @staticmethod
    def nodes_at_edge(shape):
        n_rows, n_cols = shape
        if n_rows == n_cols == 1:
            return (np.array([0]),) + (np.array([], dtype=int),) * 3
        (
            northeast,
            northwest,
            southwest,
            southeast,
        ) = HorizontalRectTriGraph.corner_nodes(shape)

        if n_rows > 1:
            south = np.arange(southwest, southeast)
        else:
            south = np.array([southwest], dtype=int)

        if n_cols > 1:
            west = np.arange(northwest, southwest, -n_cols)
        else:
            west = np.array([northwest], dtype=int)

        return (
            np.arange(southeast, northeast, n_cols),
            np.arange(northeast, northwest, -1),
            west,
            south,
        )


class VerticalRectTriGraph:
    @staticmethod
    def number_of_nodes(shape):
        n_rows, n_cols = shape
        return n_rows * n_cols

    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        n_rows, n_cols = shape
        x_spacing, y_spacing = spacing * np.sin(np.pi / 3.0), spacing

        x_of_node = np.empty((n_rows, n_cols), dtype=float)
        y_of_node = np.empty((n_rows, n_cols), dtype=float)

        x_of_node[:, (n_cols + 1) // 2 :] = (
            np.arange(n_cols // 2) * x_spacing * 2.0 + x_spacing + xy_of_lower_left[0]
        )
        x_of_node[:, : (n_cols + 1) // 2] = (
            np.arange((n_cols + 1) // 2) * x_spacing * 2.0 + xy_of_lower_left[0]
        )

        y_of_node[:, : (n_cols + 1) // 2] = (
            np.arange(n_rows) * y_spacing + xy_of_lower_left[1]
        ).reshape((n_rows, 1))
        y_of_node[:, (n_cols + 1) // 2 :] = (
            np.arange(n_rows) * y_spacing + xy_of_lower_left[1] + y_spacing * 0.5
        ).reshape((n_rows, 1))

        return x_of_node.reshape(-1), y_of_node.reshape(-1)

    @staticmethod
    def corner_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import VerticalRectTriGraph
        >>> VerticalRectTriGraph.corner_nodes((4, 3))
        (10, 9, 0, 1)
        >>> VerticalRectTriGraph.corner_nodes((4, 4))
        (15, 12, 0, 3)
        >>> VerticalRectTriGraph.corner_nodes((3, 2))
        (5, 4, 0, 1)
        >>> VerticalRectTriGraph.corner_nodes((7, 1))
        (6, 6, 0, 0)
        >>> VerticalRectTriGraph.corner_nodes((1, 3))
        (1, 0, 0, 1)
        >>> VerticalRectTriGraph.corner_nodes((2, 3))
        (4, 3, 0, 1)
        """
        n_rows, n_cols = shape
        if n_cols % 2 == 0:
            return (n_rows * n_cols - 1, n_cols * (n_rows - 1), 0, n_cols - 1)
        else:
            return (
                n_rows * n_cols - 1 - n_cols // 2,
                n_cols * (n_rows - 1),
                0,
                n_cols // 2,
            )

    @staticmethod
    def number_of_perimeter_nodes(shape):
        if 1 in shape:
            return np.prod(shape)
        return 2 * shape[1] + 2 * (shape[0] - 2)

    @staticmethod
    def perimeter_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import VerticalRectTriGraph
        >>> VerticalRectTriGraph.perimeter_nodes((3, 2))
        array([1, 3, 5, 4, 2, 0])
        >>> VerticalRectTriGraph.perimeter_nodes((2, 3))
        array([1, 4, 5, 3, 0, 2])
        >>> VerticalRectTriGraph.perimeter_nodes((2, 4))
        array([3, 7, 5, 6, 4, 0, 2, 1])
        """
        return np.concatenate(VerticalRectTriGraph.nodes_at_edge(shape))

    @staticmethod
    def nodes_at_edge(shape):
        n_rows, n_cols = shape
        if n_rows == n_cols == 1:
            return (np.array([0]),) + (np.array([], dtype=int),) * 3
        n_nodes = n_rows * n_cols
        (
            northeast,
            northwest,
            southwest,
            southeast,
        ) = VerticalRectTriGraph.corner_nodes(shape)

        if n_cols == 1:
            southwest = northwest - n_cols

        north = np.empty(n_cols - 1, dtype=int)
        north[::2] = n_nodes - n_cols // 2 + np.arange(n_cols // 2)
        north[1::2] = northwest + np.arange(1, n_cols - n_cols // 2)

        if n_rows > 1:
            south = np.empty(n_cols - 1, dtype=int)
            south[::2] = np.arange(0, n_cols // 2)
            south[1::2] = (n_cols + 1) // 2 + np.arange(n_cols - n_cols // 2 - 1)
        else:
            south = np.array([southwest], dtype=int)

        return (
            np.arange(southeast, northeast, n_cols),
            north[::-1],
            np.arange(northwest, southwest, -n_cols),
            south,
        )


class HorizontalHexTriGraph:
    @staticmethod
    def number_of_nodes(shape):
        n_rows, n_cols = shape
        return n_rows * n_cols + (n_rows // 2) ** 2

    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        n_rows, n_cols = shape
        x_spacing, y_spacing = spacing, spacing * np.sin(np.pi / 3.0)

        length_of_row = np.concatenate(
            (
                np.arange((n_rows + 2) // 2) + n_cols,
                (n_rows + 2) // 2
                + n_cols
                - 1
                - np.arange(1, n_rows - (n_rows + 2) // 2 + 1),
            )
        )
        offset_to_row = np.concatenate((np.array([0]), length_of_row)).cumsum()
        rows = [
            slice(start, end)
            for start, end in zip(offset_to_row[:-1], offset_to_row[1:])
        ]

        y_of_node = np.empty(HorizontalHexTriGraph.number_of_nodes(shape), dtype=float)
        for row, inds in enumerate(rows):
            y_of_node[inds] = row * y_spacing + xy_of_lower_left[1]

        x_of_node = np.empty(HorizontalHexTriGraph.number_of_nodes(shape), dtype=float)

        x_of_row = (
            np.abs((n_rows + 2) // 2 - 1 - np.arange(n_rows)) * x_spacing * 0.5
            + xy_of_lower_left[0]
        )
        for row, inds in enumerate(rows):
            x_of_node[inds] = x_of_row[row] + np.arange(length_of_row[row]) * x_spacing

        return x_of_node.reshape(-1), y_of_node.reshape(-1)

    @staticmethod
    def corner_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import HorizontalHexTriGraph
        >>> HorizontalHexTriGraph.corner_nodes((3, 2))
        (4, 6, 5, 2, 0, 1)
        >>> HorizontalHexTriGraph.corner_nodes((7, 1))
        (9, 15, 15, 6, 0, 0)
        >>> HorizontalHexTriGraph.corner_nodes((6, 1))
        (9, 14, 13, 6, 0, 0)
        >>> HorizontalHexTriGraph.corner_nodes((4, 2))
        (8, 11, 9, 5, 0, 1)
        """
        n_rows, n_cols = shape

        n_nodes_in_middle_row = n_rows // 2 + n_cols
        n_nodes = n_rows * n_cols + (n_rows // 2) ** 2

        east = (n_nodes_in_middle_row + n_cols) * ((n_rows // 2 + 1) // 2) - 1
        if (n_rows // 2) % 2 == 0:
            east += (n_nodes_in_middle_row + n_cols) // 2

        return (
            east,
            n_nodes - 1,
            n_nodes - (n_cols + (n_rows + 1) % 2),
            east - (n_nodes_in_middle_row - 1),
            0,
            n_cols - 1,
        )

    @staticmethod
    def number_of_perimeter_nodes(shape):
        if shape[0] == 1:
            return shape[1]
        return 2 * shape[0] + 2 * (shape[1] - 2) + (shape[0] + 1) % 2

    @staticmethod
    def perimeter_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import HorizontalHexTriGraph
        >>> HorizontalHexTriGraph.perimeter_nodes((3, 2))
        array([4, 6, 5, 2, 0, 1])
        >>> HorizontalHexTriGraph.perimeter_nodes((1, 3))
        array([2, 1, 0])
        """
        return np.concatenate(HorizontalHexTriGraph.nodes_at_edge(shape))

    @staticmethod
    def nodes_at_edge(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import HorizontalHexTriGraph
        >>> HorizontalHexTriGraph.nodes_at_edge((5, 3))
        (array([11, 15]), array([18, 17]), array([16, 12]), array([7, 3]),
         array([0, 1]), array([2, 6]))
        >>> HorizontalHexTriGraph.nodes_at_edge((4, 3))
        (array([11]), array([15, 14, 13]), array([12]), array([7, 3]), array([0, 1]),
         array([2, 6]))
        """
        n_rows, n_cols = shape
        (
            east,
            northeast,
            northwest,
            west,
            southwest,
            southeast,
        ) = HorizontalHexTriGraph.corner_nodes(shape)

        if n_rows == 1:
            nodes_at_south_edge = np.asarray([southwest], dtype=int)
        else:
            nodes_at_south_edge = np.arange(southwest, southeast)

        return (
            as_id_array(
                northeast
                - np.arange(northeast - northwest + 1, east - west + 1).cumsum()[::-1]
            ),
            np.arange(northeast, northwest, -1),
            as_id_array(
                west
                + np.arange(east - west + 1, northeast - northwest + 1, -1).cumsum()[
                    ::-1
                ]
            ),
            as_id_array(
                southwest + np.arange(n_cols, n_rows // 2 + n_cols).cumsum()[::-1]
            ),
            nodes_at_south_edge,
            as_id_array(
                east - np.arange(n_cols + n_rows // 2, n_cols, -1).cumsum()[::-1]
            ),
        )


class VerticalHexTriGraph:
    @staticmethod
    def number_of_nodes(shape):
        n_rows, n_cols = shape
        return n_rows * n_cols + (n_cols // 2) ** 2

    @staticmethod
    def xy_of_node(shape, spacing=1.0, xy_of_lower_left=(0.0, 0.0)):
        n_rows, n_cols = shape
        x_spacing, y_spacing = spacing * np.sin(np.pi / 3.0), spacing

        length_of_middle_rows = np.full(2 * n_rows - 1, n_cols // 2)
        if n_cols % 2 == 1:
            length_of_middle_rows[::2] += 1

        length_of_row = np.concatenate(
            (
                np.arange(1, n_cols // 2 + 1),
                length_of_middle_rows,
                np.arange(n_cols // 2, 0, -1),
            )
        )
        offset_to_row = np.concatenate((np.array([0]), length_of_row)).cumsum()
        rows = [
            slice(start, end)
            for start, end in zip(offset_to_row[:-1], offset_to_row[1:])
        ]

        y_of_node = np.empty(VerticalHexTriGraph.number_of_nodes(shape), dtype=float)
        for row, inds in enumerate(rows):
            y_of_node[inds] = row * y_spacing * 0.5

        x_of_node = np.empty(VerticalHexTriGraph.number_of_nodes(shape), dtype=float)

        x_of_middle_rows = np.zeros(2 * n_rows - 1)
        x_of_middle_rows[1::2] += 1.0
        x_of_row = (
            np.concatenate(
                (
                    np.arange(n_cols // 2, 0, -1),
                    x_of_middle_rows,
                    np.arange(1, n_cols // 2 + 1),
                )
            )
            * x_spacing
        )
        for row, inds in enumerate(rows):
            x_of_node[inds] = (
                x_of_row[row] + np.arange(length_of_row[row]) * 2.0 * x_spacing
            )

        x_of_node += xy_of_lower_left[0]
        y_of_node += xy_of_lower_left[1]

        return x_of_node.reshape(-1), y_of_node.reshape(-1)

    @staticmethod
    def corner_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import VerticalHexTriGraph
        >>> VerticalHexTriGraph.corner_nodes((2, 5))
        (10, 13, 8, 3, 0, 5)
        >>> VerticalHexTriGraph.corner_nodes((2, 3))
        (5, 6, 4, 1, 0, 2)
        >>> VerticalHexTriGraph.corner_nodes((2, 4))
        (10, 11, 7, 3, 0, 2)
        >>> VerticalHexTriGraph.corner_nodes((2, 2))
        (4, 4, 3, 1, 0, 0)
        >>> VerticalHexTriGraph.corner_nodes((3, 1))
        (2, 2, 2, 0, 0, 0)
        >>> VerticalHexTriGraph.corner_nodes((1, 3))
        (2, 3, 1, 1, 0, 2)
        """
        n_rows, n_cols = shape
        n_nodes = n_rows * n_cols + (n_cols // 2) ** 2

        n = n_cols // 2
        tri_nodes = (n + 1) * (n // 2)
        if n % 2 == 1:
            tri_nodes += (n + 1) // 2

        if n_cols % 2 == 0:
            southeast = tri_nodes - 1
            southwest = tri_nodes
            northwest = n_nodes - 1 - (tri_nodes - 1 + n_cols // 2)
            northeast = n_nodes - 1 - (tri_nodes - n)

        else:
            southwest = tri_nodes
            southeast = tri_nodes + n_cols // 2
            northwest = n_nodes - (tri_nodes + (n_cols + 1) // 2)
            northeast = n_nodes - 1 - tri_nodes

        south = 0
        north = n_nodes - 1

        return (northeast, north, northwest, southwest, south, southeast)

    @staticmethod
    def number_of_perimeter_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import VerticalHexTriGraph
        >>> VerticalHexTriGraph.number_of_perimeter_nodes((2, 3))
        6
        >>> VerticalHexTriGraph.number_of_perimeter_nodes((2, 2))
        5
        """
        if shape[1] == 1:
            return shape[0]
        return 2 * shape[1] + 2 * (shape[0] - 2) + (shape[1] + 1) % 2

    @staticmethod
    def perimeter_nodes(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import VerticalHexTriGraph
        >>> VerticalHexTriGraph.perimeter_nodes((3, 7))
        array([ 9, 16, 23, 26, 28, 29, 27, 24, 20, 13,  6,  3,  1,  0,  2,  5])
        >>> VerticalHexTriGraph.perimeter_nodes((2, 3))
        array([2, 5, 6, 4, 1, 0])
        >>> VerticalHexTriGraph.perimeter_nodes((2, 4))
        array([ 2, 6, 10, 11, 9, 7, 3, 1, 0])
        >>> VerticalHexTriGraph.perimeter_nodes((2, 2))
        array([0, 2, 4, 3, 1])
        >>> VerticalHexTriGraph.perimeter_nodes((3, 1))
        array([0, 1, 2])
        """
        return np.concatenate(VerticalHexTriGraph.nodes_at_edge(shape))

    @staticmethod
    def nodes_at_edge(shape):
        """
        Examples
        --------
        >>> from landlab.graph.hex.hex import VerticalHexTriGraph
        >>> VerticalHexTriGraph.nodes_at_edge((3, 7))
        (array([ 9, 16]), array([23, 26, 28]), array([29, 27, 24]), array([20, 13]),
         array([6, 3, 1]), array([0, 2, 5]))
        >>> VerticalHexTriGraph.nodes_at_edge((2, 3))
        (array([2]), array([5]), array([6]), array([4]), array([1]), array([0]))
        >>> VerticalHexTriGraph.nodes_at_edge((2, 4))
        (array([2, 6]), array([10]), array([11, 9]), array([7]), array([3, 1]), array([0]))
        """
        n_rows, n_cols = shape
        (
            northeast,
            north,
            northwest,
            southwest,
            south,
            southeast,
        ) = VerticalHexTriGraph.corner_nodes(shape)

        if shape[1] == 1:
            southwest = northwest - n_cols

        return (
            np.arange(southeast, northeast, n_cols),
            as_id_array(north - np.arange(1, (n_cols + 1) // 2).cumsum())[::-1],
            as_id_array(north - np.arange(1, (n_cols + 2) // 2).cumsum() + 1),
            np.arange(northwest, southwest, -n_cols),
            as_id_array(south + np.arange(1, (n_cols + 2) // 2).cumsum())[::-1],
            as_id_array(south + np.arange(1, (n_cols + 1) // 2).cumsum() - 1),
        )


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
        >>> graph = TriGraph((3, 4), node_layout="rect")
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
        >>> graph = TriGraph((3, 4), node_layout="rect")
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
        >>> graph = TriGraph((3, 4), node_layout="rect")
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
        >>> graph = TriGraph((3, 4), node_layout="rect")
        >>> graph.nodes_at_bottom_edge
        array([0, 1, 2, 3])
        """
        return np.arange(self.shape[1], dtype=int)

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def length_of_link(self):
        return np.full(self.number_of_links, self.spacing, dtype=float)

    @cached_property
    def orientation_of_link(self):
        """Return array of link orientation codes (one value per link).

        Orientation codes are defined by :class:`~.LinkOrientation`;
        1 = E, 2 = ENE, 4 = NNE, 8 = N, 16 = NNW, 32 = ESE (using powers
        of 2 allows for future applications that might want additive
        combinations).

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> import numpy as np
        >>> grid = HexModelGrid((3, 2))
        >>> grid.orientation_of_link
        array([ 1, 16,  4, 16,  4,  1,  1,  4, 16,  4, 16,  1], dtype=uint8)
        >>> grid = HexModelGrid((2, 3), orientation="vertical")
        >>> grid.orientation_of_link
        array([32,  2,  8,  2, 32,  8,  8, 32,  2,  8,  2, 32], dtype=uint8)
        """
        code = np.round(self.angle_of_link * 6.0 / np.pi).astype(np.uint8)
        code[code == 11] = 5
        code[:] = 2**code

        return code

    @cached_property
    def parallel_links_at_link(self):
        """Return similarly oriented links connected to each link.

        Return IDs of links of the same orientation that are connected to
        each given link's tail or head node.

        The data structure is a *numpy* array of shape ``(n_links, 2)`` containing the
        IDs of the "tail-wise" (connected to tail node) and "head-wise" (connected
        to head node) links, or -1 if the link is inactive (e.g., on the perimeter)
        or it has no attached parallel neighbor in the given direction.

        For instance, consider a 3x3 hex, in which link IDs are as shown::

               o---17--o---18--o
              / .     / .     / .
             11 12   13  14  15  16
            /     . /     . /     .
           o---8---o---9---o---10--o
            .     / .     / .     /
             2   3   4   5   6   7
              . /     . /     . /
               o---0---o---1---o

        Here's a mapping of the tail-wise and head-wise links, where
        there are valid parallel links::

               o-------o-------o
              / .     / .     / .
             /   .   /   .   /   .
            /     4 3     6 5     .
           o-----9-o-8--10-o-9-----o
            .    13 12   15 14    /
             .   /   .   /   .   /
              . /     . /     . /
               o-------o-------o

        The corresponding data structure would be mostly filled with -1, but
        for the 11 active links, it would look like::

            3: [[-1, 13],
            4:  [-1, 12],
            5:  [-1, 15],
            6:  [-1, 14],
            8:  [-1,  9],
            9:  [ 8, 10],
            10: [ 9, -1],
            12: [ 4, -1],
            13: [ 3, -1],
            14: [ 6, -1],
            15: [ 5, -1]]

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid((3, 3))
        >>> pll = grid.parallel_links_at_link
        >>> pll[3:16]
        array([[-1, 13],
               [-1, 12],
               [-1, 15],
               [-1, 14],
               [-1, -1],
               [-1,  9],
               [ 8, 10],
               [ 9, -1],
               [-1, -1],
               [ 4, -1],
               [ 3, -1],
               [ 6, -1],
               [ 5, -1]])
        """
        if self.orientation == "horizontal":
            orientations = (LinkOrientation.E, LinkOrientation.NNE, LinkOrientation.NNW)
        else:
            orientations = (LinkOrientation.ENE, LinkOrientation.N, LinkOrientation.ESE)

        links_at_node = self._oriented_links_at_node()

        pll = np.full((self.number_of_links, 2), -1, dtype=int)

        for col, orientation in enumerate(orientations):
            links = self.orientation_of_link == orientation

            if orientation == LinkOrientation.ESE:
                pll[links, 0] = links_at_node[self.node_at_link_tail[links]][:, col]
                pll[links, 1] = links_at_node[self.node_at_link_head[links]][:, col + 3]
            else:
                pll[links, 0] = links_at_node[self.node_at_link_tail[links]][:, col + 3]
                pll[links, 1] = links_at_node[self.node_at_link_head[links]][:, col]

        return pll

    def _oriented_links_at_node(self):
        if self.orientation == "horizontal":
            orientations = (LinkOrientation.E, LinkOrientation.NNE, LinkOrientation.NNW)
        else:
            orientations = (LinkOrientation.ENE, LinkOrientation.N, LinkOrientation.ESE)

        links_at_node = np.full((self.number_of_nodes, 6), -1, dtype=int)

        for col, orientation in enumerate(orientations):
            links = np.where(self.orientation_of_link == orientation)[0]

            if orientation == LinkOrientation.ESE:
                links_at_node[self.node_at_link_tail[links], col + 3] = links
                links_at_node[self.node_at_link_head[links], col] = links
            else:
                links_at_node[self.node_at_link_tail[links], col] = links
                links_at_node[self.node_at_link_head[links], col + 3] = links

        return links_at_node


class TriGraph(HexGraphExtras, DelaunayGraph):
    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import TriGraph

    >>> graph = TriGraph((3, 2))
    >>> graph.number_of_nodes == 6
    True
    >>> np.round(graph.y_of_node * 2.0 / np.sqrt(3))
    array([0., 0., 1., 1., 2., 2.])
    >>> graph.x_of_node
    array([0. , 1. , 0.5, 1.5, 0. , 1. ])
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

        if 1 in shape:
            Graph.__init__(
                self,
                (y_of_node, x_of_node),
                links=list(
                    zip(np.arange(len(y_of_node) - 1), np.arange(1, len(y_of_node)))
                ),
                sort=False,
            )
        else:
            DelaunayGraph.__init__(
                self,
                (y_of_node, x_of_node),
                perimeter_links=perimeter_links,
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

    @cached_property
    @make_return_array_immutable
    def perimeter_nodes(self):
        return self._perimeter_nodes
