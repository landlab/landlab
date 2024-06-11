""" Implementation of the FramedVoronoiGraph and its static layout:
HorizontalRectVoronoiGraph. This pattern is inspired from the developments of the HexModelGrid

.. codeauthor:: sebastien lenard
"""

from functools import cached_property

import numpy as np

from ...utils.decorators import make_return_array_immutable
from ..graph import Graph
from ..voronoi.voronoi import DelaunayGraph


class HorizontalRectVoronoiGraph:
    """The horizontal rectangular frame for the FramedVoronoiGraph."""

    @staticmethod
    def number_of_nodes(shape):
        """
        Parameters
        ----------
        shape: tuple of int
            Number of rows and number of columns

        Returns
        -------
        int
            Number of nodes

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import (
        ...     HorizontalRectVoronoiGraph,
        ... )
        >>> HorizontalRectVoronoiGraph.number_of_nodes((3, 2))
        6
        """
        n_rows, n_cols = shape
        return n_rows * n_cols

    @staticmethod
    def xy_of_node(
        shape,
        xy_spacing=(1.0, 1.0),
        xy_of_lower_left=(0.0, 0.0),
        xy_min_spacing=(0.5, 0.5),
        seed=200,
    ):
        """The x and y coordinates of the graph's nodes.

        Calculation of the x-y coordinates is done following these steps:

        1. Generate a rectangular, regular meshgrid.
        2. Move the coordinates of the core nodes over a random distance around their
           initial position, within a threshold calculated from *xy_spacing* and
           *xy_min_spacing*.
        3. Rectify the y-coordinates of the nodes of the left and right to ensure
           that the leftmost node of a row has a lower y than the rightmost node.
           This ensures that the ids of these nodes are not modified by subsequent
           sorting operations on the graph and make it possible to get the
           perimeter nodes in simple way.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of nodes.
        xy_spacing : float or tuple of float, optional
            Node spacing along x and y coordinates. If ``float``, same spacing at *x* and *y*.
        xy_of_lower_left : tuple, optional
            Minimum *x*-of-node and *y*-of-node values. Depending on the grid,
            a node may not be present at this location.
        xy_min_spacing: float or tuple of float, optional
            Final minimal spacing between nodes. Random moves of the core nodes
            around their initial positions cannot be above this threshold:
            ``(xy_spacing - xy_min_spacing) / 2``.  If ``float``, same minimal
            spacing for *x* and *y*.
        seed: int, optional
            Seed used to generate the random *x* and *y* moves. When set, controls a
            pseudo-randomness of moves to ensure reproducibility.
            When ``None``, seed is random and the moves of coordinates are
            completely random.

        Returns
        -------
        x_of_node, y_of_node : ndarray of float
            The arrays of *x* and *y* coordinates.

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import (
        ...     HorizontalRectVoronoiGraph,
        ... )

        >>> x_of_node, y_of_node = HorizontalRectVoronoiGraph.xy_of_node(
        ...     (3, 3), seed=200
        ... )

        Coordinates of the lower left node,

        >>> x_of_node[0], y_of_node[0]
        (0.0, 0.0)

        *x* coordinates of the left and right edges,

        >>> x_of_node[3], x_of_node[5]
        (0.0, 2.0)

        *y* coordinate of the middle row of the left edge,

        >>> y_of_node[3]
        0.749
        """

        n_rows, n_cols = shape
        max_move = (
            (xy_spacing[0] - xy_min_spacing[0]) / 2,
            (xy_spacing[1] - xy_min_spacing[1]) / 2,
        )

        if max_move[0] < 0.0 or max_move[1] < 0.0:
            raise ValueError("minimum spacing must be greater than node spacing")
        if np.allclose(max_move, 0.0):
            raise ValueError("at least one of x and y moves must be greater than zero")

        # Generation of a rectangular grid, coordinates must be float
        x_of_node, y_of_node = np.meshgrid(
            np.arange(n_cols, dtype=float) * xy_spacing[0] + xy_of_lower_left[0],
            np.arange(n_rows, dtype=float) * xy_spacing[1] + xy_of_lower_left[1],
        )
        # Randomly move the coordinates of the core nodes of the grid. Move
        # below +/- (spacing - min_spacing) / 2
        xy_random_generator = np.random.default_rng(seed=seed)

        x_moves = xy_random_generator.uniform(-max_move[0], max_move[0], shape)
        y_moves = xy_random_generator.uniform(-max_move[1], max_move[1], shape)

        x_of_node[1:-1, 1:-1] += x_moves[1:-1, 1:-1]
        y_of_node[1:-1, 1:-1] += y_moves[1:-1, 1:-1]
        # Control the node id attribution for left and right edge. For instance, for a 3x3 grid,
        # make sure that node 3 is at the left of the 2nd row and node 5 at the right.
        # For this, for each core row, set y of the leftmost node as the minimal y of the row
        # and set y of the rightmost node as the maximal y of the row
        for i in range(1, n_rows - 1):
            y_of_node[i, 0] -= max_move[1] + 1.0e-3
            y_of_node[i, n_cols - 1] += max_move[1] + 1.0e-3

        return x_of_node.reshape(-1), y_of_node.reshape(-1)

    @staticmethod
    def corner_nodes(shape):
        """
        Parameters
        ----------
        shape: tuple of int
            Number of rows and number of columns

        Returns
        -------
        ndarray of int
            Ids of the corner nodes

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import (
        ...     HorizontalRectVoronoiGraph,
        ... )
        >>> HorizontalRectVoronoiGraph.corner_nodes((3, 4))
        (11, 8, 0, 3)
        """
        n_rows, n_cols = shape
        return (n_rows * n_cols - 1, n_cols * (n_rows - 1), 0, n_cols - 1)

    @staticmethod
    def number_of_perimeter_nodes(shape):
        """
        Parameters
        ----------
        shape: tuple of int
            Number of rows and number of columns

        Returns
        -------
        int
            Number of perimeter nodes

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import (
        ...     HorizontalRectVoronoiGraph,
        ... )
        >>> HorizontalRectVoronoiGraph.number_of_perimeter_nodes((3, 4))
        10
        """
        if 1 in shape:
            return np.prod(shape)
        return 2 * shape[0] + 2 * (shape[1] - 2)

    @staticmethod
    def perimeter_nodes(shape):
        """
        Parameters
        ----------
        shape: tuple of int
            Number of rows and number of columns

        Returns
        -------
        ndarray of int
            Ids of the perimeter nodes

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import (
        ...     HorizontalRectVoronoiGraph,
        ... )
        >>> HorizontalRectVoronoiGraph.perimeter_nodes((3, 3))
        array([2, 5, 8, 7, 6, 3, 0, 1])
        """
        return np.concatenate(HorizontalRectVoronoiGraph.nodes_at_edge(shape))

    @staticmethod
    def nodes_at_edge(shape):
        """
        Parameters
        ----------
        shape: tuple of int
            Number of rows and number of columns

        Returns
        -------
        right, top, left, bottom : ndarray of int
            For each edge give the ids of the nodes present at the edge

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import (
        ...     HorizontalRectVoronoiGraph,
        ... )
        >>> HorizontalRectVoronoiGraph.nodes_at_edge((3, 3))
        (array([2, 5]), array([8, 7]), array([6, 3]), array([0, 1]))
        """
        n_rows, n_cols = shape
        if n_rows == n_cols == 1:
            return (np.array([0]),) + (np.array([], dtype=int),) * 3
        (
            northeast,
            northwest,
            southwest,
            southeast,
        ) = HorizontalRectVoronoiGraph.corner_nodes(shape)

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


class FramedVoronoiGraph(DelaunayGraph):
    """VoronoiDelaunay graph based on a fixed lattice.

    Graph of an unstructured grid of Voronoi Delaunay cells and
    irregular patches. It is a special type of :class`~.VoronoiDelaunayGraph` in which
    the initial set of points is arranged in a fixed lattice (e.g. like a
    :class:`~.RasterModelGrid`) named here "layout" and the core points are
    then moved from their initial position by a random distance, lower than a
    certain threshold.

    Examples
    --------
    >>> from landlab.graph import FramedVoronoiGraph

    >>> graph = FramedVoronoiGraph((3, 3), seed=200)
    >>> graph.number_of_nodes
    9

    >>> graph.x_of_node[2:4]
    array([2., 0.])
    >>> graph.y_of_node[2:4]
    array([0.   , 0.749])
    >>> graph.y_of_node[5]
    1.251

    >>> graph.number_of_links
    16
    >>> graph.number_of_patches
    8
    """

    def __init__(
        self,
        shape,
        xy_spacing=(1.0, 1.0),
        xy_of_lower_left=(0.0, 0.0),
        sort=False,
        xy_min_spacing=(0.5, 0.5),
        seed=200,
    ):
        """Create the graph.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of nodes.
        xy_spacing : float or tuple of float, optional
            Node spacing along *x* and *y* coordinates. If ``float``, same spacing *x* and *y*
            spacing.
        xy_of_lower_left : tuple, optional
            Minimum *x*-of-node and *y*-of-node values. Depending on the grid,
            a node may not be present at this location.
        sort: bool
            If ``True``, nodes, links and patches are re-numbered according
            certain their positions.  Currently not used.
        xy_min_spacing: float or tuple of float, optional
            Final minimal spacing between nodes. Random moves of the core nodes
            around their position cannot be above this threshold:
            ``(xy_spacing - xy_min_spacing) / 2``
            If ``float``, same minimal spacing for *x* and *y*.
        seed: int, optional
            Seed used to generate the random *x* and *y* moves.
            When set, controls a pseudo-randomness of moves to ensure
            reproducibility. When ``None``, seed is random and the moves of coordinates are
            completely random.

        Returns
        -------
        FramedVoronoiGraph
            A newly-created graph.

        Examples
        --------
        Create a grid with 3 rows and 2 columns of nodes.

        >>> from landlab.graph import FramedVoronoiGraph
        >>> graph = FramedVoronoiGraph((3, 2))
        >>> graph.number_of_nodes
        6
        """
        # 1. Check and format input arguments
        #####################################
        self._shape = shape
        self._seed = seed

        try:
            xy_spacing = np.asfarray(np.broadcast_to(xy_spacing, 2))
        except TypeError as exc:
            raise TypeError("spacing must be a float or a tuple of floats") from exc
        else:
            self._xy_spacing = xy_spacing[0], xy_spacing[1]

        try:
            xy_of_lower_left = np.asfarray(np.broadcast_to(xy_of_lower_left, 2))
        except TypeError as exc:
            raise TypeError(
                "xy of lower left must be a float or a tuple of floats"
            ) from exc
        else:
            self._xy_of_lower_left = xy_of_lower_left[0], xy_of_lower_left[1]

        node_layout = self._node_layout = "rect"
        orientation = self._orientation = "horizontal"

        layouts = {
            "horizontal_rect": HorizontalRectVoronoiGraph,
        }
        layout = layouts["_".join([orientation, node_layout])]

        try:
            xy_min_spacing = np.asfarray(np.broadcast_to(xy_min_spacing, 2))
        except TypeError as exc:
            raise TypeError(
                "minimal spacing must be a float or a tuple of floats"
            ) from exc
        else:
            self._xy_min_spacing = xy_min_spacing[0], xy_min_spacing[1]

        # 2. Construction of the layout and the x-y coordinates of nodes
        ################################################################
        x_of_node, y_of_node = layout.xy_of_node(
            self._shape,
            xy_spacing=self._xy_spacing,
            xy_of_lower_left=self._xy_of_lower_left,
            xy_min_spacing=self._xy_min_spacing,
            seed=self._seed,
        )

        # 3. Determination of the perimeter and edge nodes
        #########################################
        self._perimeter_nodes = layout.perimeter_nodes(self._shape)

        (right, top, left, bottom) = layout.nodes_at_edge(self._shape)
        self._nodes_at_right_edge = np.sort(np.append(right, top[0]))
        self._nodes_at_top_edge = np.sort(np.append(top, left[0]))
        self._nodes_at_left_edge = np.sort(np.append(left, bottom[0]))
        self._nodes_at_bottom_edge = np.sort(np.append(bottom, right[0]))

        perimeter_links = np.empty((len(self._perimeter_nodes), 2), dtype=int)
        perimeter_links[:, 0] = self._perimeter_nodes
        perimeter_links[:-1, 1] = self._perimeter_nodes[1:]
        perimeter_links[-1, 1] = self._perimeter_nodes[0]

        self._x_of_node = x_of_node
        self._y_of_node = y_of_node
        self._perimeter_links = perimeter_links

        # 3. Instantiation of the parent class
        ######################################
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

    @property
    def shape(self):
        return self._shape

    @property
    def xy_spacing(self):
        return self._xy_spacing

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

    @property
    @make_return_array_immutable
    def nodes_at_right_edge(self):
        return self._nodes_at_right_edge

    @property
    @make_return_array_immutable
    def nodes_at_top_edge(self):
        return self._nodes_at_top_edge

    @property
    @make_return_array_immutable
    def nodes_at_left_edge(self):
        return self._nodes_at_left_edge

    @property
    @make_return_array_immutable
    def nodes_at_bottom_edge(self):
        return self._nodes_at_bottom_edge
