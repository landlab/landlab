""" Implementation of the FramedVoronoiGraph and its static layout:
HorizontalRectVoronoiGraph. This pattern is inspired from the developments of the HexModelGrid

@author sebastien lenard
@date 2022, Aug
"""

from functools import lru_cache

import numpy as np

from ...utils.decorators import make_return_array_immutable
from ..graph import Graph
from ..voronoi.voronoi import DelaunayGraph


class HorizontalRectVoronoiGraph:
    """This static class implements the horizontal rectangular frame for
    the FramedVoronoiGraph.
    """

    @staticmethod
    def number_of_nodes(shape):
        """
        Parameters
        ----------
        shape: tuple of int
            Number of rows and number of columns

        Returns
        -------
        Int
            Number of nodes

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph
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
        seed=(200, 500),
    ):
        """Calculation of the x-y coordinates this way:

        1. Generate a rectangular, regular meshgrid.

        2. Move the coordinates of the core nodes over a random distance around their
        initial position, within a threshold calculated from xy_spacing and
        xy_min_spacing.

        3. Rectify the y-coordinates of the nodes of the left and right to ensure
        that the leftmost node of a row has a lower y than the rightmost node.
        This ensures that the ids of these nodes are not modified by subsequent
        sorting operations on the graph and make it possible to get the
        perimeter nodes in simple way.

        Parameters
        ----------
        shape : tuple of int.
            Number of rows and columns of nodes.
        xy_spacing : float or tuple of float, optional.
            Node spacing along x and y coordinates. If float, same spacing at x and y.
        xy_of_lower_left : tuple, optional.
            Minimum x-of-node and y-of-node values. Depending on the grid.
            no node may be present at this coordinate. Default is (0., 0.).
        xy_min_spacing: float or tuple of float, optional.
            Final minimal spacing between nodes. Random moves of the core nodes
            around their position cannot be above this threshold:
            (xy_spacing - xy_min_spacing) /2
            If float, same minimal spacing for x and y.
        seed: tuple of int, optional
            Seeds used to generate the random x and y moves.
            When set, controls a pseudo-randomness of moves to ensure
            reproducibility.
            When None, seed is random and the moves of coordinates are
            completely random.

        Returns
        -------
        tuple of ndarray(int).
            The arrays of x and y coordinates.

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph
        >>> HorizontalRectVoronoiGraph.xy_of_node((3, 3), seed=(200, 500))[0][3]       # doctest: +NORMALIZE_WHITESPACE
        0.0
        >>> HorizontalRectVoronoiGraph.xy_of_node((3, 3), seed=(200, 500))[0][5]       # doctest: +NORMALIZE_WHITESPACE
        2.0
        >>> HorizontalRectVoronoiGraph.xy_of_node((3, 3), seed=(200, 500))[1][3]       # doctest: +NORMALIZE_WHITESPACE
        0.749
        """

        n_rows, n_cols = shape
        max_move = [
            (xy_spacing[0] - xy_min_spacing[0]) / 2,
            (xy_spacing[1] - xy_min_spacing[1]) / 2,
        ]
        n_nodes = HorizontalRectVoronoiGraph.number_of_nodes(shape)
        n_core_nodes = n_nodes - HorizontalRectVoronoiGraph.number_of_perimeter_nodes(
            shape
        )

        # Generation of a rectangular grid, coordinates must be float
        x_of_node, y_of_node = np.meshgrid(
            np.arange(n_cols, dtype=float) * xy_spacing[0] + xy_of_lower_left[0],
            np.arange(n_rows, dtype=float) * xy_spacing[1] + xy_of_lower_left[1],
        )
        # Randomly move the coordinates of the core nodes of the grid. Move below +/- (spacing - min_spacing)/2
        if seed is None:
            xy_random_generator = (np.random.default_rng(), np.random.default_rng())
        else:
            xy_random_generator = (
                np.random.Generator(np.random.PCG64(seed[0])),
                np.random.Generator(np.random.PCG64(seed[1])),
            )

        x_y_random_move = (
            xy_random_generator[0].random(n_core_nodes) * max_move[0] * 2 - max_move[0],
            xy_random_generator[1].random(n_core_nodes) * max_move[1] * 2 - max_move[1],
        )

        k = 0
        for i in range(1, n_rows - 1):
            for j in range(1, n_cols - 1):
                x_of_node[i, j] += x_y_random_move[0][k]
                y_of_node[i, j] += x_y_random_move[1][k]
                k += 1
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
        >>> from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph
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
        Int
            Number of perimeter nodes

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph
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
        ndarray of Int
            Ids of the perimeter nodes

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph
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
        tuple(ndarray of Int)
            For each edge give the ids of the nodes present at the edge

        Examples
        --------
        >>> from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph
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
    """Graph of a unstructured grid of Voronoi Delaunay cells and
    irregular patches. It is a special type of VoronoiDelaunay graph in which
    the initial set of points is arranged in a fixed lattice (e.g. like a rectangular
    raster grid) named here "layout" and the core points are then moved aroung their
    initial position by a random distance, lower than a certain threshold.

    Examples
    --------
    >>> from landlab.graph import FramedVoronoiGraph

    >>> graph = FramedVoronoiGraph((3, 3), seed=(200, 500))
    >>> graph.number_of_nodes
    9

    >>> graph.x_of_node[2:4]    # doctest: +NORMALIZE_WHITESPACE
    array([ 2.,  0.])
    >>> graph.y_of_node[2:4]    # doctest: +NORMALIZE_WHITESPACE
    array([ 0.   ,  0.749])
    >>> graph.y_of_node[5]    # doctest: +NORMALIZE_WHITESPACE
    1.2509999999999999

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
        seed=(200, 500),
    ):
        """Create the graph.

        Parameters
        ----------
        shape : int or tuple of int
            For a rectangular layout, number of rows and columns of nodes.
            If int, rows number = columns number = value
        xy_spacing : float or tuple of float, optional
            Node spacing along x and y coordinates. If float, same spacing at x and y.
        xy_of_lower_left : tuple, optional
            Minimum x-of-node and y-of-node values. Depending on the grid
            no node may be present at this coordinate. Default is (0., 0.).
        sort: bool
            If true, nodes, links and patches are re-numbered according certain criterias of position.
            Currently not used.
        xy_min_spacing: float or tuple of float, optional
            Final minimal spacing between nodes. Random moves of the core nodes
            around their position cannot be above this threshold:
            (xy_spacing - xy_min_spacing) /2
            If float, same minimal spacing for x and y.
        seed: tuple of int, optional
            Seeds used to generate the random x and y moves.
            When set, controls a pseudo-randomness of moves to ensure
            reproducibility.
            When None, seed is random and the moves of coordinates are
            completely random.

        Returns
        -------
        FramedVoronoiGraph
            A newly-created graph.

        Examples
        --------
        Create a grid with 2 rows and 3 columns of nodes.

        >>> from landlab.graph import FramedVoronoiGraph
        >>> graph = FramedVoronoiGraph((3, 2))
        >>> graph.number_of_nodes
        6
        """
        # 1. Check and format input arguments
        #####################################
        try:
            shape_ = np.asarray(np.broadcast_to(shape, 2))
            self._shape = (int(shape_[0]), int(shape_[1]))
        except TypeError:
            raise TypeError("shape must be a tuple of ints")
        try:
            xy_spacing_ = np.asfarray(np.broadcast_to(xy_spacing, 2))
            self._xy_spacing = (float(xy_spacing_[0]), float(xy_spacing_[1]))
        except TypeError:
            raise TypeError("spacing must be a float or a tuple of floats")
        try:
            xy_of_lower_left_ = np.asfarray(np.broadcast_to(xy_of_lower_left, 2))
            self._xy_of_lower_left = (
                float(xy_of_lower_left_[0]),
                float(xy_of_lower_left_[1]),
            )
        except TypeError:
            raise TypeError("xy of lower left must be a float or a tuple of floats")

        node_layout = self._node_layout = "rect"
        orientation = self._orientation = "horizontal"

        layouts = {
            "horizontal_rect": HorizontalRectVoronoiGraph,
        }
        layout = layouts["_".join([orientation, node_layout])]

        try:
            xy_min_spacing_ = np.asfarray(np.broadcast_to(xy_min_spacing, 2))
            self._xy_min_spacing = (
                float(xy_min_spacing_[0]),
                float(xy_min_spacing_[1]),
            )
        except TypeError:
            raise TypeError("minimal spacing must be a float or a tuple of floats")

        self._seed = seed
        if seed is not None:
            try:
                seed_ = np.asarray(np.broadcast_to(seed, 2))
                self._seed = (int(seed_[0]), int(seed_[1]))
            except TypeError:
                raise TypeError("seed must be None or a tuple of ints")

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

        # 3. Instanciation of the parent class
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

    @property
    @lru_cache()
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
