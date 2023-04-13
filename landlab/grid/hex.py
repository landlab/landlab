#! /usr/env/python
"""Python implementation of HexModelGrid, a grid class used to create and
manage structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a
semi- automated fashion. To modify the text seen on the web, edit the
files `docs/text_for_[gridfile].py.txt`.
"""
from functools import cached_property

import numpy
import xarray as xr

from ..core.utils import as_id_array
from ..graph import DualHexGraph
from .base import ModelGrid
from .linkorientation import LinkOrientation


class HexModelGrid(DualHexGraph, ModelGrid):
    """A grid of hexagonal cells.

    This inherited class implements a regular 2D grid with hexagonal cells and
    triangular patches. It is a special type of VoronoiDelaunay grid in which
    the initial set of points is arranged in a triangular/hexagonal lattice.

    Examples
    --------
    Create a hex grid with 2 rows of nodes. The first and third rows will
    have 2 nodes, and the second nodes.

    >>> from landlab import HexModelGrid
    >>> grid = HexModelGrid((3, 2), spacing=1.0)
    >>> grid.number_of_nodes
    7

    >>> grid = HexModelGrid((3, 3), node_layout="rect", spacing=2.0)
    >>> grid.status_at_node
    array([1, 1, 1, 1, 0, 1, 1, 1, 1], dtype=uint8)
    >>> grid = HexModelGrid((3, 3), node_layout="rect", orientation="vertical")
    >>> grid.status_at_node
    array([1, 1, 1, 1, 1, 0, 1, 1, 1], dtype=uint8)
    >>> grid = HexModelGrid((4, 4), node_layout='rect', orientation="vertical")
    >>> grid.status_at_node
    array([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
    >>> grid.boundary_nodes
    array([ 0,  1,  2,  3,  4,  7,  8, 11, 12, 13, 14, 15])
    >>> grid = HexModelGrid((3, 4), node_layout="rect")
    >>> grid.status_at_node
    array([1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=uint8)
    """

    def __init__(
        self,
        shape,
        spacing=1.0,
        xy_of_lower_left=(0.0, 0.0),
        orientation="horizontal",
        node_layout="hex",
        reorient_links=True,
        xy_of_reference=(0.0, 0.0),
        xy_axis_name=("x", "y"),
        xy_axis_units="-",
    ):
        """Create a grid of hexagonal cells.

        Create a regular 2D grid with hexagonal cells and triangular patches.
        It is a special type of :class:`~.VoronoiModelGrid` in which the initial set
        of points is arranged in a triangular/hexagonal lattice.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of nodes.
        spacing : float, optional
            Node spacing.
        xy_of_lower_left : tuple of float, optional
            Minimum x-of-node and y-of-node values. Depending on the grid
            no node may be present at this coordinate. Default is (0., 0.).
        xy_of_reference : tuple of float, optional
            Coordinate value in projected space of the reference point,
            `xy_of_lower_left`. Default is (0., 0.)
        orientation : str, optional
            One of the 3 cardinal directions in the grid, either 'horizontal'
            (default) or 'vertical'
        node_layout : {"hex", "rect"}
            The grid layout of nodes.
        reorient_links : bool, optional
            Whether or not to re-orient all links to point between -45 deg
            and +135 deg clockwise from "north" (i.e., along y axis). default
            is True.

        Returns
        -------
        HexModelGrid
            A newly-created grid.

        Examples
        --------
        Create a hex grid with 2 rows of nodes. The first and third rows will
        have 2 nodes, and the second nodes.

        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid((3, 2), spacing=1.0)
        >>> hmg.number_of_nodes
        7
        """
        self._xy_of_lower_left = tuple(numpy.asfarray(xy_of_lower_left))

        DualHexGraph.__init__(
            self,
            shape,
            spacing=spacing,
            xy_of_lower_left=self.xy_of_lower_left,
            orientation=orientation,
            node_layout=node_layout,
            sort=True,
        )
        ModelGrid.__init__(
            self,
            xy_axis_name=xy_axis_name,
            xy_axis_units=xy_axis_units,
            xy_of_reference=xy_of_reference,
        )

        self._node_status = numpy.full(
            self.number_of_nodes, self.BC_NODE_IS_CORE, dtype=numpy.uint8
        )
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE

    @classmethod
    def from_dict(cls, kwds):
        args = (kwds.pop("shape"),)
        return cls(*args, **kwds)

    @classmethod
    def from_dataset(cls, dataset):
        return cls(
            tuple(dataset["shape"].values),
            spacing=dataset["spacing"],
            xy_of_lower_left=dataset["xy_of_lower_left"],
            orientation=dataset.attrs["orientation"],
            node_layout=dataset.attrs["node_layout"],
        )

    def as_dataset(self, include="*", exclude=None, time=None):
        dataset = xr.Dataset(
            {
                "shape": (("dim",), list(self.shape)),
                "spacing": self.spacing,
                "xy_of_lower_left": (("dim",), list(self.xy_of_lower_left)),
            },
            attrs={
                "grid_type": "triangular",
                "node_layout": self.node_layout,
                "orientation": self.orientation,
            },
        )
        return dataset.update(
            super().as_dataset(include=include, exclude=exclude, time=None)
        )

    @property
    def xy_of_lower_left(self):
        """Return (x, y) of the reference point."""
        return self._xy_of_lower_left

    @xy_of_lower_left.setter
    def xy_of_lower_left(self, xy_of_lower_left):
        """Set a new value for the xy_of_lower_left."""
        dx = self.xy_of_lower_left[0] - xy_of_lower_left[0]
        dy = self.xy_of_lower_left[1] - xy_of_lower_left[1]
        # self._xy_of_node -= (dx, dy)
        with self.thawed():
            self.x_of_node[:] -= dx
            self.y_of_node[:] -= dy

        self._xy_of_lower_left = tuple(xy_of_lower_left)

    @property
    def number_of_node_columns(self):
        """Number of node columns hex grid.

        Number of node columns in a rectangular-shaped and/or
        vertically oriented hex grid.

        Returns the number of columns, including boundaries.

        Notes
        -----
        Will generate an error if called with a hex-shaped, horizontally
        aligned grid.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid((5, 5), node_layout="rect")
        >>> grid.number_of_node_columns
        5

        :meta landlab: info-grid, info-node
        """
        return self.shape[1]

    @property
    def number_of_node_rows(self):
        """Number of node rows in a rectangular-shaped and/or horizontally
        oriented hex grid.

        Returns the number of rows, including boundaries.

        Notes
        -----
        Will generate an error if called with a hex-shaped, vertically
        aligned grid.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid((5, 5), node_layout="rect")
        >>> grid.number_of_node_rows
        5

        :meta landlab: info-grid, info-node
        """
        return self._shape[0]

    def node_row_and_column(self, node_id):
        """Row and column from node ID, FOR VERT RECT CONFIGURATION ONLY.

        Examples
        --------
        >>> from landlab import HexModelGrid
        >>> grid = HexModelGrid((3, 4), node_layout='rect', orientation="vertical")
        >>> grid.node_row_and_column(5)
        (1, 2)
        >>> grid = HexModelGrid((3, 5), node_layout='rect', orientation="vertical")
        >>> grid.node_row_and_column(13)
        (2, 1)
        """
        assert self.orientation[0] == "v", "grid orientation must be vertical"
        try:
            (nr, nc) = self._shape
        except AttributeError as exc:
            raise AttributeError(
                "Only rectangular Hex grids have defined rows and columns."
            ) from exc

        row = node_id // nc
        n_mod_nc = node_id % nc
        half_nc = (nc + 1) // 2
        col = 2 * (n_mod_nc % half_nc) + n_mod_nc // half_nc
        return (row, col)

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
        array([ 1, 16,  4, 16,  4,  1,  1,  4, 16,  4, 16,  1])
        >>> grid = HexModelGrid((2, 3), orientation="vertical")
        >>> grid.orientation_of_link
        array([32,  2,  8,  2, 32,  8,  8, 32,  2,  8,  2, 32])
        """
        code = numpy.round(self.angle_of_link * 6.0 / numpy.pi).astype(numpy.uint8)
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

        pll = numpy.full((self.number_of_links, 2), -1, dtype=int)

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

        links_at_node = numpy.full((self.number_of_nodes, 6), -1, dtype=int)

        for col, orientation in enumerate(orientations):
            links = numpy.where(self.orientation_of_link == orientation)[0]

            if orientation == LinkOrientation.ESE:
                links_at_node[self.node_at_link_tail[links], col + 3] = links
                links_at_node[self.node_at_link_head[links], col] = links
            else:
                links_at_node[self.node_at_link_tail[links], col] = links
                links_at_node[self.node_at_link_head[links], col + 3] = links

        return links_at_node

    def _configure_hexplot(self, data, data_label=None, color_map=None):
        """Sets up necessary information for making plots of the hexagonal grid
        colored by a given data element.

        Parameters
        ----------
        data : str or (n_link,) ndarray
            Data field to be colored
        data_label : str, optional
            Label for colorbar
        color_map : matplotlib colormap object, optional
            Color map to apply (defaults to "jet")

        Notes
        -----
        Creates and stores a PatchCollection representing the hexagons. Also
        stores a handle to the current plotting axis. Both of these are then
        used by hexplot().
        """
        import matplotlib
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Polygon
        from numpy import array, sqrt, zeros

        # color
        if color_map is None:
            color_map = matplotlib.cm.jet

        # geometry
        apothem = self.spacing / 2.0
        # distance from node to each hexagon cell vertex
        radius = 2.0 * apothem / sqrt(3.0)

        # offsets from node x,y position
        offsets = zeros((6, 2))
        poly_verts = zeros((6, 2))

        # Figure out whether the orientation is horizontal or vertical
        if self.orientation[0] == "h":  # horizontal
            offsets[:, 0] = array([0.0, apothem, apothem, 0.0, -apothem, -apothem])
            offsets[:, 1] = array(
                [
                    radius,
                    radius / 2.0,
                    -radius / 2.0,
                    -radius,
                    -radius / 2.0,
                    radius / 2.0,
                ]
            )
        else:  # vertical
            offsets[:, 0] = array(
                [
                    radius / 2.0,
                    radius,
                    radius / 2.0,
                    -radius / 2.0,
                    -radius,
                    -radius / 2.0,
                ]
            )
            offsets[:, 1] = array([apothem, 0.0, -apothem, -apothem, 0.0, apothem])

        patches = []
        for i in range(self.number_of_nodes):
            poly_verts[:, 0] = self.node_x[i] + offsets[:, 0]
            poly_verts[:, 1] = self.node_y[i] + offsets[:, 1]
            p = Polygon(poly_verts, True)
            patches.append(p)

        self._hexplot_pc = PatchCollection(
            patches, cmap=color_map, edgecolor="none", linewidth=0.0
        )

        self._hexplot_configured = True

    def hexplot(self, data, data_label=None, color_map=None):
        """Create a plot of the grid elements.

        Creates a plot of the grid and one node-data field, showing hexagonal
        cells colored by values in the field.

        Parameters
        ----------
        data : str or (n_nodes,) ndarray
            Data field to be colored.
        data_label : str, optional
            Label for colorbar.
        color_map : matplotlib colormap object, None
            Color map to apply (defaults to "jet")

        See also
        --------
        :func:`~.plot.imshow_grid`
            Another Landlab function capable of producing hexplots, with a
            fuller-featured set of options.


        :meta landlab: info-grid
        """
        import copy

        import matplotlib.pyplot as plt
        from numpy import amax, amin, array

        try:
            self._hexplot_configured
        except AttributeError:
            self._configure_hexplot(data, data_label, color_map)
        else:
            if self._hexplot_pc.cmap != color_map:
                self._configure_hexplot(data, data_label, color_map)

        # Handle *data*: if it's a numpy array, then we consider it the
        # data to be plotted. If it's a string, we consider it the name of the
        # node-field to plot, and we fetch it.
        if type(data) is str:
            data_label = data
            data = self.at_node[data]

        ax = plt.gca()
        self._hexplot_pc.set_array(array(data))
        copy_of_pc = copy.copy(self._hexplot_pc)
        ax.add_collection(copy_of_pc)
        plt.xlim([amin(self.node_x) - self.spacing, amax(self.node_x) + self.spacing])
        plt.ylim([amin(self.node_y) - self.spacing, amax(self.node_y) + self.spacing])

        return ax

    def set_watershed_boundary_condition_outlet_id(
        self, outlet_id, node_data, nodata_value=-9999.0
    ):
        """Set the boundary conditions for a watershed.

        All nodes with ``nodata_value`` are set to :attr:`~.NodeStatus.CLOSED`.
        All nodes with data values are set to :attr:`~.NodeStatus.CORE`, with the
        exception that the outlet node is set to a :attr:`~.NodeStatus.FIXED_VALUE`.

        Note that the outer ring of the HexModelGrid is set to
        :attr:`~.NodeStatus.CLOSED`, even if there are nodes that have values.
        The only exception to this would be if the outlet node is on the boundary,
        which is acceptable.

        Assumes that the id of the outlet is already known.

        This assumes that the grid has a single watershed.  If this is not
        the case this will not work.

        Parameters
        ----------
        outlet_id : int
            id of the outlet node
        node_data : str or (n_nodes,) ndarray
            At-node field name or at-node data values to use for identifying
            watershed location.
        nodata_value : float, optional
            Value that indicates an invalid value.

        Examples
        --------
        The example will use a *HexModelGrid* with node data values
        as illustrated::

                1. ,  2. ,  3. ,  4. ,
            0.5,  1.5,  2.5,  3.5,  4.5,
          0. ,  1. ,  2. ,  3. ,  4. ,  5.,
            0.5,  1.5,  2.5,  3.5,  4.5,
                1. ,  2. ,  3. ,  4.

        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid((5, 4))
        >>> z = hmg.add_zeros("topographic__elevation", at="node")
        >>> z += hmg.x_of_node + 1.0

        >>> hmg.status_at_node
        array([1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1,
           1], dtype=uint8)

        >>> outlet = hmg.set_watershed_boundary_condition_outlet_id(9, z, -9999.)
        >>> hmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 0, 4, 1, 0, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4,
           4], dtype=uint8)

        :meta landlab: boundary-condition
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # make ring of no data nodes
        self.status_at_node[self.boundary_nodes] = self.BC_NODE_IS_CLOSED

        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # set the boundary condition (fixed value) at the outlet_node
        self.status_at_node[outlet_id] = self.BC_NODE_IS_FIXED_VALUE

    def set_watershed_boundary_condition(
        self, node_data, nodata_value=-9999.0, return_outlet_id=False
    ):
        """Find the node adjacent to a boundary node with the smallest value.

        This node is set as the outlet.  The outlet node must have a data
        value.  Can return the outlet id as a one element numpy array if
        ``return_outlet_id`` is set to `True`.

        All nodes with ``nodata_value`` are set to :attr:`~.NodeStatus.CLOSED`
        (grid.status_at_node == 4). All nodes with data values are set to
        :attr:`~.NodeStatus.CORE` (grid.status_at_node == 0), with the exception
        that the outlet node is set to a :attr:`~.NodeStatus.FIXED_VALUE`
        (grid.status_at_node == 1).

        Note that the outer ring (perimeter) of the grid is set to
        :attr:`~.NodeStatus.CLOSED`, even if there are nodes that have values. The only
        exception to this would be if the outlet node is on the perimeter, which
        is acceptable.

        This routine assumes that all of the ``nodata_value`` are on the outside of
        the data values. In other words, there are no islands of ``nodata_value``
        surrounded by nodes with data.

        This also assumes that the grid has a single watershed (that is a single
        outlet node).

        Parameters
        ----------
        node_data : str or (n_nodes,) ndarray
            At-node field name or at-node data values to use for identifying
            watershed location.
        nodata_value : float, optional
            Value that indicates an invalid value.
        return_outlet_id : bool, optional
            Indicates whether or not to return the id of the found outlet

        Examples
        --------
        The example will use a :class:`~.HexModelGrid` with node data values
        as illustrated::

                1. ,  2. ,  3. ,  4. ,
            0.5,  1.5,  2.5,  3.5,  4.5,
          0. ,  1. ,  2. ,  3. ,  4. ,  5.,
            0.5,  1.5,  2.5,  3.5,  4.5,
                1. ,  2. ,  3. ,  4.

        >>> from landlab import HexModelGrid
        >>> hmg = HexModelGrid((5, 4))
        >>> z = hmg.add_zeros("topographic__elevation", at="node")
        >>> z += hmg.x_of_node + 1.0
        >>> out_id = hmg.set_watershed_boundary_condition(z, -9999., True)
        >>> out_id
        array([9])
        >>> hmg.status_at_node
        array([4, 4, 4, 4, 4, 0, 0, 0, 4, 1, 0, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4,
           4], dtype=uint8)

        :meta landlab: boundary-condition
        """
        # get node_data if a field name
        node_data = self.return_array_or_field_values("node", node_data)

        # make ring of no data nodes
        self.status_at_node[self.boundary_nodes] = self.BC_NODE_IS_CLOSED

        # set no data nodes to inactive boundaries
        self.set_nodata_nodes_to_closed(node_data, nodata_value)

        # locs is a list that contains locations where
        # node data is not equal to the nodata value
        locs = numpy.where(node_data != nodata_value)
        if len(locs) < 1:
            raise ValueError("All data values are no_data values")

        # now find minimum of the data values
        min_val = numpy.min(node_data[locs])

        # now find where minimum values are
        min_locs = numpy.where(node_data == min_val)[0]

        # check all the locations with the minimum value to see if one
        not_found = True
        while not_found:
            # now check the min locations to see if any are next to
            # a boundary node
            local_not_found = True
            # next_to_boundary = []

            # check all nodes rather than selecting the first node that meets
            # the criteria
            # for i in range(len(min_locs)):
            #     next_to_boundary.append(self.node_has_boundary_neighbor()[min_locs[i])]
            next_to_boundary = self.node_has_boundary_neighbor()[(min_locs,)]
            # if any of those nodes were adjacent to the boundary, check
            # that  there is only one. If only one, set as outlet loc, else,
            # raise a value error
            if numpy.any(next_to_boundary):
                local_not_found = False
                if sum(next_to_boundary) > 1:
                    potential_locs = min_locs[
                        numpy.where(numpy.asarray(next_to_boundary))[0]
                    ]
                    raise ValueError(
                        "Grid has two potential outlet nodes."
                        "They have the following node IDs: \n"
                        + str(potential_locs)
                        + "\nUse the method set_watershed_boundary_condition_outlet_id "
                        "to explicitly select one of these "
                        "IDs as the outlet node."
                    )
                else:
                    outlet_loc = min_locs[numpy.where(next_to_boundary)[0][0]]

            # checked all of the min vals, (so done with inner while)
            # and none of the min values were outlet candidates
            if local_not_found:
                # need to find the next largest minimum value
                # first find the locations of all values greater
                # than the old minimum
                # not done with outer while
                locs = numpy.where((node_data > min_val) & (node_data != nodata_value))
                # now find new minimum of these values
                min_val = numpy.min(node_data[locs])
                min_locs = numpy.where(node_data == min_val)[0]
            else:
                # if locally found, it is also globally found
                # so done with outer while
                not_found = False

        # set outlet boundary condition
        self.status_at_node[outlet_loc] = self.BC_NODE_IS_FIXED_VALUE

        if return_outlet_id:
            return as_id_array(numpy.array([outlet_loc]))
