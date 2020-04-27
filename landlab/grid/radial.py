#! /usr/env/python
"""Python implementation of RadialModelGrid, a grid class used to create and
manage structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a
semi- automated fashion. To modify the text seen on the web, edit the
files `docs/text_for_[gridfile].py.txt`.
"""

import numpy as np
import xarray as xr

from ..graph import DualRadialGraph
from .base import ModelGrid


class RadialModelGrid(DualRadialGraph, ModelGrid):

    """Grid of concentric circles.

    This inherited class implements a circular grid in which grid nodes are
    placed at regular radial and semi-regular arc-wise intervals. That is,
    if the radial spacing between *shells* is *dr*, the nodes are placed around
    the circular shell at regular intervals that get as close as possible to
    *dr*. The points are then arranged in a Delaunay triangulation with Voronoi
    cells. Within each ring, nodes are numbered according to Landlab
    convention, from the first node counterclockwise of east. Numbering
    begins at the centermost node and works outwards through the rings.
    """

    def __init__(
        self,
        n_rings=0,
        nodes_in_first_ring=6,
        spacing=1.0,
        xy_of_center=(0.0, 0.0),
        xy_of_reference=(0.0, 0.0),
        xy_axis_name=("x", "y"),
        xy_axis_units="-",
    ):
        """Create a circular grid.

        Create a circular grid in which grid nodes are placed at regular
        radial and semi-regular arc-wise intervals. That is, if the radial
        spacing between *shells* is *dr*, the nodes are placed around the
        circular shell at regular intervals that get as close as possible to
        *dr*.  The points are then arranged in a Delaunay triangulation with
        Voronoi cells.


        Parameters
        ----------
        n_rings : int
            Number of rings in the grid.
        nodes_in_first_ring : int, optional
            Number of nodes in the first ring.
        spacing : float, optional
            Distance between rings.
        xy_of_center : tuple, optional
            (x, y) coordinate of center point. Default
            is (0., 0.)
        xy_of_reference : tuple, optional
            Coordinate value in projected space of the reference point,
            `xy_of_lower_left`. Default is (0., 0.)

        Returns
        -------
        RadialModelGrid
            A newly-created grid.

        Examples
        --------
        A grid with just one ring will have a node at the origin surrounded
        by six other nodes by default. This can be changed by providing the
        keyword argument `nodes_in_first_ring`.

        >>> import numpy as np
        >>> from landlab import RadialModelGrid
        >>> omg = RadialModelGrid(
        ...     n_rings=1, nodes_in_first_ring=8, xy_of_center=(0., 0.)
        ... )
        >>> omg.number_of_nodes
        9
        >>> omg.number_of_cells
        1

        A second rings will have 16 nodes (1 + 8 + 16 = 25).

        >>> omg = RadialModelGrid(2, nodes_in_first_ring=8)
        >>> omg.number_of_nodes
        25

        >>> np.round(omg.radius_at_node)
        array([ 2.,  2.,  2.,  2.,  2.,  1.,  2.,  2.,  1.,  1.,  2.,  1.,  0.,
                1.,  2.,  1.,  1.,  2.,  2.,  1.,  2.,  2.,  2.,  2.,  2.])
        """
        xy_of_center = tuple(xy_of_center)

        DualRadialGraph.__init__(
            self,
            (n_rings, nodes_in_first_ring),
            spacing=spacing,
            xy_of_center=xy_of_center,
            sort=True,
        )
        ModelGrid.__init__(
            self,
            xy_of_reference=xy_of_reference,
            xy_axis_name=xy_axis_name,
            xy_axis_units=xy_axis_units,
        )

        self._node_status = np.full(
            self.number_of_nodes, self.BC_NODE_IS_CORE, dtype=np.uint8
        )
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE
        self._xy_of_center = tuple(xy_of_center)

    @classmethod
    def from_dict(cls, kwds):
        args = ()
        return cls(*args, **kwds)

    @classmethod
    def from_dataset(cls, dataset):
        return cls(
            n_rings=int(dataset["n_rings"]),
            nodes_in_first_ring=int(dataset["nodes_in_first_ring"]),
            spacing=float(dataset["spacing"]),
            xy_of_center=dataset["xy_of_center"],
        )

    def as_dataset(self, include="*", exclude=None):
        dataset = xr.Dataset(
            {
                "n_rings": self.number_of_rings,
                "nodes_in_first_ring": self.number_of_nodes_in_ring[0],
                "spacing": self.spacing_of_rings,
                "xy_of_center": (("dim",), list(self.xy_of_center)),
            },
            attrs={"grid_type": "radial"},
        )
        return dataset.update(
            super(RadialModelGrid, self).as_dataset(include=include, exclude=exclude)
        )

    @property
    def xy_of_center(self):
        """Return (x, y) of the reference point."""
        return self._xy_of_center

    @xy_of_center.setter
    def xy_of_center(self, xy_of_center):
        """Set a new value for the xy_of_lower_left."""
        dx = self.xy_of_center[0] - xy_of_center[0]
        dy = self.xy_of_center[1] - xy_of_center[1]
        with self.thawed():
            self.x_of_node[:] -= dx
            self.y_of_node[:] -= dy
        self._xy_of_center = xy_of_center
