"""Python implementation of TriangleMeshGrid, a class used to create and
manage unstructured, irregular grids for 2D numerical models.
"""

import numpy as np
import pathlib

from ..graph.triangle import TriangleGraph
from .base import ModelGrid

class TriangleMeshGrid(TriangleGraph, ModelGrid):

    """This inherited class implements an unstructured grid from dual
    Delaunay and Voronoi graphs. By convention, nodes, links, and patches
    compose a Delaunay triangulation, while corners, faces, and cells
    compose the corresponding Voronoi tesselation. Uses the Triangle
    software package to build the mesh.

    Create an unstructured grid from points whose coordinates are given
    by the arrays *x*, *y*.

    Returns
    -------
    TriangleMeshGrid
        A newly-created grid.

    See also
    --------
    TriangleMeshGrid.from_shapefile
        Constructs the grid from a shapefile, geojson, geopackage, etc.

    Examples
    --------
    """

    def __init__(
        self,
        x=None,
        y=None,
        holes=None,
        triangle_opts="",
        timeout=10,
        reorient_links=False,
        xy_of_reference=(0.0, 0.0),
        xy_axis_name=("x", "y"),
        xy_axis_units="-"
    ):
        """Create a TriangleMeshGrid from a set of points.

        Create an unstructured grid from points whose coordinates are given
        by the arrays *x*, *y*.

        Parameters
        ----------
        x : array_like
            x-coordinate of points
        y : array_like
            y-coordinate of points
        holes : array_like
            (N, 2) shaped array with coordinates of any holes in the domain
        triangle_opts : str
            command-line options for the Triangle meshing software
        timeout : float
            how many seconds to allow Triangle to run before terminating
        reorient_links (optional) : bool
            whether to point all links to the upper-right quadrant
        xy_of_reference : tuple, optional
            Coordinate value in projected space of (0., 0.)
            Default is (0., 0.)

        Returns
        -------
        TriangleMeshGrid
            A newly-created grid.

        Examples
        --------
        """
        TriangleGraph.__init__(
            self, 
            (y, x), 
            holes=holes, 
            triangle_opts=triangle_opts,
            timeout=timeout
        )
        ModelGrid.__init__(
            self,
            xy_axis_name=xy_axis_name,
            xy_axis_units=xy_axis_units,
            xy_of_reference=xy_of_reference
        )

        self._node_status = np.full(
            self.number_of_nodes, self.BC_NODE_IS_CORE, dtype=np.uint8
        )
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE

