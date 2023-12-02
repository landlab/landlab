"""Python implementation of TriangleMeshGrid, a class used to create and
manage unstructured, irregular grids for 2D numerical models.
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np

from landlab.graph.triangle.dual_triangle import DualTriGraph
from landlab.grid.base import ModelGrid


class TriangleMeshGrid(DualTriGraph, ModelGrid):

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
    TriGraph.from_shapefile
        Constructs the grid from a shapefile, geojson, geopackage, etc.

    Examples
    --------
    """

    def __init__(
        self,
        exterior_y_and_x: tuple[np.ndarray, np.ndarray],
        holes=None,
        triangle_opts="pqDevjz",
        timeout=10,
        reorient_links=False,
        xy_of_reference=(0.0, 0.0),
        xy_axis_name=("x", "y"),
        xy_axis_units="-",
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
        DualTriGraph.__init__(
            self,
            exterior_y_and_x,
            holes=holes,
            triangle_opts=triangle_opts,
            timeout=timeout,
        )
        ModelGrid.__init__(
            self,
            xy_axis_name=xy_axis_name,
            xy_axis_units=xy_axis_units,
            xy_of_reference=xy_of_reference,
        )

        self._node_status = np.full(
            self.number_of_nodes, self.BC_NODE_IS_CORE, dtype=np.uint8
        )
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE

    @classmethod
    def from_dict(cls, kwds):
        """Initialize a new TriangleMeshGrid from a dict with "x" and "y" keys."""
        args = (kwds.pop("x"), kwds.pop("y"))
        return cls(*args, **kwds)

    def plot_nodes_and_links(
        self,
        nodes_args: dict = None,
        links_args: dict = None,
        subplots_args: dict = None,
    ):
        """Produce a plot of nodes and links."""
        if nodes_args is None:
            nodes_args = {}
        if links_args is None:
            links_args = {}
        if subplots_args is None:
            subplots_args = {}

        fig, ax = plt.subplots(**subplots_args)

        for link in np.arange(self.number_of_links):
            head, tail = self.nodes_at_link[link]
            ax.plot(
                [self.x_of_node[head], self.x_of_node[tail]],
                [self.y_of_node[head], self.y_of_node[tail]],
                **links_args,
            )

        ax.scatter(self.x_of_node, self.y_of_node, **nodes_args)

        return fig

    def save(self, path, clobber=False):
        """Save a grid and fields.

        This method uses pickle to save the grid as a pickle file.
        At the time of coding, this is the only convenient output format
        for unstructured grids, but support for netCDF is likely coming.

        All fields will be saved, along with the grid.

        The recommended suffix for the save file is '.grid'. This will
        be added to your save if you don't include it.

        This method is equivalent to
        :py:func:`~landlab.io.native_landlab.save_grid`, and
        :py:func:`~landlab.io.native_landlab.load_grid` can be used to
        load these files.

        Caution: Pickling can be slow, and can produce very large files.
        Caution 2: Future updates to Landlab could potentially render old
        saves unloadable.

        Parameters
        ----------
        path : str
            Path to output file.
        clobber : bool (defaults to false)
            Set to true to allow overwriting

        Returns
        -------
        str
            The name of the saved file (with the ".grid" extension).

        Examples
        --------
        """
        import pickle

        path = pathlib.Path(path)

        if path.suffix != ".grid":
            path = path.with_suffix(path.suffix + ".grid")

        if path.exists() and not clobber:
            raise ValueError(
                f"File exists: {str(path)!r}. "
                "Either remove this file and try again or set the "
                "'clobber' keyword to True"
            )

        with open(path, "wb") as fp:
            pickle.dump(self, fp)

        return str(path)
