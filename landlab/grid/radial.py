#! /usr/env/python
"""
Python implementation of RadialModelGrid, a grid class used to create and
manage structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""

import numpy as np

from landlab.utils.decorators import deprecated

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

    Parameters
    ----------
    num_shells : int
        Number of rings in the grid.
    dr : float, optional
        Radial interval for rings.
    origin_x : float, optional
        x-coordinate of origin node.
    origin_y : float, optional
        y-coordinate of origin node.

    Returns
    -------
    RadialModelGrid
        A newly-created grid.

    Examples
    --------
    A grid with just one ring will have a node at the origin surrounded
    by six other nodes.

    >>> from landlab import RadialModelGrid
    >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0., origin_y=0.)
    >>> omg.number_of_nodes
    7
    >>> omg.number_of_cells
    1

    A second rings will have 13 nodes.

    >>> omg = RadialModelGrid(2)
    >>> omg.number_of_nodes
    19
    """

    def __init__(self, num_shells=0, dr=1.0, origin_x=0.0, origin_y=0.0, **kwds):
        """Create a circular grid.

        Create a circular grid in which grid nodes are placed at regular
        radial and semi-regular arc-wise intervals. That is, if the radial
        spacing between *shells* is *dr*, the nodes are placed around the
        circular shell at regular intervals that get as close as possible to
        *dr*.  The points are then arranged in a Delaunay triangulation with
        Voronoi cells.

        Parameters
        ----------
        num_shells : int
            Number of rings in the grid.
        dr : float, optional
            Radial interval for rings.
        origin_x : float, optional
            x-coordinate of origin node.
        origin_y : float, optional
            y-coordinate of origin node.

        Returns
        -------
        RadialModelGrid
            A newly-created grid.

        Examples
        --------
        A grid with just one ring will have a node at the origin surrounded
        by six other nodes.

        >>> from landlab import RadialModelGrid
        >>> omg = RadialModelGrid(num_shells=1, dr=1.)
        >>> omg.number_of_nodes
        7
        >>> omg.number_of_cells
        1

        A second rings will have 13 nodes.

        >>> omg = RadialModelGrid(2)
        >>> omg.number_of_nodes
        19
        """
        if "xy_of_center" not in kwds:
            DualRadialGraph.__init__(
                self,
                (num_shells, int(2. * np.pi / dr)),
                **kwds,
                xy_of_center=(origin_x, origin_y)
            )
        else:
            DualRadialGraph.__init__(self, (num_shells, int(2. * np.pi / dr)), **kwds)
        ModelGrid.__init__(self, **kwds)

    @classmethod
    def from_dict(cls, params):
        """
        LLCATS: GINF
        """
        num_shells = params["num_shells"]
        dr = params.get("dr", 1.)
        origin = params.get("origin", (0., 0.))

        return cls(num_shells=num_shells, dr=dr, origin_x=origin[0], origin_y=origin[1])

    def _initialize(self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        [pts, npts] = self._create_radial_points(num_shells, dr)
        self._n_shells = int(num_shells)
        self._dr = dr
        super(RadialModelGrid, self)._initialize(pts[:, 0], pts[:, 1])

    def _create_radial_points(self, num_shells, dr, origin_x=0.0, origin_y=0.0):
        """Create a set of points on concentric circles.

        Creates and returns a set of (x,y) points placed in a series of
        concentric circles around the origin.
        """
        shells = np.arange(0, num_shells) + 1
        twopi = 2 * np.pi
        # number of points in each shell
        n_pts_in_shell = np.round(twopi * shells)
        dtheta = twopi / n_pts_in_shell
        npts = int(sum(n_pts_in_shell) + 1)
        pts = np.zeros((npts, 2))
        r = shells * dr
        startpt = 1
        for i in np.arange(0, num_shells):
            theta = dtheta[i] * np.arange(0, n_pts_in_shell[i]) + dtheta[i] / (i + 1)
            ycoord = r[i] * np.sin(theta)
            if np.isclose(ycoord[-1], 0.):
                # this modification necessary to force the first ring to
                # follow our new CCW from E numbering convention (DEJH, Nov15)
                ycoord[-1] = 0.
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 0] = np.roll(
                    r[i] * np.cos(theta), 1
                )
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 1] = np.roll(
                    ycoord, 1
                )
            else:
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 0] = r[i] * np.cos(
                    theta
                )
                pts[startpt : (startpt + int(n_pts_in_shell[i])), 1] = ycoord
            startpt += int(n_pts_in_shell[i])
        pts[:, 0] += origin_x
        pts[:, 1] += origin_y

        return pts, npts

    @property
    @deprecated(use="spacing_of_shells", version=1.0)
    def shell_spacing(self):
        """Fixed distance between shells.

        LLCATS: DEPR GINF MEAS
        """
        return self.spacing_of_shells
