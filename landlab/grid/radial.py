#! /usr/env/python
"""
Python implementation of RadialModelGrid, a grid class used to create and
manage structured Voronoi-Delaunay grids for 2D numerical models.

Do NOT add new documentation here. Grid documentation is now built in a semi-
automated fashion. To modify the text seen on the web, edit the files
`docs/text_for_[gridfile].py.txt`.
"""

import numpy as np

from .base import ModelGrid
from ..graph import DualRadialGraph


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
    20

    >>> omg.radius_at_node
    array([ 2.,  2.,  2.,  2.,  1.,  2.,  1.,  1.,  2.,  2.,  0.,  2.,  1.,
            1.,  2.,  1.,  2.,  2.,  2.,  2.])
    """

    def __init__(self, num_shells=0, dr=1.0, origin_x=0.0, origin_y=0.0,
                 **kwds):
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
        >>> omg = RadialModelGrid(num_shells=1, dr=1., origin_x=0.,
        ...                       origin_y=0.)
        >>> omg.number_of_nodes
        7
        >>> omg.number_of_cells
        1

        A second rings will have 13 nodes.

        >>> omg = RadialModelGrid(2)
        >>> omg.number_of_nodes
        20
        """
        shape = (num_shells + 1, 6)
        spacing = dr
        origin = origin_y, origin_x

        DualRadialGraph.__init__(self, shape, spacing=spacing, origin=origin)
        ModelGrid.__init__(self, **kwds)

        self._node_status = np.full(self.number_of_nodes,
                                    self.BC_NODE_IS_CORE, dtype=np.uint8)
        self._node_status[self.perimeter_nodes] = self.BC_NODE_IS_FIXED_VALUE

    @classmethod
    def from_dict(cls, params):
        """
        LLCATS: GINF
        """
        num_shells = params['num_shells']
        dr = params.get('dr', 1.)
        origin = params.get('origin', (0., 0.))

        return cls(num_shells=num_shells, dr=dr, origin_x=origin[0],
                   origin_y=origin[1])
