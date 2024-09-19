#! /usr/env/python

"""Create 2D grid with randomly generated fractures.

Created: September 2013 by Greg Tucker
Last significant modification: conversion to proper component 7/2019 GT
"""

import numpy as np

from landlab import Component
from landlab import HexModelGrid
from landlab import RasterModelGrid


def _calc_fracture_starting_position_raster(shape):
    """Choose a random starting position along one of the sides of the grid.

    Parameters
    ----------
    shape : tuple of int
        Number of rows and columns in the grid

    Returns
    -------
    (c, r) : tuple of int
        Fracture starting coordinates (column and row IDs)
    """
    grid_side = np.random.randint(0, 3)  # east, north, west, south

    if (grid_side % 2) == 0:  # east or west
        c = (1 - grid_side // 2) * (shape[1] - 1)
        r = np.random.randint(0, shape[0] - 1)
    else:
        c = np.random.randint(0, shape[1] - 1)
        r = (1 - grid_side // 2) * (shape[0] - 1)

    return (c, r)


def _calc_fracture_starting_position_and_angle_hex(shape, is_horiz, spacing):
    """Choose a random starting position along one of the sides of the grid.

    Parameters
    ----------
    shape : tuple of int
        Number of rows and columns in the grid

    Returns
    -------
    (x, y, ang) : tuple of float
        Fracture starting coordinates and angle (radians)
    """
    grid_side = np.random.randint(0, 3)  # east, north, west, south
    ang = np.pi * np.random.rand()
    if is_horiz:
        row_spacing = 0.5 * 3.0**0.5 * spacing
        col_spacing = spacing
    else:
        row_spacing = spacing
        col_spacing = 0.5 * 3.0**0.5 * spacing

    if (grid_side % 2) == 0:  # east or west
        c = (1 - grid_side // 2) * (shape[1] - 1)
        r = np.random.randint(0, shape[0] - 1)
        if grid_side == 0:  # east
            ang += np.pi / 2
        else:  # west
            ang += 1.5 * np.pi
    else:
        c = np.random.randint(0, shape[1] - 1)
        r = (1 - grid_side // 2) * (shape[0] - 1)
        if grid_side == 1:  # north
            ang += np.pi

    x = c * col_spacing
    y = r * row_spacing

    epsilon = 0.001 * spacing  # tiny offset to ensure points start inside grid
    if x > epsilon:
        x -= epsilon
    if y > epsilon:
        y -= epsilon

    return (x, y, ang)


def _calc_fracture_orientation(coords, shape):
    """Choose a random orientation for the fracture.

    Parameters
    ----------
    coords : tuple of int
        Starting coordinates (one of which should be zero) as *x*, *y*.
    shape : tuple of int
        Number of rows and columns

    Returns
    -------
    ang : float
        Fracture angle relative to horizontal

    Notes
    -----
    The angle depends on which side of the grid the fracture starts on:
        - east: 90-270
        - north: 180-360
        - west: 270-450 (mod 360)
        - south: 0-180
    """
    x, y = coords
    ang = np.pi * np.random.rand()

    if x == shape[1] - 1:  # east
        ang += np.pi / 2
    elif y == shape[0] - 1:  # north
        ang += np.pi
    elif x == 0:  # west
        ang += 1.5 * np.pi

    return ang


def _calc_fracture_step_sizes(ang):
    """Calculate the sizes of steps dx and dy to be used when "drawing" the
    fracture onto the grid.

    Parameters
    ----------
    start_yx : tuple of int
        Starting grid coordinates
    ang : float
        Fracture angle relative to horizontal (radians)

    Returns
    -------
    (dy, dx) : tuple of float
        Step sizes in y and x directions. One will always be unity, and the
        other will always be <1.

    Examples
    --------
    >>> np.round(_calc_fracture_step_sizes(0 * np.pi / 6), 3)
    array([1., 0.])
    >>> np.round(_calc_fracture_step_sizes(1 * np.pi / 6), 3)
    array([1.   , 0.577])
    >>> np.round(_calc_fracture_step_sizes(2 * np.pi / 6), 3)
    array([0.577, 1.   ])
    >>> np.round(_calc_fracture_step_sizes(3 * np.pi / 6), 3)
    array([0., 1.])
    >>> np.round(_calc_fracture_step_sizes(4 * np.pi / 6), 3)
    array([-0.577,  1.   ])
    >>> np.round(_calc_fracture_step_sizes(5 * np.pi / 6), 3)
    array([-1.   ,  0.577])
    >>> np.round(_calc_fracture_step_sizes(6 * np.pi / 6), 3)
    array([-1.,  0.])
    >>> np.round(_calc_fracture_step_sizes(7 * np.pi / 6), 3)
    array([-1.   , -0.577])
    >>> np.round(_calc_fracture_step_sizes(8 * np.pi / 6), 3)
    array([-0.577, -1.   ])
    >>> np.round(_calc_fracture_step_sizes(9 * np.pi / 6), 3)
    array([-0., -1.])
    >>> np.round(_calc_fracture_step_sizes(10 * np.pi / 6), 3)
    array([ 0.577, -1.   ])
    >>> np.round(_calc_fracture_step_sizes(11 * np.pi / 6), 3)
    array([ 1.   , -0.577])
    >>> np.round(_calc_fracture_step_sizes(12 * np.pi / 6), 3)
    array([ 1., -0.])
    """
    dx = np.cos(ang)
    dy = np.sin(ang)
    multiplier = 1.0 / max(np.abs(dx), np.abs(dy))
    dx *= multiplier
    dy *= multiplier

    return dx, dy


def _trace_fracture_through_grid_raster(m, start_xy, spacing):
    """Create a 2D fracture in a grid.

    Creates a "fracture" in a 2D grid, m, by setting cell values to unity
    along the trace of the fracture (i.e., "drawing" a line throuh the
    grid).

    Parameters
    ----------
    m : 2D Numpy array
        Array that represents the grid
    start_xy : tuple of int
        Starting grid coordinates (col, row) for fracture
    spacing : tuple of float
        Step sizes in x and y directions

    Returns
    -------
    None, but changes contents of m
    """
    x0, y0 = start_xy
    dx, dy = spacing

    x = x0
    y = y0

    while (
        round(x) < np.size(m, 1)
        and round(y) < np.size(m, 0)
        and round(x) >= 0
        and round(y) >= 0
    ):
        m[int(y + 0.5)][int(x + 0.5)] = 1
        x += dx
        y += dy


class FractureGridGenerator(Component):
    """Create a 2D grid with randomly generated fractures.

    The grid contains the value 1 where fractures (one cell wide) exist, and
    0 elsewhere. The idea is to use this for simulations based on weathering
    and erosion of, and/or flow within, fracture networks.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((5, 5))
    >>> fg = FractureGridGenerator(grid=grid, frac_spacing=3)
    >>> fg.run_one_step()
    >>> grid.at_node["fracture_at_node"].reshape((5, 5))
    array([[1, 0, 0, 1, 0],
           [0, 1, 1, 1, 1],
           [0, 0, 0, 1, 1],
           [0, 0, 0, 1, 1],
           [0, 0, 0, 1, 0]], dtype=int8)

    Notes
    -----
    Potential improvements:

    - Fractures could be defined by links rather than nodes (i.e., return a
        link array with a code indicating whether the link crosses a fracture
        or not)
    - Fractures could have a finite length rather than extending all the way
        across the grid
    - Use of starting position along either x or y axis makes fracture net
        somewhat asymmetric. One would need a different algorithm to make it
        fully (statistically) symmetric.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed

    """

    _name = "FractureGridGenerator"

    _unit_agnostic = True

    _info = {
        "fracture_at_node": {
            "dtype": np.int8,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "presence (1) or absence (0) of fracture",
        }
    }

    def __init__(self, grid, frac_spacing=10.0, seed=0):
        """Initialize the FractureGridGenerator.

        Parameters
        ----------
        frac_spacing : int, optional
            Average spacing of fractures (in grid cells) (default = 10)
        seed : int, optional
            Seed used for random number generator (default = 0)

        """

        self._frac_spacing = frac_spacing
        self._seed = seed
        super().__init__(grid)

        if isinstance(grid, RasterModelGrid):
            self._make_frac_grid = self.make_frac_grid_raster
        elif isinstance(grid, HexModelGrid):
            self._make_frac_grid = self.make_frac_grid_hex
        else:
            raise TypeError("grid must be RasterModelGrid or HexModelGrid")

        if "fracture_at_node" not in grid.at_node:
            grid.add_zeros("node", "fracture_at_node", dtype=np.int8)

        np.random.seed(seed)

    def run_one_step(self):
        """Run FractureGridGenerator and create a random fracture grid."""
        self._make_frac_grid(self._frac_spacing)

    def make_frac_grid_raster(self, frac_spacing):
        """Create a raster grid that contains a network of random fractures.

        Creates a grid containing a network of random fractures, which are
        represented as 1's embedded in a grid of 0's. The grid is stored in
        the "fracture_at_node" field.

        Parameters
        ----------
        frac_spacing : int
            Average spacing of fractures (in grid cells)
        """
        # Make an initial grid of all zeros. If user specified a model grid,
        # use that. Otherwise, use the given dimensions.
        nr = self._grid.number_of_node_rows
        nc = self._grid.number_of_node_columns
        m = self._grid.at_node["fracture_at_node"].reshape((nr, nc))

        # Add fractures to grid
        nfracs = (nr + nc) // frac_spacing
        for _ in range(nfracs):
            (c, r) = _calc_fracture_starting_position_raster((nr, nc))
            ang = _calc_fracture_orientation((c, r), (nr, nc))
            (dx, dy) = _calc_fracture_step_sizes(ang)

            _trace_fracture_through_grid_raster(m, (c, r), (dx, dy))

    def make_frac_grid_hex(self, frac_spacing):
        """Create a hex grid that contains a network of random fractures.

        Creates a grid containing a network of random fractures, which are
        represented as 1's embedded in a grid of 0's. The grid is stored in
        the "fracture_at_node" field.

        Parameters
        ----------
        frac_spacing : int
            Average spacing of fractures (in # of grid-cell widths)
        """
        # Make an initial grid of all zeros
        nr = self._grid.number_of_node_rows
        nc = self._grid.number_of_node_columns
        m = self._grid.at_node["fracture_at_node"]  # .reshape((nr, nc))

        # Add fractures to grid
        nfracs = (nr + nc) // frac_spacing
        for _ in range(nfracs):
            (x, y, ang) = _calc_fracture_starting_position_and_angle_hex(
                (nr, nc),
                is_horiz=(self._grid.orientation[0] == "h"),
                spacing=self._grid.spacing,
            )

            dx = 0.5 * np.cos(ang)
            dy = 0.5 * np.sin(ang)

            # the following is a TERRIBLE brute-force, algorithm, with its only
            # redeeming feature being that it was quick to code
            xmax = np.amax(self._grid.x_of_node)
            ymax = np.amax(self._grid.y_of_node)
            while x >= 0 and x <= xmax and y >= 0 and y <= ymax:
                distx2 = (self._grid.x_of_node - x) ** 2
                disty2 = (self._grid.y_of_node - y) ** 2
                dist = np.sqrt(distx2 + disty2)
                closest_node = np.argmin(dist)
                x += dx
                y += dy
                m[closest_node] = 1
