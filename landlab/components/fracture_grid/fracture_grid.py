#! /usr/env/python

"""Create 2D grid with randomly generated fractures.

Created: September 2013 by Greg Tucker
Last significant modification: conversion to proper component 7/2019 GT
"""

import numpy as np

from landlab import Component, RasterModelGrid, HexModelGrid


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
    grid_side = np.random.randint(0, 3) # east, north, west, south

    if (grid_side % 2) == 0:  # east or west
        c = (1 - grid_side // 2) * (shape[1] - 1)
        r = np.random.randint(0, shape[0] - 1)
    else:
        c = np.random.randint(0, shape[1] - 1)
        r = (1 - grid_side // 2) * (shape[0] - 1)
    print('starting on side ' + str(grid_side) + ' at ' + str((c,r)))
    return (c, r)

def _calc_fracture_starting_position_hex(shape, is_horiz):
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
    grid_side = np.random.randint(0, 3) # east, north, west, south

    if (grid_side % 2) == 0:  # east or west
        c = (1 - grid_side // 2) * (shape[1] - 1)
        r = np.random.randint(0, shape[0] - 1)
    else:
        c = np.random.randint(0, shape[1] - 1)
        r = (1 - grid_side // 2) * (shape[0] - 1)
    print('starting on side ' + str(grid_side) + ' at ' + str((c,r)))
    return (c, r)

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
    print('in here x = ' + str(x))
    ang = np.pi * np.random.rand()
    print('raw ang ' + str(np.degrees(ang)))
    print(coords)
    if x == shape[1] - 1:  # east
        ang += np.pi / 2
        print('east ' + str(np.degrees(ang)))
    elif y == shape[0] - 1:  # north
        ang += np.pi
        print('n ' + str(np.degrees(ang)))
    elif x == 0:  # west
        ang += 1.5 * np.pi
        print('w ' + str(np.degrees(ang)))
    else:
        print('s ' + str(np.degrees(ang)))
    print('ang ' + str(np.degrees(ang)))
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
    array([ 1.,  0.])
    >>> np.round(_calc_fracture_step_sizes(1 * np.pi / 6), 3)
    array([ 1.   , 0.577])
    >>> np.round(_calc_fracture_step_sizes(2 * np.pi / 6), 3)
    array([ 0.577,  1.   ])
    >>> np.round(_calc_fracture_step_sizes(3 * np.pi / 6), 3)
    array([ 0., 1.])
    >>> np.round(_calc_fracture_step_sizes(4 * np.pi / 6), 3)
    array([-0.577, 1.   ])
    >>> np.round(_calc_fracture_step_sizes(5 * np.pi / 6), 3)
    array([-1.   , 0.577])
    >>> np.round(_calc_fracture_step_sizes(6 * np.pi / 6), 3)
    array([-1., 0.])
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
    print('dy=' + str(dy) + ', dx=' + str(dx))
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


def _coords_to_row_col_hex(x, y, is_horiz):
    """Convert (x,y) coordinates to (row,col) in a rect-layout hex grid."""
    if is_horiz:
        row = int((2 / 3.0**0.5) * y + 0.5)
        col = int(x + 0.5)  # col number is rounded x coordinate
    else:
        row = int(y + 0.5)  # row number is rounded y coordinate
        col = int((2 / 3.0**0.5) * x + 0.5)


def _trace_fracture_through_grid_hex(m, start_xy, spacing, is_horiz=True):
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
    root3over2 = 3.0**0.5 / 2
    dx *= 0.5
    dy *= 0.5

    x = x0
    y = y0

    if is_horiz:
        while (
            round(x) < np.size(m, 1)
            and round(y) < np.size(m, 0)
            and round(x) >= 0
            and round(y) >= 0
        ):
            xshift = 0.5 - 0.5 * (round(y) % 2) # shift 0.5,1 for even,odd col
            m[int(y + 0.5)][int(x + xshift)] = 1
            x += dx
            y += dy

    else:
        while (
            round(x) < np.size(m, 1)
            and round(y) < np.size(m, 0)
            and round(x) >= 0
            and round(y) >= 0
        ):
            yshift = 0.5 - 0.5 * (round(x) % 2) # shift 0.5,1 for even,odd col
            m[int(y + yshift)][int(x + 0.5)] = 1
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
    >>> grid.at_node['fracture_at_node']
    array([1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0], dtype=int8)

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
        for i in range(nfracs):

            (c, r) = _calc_fracture_starting_position((nr, nc))
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
        # Make an initial grid of all zeros. If user specified a model grid,
        # use that. Otherwise, use the given dimensions.
        nr = self._grid.number_of_node_rows
        nc = self._grid.number_of_node_columns
        m = self._grid.at_node["fracture_at_node"].reshape((nr, nc))

        # Add fractures to grid
        nfracs = (nr + nc) // frac_spacing
        for i in range(nfracs):

            (c, r) = _calc_fracture_starting_position((nr, nc))
            ang = _calc_fracture_orientation((c, r), (nr, nc))
            (dx, dy) = _calc_fracture_step_sizes(ang)
            x, y =
            dx /= 2.0
            dy /= 2.0

            _trace_fracture_through_grid_hex(m, (c, r), (dx, dy))

    def _calc_fracture_starting_position(self):
        """Choose a random starting position along one of the sides of the grid.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns in the grid

        Returns
        -------
        (x, y) : tuple of float
            Fracture starting coordinates (normalized? x and y coords)
        """
        #
        grid_height = np.amax(self.grid.y_of_node)
        grid_width = np.amax(self.grid.x_of_node)

        try:  # raster
            grid_dx, grid_dy = self.grid.spacing
        except TypeError:  # hex
            if orientation == 'horizontal':
                grid_dx = self.grid.spacing
                grid_dy = (3.0**0.5 / 2.0) * self.grid.spacing
            else:
                grid_dy = self.grid.spacing
                grid_dx = (3.0**0.5 / 2.0) * self.grid.spacing

        if np.random.rand() < grid_height / (grid_height + grid_width): # e or w
            grid_side = 0
        else: # n or s
            grid_side = 1
        coin_flip = np.random.randint(0, 1)
        grid_side += 2 * coin_flip  # 0 = e or n, 1 = w or s

        if (grid_side % 2) == 0:  # east or west
            x = (1 - grid_side // 2) * grid_width
            y = np.random.rand() * grid_height / grid_dy
        else:
            x = np.random.rand() * grid_width / grid_dx
            y = (1 - grid_side // 2) * (shape[0] - 1)
        print('starting on side ' + str(grid_side) + ' at ' + str((x,y)))
        return (x, y)
