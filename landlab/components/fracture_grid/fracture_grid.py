#! /usr/env/python

"""Create 2D grid with randomly generated fractures.

Created: September 2013 by Greg Tucker
Last significant modification: conversion to proper component 7/2019 GT
"""

import numpy as np

from landlab import Component


def _calc_fracture_starting_position(shape, seed):
    """Choose a random starting position along the x or y axis (random choice).

    Parameters
    ----------
    shape : tuple of int
        Number of rows and columns in the grid
    seed : int
        Seeds the random number generator, so that a particular random
        sequence can be recreated.

    Returns
    -------
    (y, x) : tuple of int
        Fracture starting coordinates
    """
    np.random.seed(seed)

    if np.random.randint(0, 1) == 0:
        x = 0
        y = np.random.randint(0, shape[0] - 1)
    else:
        x = np.random.randint(0, shape[1] - 1)
        y = 0
    return (y, x)


def _calc_fracture_orientation(coords, seed):
    """Choose a random orientation for the fracture.

    Parameters
    ----------
    coords : tuple of int
        Starting coordinates (one of which should be zero) as *y*, *x*.
    seed : int
        Seed value for random number generator

    Returns
    -------
    ang : float
        Fracture angle relative to horizontal

    Notes
    -----
    If the fracture starts along the bottom of the grid (y=0), then the
    angle will be between 45 and 135 degrees from horizontal
    (counter-clockwise). Otherwise, it will be between -45 and 45 degrees.
    """
    y, x = coords

    np.random.seed(seed)
    ang = (np.pi / 2) * np.random.rand()
    if y == 0:
        ang += np.pi / 4
    else:
        ang -= np.pi / 4

    return ang


def _calc_fracture_step_sizes(start_yx, ang):
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
    """
    starty, startx = start_yx
    if startx == 0:  # frac starts on left side
        dx = 1
        dy = np.tan(ang)
    else:  # frac starts on bottom side
        dy = 1
        dx = -np.tan(ang - np.pi / 2)

    return (dy, dx)


def _trace_fracture_through_grid(m, start_yx, spacing):
    """Create a 2D fracture in a grid.

    Creates a "fracture" in a 2D grid, m, by setting cell values to unity
    along the trace of the fracture (i.e., "drawing" a line throuh the
    grid).

    Parameters
    ----------
    m : 2D Numpy array
        Array that represents the grid
    start_yx : tuple of int
        Starting grid coordinates for fracture
    spacing : tuple of float
        Step sizes in y and x directions

    Returns
    -------
    None, but changes contents of m
    """
    y0, x0 = start_yx
    dy, dx = spacing

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

    def run_one_step(self):
        """Run FractureGridGenerator and create a random fracture grid."""
        self._make_frac_grid(self._frac_spacing, self._seed)

    def _make_frac_grid(self, frac_spacing, seed):
        """Create a grid that contains a network of random fractures.

        Creates a grid containing a network of random fractures, which are
        represented as 1's embedded in a grid of 0's. The grid is stored in
        the "fracture_at_node" field.

        Parameters
        ----------
        frac_spacing : int
            Average spacing of fractures (in grid cells)
        seed : int
            Seed used for random number generator
        """
        # Make an initial grid of all zeros. If user specified a model grid,
        # use that. Otherwise, use the given dimensions.
        nr = self._grid.number_of_node_rows
        nc = self._grid.number_of_node_columns
        m = self._grid.at_node["fracture_at_node"].reshape((nr, nc))

        # Add fractures to grid
        nfracs = (nr + nc) // frac_spacing
        for i in range(nfracs):

            (y, x) = _calc_fracture_starting_position((nr, nc), seed + i)
            ang = _calc_fracture_orientation((y, x), seed + i)
            (dy, dx) = _calc_fracture_step_sizes((y, x), ang)

            _trace_fracture_through_grid(m, (y, x), (dy, dx))
