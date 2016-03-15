#! /usr/env/python

"""Create 2D grid with randomly generated fractures.

fracture_grid: creates and returns a 2D grid with randomly generated fractures.
The grid contains the value 1 where fractures (one cell wide) exist, and
0 elsewhere. The idea is to use this for simulations based on weathering and
erosion of, and/or flow within, fracture networks.

The entry point is the function::

    make_frac_grid(frac_spacing, numrows=50, numcols=50, model_grid=None)

If called with a Landlab RasterModelGrid, the function returns the fracture
grid as a node array. Otherwise, it returns a numrows x numcols Numpy array.

Potential improvements:
    - Add doctests
    - Fractures could be defined by links rather than nodes (i.e., return a
        link array with a code indicating whether the link crosses a fracture
        or not)
    - Fractures could have a finite length rather than extending all the way
        across the grid
    - Use of starting position along either x or y axis makes fracture net
        somewhat asymmetric. One would need a different algorithm to make it
        fully (statistically) symmetric.

Created: September 2013 by Greg Tucker
Last significant modification: cleanup and unit tests Oct 2015 GT
"""

from numpy import tan, pi, size, zeros
import numpy as np


def calculate_fracture_starting_position(shape, seed):
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


def calculate_fracture_orientation(coords, seed):
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
    If the fracture starts along the bottom of the grid (y=0), then the angle
    will be between 45 and 135 degrees from horizontal (counter-clockwise).
    Otherwise, it will be between -45 and 45 degrees.
    """
    y, x = coords

    np.random.seed(seed)
    ang = (pi/2)*np.random.rand()
    if y == 0:
        ang += pi/4
    else:
        ang -= pi/4

    return ang


def calculate_fracture_step_sizes(start_yx, ang):
    """
    Calculate the sizes of steps dx and dy to be used when "drawing" the
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
    if startx==0:  # frac starts on left side
        dx = 1
        dy = tan(ang)
    else:  # frac starts on bottom side
        dy = 1
        dx = -tan(ang-pi/2)

    return (dy, dx)


def trace_fracture_through_grid(m, start_yx, spacing):
    """Create a 2D fracture in a grid.

    Creates a "fracture" in a 2D grid, m, by setting cell values to unity along
    the trace of the fracture (i.e., "drawing" a line throuh the grid).

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
    
    while round(x) < size(m, 1) and round(y) < size(m, 0) \
            and round(x) >= 0 and round(y) >= 0:
        m[int(y + 0.5)][int(x + 0.5)] = 1
        x += dx
        y += dy


def make_frac_grid(frac_spacing, numrows=50, numcols=50, model_grid=None,
                   seed=0):
    """Create a grid that contains a network of random fractures.

    Creates and returns a grid containing a network of random fractures, which
    are represented as 1's embedded in a grid of 0's.

    Parameters
    ----------
    frac_spacing : int
        Average spacing of fractures (in grid cells)
    numrows : int, optional
        Number of rows in grid (if model_grid parameter is given,
        uses values from the model grid instead)
    numcols : int, optional
        Number of columns in grid (if model_grid parameter is given,
        uses values from the model grid instead)
    model_grid : Landlab RasterModelGrid object, optiona
        RasterModelGrid to use for grid size
    seed : int, optional
        Seed used for random number generator

    Returns
    -------
    m : Numpy array
        Array containing fracture grid, represented as 0's (matrix) and 1's
        (fractures). If model_grid parameter is given, returns a 1D array
        corresponding to a node-based array in the model grid. Otherwise,
        returns a 2D array with dimensions given by numrows, numcols.
    """
    # Make an initial grid of all zeros. If user specified a model grid,
    # use that. Otherwise, use the given dimensions.
    if model_grid is not None:
        numrows = model_grid.number_of_node_rows
        numcols = model_grid.number_of_node_columns
    m = zeros((numrows,numcols), dtype=int)

    # Add fractures to grid
    nfracs = (numrows + numcols) // frac_spacing
    for i in range(nfracs):

        (y, x) = calculate_fracture_starting_position((numrows, numcols), seed+i)
        ang = calculate_fracture_orientation((y, x), seed+i)
        (dy, dx) = calculate_fracture_step_sizes((y, x), ang)

        trace_fracture_through_grid(m, (y, x), (dy, dx))

    # If we have a model_grid, flatten the frac grid so it's equivalent to
    # a node array.
    if model_grid is not None:
        m.shape = (m.shape[0]*m.shape[1])

    return m
