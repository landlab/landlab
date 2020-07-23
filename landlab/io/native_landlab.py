#! /usr/bin/env python
"""Read data from a pickled Landlab grid file into a RasterModelGrid.

Read Landlab native
+++++++++++++++++++

.. autosummary::

    ~landlab.io.native_landlab.load_grid
    ~landlab.io.native_landlab.save_grid
"""

import os
import pickle

from landlab import ModelGrid


def save_grid(grid, path, clobber=False):
    """Save a grid and fields to a Landlab "native" format.

    This method uses pickle to save a grid as a pickle file.
    All fields will be saved, along with the grid.

    The recommended suffix for the save file is '.grid'. This will
    be added to your save if you don't include it.

    Caution: Pickling can be slow, and can produce very large files.
    Caution 2: Future updates to Landlab could potentially render old
    saves unloadable.

    Parameters
    ----------
    grid : object of subclass ModelGrid
        Grid object to save
    path : str
        Path to output file, either without suffix, or '.grid'
    clobber : bool (default False)
        Set to True to allow overwrites of existing files

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.io.native_landlab import save_grid
    >>> import tempfile
    >>> grid_out = RasterModelGrid((4, 5), xy_spacing=2.)
    >>> with tempfile.TemporaryDirectory() as tmpdirname:
    ...     fname = os.path.join(tmpdirname, 'testsavedgrid.grid')
    ...     save_grid(grid_out, fname, clobber=True)
    """
    if os.path.exists(path) and not clobber:
        raise ValueError("file exists")

    # test it's a grid
    assert issubclass(type(grid), ModelGrid)

    (base, ext) = os.path.splitext(path)
    if ext != ".grid":
        ext = ext + ".grid"
    path = base + ext

    with open(path, "wb") as file_like:
        pickle.dump(grid, file_like)


def load_grid(path):
    """Load a grid and its fields from a Landlab "native" format.

    This method uses pickle to load a saved grid.
    It assumes you saved using vmg.save() or save_grid, i.e., that the
    pickle file is a .grid file.

    Caution: Pickling can be slow, and can produce very large files.
    Caution 2: Future updates to Landlab could potentially render old
    saves unloadable.

    Parameters
    ----------
    path : str
        Path to output file, either without suffix, or '.grid'

    Examples
    --------
    >>> from landlab import VoronoiDelaunayGrid
    >>> from landlab.io.native_landlab import load_grid, save_grid
    >>> import numpy as np
    >>> import tempfile
    >>> x = np.random.rand(20)
    >>> y = np.random.rand(20)
    >>> grid_out = VoronoiDelaunayGrid(x, y)
    >>> with tempfile.TemporaryDirectory() as tmpdirname:
    ...     fname = os.path.join(tmpdirname, 'testsavedgrid.grid')
    ...     save_grid(grid_out, fname, clobber=True)
    ...     grid_in = load_grid(fname)
    """
    (base, ext) = os.path.splitext(path)
    if ext != ".grid":
        ext = ext + ".grid"
    path = base + ext
    with open(path, "rb") as file_like:
        loaded_grid = pickle.load(file_like)
    assert issubclass(type(loaded_grid), ModelGrid)
    return loaded_grid
