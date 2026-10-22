# =============================================================================
#
#   CRISIS — Coupled Regional Rainfall-Induced and Seismic Slope Instability
#            Simulations
#
#   Module  : _crisis_core.py
#   Purpose : CRISIS core — physics and input/output functions used by the CRISIS
#             Landlab component (crisis_landlab_component.py)
#
#   Authors : Mirna Kassem
#             University of California, Berkeley
#             kassem_mirna@berkeley.edu
#
#             Dimitrios Zekkos
#             University of California, Berkeley
#             zekkos@berkeley.edu
#
#   This module has no driver of its own: CrisisComponent supplies the grid
#   arrays (elevation, slope, flow direction, drainage area) from Landlab
#   fields and calls these functions once per time step.
#
# =============================================================================


# =============================================================================
# Standard-library and third-party imports
# =============================================================================

import math
from contextlib import nullcontext
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import geopandas as gpd
from shapely.geometry import box as _shapely_box
from shapely.ops import unary_union as _unary_union


# =============================================================================
# Module-level constants
# =============================================================================

# Pre-allocated empty integer array used as a fast stand-in when a cell has
# no upslope neighbors, avoiding repeated np.array([]) allocations in the loop
_EMPTY_ARR = np.array([], dtype=np.int64)


# =============================================================================
# Utility helpers
# =============================================================================

def load_h5(file_path):
    """Load the 'data' dataset from an HDF5 file and return it as a float32 array.

    The array is flipped vertically (``np.flipud``) on load: the files are
    north-up (row 0 = northern edge), while CRISIS uses Landlab's row order
    (row 0 = southern edge).

    Parameters
    ----------
    file_path : str, Path or file-like
        Path to an HDF5 file containing a dataset named ``'data'``, or an
        in-memory HDF5 file (e.g. ``io.BytesIO``).

    Returns
    -------
    ndarray of float32
        The dataset, vertically flipped.

    Examples
    --------
    A small HDF5 file, created in memory (a file path works the same way):

    >>> import io, h5py, numpy as np
    >>> from landlab.components.Crisis._crisis_core import load_h5
    >>> h5 = io.BytesIO()
    >>> with h5py.File(h5, "w") as f:
    ...     _ = f.create_dataset("data", data=np.array([[1.0, 2.0],
    ...                                                  [3.0, 4.0]]))
    >>> load_h5(h5)
    array([[3., 4.],
           [1., 2.]], dtype=float32)
    """
    with h5py.File(file_path, 'r') as f:
        return np.flipud(f['data'][:].astype(np.float32))  # north-up file -> Landlab order (row 0 = south)

# -----------------------------------------------------------------------------------------------------
def pad_rows(data):
    """Pad a list-of-lists to uniform row length using ``None`` as filler.

    Required before writing ragged (variable-length) rows to Excel, since
    openpyxl needs rectangular data.

    Parameters
    ----------
    data : list of list
        Ragged rows of arbitrary (possibly unequal) length.

    Returns
    -------
    list of list
        Same rows, each padded with ``None`` at the end to match the length
        of the longest row. Returns ``[]`` unchanged if `data` is empty.

    Examples
    --------
    >>> from landlab.components.Crisis._crisis_core import pad_rows
    >>> pad_rows([[1, 2, 3], [4, 5]])
    [[1, 2, 3], [4, 5, None]]
    >>> pad_rows([])
    []
    """
    if data:
        max_len = max(len(row) for row in data)
        return [list(row) + [None] * (max_len - len(row)) for row in data]
    return []

# -----------------------------------------------------------------------------------------------------
def read_ascii_header(file_path):
    """Read the header of an Esri ASCII grid (``.asc``) DEM.

    Header lines are ``key value`` pairs at the top of the file; reading
    stops at the first line that is not one (the start of the data block).
    Keys are returned in lower case, with ``ncols`` and ``nrows`` as int and
    all other values as float.

    Parameters
    ----------
    file_path : str, Path or file-like
        Path to an Esri ASCII grid file, or an open text stream (e.g.
        ``io.StringIO``) holding its content.

    Returns
    -------
    dict
        Header values. Always contains ``ncols``, ``nrows``, ``cellsize``,
        one of ``xllcorner`` / ``xllcenter``, and one of ``yllcorner`` /
        ``yllcenter``; may also contain ``nodata_value``.

    Raises
    ------
    ValueError
        If a required key is missing, or both the corner and center form
        of the lower-left coordinate are given.

    Examples
    --------
    The content of a small DEM, given in memory (a file path works the same way):

    >>> import io
    >>> from landlab.components.Crisis._crisis_core import read_ascii_header
    >>> dem = io.StringIO("ncols 3\\nnrows 2\\nxllcorner 1000.0\\nyllcorner 5000.0\\n"
    ...                   "cellsize 10.0\\nNODATA_value -9999\\n1 2 3\\n4 5 6\\n")
    >>> header = read_ascii_header(dem)
    >>> header["ncols"], header["nrows"], header["xllcorner"], header["cellsize"]
    (3, 2, 1000.0, 10.0)
    >>> header["nodata_value"]
    -9999.0
    """
    header = {}
    # a path is opened here; an already open text stream is read from its start
    if hasattr(file_path, "read"):
        file_path.seek(0)
    with (open(file_path, "r") if not hasattr(file_path, "read") else nullcontext(file_path)) as f:
        for line in f:
            parts = line.split()
            # header lines are exactly "key value" with a non-numeric key
            if len(parts) != 2:
                break
            try:
                float(parts[0])
                break          # numeric first token: the data block has started
            except ValueError:
                pass
            header[parts[0].lower()] = float(parts[1])

    # check the required keys
    for key in ("ncols", "nrows", "cellsize"):
        if key not in header:
            raise ValueError(f"ASCII header of '{file_path}' is missing '{key}'.")
    for corner, center in (("xllcorner", "xllcenter"), ("yllcorner", "yllcenter")):
        if (corner in header) == (center in header):
            raise ValueError(
                f"ASCII header of '{file_path}' must contain exactly one of "
                f"'{corner}' or '{center}'."
            )

    header["ncols"] = int(header["ncols"])
    header["nrows"] = int(header["nrows"])
    return header

# -----------------------------------------------------------------------------------------------------
def lower_left_cell_center(header):
    """Return the (x, y) coordinates of the center of the lower-left cell.

    An Esri ASCII header gives the lower-left position either as the outer
    corner of the lower-left cell (``xllcorner`` / ``yllcorner``) or as its
    center (``xllcenter`` / ``yllcenter``). A corner is shifted by half a
    cell to obtain the center.

    Parameters
    ----------
    header : dict
        Header as returned by `read_ascii_header`.

    Returns
    -------
    xmin, ymin : float
        Easting / northing of the center of the lower-left cell.

    Examples
    --------
    >>> from landlab.components.Crisis._crisis_core import lower_left_cell_center
    >>> lower_left_cell_center({"xllcorner": 1000.0, "yllcorner": 5000.0, "cellsize": 10.0})
    (1005.0, 5005.0)
    >>> lower_left_cell_center({"xllcenter": 1005.0, "yllcenter": 5005.0, "cellsize": 10.0})
    (1005.0, 5005.0)
    """
    half = 0.5 * header["cellsize"]
    xmin = header["xllcenter"] if "xllcenter" in header else header["xllcorner"] + half
    ymin = header["yllcenter"] if "yllcenter" in header else header["yllcorner"] + half
    return float(xmin), float(ymin)

# -----------------------------------------------------------------------------------------------------
def compute_coordinate_vectors(rows, columns, xmin, ymin, Cell_size):
    """Build the x (easting) and y (northing) cell-center coordinate vectors.

    The grid is fully determined by its size, its cell size, and the center
    of its lower-left cell, so no separate coordinate files are needed.

    Convention (Landlab orientation, as used by the CRISIS Landlab
    component): row 0 is the southern-most (bottom) row and column 0 the
    western-most column, so (xmin, ymin) is the center of cell (row 0,
    col 0). x increases eastward with the column index; y increases
    northward with the row index.

    Parameters
    ----------
    rows, columns : int
        Grid dimensions.
    xmin, ymin : float
        Easting / northing of the center of the lower-left cell (m).
    Cell_size : float
        Grid cell size (m).

    Returns
    -------
    x : ndarray of float64, shape (columns,)
        Easting of the center of each column.
    y : ndarray of float64, shape (rows,)
        Northing of the center of each row (row 0 = southern-most row).

    Examples
    --------
    >>> from landlab.components.Crisis._crisis_core import compute_coordinate_vectors
    >>> x, y = compute_coordinate_vectors(2, 3, 1005.0, 5005.0, 10.0)
    >>> x
    array([1005., 1015., 1025.])
    >>> y
    array([5005., 5015.])
    """
    x = xmin + np.arange(columns, dtype=np.float64) * Cell_size
    y = ymin + np.arange(rows, dtype=np.float64) * Cell_size
    return x, y


# =============================================================================
# Function 1 — Upslope grid cells (D8 flow routing)
# =============================================================================

def T_UpGrid(FD_Mat, indID):
    """
    Find all upslope (inflow) cells for a target cell using the D8 convention.

    Parameters
    ----------
    FD_Mat : ndarray (m x n)   Flow direction matrix.
    indID  : int                Linear index of the target cell.

    Returns
    -------
    UpGrid_Cells : ndarray      Linear indices of upslope cells (length <= 8).

    Examples
    --------
    3x3 grid where all four orthogonal neighbors of the center cell (index 4)
    drain into it (codes 4 = row +1, 1 = col +1, 16 = col -1, 64 = row -1,
    each pointing toward the center):

    >>> import numpy as np
    >>> from landlab.components.Crisis._crisis_core import T_UpGrid
    >>> FD_Mat = np.array([
    ...     [0,  4,  0],
    ...     [1,  0, 16],
    ...     [0, 64,  0],
    ... ])
    >>> T_UpGrid(FD_Mat, 4)
    array([5, 7, 3, 1])
    """
    # Get the grid dimensions and convert the linear index to row/col coordinates
    nrow, ncol = FD_Mat.shape
    T_row, T_col = np.unravel_index(indID, (nrow, ncol))

    # D8 direction codes (ArcGIS numbering 1-128) mapped to (row, col) offsets.
    # Each key is a D8 code; value is the (row, col) offset of the neighbor
    # that code points to. The codes are used purely as offset labels.
    directions = {
        1:   (0,  1),
        2:   (1,  1),
        4:   (1,  0),
        8:   (1, -1),
        16:  (0, -1),
        32:  (-1, -1),
        64:  (-1, 0),
        128: (-1, 1),
    }
    # Opposite direction codes: a cell flows INTO the target if its flow
    # direction points back toward the target cell
    opposite = {1: 16, 2: 32, 4: 64, 8: 128, 16: 1, 32: 2, 64: 4, 128: 8}

    # Loop over all 8 neighbors and check whether they drain into the target cell
    UpGrid_Cells = []
    for code, (dr, dc) in directions.items():
        r, c = T_row + dr, T_col + dc
        # Only accept neighbors that are within grid bounds and whose flow
        # direction points back to the target cell
        if 0 <= r < nrow and 0 <= c < ncol and FD_Mat[r, c] == opposite[code]:
            UpGrid_Cells.append(np.ravel_multi_index((r, c), (nrow, ncol)))

    return np.array(UpGrid_Cells, dtype=int)


def _build_upslope_map(FD_Mat):
    """
    Precompute upslope (inflow) neighbors for EVERY cell in one vectorized pass.
    This replaces repeated calls to T_UpGrid during landslide delineation,
    reducing the complexity from O(cells × BFS steps) to O(cells) pre-processing.

    Parameters
    ----------
    FD_Mat : ndarray (m x n)
        D8 flow direction matrix (codes 1-128, see `T_UpGrid`).

    Returns
    -------
    dict : linear_index (int) -> ndarray of upslope linear indices (int64)

    Examples
    --------
    Same 3x3 grid as :func:`T_UpGrid`, where every orthogonal neighbor of
    the center cell (index 4) drains into it:

    >>> import numpy as np
    >>> from landlab.components.Crisis._crisis_core import _build_upslope_map
    >>> FD_Mat = np.array([
    ...     [0,  4,  0],
    ...     [1,  0, 16],
    ...     [0, 64,  0],
    ... ])
    >>> upslope = _build_upslope_map(FD_Mat)
    >>> upslope[4].tolist()
    [5, 7, 3, 1]
    """
    nrow, ncol = FD_Mat.shape

    # D8 direction codes and their (row, col) offsets
    directions = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1),
                  16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}
    # Opposite direction lookup for determining inflow relationships
    opposite   = {1: 16, 2: 32, 4: 64, 8: 128, 16: 1, 32: 2, 64: 4, 128: 8}

    # Accumulators for all (source → target) upslope pairs across the whole grid
    all_src_list, all_tgt_list = [], []

    for code, (dr, dc) in directions.items():
        opp = opposite[code]

        # Find all cells whose flow direction equals the opposite code:
        # these cells flow in the direction opposite to 'code', i.e. they
        # drain into the cell located (-dr, -dc) steps away (the target below)
        r_src, c_src = np.where(FD_Mat == opp)
        if len(r_src) == 0:
            continue

        # Compute the coordinates of the target (downstream) cell for each source
        r_tgt = r_src - dr
        c_tgt = c_src - dc

        # Keep only pairs where the target cell is still within grid bounds
        valid = (r_tgt >= 0) & (r_tgt < nrow) & (c_tgt >= 0) & (c_tgt < ncol)
        if not np.any(valid):
            continue

        # Convert valid (row, col) pairs to flat linear indices and store them
        all_src_list.append(np.ravel_multi_index((r_src[valid], c_src[valid]), (nrow, ncol)))
        all_tgt_list.append(np.ravel_multi_index((r_tgt[valid], c_tgt[valid]), (nrow, ncol)))

    # If no upslope relationships were found (e.g. flat DEM), return empty dict
    if not all_src_list:
        return {}

    # Concatenate all pairs into flat arrays
    all_src = np.concatenate(all_src_list)
    all_tgt = np.concatenate(all_tgt_list)

    # Sort pairs by target index so we can group upslope sources per target cell
    sort_idx = np.argsort(all_tgt, kind='stable')
    all_src  = all_src[sort_idx]
    all_tgt  = all_tgt[sort_idx]

    # Group all source cells by their common target cell using np.unique + np.split
    uniq_tgt, counts = np.unique(all_tgt, return_counts=True)
    splits = np.split(all_src, np.cumsum(counts[:-1]))

    # Return a dictionary: target_cell_index → array of upslope source indices
    return dict(zip(uniq_tgt.tolist(), splits))


# =============================================================================
# Function 2 — Pseudo-3D landslide delineation
# =============================================================================

def Determine_Pseudo_3D_landslides(Elev_Mat, Trig_Cells, Area_Mat, Slope_Mat,
                                   FD_Mat, Cell_size, Trig_Depths, Trig_sources,
                                   Fail_Index_all_times, FID_all_times, Min_LD_area,
                                   upslope_map=None):
    """
    Form pseudo-3D landslide bodies from triggering cells using upslope BFS expansion.
    Based on Gong et al. (2023) and Gallen et al. (2015) with modifications.

    Starting from each triggering cell, the algorithm walks upslope through the
    flow network and accumulates all cells that satisfy the geometric failure
    criterion (positive failure depth). The resulting set of cells defines one
    landslide body. Triggering cells with no contributing area, or that already
    failed, are skipped; landslides with 2 cells or fewer, or smaller than
    `Min_LD_area`, are rejected.

    Parameters
    ----------
    Elev_Mat : ndarray (rows, clms)
        Surface elevation (m).
    Trig_Cells : ndarray of int
        Linear indices of candidate triggering cells for this time step
        (e.g. from `_find_triggering_cells_vectorized`).
    Area_Mat : ndarray (rows, clms)
        Upstream contributing area (m^2); used to process large-catchment
        triggers before nested smaller ones.
    Slope_Mat : ndarray (rows, clms)
        Slope angle (degrees).
    FD_Mat : ndarray (rows, clms)
        D8 flow direction matrix (codes 1-128, see `T_UpGrid`). Only used
        when `upslope_map` is not supplied.
    Cell_size : float
        Grid cell size (m).
    Trig_Depths : array-like of float
        Failure depth (m) at each entry of `Trig_Cells`.
    Trig_sources : array-like of int
        Trigger source code at each entry of `Trig_Cells`
        (1=rain, 2=rain+eq, 3=eq only).
    Fail_Index_all_times : ndarray of bool (rows, clms)
        Cells that have already failed in any previous time step (mutated
        in place: cells belonging to newly accepted landslides are set to
        `True`).
    FID_all_times : ndarray of float (rows, clms)
        Cumulative landslide-ID raster (mutated in place: cells belonging
        to newly accepted landslides are stamped with their new ID).
    Min_LD_area : float
        Minimum accepted landslide planar area (m^2); landslides with fewer
        than `ceil(Min_LD_area / Cell_size**2)` cells are rejected.
    upslope_map : dict or None, optional
        Precomputed output of `_build_upslope_map(FD_Mat)`. When `None`,
        upslope neighbors are computed on the fly with `T_UpGrid` (slower).

    Returns
    -------
    Vol : ndarray of float32
        Volume (m^3) of each accepted landslide.
    Area : ndarray of float32
        Planar area (m^2) of each accepted landslide.
    NL : int
        Number of landslides accepted at this time step.
    Trig_Index_Accept : ndarray of int64
        Linear index of the triggering cell of each accepted landslide.
    Slope_Trig_Index_Accept : ndarray of float32
        Slope angle (degrees) at each accepted landslide's triggering cell.
    Fail_Index_all_times : ndarray of bool
        Updated failure mask (same array passed in, mutated in place).
    Trig_Depth_Accept : list of float
        Failure depth at each accepted landslide's triggering cell.
    Trig_source_Accept : list of int
        Trigger source code at each accepted landslide's triggering cell.
    FID_all_times : ndarray of float32
        Updated cumulative landslide-ID raster (same array passed in,
        mutated in place).
    LD_Z_Accept : ndarray of float
        Failure depth of every individual cell across all accepted
        landslides at this time step (concatenated, not grouped by
        landslide). An empty list if no landslide is accepted.
    LD_Index_Accept : ndarray of int64
        Linear index of every individual cell across all accepted
        landslides at this time step (same order/length as `LD_Z_Accept`).
        An empty list if no landslide is accepted.
    FID_Index : ndarray (same shape as `Area_Mat`)
        Landslide ID *local to this time step* (1..NL; 0 elsewhere), as
        opposed to `FID_all_times` which carries IDs cumulative across the
        whole run.

    Examples
    --------
    A 1x5 upslope chain (column 4 = headwater, column 0 = outlet), all
    flowing West, with a 45-degree slope and one triggering cell at the
    outlet (column 0) with a 1 m failure depth. The BFS walks upslope and
    accepts columns 0-3; column 4 is rejected because refining the failure
    plane geometry along the way makes its computed depth negative.

    >>> import numpy as np
    >>> from landlab.components.Crisis._crisis_core import Determine_Pseudo_3D_landslides
    >>> Elev_Mat  = np.array([[0.0, 2.0, 4.0, 6.0, 8.0]])
    >>> Slope_Mat = np.full((1, 5), 45.0)
    >>> Area_Mat  = np.array([[500.0, 400.0, 300.0, 200.0, 100.0]])
    >>> FD_Mat    = np.full((1, 5), 16, dtype=int)   # every cell flows West
    >>> Trig_Cells   = np.array([0])
    >>> Trig_Depths  = np.array([1.0])
    >>> Trig_sources = np.array([1])
    >>> Fail_Index_all_times = np.zeros((1, 5), dtype=bool)
    >>> FID_all_times = np.full((1, 5), np.nan, dtype=np.float32)
    >>> result = Determine_Pseudo_3D_landslides(
    ...     Elev_Mat, Trig_Cells, Area_Mat, Slope_Mat, FD_Mat, Cell_size=2.0,
    ...     Trig_Depths=Trig_Depths, Trig_sources=Trig_sources,
    ...     Fail_Index_all_times=Fail_Index_all_times, FID_all_times=FID_all_times,
    ...     Min_LD_area=8.0, upslope_map=None)
    >>> Vol, Area, NL, Trig_Index_Accept = result[0], result[1], result[2], result[3]
    >>> Vol
    array([11.903158], dtype=float32)
    >>> Area
    array([16.], dtype=float32)
    >>> NL
    1
    >>> result[5]  # updated Fail_Index_all_times
    array([[ True,  True,  True,  True, False]])
    >>> result[8]  # updated FID_all_times: columns 0-3 belong to landslide 1
    array([[ 1.,  1.,  1.,  1., nan]], dtype=float32)
    """
    # Retrieve the highest landslide ID assigned in previous time steps so new
    # IDs continue the sequence rather than restarting from 1
    NL_previous = int(np.nanmax(FID_all_times)) if not np.isnan(FID_all_times).all() else 0
    # Counter for landslides accepted at this time step
    NL = 0

    # Sort triggering cells by drainage area in descending order.
    # Processing large-drainage-area cells first ensures that their upslope
    # domain is claimed before smaller nested cells can steal those cells.
    # The stable sort keeps cells with equal drainage area in their original
    # order (lower linear index first), so results do not depend on the NumPy version.
    Idx = np.argsort(-np.take(Area_Mat, Trig_Cells), kind="stable")
    Trig_Cells_Desc = np.take(Trig_Cells, Idx)
    NTC = len(Trig_Cells_Desc)

    # Output accumulators for this time step
    Vol, Area_out = [], []
    Trig_Index_Accept, Slope_Trig_Index_Accept = [], []
    Trig_Depth_Accept, Trig_source_Accept = [], []

    m, n = Elev_Mat.shape
    # Minimum number of grid cells a landslide must contain to be accepted
    min_required_cells = math.ceil(Min_LD_area / Cell_size**2)

    LD_Index_Accept = []
    LD_Z_Accept = []
    FID_Index  =np.zeros_like(Area_Mat)
    # ── Process each triggering cell ──────────────────────────────────────────
    for j in range(NTC):
        Trig_Index = Trig_Cells_Desc[j]

        # Recover original position in the unsorted array to get depth and source
        orig_pos = np.where(Trig_Cells == Trig_Index)[0][0]
        Trig_Depth_cell  = Trig_Depths[orig_pos]
        Trig_source_cell = Trig_sources[orig_pos]

        # Skip this cell if it has zero contributing area (sink/flat) or if it
        # has already been assigned to a previously delineated landslide body
        if np.take(Area_Mat, Trig_Index) <= 0 or np.take(Fail_Index_all_times, Trig_Index) != 0:
            continue

        # Extract topographic properties at the triggering cell
        Elev_Trig  = np.take(Elev_Mat, Trig_Index)   # surface elevation (m)
        Trig_Slope = np.take(Slope_Mat, Trig_Index)   # slope angle (degrees)
        # Elevation of the failure plane at the triggering cell
        Hf_Trig    = Elev_Trig - Trig_Depth_cell
        Trig_row, Trig_col = np.unravel_index(Trig_Index, (m, n))

        # Initialize the BFS (breadth-first search) queue with the triggering cell.
        # A set is used alongside the list for O(1) duplicate checking.
        LD_set  = {int(Trig_Index)}
        LD_list = [int(Trig_Index)]
        Counter = 0


        # ── BFS: expand landslide body upslope ────────────────────────────────
        while Counter < len(LD_list):
            index = LD_list[Counter]

            # Retrieve upslope neighbors from the pre-computed map (O(1) lookup)
            # or compute them on-the-fly with T_UpGrid if the map is unavailable
            if upslope_map is not None:
                UpGrid_Cells = upslope_map.get(index, _EMPTY_ARR)
            else:
                UpGrid_Cells = T_UpGrid(FD_Mat, index)

            # Condition 1: exclude neighbor cells that already belong to a
            # previously delineated landslide (they cannot be re-assigned)
            if len(UpGrid_Cells) > 0:
                UpGrid_Cells = UpGrid_Cells[
                    Fail_Index_all_times.flat[UpGrid_Cells] == 0
                ]
            # Condition 2: exclude neighbor cells whose contributing area is
            # larger than the current cell — upslope cells should always drain
            # through the current cell, so larger area indicates a different branch
            if len(UpGrid_Cells) > 0:
                Acc_n = np.take(Area_Mat, index)
                UpGrid_Cells = UpGrid_Cells[np.take(Area_Mat, UpGrid_Cells) <= Acc_n]
            # Condition 3: exclude neighbor cells that are lower in elevation
            # than the triggering cell — the landslide body must extend uphill
            if len(UpGrid_Cells) > 0:
                UpGrid_Cells = UpGrid_Cells[np.take(Elev_Mat, UpGrid_Cells) >= Elev_Trig]
            # Condition 4: exclude neighbor cells where the projected failure
            # plane depth would be zero or negative (no room for a slide body)
            if len(UpGrid_Cells) > 0:
                UpGrid_row, UpGrid_col = np.unravel_index(UpGrid_Cells, (m, n))
                # Horizontal distance from the triggering cell to each candidate neighbor
                UpGrid_dist = np.sqrt(
                    (UpGrid_row - Trig_row)**2 + (UpGrid_col - Trig_col)**2
                ) * Cell_size
                # Failure plane elevation at each neighbor, projected from the
                # triggering cell using its slope angle
                Hf = Hf_Trig + UpGrid_dist * np.tan(np.radians(Trig_Slope))
                # Failure depth = surface elevation minus failure plane elevation
                UpGrid_Z = np.take(Elev_Mat, UpGrid_Cells) - Hf
                UpGrid_Cells = UpGrid_Cells[UpGrid_Z > 0]

            # Add surviving neighbors to the landslide body (no duplicates)
            for cell in UpGrid_Cells.tolist():
                if cell not in LD_set:
                    LD_set.add(cell)
                    LD_list.append(cell)

            Counter += 1

        # Convert the landslide cell list to a numpy array for vectorized ops
        LD_Index = np.array(LD_list, dtype=np.int64)

        # Reject landslides with only 1 or 2 cells (geometrically degenerate)
        if len(LD_Index) <= 2:
            continue

        # ── Refine failure plane geometry ─────────────────────────────────────
        # Re-estimate the slope of the failure surface using the maximum surface
        # elevation within the landslide body, then recompute failure depths

        Surf_Elev = np.take(Elev_Mat, LD_Index)
        LD_row, LD_col = np.unravel_index(LD_Index, (m, n))
        # Horizontal distance from the triggering cell to every landslide cell
        LD_dist = np.sqrt((LD_row - Trig_row)**2 + (LD_col - Trig_col)**2) * Cell_size

        # Find the highest elevation in the landslide and its distance from the trigger
        Max_elev    = np.max(Surf_Elev)
        Min_LD_dist = np.min(LD_dist[Surf_Elev == Max_elev])
        Elev_diff   = Max_elev - Hf_Trig
        if Min_LD_dist > 0:
            # Refined mean slope angle from triggering cell to the highest point
            M_slope  = math.ceil(math.degrees(math.atan(Elev_diff / Min_LD_dist)))
            # Rate of slope change along the landslide (slope gradient)
            dx_slope = (M_slope - Trig_Slope) / Min_LD_dist
        else:
            # the highest point is the triggering cell itself: keep the constant slope
            dx_slope = 0
        # Set negative slope changes to zero
        if dx_slope < 0:
            dx_slope = 0

        # Recompute failure plane depths using the refined slope profile
        Hf    = Hf_Trig + LD_dist * np.tan(np.radians(Trig_Slope + LD_dist * dx_slope))
        LD_Z  = Surf_Elev - Hf
        # Keep only cells where the refined failure depth is still positive
        LD_Index = np.array(LD_Index)[LD_Z > 0]
        LD_Z     = LD_Z[LD_Z > 0]

        # Reject landslides that fall below the minimum area threshold after refinement
        if len(LD_Index) < min_required_cells:
            continue

        # ── Accept this landslide ─────────────────────────────────────────────
        if len(LD_Index_Accept)>0:
            LD_Index_Accept = np.concatenate((LD_Index_Accept,LD_Index))
            LD_Z_Accept = np.concatenate((LD_Z_Accept, LD_Z))
        else:
            LD_Index_Accept = np.copy(LD_Index)
            LD_Z_Accept = np.copy(LD_Z)

        NL += 1
        Trig_Index_Accept.append(Trig_Index)
        Slope_Trig_Index_Accept.append(Trig_Slope)
        Trig_Depth_Accept.append(Trig_Depth_cell)
        Trig_source_Accept.append(Trig_source_cell)
        # Landslide planar area = number of cells × cell area
        Area_out.append(len(LD_Index) * Cell_size**2)
        # Landslide volume = sum of failure depths × cell area (prism volumes)
        Vol.append(np.sum(LD_Z * Cell_size**2))

        # Mark all cells of this landslide in the global failure and FID arrays
        FID_all_times.flat[LD_Index]       = NL_previous + NL
        FID_Index.flat[LD_Index]                = NL
        Fail_Index_all_times.flat[LD_Index] = 1

    # Return results as numpy arrays (trim to accepted count NL)
    return (
        np.array(Vol[:NL],                       dtype=np.float32),
        np.array(Area_out[:NL],                  dtype=np.float32),
        NL,
        np.array(Trig_Index_Accept[:NL],         dtype=np.int64),
        np.array(Slope_Trig_Index_Accept[:NL],   dtype=np.float32),
        Fail_Index_all_times,
        Trig_Depth_Accept,
        Trig_source_Accept,
        FID_all_times,
        LD_Z_Accept,
        LD_Index_Accept,
        FID_Index
    )


# =============================================================================
# Function 3 — Closest-distance to fault rupture plane
# =============================================================================

def Rupture_Dist(x_2d, y_2d, z_2d, corners):
    """
    Compute closest 3-D distance (km) from every DEM cell to a planar
    quadrilateral fault rupture plane.

    The algorithm projects each site onto the fault plane, then returns
    the perpendicular distance to the plane when the projection falls inside
    the quadrilateral, or the shortest distance to the nearest edge otherwise.

    Parameters
    ----------
    x_2d, y_2d : ndarray (rows, clms)
        Projected easting / northing of each DEM cell (meters, same CRS
        as the fault corners).
    z_2d : ndarray (rows, clms)
        Surface elevation of each DEM cell (meters, positive upward).
    corners : array-like, shape (4, 3)
        [X, Y, Z] of the four fault-plane corners in meters.
        X, Y: easting / northing in the same CRS as the DEM.
        Z: depth positive downward (surface = 0, e.g. 15 000 for 15 km depth).
        Winding order: top-left, top-right, bottom-right, bottom-left
        (top = shallow edge, bottom = deep edge).

    Returns
    -------
    R_2D : ndarray (rows, clms), float64
        Closest-distance-to-rupture in km.

    Examples
    --------
    A vertical fault plane along x=0, striking north-south from y=0 to
    y=1000 m, extending from the surface (0 m) to 1000 m depth. One site is
    500 m east of the fault (inside the quad, so the answer is the
    perpendicular distance); another is 300 m west of it.

    >>> import numpy as np
    >>> from landlab.components.Crisis._crisis_core import Rupture_Dist
    >>> x_2d = np.array([[500.0, -300.0]])
    >>> y_2d = np.array([[500.0, 500.0]])
    >>> z_2d = np.array([[0.0, 0.0]])
    >>> corners = [
    ...     [0.0,    0.0,    0.0],    # top-left    (shallow)
    ...     [0.0, 1000.0,    0.0],    # top-right   (shallow)
    ...     [0.0, 1000.0, 1000.0],    # bottom-right (deep)
    ...     [0.0,    0.0, 1000.0],    # bottom-left  (deep)
    ... ]
    >>> Rupture_Dist(x_2d, y_2d, z_2d, corners)
    array([[0.5, 0.3]])
    """
    # Convert corners to a numpy array and preserve the original grid shape
    corners = np.asarray(corners, dtype=np.float64)
    shape   = x_2d.shape

    # The user supplies fault depth as positive-downward; negate the Z column
    # so the fault sits in the same right-handed coordinate system as the DEM
    # (where elevation is positive upward)
    corners = corners.copy()
    corners[:, 2] = -corners[:, 2]

    # Flatten all DEM cells to an (N, 3) array of 3-D site coordinates
    N = x_2d.size
    P = np.stack([x_2d.ravel(), y_2d.ravel(), z_2d.ravel()], axis=1)   # (N, 3)

    # ── Define the fault plane via its unit normal vector ─────────────────────
    # Two edge vectors spanning the fault plane from corner 0
    v1    = corners[1] - corners[0]
    v2    = corners[3] - corners[0]
    # Cross product gives the plane normal; normalize to unit length
    n     = np.cross(v1, v2)
    n_len = np.linalg.norm(n)
    if n_len < 1e-10:
        raise ValueError("Fault corners are collinear — cannot define a plane.")
    n_hat = n / n_len   # unit normal vector

    # ── Compute perpendicular distance from each site to the fault plane ──────
    # Vector from fault corner 0 to each site
    AP             = P - corners[0]
    # Signed scalar projection onto the unit normal = perpendicular distance
    dist_plane     = AP @ n_hat
    abs_dist_plane = np.abs(dist_plane)
    # Orthogonal projection of each site onto the fault plane
    P_proj         = P - dist_plane[:, np.newaxis] * n_hat

    # ── Test whether each projected point lies inside the fault quadrilateral ─
    # For a convex quad, a point is inside if it is on the same side of all edges.
    # edge_sign returns the signed cross-product component along the normal.
    def edge_sign(Q, A, B):
        e  = B - A          # edge direction vector
        qA = Q - A          # vectors from edge start to each query point
        c  = np.cross(e, qA)
        return c @ n_hat    # scalar: positive = left of edge, negative = right

    # Evaluate sign for all four edges of the quadrilateral
    s0 = edge_sign(P_proj, corners[0], corners[1])
    s1 = edge_sign(P_proj, corners[1], corners[2])
    s2 = edge_sign(P_proj, corners[2], corners[3])
    s3 = edge_sign(P_proj, corners[3], corners[0])

    # A point is inside if all four signs agree (all positive OR all negative)
    inside = (
        ((s0 >= 0) & (s1 >= 0) & (s2 >= 0) & (s3 >= 0)) |
        ((s0 <= 0) & (s1 <= 0) & (s2 <= 0) & (s3 <= 0))
    )

    # ── Compute distance to each edge (used for points outside the quad) ──────
    def dist_to_segment(Q, A, B):
        """Shortest 3-D distance from each point in Q (N×3) to segment A–B."""
        AB  = B - A
        AB2 = AB @ AB
        # Scalar projection onto the segment direction, clamped to [0, 1]
        t   = np.clip((Q - A) @ AB / (AB2 + 1e-30), 0.0, 1.0)
        # Closest point on segment for each query point
        C   = A + t[:, np.newaxis] * AB
        return np.linalg.norm(Q - C, axis=1)

    # Find the minimum distance across all four edges for each site
    d_edge = np.full(N, np.inf)
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 0)]:
        d_edge = np.minimum(d_edge, dist_to_segment(P, corners[i], corners[j]))

    # ── Choose the correct distance ───────────────────────────────────────────
    # For inside points: perpendicular (plane) distance
    # For outside points: distance to the nearest edge
    dist_m = np.where(inside, abs_dist_plane, d_edge)

    # Convert from meters to kilometers and reshape back to the DEM grid
    return (dist_m / 1000.0).reshape(shape)


# =============================================================================
# Function 4 — Spectral acceleration interpolation (Abrahamson & Silva 2008)
# =============================================================================

# Look-up table from Abrahamson & Silva (2008), Table 5.
# Each row defines the GMM coefficients for a period bin [T1, T2).
# Used by the Bray & Macedo (2019) displacement model to compute Sa at 1.3·Ts.
# Columns: T1, T2, then the coefficients at period T1 (a1_F, a2_F, a8_F)
# and at period T2 (a1_S, a2_S, a8_S)
_SA_TABLE = [
    (0.000, 0.010,  0.804, -0.9679, -0.0372,  0.811, -0.9676, -0.0372),
    (0.010, 0.020,  0.811, -0.9676, -0.0372,  0.855, -0.9774, -0.0372),
    (0.020, 0.030,  0.855, -0.9774, -0.0372,  0.962, -1.0024, -0.0372),
    (0.030, 0.040,  0.962, -1.0024, -0.0372,  1.037, -1.0289, -0.0315),
    (0.040, 0.050,  1.037, -1.0289, -0.0315,  1.133, -1.0508, -0.0271),
    (0.050, 0.075,  1.133, -1.0508, -0.0271,  1.375, -1.0810, -0.0191),
    (0.075, 0.100,  1.375, -1.0810, -0.0191,  1.563, -1.0833, -0.0166),
    (0.100, 0.150,  1.563, -1.0833, -0.0166,  1.716, -1.0375, -0.0254),
    (0.150, 0.200,  1.716, -1.0375, -0.0254,  1.687, -0.9700, -0.0396),
    (0.200, 0.250,  1.687, -0.9700, -0.0396,  1.646, -0.9202, -0.0539),
    (0.250, 0.300,  1.646, -0.9202, -0.0539,  1.601, -0.8974, -0.0656),
    (0.300, 0.400,  1.601, -0.8974, -0.0656,  1.511, -0.8677, -0.0807),
    (0.400, 0.500,  1.511, -0.8677, -0.0807,  1.397, -0.8475, -0.0924),
    (0.500, 0.750,  1.397, -0.8475, -0.0924,  1.137, -0.8206, -0.1137),
    (0.750, 1.000,  1.137, -0.8206, -0.1137,  0.915, -0.8088, -0.1289),
    (1.000, 1.500,  0.915, -0.8088, -0.1289,  0.510, -0.7995, -0.1534),
    (1.500, 2.000,  0.510, -0.7995, -0.1534,  0.192, -0.7960, -0.1708),
    (2.000, 3.000,  0.192, -0.7960, -0.1708, -0.280, -0.7960, -0.1954),
    (3.000, 4.000, -0.280, -0.7960, -0.1954, -0.639, -0.7960, -0.2128),
    (4.000, 5.000, -0.639, -0.7960, -0.2128, -0.936, -0.7960, -0.2263),
]


def _compute_Sa_vec(Ts_13, M, R_2D):
    """Compute the spectral acceleration Sa (g) at period Ts_13 for every DEM cell.

    Interpolates within the Abrahamson & Silva (2008) Table 5 using the
    magnitude- and distance-dependent GMM equations, then linearly interpolates
    between period bin boundaries (vectorized over the whole grid).

    Parameters
    ----------
    Ts_13 : float
        Target period = 1.3 * Ts (s).
    M : float
        Moment magnitude.
    R_2D : ndarray
        2-D array of closest rupture distances (km).

    Returns
    -------
    Sa_2D : ndarray
        Same shape as R_2D, Sa in g. NaN everywhere if Ts_13 is outside the
        table (Ts_13 >= 5 s, i.e. Ts above about 3.85 s).

    Examples
    --------
    Sa at period 0.13 s (falls in the [0.100, 0.150) table bin) for an
    M7.0 earthquake, at rupture distances of 10 and 20 km:

    >>> import numpy as np
    >>> from landlab.components.Crisis._crisis_core import _compute_Sa_vec
    >>> R_2D = np.array([[10.0, 20.0]])
    >>> np.round(_compute_Sa_vec(0.13, 7.0, R_2D), 5)
    array([[0.46497, 0.23461]])
    """
    # Fixed GMM coefficients common to all period bins
    a3, a4, a5, c1 = 0.265, -0.231, -0.398, 6.75

    R_2D = np.asarray(R_2D, dtype=float)

    # Select magnitude scaling coefficient based on whether M is above or below c1
    a_mag = a4 if M <= c1 else a5

    # Natural log of rupture distance (floor at 1e-6 km to avoid log(0))
    logR = np.log(np.maximum(R_2D, 1e-6))

    # Search the table for the period bin containing Ts_13
    for entry in _SA_TABLE:
        T1, T2 = entry[0], entry[1]
        if T1 <= Ts_13 < T2:
            # Unpack all coefficients for this bin
            T1, T2, a1_F, a2_F, a8_F, a1_S, a2_S, a8_S = entry

            # Compute log(Sa) at the lower (T1, "_F") and upper (T2, "_S") period
            # of the bin
            f1_F = a1_F + a_mag*(M - c1) + a8_F*(8.5 - M)**2 + (a2_F + a3*(M - c1))*logR
            f1_S = a1_S + a_mag*(M - c1) + a8_S*(8.5 - M)**2 + (a2_S + a3*(M - c1))*logR

            # Exponentiate to get Sa in g at each bin boundary period
            Sa_F = np.exp(f1_F)
            Sa_S = np.exp(f1_S)

            # Linearly interpolate between the two period bin boundaries
            return (Sa_S - Sa_F) / (T2 - T1) * (Ts_13 - T1) + Sa_F

    # Return NaN if Ts_13 falls outside all defined period bins
    return np.full(R_2D.shape, np.nan)


# =============================================================================
# Function 5 — Find triggering cells (initiation points of landslides)
# =============================================================================

def _find_triggering_cells_vectorized(
        rows, clms, Elev_Mat, Slope_Mat, Depths_vector,
        Hydro_indicator, Pressures_3D, Theta_3D, Ru,
        Triggering_Indicator, Strength_indicator, c, Phi,
        Theta_s, Theta_r, Cohesion_3D_array, Phi_3D_array,
        Gamma_W, Gamma, t, Earthquake_time_point, PGA_2D_array,
        Seismic_displacement_threshold, M,
        Strength_variation_with_theta, top_bot,
        Seismic_displacement_model, PGV_2D_array,
        Ts, Rathje_method, Rupture_2D_array):
    """
    Identify all triggering cells (landslide initiation points) for a given
    time step using fully vectorized numpy operations over the entire grid.

    Instead of looping over each cell individually, FS and Newmark displacement
    are computed simultaneously for all (rows × cols × depth_layers) combinations
    using broadcasting. Cells that satisfy the failure criterion in any depth
    layer are flagged as triggering cells.

    Parameters
    ----------
    rows, clms : int
        Grid dimensions.
    Elev_Mat, Slope_Mat : ndarray (rows, clms)
        Surface elevation (m) and slope angle (degrees).
    Depths_vector : ndarray, shape (n_dep,)
        Depths at which FS/displacement are evaluated (m). Every layer is
        checked for failure, except a layer at depth 0 (the ground surface).
    Hydro_indicator : {1, 2, 3}
        1 = constant Ru pore pressure coefficient (`Ru` used, `Pressures_3D`
        ignored); 2 or 3 = time-varying pore pressure supplied via
        `Pressures_3D` (2: from HDF5 files, 3: from a Landlab grid field;
        the core treats both the same).
    Pressures_3D : ndarray (rows, clms, n_dep) or None
        Pore pressure head (m) per cell and depth, rows in Landlab's order
        (row 0 = southern edge). Required when `Hydro_indicator` is 2 or 3.
    Theta_3D : ndarray (rows, clms, n_dep) or None
        Volumetric water content per cell and depth, rows in Landlab's
        order. Used only when `Strength_variation_with_theta` is 2 or 3 and
        `Hydro_indicator` is 2 or 3.
    Ru : float
        Pore pressure coefficient. Used only when `Hydro_indicator == 1`.
    Triggering_Indicator : {1, 2}
        1 = rainfall only; 2 = rainfall + seismic (seismic path only runs
        when `t == Earthquake_time_point`).
    Strength_indicator : {1, 2}
        1 = uniform strength (`c`, `Phi`); 2 = spatially variable
        (`Cohesion_3D_array`, `Phi_3D_array`).
    c, Phi : float
        Uniform cohesion (kPa) and friction angle (degrees). Used only
        when `Strength_indicator == 1`.
    Theta_s, Theta_r : float
        Saturated / residual volumetric water content (-).
    Cohesion_3D_array, Phi_3D_array : ndarray (rows, clms, n_dep) or None
        Spatially variable cohesion (kPa) / friction angle (degrees).
        Used only when `Strength_indicator == 2`.
    Gamma_W, Gamma : float
        Unit weight of water and soil (kN/m^3).
    t : int
        Current time step index.
    Earthquake_time_point : int or None
        Time-step index at which the earthquake occurs. The seismic path
        only runs when `t == Earthquake_time_point`.
    PGA_2D_array : ndarray (rows, clms) or None
        Peak Ground Acceleration (g). Required for the seismic path.
    Seismic_displacement_threshold : float or None
        Newmark displacement (cm) at or above which a cell is seismically
        triggered.
    M : float or None
        Earthquake moment magnitude.
    Strength_variation_with_theta : {1, 2, 3}
        1 = fully saturated (Bishop chi = 1); 2 or 3 = chi = (theta - theta_r)
        / (theta_s - theta_r), clipped to [0, 1], with theta from `Theta_3D`
        (2: from HDF5 files, 3: from a Landlab grid field; the core treats
        both the same).
    top_bot : {1, 2}
        1 = report the shallowest failing depth per cell; 2 = report the
        deepest (`top_or_bottom_indicator`).
    Seismic_displacement_model : {1, 2, 3, 4} or None
        1 = Jibson (2007) rigid block; 2 = Rathje et al. (2014) rigid
        block; 3 = Rathje & Antonakos (2011) flexible/decoupled; 4 = Bray
        & Macedo (2019). None when `Triggering_Indicator` is 1.
    PGV_2D_array : ndarray (rows, clms) or None
        Peak Ground Velocity (cm/s). Required for the PGAV variant of
        models 2/3 and for model 4; None otherwise.
    Ts : float or None
        Sliding-mass natural period (s). Required for models 3 and 4;
        None otherwise.
    Rathje_method : {1, 2} or None
        1 = PGAM (uses PGA only); 2 = PGAV (uses PGA and PGV). Required for
        models 2 and 3; None otherwise.
    Rupture_2D_array : ndarray (rows, clms) or None
        Closest-distance-to-rupture (km). Required for models 3 and 4;
        None otherwise.

    Returns
    -------
    trig_cells   : int64 ndarray   Flat linear indices of triggering cells.
    trig_depths  : float32 ndarray Failure depth at each triggering cell (m).
    trig_sources : int64 ndarray   Trigger source: 1=rain, 2=rain+eq, 3=eq only.

    Notes
    -----
    A depth layer at depth 0 (the ground surface, e.g. ``z_min = 0``)
    makes `denom` exactly zero when computing FS, which raises harmless
    numpy `RuntimeWarning`s ("divide by zero", "invalid value"). The
    resulting NaN/inf values at that layer are overwritten with FS = 2.0
    (stable) immediately afterward, so they never affect the returned
    result. Every layer below the ground surface, including the first one
    when ``z_min > 0``, is checked for failure.

    Examples
    --------
    Rainfall-only trigger check (`Triggering_Indicator=1`) on a 2x2 grid
    using a constant Ru pore pressure coefficient (`Hydro_indicator=1`).
    Cell (0, 0) has a 45-degree slope and triggers at 1 m depth; cell
    (1, 0) has a gentler 10-degree slope and only triggers at the deeper
    2 m layer; cell (0, 1) is excluded for being essentially flat, and
    cell (1, 1) is excluded for having a NaN (no-data) elevation.

    >>> import numpy as np
    >>> from landlab.components.Crisis._crisis_core import _find_triggering_cells_vectorized
    >>> Elev_Mat  = np.array([[100.0, 90.0], [80.0, np.nan]])
    >>> Slope_Mat = np.array([[45.0, 0.00005], [10.0, 20.0]])
    >>> Depths_vector = np.array([0.0, 1.0, 2.0], dtype=np.float32)
    >>> trig_cells, trig_depths, trig_sources = _find_triggering_cells_vectorized(
    ...     rows=2, clms=2, Elev_Mat=Elev_Mat, Slope_Mat=Slope_Mat,
    ...     Depths_vector=Depths_vector,
    ...     Hydro_indicator=1, Pressures_3D=None, Theta_3D=None, Ru=0.8,
    ...     Triggering_Indicator=1, Strength_indicator=1, c=2.0, Phi=25.0,
    ...     Theta_s=0.45, Theta_r=0.05, Cohesion_3D_array=None, Phi_3D_array=None,
    ...     Gamma_W=9.81, Gamma=18.0, t=0, Earthquake_time_point=None,
    ...     PGA_2D_array=None, Seismic_displacement_threshold=None, M=None,
    ...     Strength_variation_with_theta=1, top_bot=1,
    ...     Seismic_displacement_model=None, PGV_2D_array=None, Ts=None,
    ...     Rathje_method=None, Rupture_2D_array=None)
    >>> trig_cells.tolist()
    [0, 2]
    >>> trig_depths
    array([1., 2.], dtype=float32)
    >>> trig_sources.tolist()
    [1, 1]

    Above, the first layer is at depth 0 (the ground surface), so it is never
    a failure depth. Every layer below the ground surface is checked,
    including the first one: with the first layer at 0.5 m instead, cell
    (0, 0) now triggers at that first layer.

    >>> Depths_vector = np.array([0.5, 1.0, 2.0], dtype=np.float32)
    >>> trig_cells, trig_depths, trig_sources = _find_triggering_cells_vectorized(
    ...     rows=2, clms=2, Elev_Mat=Elev_Mat, Slope_Mat=Slope_Mat,
    ...     Depths_vector=Depths_vector,
    ...     Hydro_indicator=1, Pressures_3D=None, Theta_3D=None, Ru=0.8,
    ...     Triggering_Indicator=1, Strength_indicator=1, c=2.0, Phi=25.0,
    ...     Theta_s=0.45, Theta_r=0.05, Cohesion_3D_array=None, Phi_3D_array=None,
    ...     Gamma_W=9.81, Gamma=18.0, t=0, Earthquake_time_point=None,
    ...     PGA_2D_array=None, Seismic_displacement_threshold=None, M=None,
    ...     Strength_variation_with_theta=1, top_bot=1,
    ...     Seismic_displacement_model=None, PGV_2D_array=None, Ts=None,
    ...     Rathje_method=None, Rupture_2D_array=None)
    >>> trig_cells.tolist(), trig_depths.tolist()
    ([0, 2], [0.5, 2.0])
    """
    n_dep = len(Depths_vector)
    # Small epsilon to prevent division-by-zero throughout the function
    eps   = 1e-10

    # Accept flat (number_of_nodes,) slope / elevation arrays: reshape them to (rows, clms)
    if np.ndim(Slope_Mat)<2:
        Slope_Mat = Slope_Mat.reshape(rows, clms)

    if np.ndim(Elev_Mat) < 2:
        Elev_Mat = Elev_Mat.reshape(rows, clms)

    # ── Step 1: Build candidate cell mask ────────────────────────────────────
    # Valid cells are those with a non-NaN DEM elevation
    valid_mask = ~np.isnan(Elev_Mat)
    # Cells with effectively zero slope are excluded (infinite slope FS is undefined)
    slope_ok   = Slope_Mat > 0.0001
    # A candidate cell must be valid AND have a non-trivial slope
    candidate  = valid_mask & slope_ok

    # ── Step 2: Build the 3D pore pressure array (rows × cols × n_dep) ───────
    if Hydro_indicator == 1:
        # Constant Ru model: pore pressure head = Ru × (γ/γ_w) × z
        # broadcast_to avoids allocating a full 3D array in memory
        P_3D = np.broadcast_to(
            (Ru * Gamma / Gamma_W) * Depths_vector,
            (rows, clms, n_dep)
        )
    else:
        # Time-varying pore pressure heads: the 3D array from the HDF5 files
        # (indicator 2) or from the grid field (indicator 3)
        P_3D = Pressures_3D

    # ── Step 3: Compute Bishop's unsaturated coefficient χ (3D array) ────────
    if Strength_variation_with_theta == 1 or Hydro_indicator == 1:
        # Fully saturated: χ = 1 everywhere (scalar broadcast)
        Bishop_3D = 1.0
    else:
        # χ scales with degree of saturation: (θ - θ_r) / (θ_s - θ_r)
        Bishop_3D = np.clip(
            (Theta_3D - Theta_r) / (Theta_s - Theta_r + eps), 0.0, 1.0
        )

    # ── Step 4: Compute slope geometry, broadcast across depth layers ─────────
    # Add depth axis so that (rows, clms) arrays broadcast with (n_dep,) arrays
    slope_rad = np.radians(Slope_Mat)[:, :, np.newaxis]   # shape: (rows, clms, 1)
    sin_s     = np.sin(slope_rad)
    cos_s     = np.cos(slope_rad)
    # Depth vector reshaped for broadcasting: shape (1, 1, n_dep)
    D         = Depths_vector[np.newaxis, np.newaxis, :]

    # ── Step 5: Compute Factor of Safety for all cells and depths at once ─────
    # Denominator of the infinite slope FS: γ · z · sin(β) · cos(β)
    denom = Gamma * D * (sin_s * cos_s + eps)
    # Effective normal stress coefficient a1 (accounts for pore pressure and Bishop χ)
    a1    = (Gamma * D * cos_s**2 - Gamma_W * P_3D * Bishop_3D) / denom
    # Cohesion coefficient a2
    a2    = 1.0 / denom

    if Strength_indicator == 1:
        # Uniform c and Phi: scalar operations broadcast over the full 3D array
        FS_3D    = a1 * np.tan(np.radians(Phi)) + a2 * c
        mask_neg = a1 <= 0
        # When normal effective stress is negative, cohesion alone provides resistance
        FS_3D[mask_neg] = (a2 * c)[mask_neg]
    else:
        # Spatially variable strength: Phi_3D and Cohesion_3D arrays shape (rows, clms, n_dep)
        FS_3D    = a1 * np.tan(np.radians(Phi_3D_array)) + a2 * Cohesion_3D_array
        mask_neg = a1 <= 0
        FS_3D[mask_neg] = (a2 * Cohesion_3D_array)[mask_neg]

    # A layer at depth 0 is the ground surface, never a valid failure depth: set FS = 2
    # (every layer below the surface, including the first one when z_min > 0, is checked)
    surface = Depths_vector <= 0
    FS_3D[:, :, surface] = 2.0
    # Invalid and flat cells are always stable — set FS = 2
    FS_3D[~candidate] = 2.0

    # Rainfall failure flag: True where FS ≤ 1 (failure criterion)
    rain_fail_3D  = FS_3D <= 1.0
    # Collapse depth dimension: True if any depth layer fails at this cell
    rain_fail_any = np.any(rain_fail_3D, axis=2)

    # ── Step 6: Rainfall-only path (or any time step other than the earthquake) ─
    # If triggering is rainfall-only, OR if this is not the earthquake time
    # step, return only the rainfall-triggered cells
    if Triggering_Indicator == 1 or t != Earthquake_time_point:
        fail_rows, fail_cols = np.where(rain_fail_any)
        # Return empty arrays if no cell fails at this timestep
        if len(fail_rows) == 0:
            return _EMPTY_ARR.copy(), np.array([], dtype=np.float32), _EMPTY_ARR.copy()

        # Extract the depth-failure profile for all failing cells
        fail_2D = rain_fail_3D[fail_rows, fail_cols, :]
        if top_bot == 1:
            # top_bot = 1: report the shallowest failing depth (first True in depth axis)
            depth_idx = np.argmax(fail_2D, axis=1)
        else:
            # top_bot = 2: report the deepest failing depth (last True in depth axis)
            depth_idx = (n_dep - 1) - np.argmax(fail_2D[:, ::-1], axis=1)

        # Convert (row, col) pairs to flat linear indices
        trig_cells   = np.ravel_multi_index((fail_rows, fail_cols), (rows, clms))
        trig_depths  = Depths_vector[depth_idx].astype(np.float32)
        # All triggers at this point are rainfall-induced (source code = 1)
        trig_sources = np.ones(len(trig_cells), dtype=np.int64)
        return trig_cells, trig_depths, trig_sources

    # ── Step 7: Earthquake path (at the earthquake timestep only) ─────────────
    # Compute critical acceleration ac: the horizontal acceleration required to
    # just bring the slope to limiting equilibrium (FS = 1)
    cos2      = cos_s ** 2
    tan_slope = np.tan(slope_rad)
    # Broadcast PGA across the depth axis for element-wise operations
    PGA_3D    = PGA_2D_array[:, :, np.newaxis]

    # Select strength parameters for ac computation based on Strength_indicator
    if Strength_indicator == 1:
        tan_phi_ac = np.tan(np.radians(Phi))
        c_ac       = c
    else:
        tan_phi_ac = np.tan(np.radians(Phi_3D_array))
        c_ac       = Cohesion_3D_array

    # Numerator of the ac / g expression (dimensionless)
    num_ac   = (c_ac / (Gamma * D * cos2 + eps)
                + tan_phi_ac * (1.0 - (P_3D * Gamma_W) / (Gamma * D * cos2 + eps))
                - tan_slope)
    # Denominator accounts for the effect of slope angle on the yield coefficient
    denom_ac = 1.0 + tan_phi_ac * tan_slope
    ky       = num_ac / denom_ac
    # Floor ac at 0.001 g to avoid numerical issues when ky ≤ 0
    ac       = np.maximum(ky, 0.001)

    # Valid cells for seismic displacement: candidate cells below the ground surface
    eq_valid_3D = np.zeros((rows, clms, n_dep), dtype=bool)
    eq_valid_3D[candidate] = True
    eq_valid_3D[:, :, surface] = False   # a layer at depth 0 (ground surface) is excluded

    # ── Step 8: Compute Newmark displacement using the selected model ─────────

    if Seismic_displacement_model == 1:
        # ── Jibson (2007) — Rigid block ───────────────────────────────────────
        # Displacement is computed only where ac < PGA (i.e. ac/PGA < 1)
        ac2PGA       = ac / np.maximum(PGA_3D, eps)
        eq_fail_mask = eq_valid_3D & (ac2PGA < 1)
        # Clip ratio to valid range for log computation; use safe values elsewhere
        safe_NM  = np.where(eq_fail_mask, np.clip(ac2PGA, eps, 0.99), 0.5)
        safe_1mN = np.where(eq_fail_mask, np.maximum(1.0 - safe_NM, eps), 0.5)
        with np.errstate(divide='ignore', invalid='ignore'):
            # Log term from Jibson (2007) Eq. 7
            log_term = np.where(
                eq_fail_mask,
                np.log10(safe_1mN ** 2.335 * safe_NM ** (-1.478)),
                0.0
            )
        # Newmark displacement in cm, capped at 50 cm
        Dn_3D = np.where(eq_fail_mask,
                         np.minimum(10 ** (-2.710 + log_term + 0.424 * M), 50.0),
                         np.nan)

    elif Seismic_displacement_model == 2:
        # ── Rathje et al. (2014) — Rigid block ───────────────────────────────
        ac2PGA       = ac / np.maximum(PGA_3D, eps)
        eq_fail_mask = eq_valid_3D & (ac2PGA < 1)
        safe_r       = np.where(eq_fail_mask, np.clip(ac2PGA, eps, 0.999), 0.5)
        PGA_sc       = np.maximum(PGA_3D, eps)
        if Rathje_method == 1:
            # PGAM variant: uses PGA and magnitude M
            a1, a2, a3, a4, a5, a6, a7 = 4.89, -4.85, -19.64, 42.49, -29.06, 0.72, 0.89
            Dn_raw = np.exp(
                a1 + a2*safe_r + a3*safe_r**2 + a4*safe_r**3 + a5*safe_r**4
                + a6*np.log(PGA_sc) + a7*(M - 6)
            )
        else:
            # PGAV variant: uses both PGA and PGV as intensity measures
            PGV_3D = PGV_2D_array[:, :, np.newaxis]
            PGV_sc = np.maximum(PGV_3D, eps)
            a1, a2, a3, a4, a5, a6, a7 = -1.56, -4.58, -20.84, 44.75, -30.50, -0.64, 1.55
            Dn_raw = np.exp(
                a1 + a2*safe_r + a3*safe_r**2 + a4*safe_r**3 + a5*safe_r**4
                + a6*np.log(PGA_sc) + a7*np.log(PGV_sc)
            )
        # Cap displacement at 50 cm; NaN where no displacement is computed
        Dn_3D = np.where(eq_fail_mask, np.minimum(Dn_raw, 50.0), np.nan)

    elif Seismic_displacement_model == 3:
        # ── Rathje & Antonakos (2011) — Flexible / Decoupled ─────────────────
        # This model accounts for the dynamic response of the sliding mass by
        # computing a maximum seismic coefficient kmax that amplifies or
        # de-amplifies PGA based on the ratio of Ts to the ground motion period Tm

        # Closest distance from every cell to the fault rupture (no default)
        if Rupture_2D_array is None:
            raise ValueError("Seismic_displacement_model = 3 requires Rupture_2D_array.")
        R_2D = Rupture_2D_array[:, :, np.newaxis]

        # Estimate the predominant period of ground motion Tm
        Tm_2D = np.exp(-1.00 + 0.18*(M - 6) + 0.0038*R_2D + 0.078*1.0)

        # Ratio of sliding-mass period to ground motion period
        ratio        = Ts / np.maximum(Tm_2D, eps)
        log_ratio_01 = np.log(np.maximum(ratio / 0.1, 1e-6))
        log_ratio_02 = np.log(np.maximum(ratio / 0.2, 1e-6))

        # Maximum seismic coefficient kmax: amplified PGA for Ts/Tm ≥ 0.1,
        # equal to PGA otherwise
        kmax = np.where(
            ratio >= 0.1,
            PGA_3D * np.exp(
                (0.459 - 0.702*PGA_3D)*log_ratio_01
                + (-0.228 + 0.076*PGA_3D)*log_ratio_01**2
            ),
            PGA_3D
        )

        # Displacement is computed only where ac < kmax
        ac_to_kmax   = ac / np.maximum(kmax, eps)
        eq_fail_mask = eq_valid_3D & (ac_to_kmax < 1)

        # Period-dependent correction factors for the displacement equation
        f1_Ts = 3.69*Ts - 1.22*Ts**2 if Ts <= 1.5 else 2.78
        f2_Ts = 1.42*Ts if Ts <= 0.5 else 0.71

        safe_r = np.where(eq_fail_mask, np.clip(ac_to_kmax, eps, 0.999), 0.5)
        kmax_s = np.where(eq_fail_mask, np.maximum(kmax, eps), 1.0)

        if Rathje_method == 1:
            # PGAM variant
            a1, a2, a3, a4, a5, a6, a7 = 4.89, -4.85, -19.64, 42.49, -29.06, 0.72, 0.89
            Dn_raw = np.exp(
                a1 + a2*safe_r + a3*safe_r**2 + a4*safe_r**3 + a5*safe_r**4
                + a6*np.log(kmax_s) + a7*(M - 6) + f1_Ts
            )
        else:
            # PGAV variant: uses maximum velocity coefficient kvel_max
            kvel_max = np.where(
                ratio >= 0.2,
                PGV_2D_array[:, :, np.newaxis] * np.exp(
                    0.240*log_ratio_02
                    + (-0.091 - 0.171*PGA_3D)*log_ratio_02**2
                ),
                PGV_2D_array[:, :, np.newaxis]
            )
            kvel_s = np.where(eq_fail_mask, np.maximum(kvel_max, eps), 1.0)
            a1, a2, a3, a4, a5, a6, a7 = -1.56, -4.58, -20.84, 44.75, -30.50, -0.64, 1.55
            Dn_raw = np.exp(
                a1 + a2*safe_r + a3*safe_r**2 + a4*safe_r**3 + a5*safe_r**4
                + a6*np.log(kmax_s) + a7*np.log(kvel_s) + f2_Ts
            )
        Dn_3D = np.where(eq_fail_mask, np.minimum(Dn_raw, 50.0), np.nan)

    else:
        # ── Bray & Macedo (2019) ──────────────────────────────────────────────
        # Uses spectral acceleration Sa(1.3Ts) and PGV as intensity measures.
        # Sa is computed from the rupture distance (no default).
        if Rupture_2D_array is None or PGV_2D_array is None:
            raise ValueError("Seismic_displacement_model = 4 requires Rupture_2D_array and PGV_2D_array.")

        # Compute effective rupture distance including 4.5 km pseudo-depth
        R_2D = np.sqrt(Rupture_2D_array**2 + 4.5**2)

        # Compute spectral acceleration Sa at period 1.3·Ts for every cell
        Sa_2D = _compute_Sa_vec(1.3 * Ts, M, R_2D)
        # Add depth axis for broadcasting
        Sa_3D = Sa_2D[:, :, np.newaxis]

        # Log terms used in the displacement equation (and in PX below)
        ln_ac = np.log(np.maximum(ac, eps))
        ln_Sa = np.log(np.maximum(Sa_3D, eps))

        # Probability term PX (probability of non-zero displacement), two
        # equations depending on Ts. Kept for reference only: CRISIS uses the
        # median displacement, not PX, in the failure decision below.
        if Ts <= 0.7:
            PX = (-2.46 - 2.98*ln_ac - 0.12*ln_ac**2 - 0.71*Ts*ln_ac + 1.69*Ts + 2.76*ln_Sa)
        else:
            PX = (-3.40 - 4.95*ln_ac - 0.30*ln_ac**2 - 0.33*Ts*ln_ac - 0.62*Ts + 2.85*ln_Sa)

        # PGV array (cm/s), with a depth axis for broadcasting
        PGV_3D = PGV_2D_array[:, :, np.newaxis]

        # Select period-dependent displacement model coefficients
        if Ts >= 0.1:
            a1_D, a2_D, a3_D = -5.894, 3.152, -0.910
        else:
            a1_D, a2_D, a3_D = -4.551, -9.690, 0.0

        # Additional PGV-dependent correction when PGV > 115 cm/s
        pgv_high = PGV_3D > 115.0
        a4_D     = np.where(pgv_high, 1.0, 0.0)
        a5_D     = np.where(pgv_high, -4.75, 0.0)

        # The displacement equation is applied where ac > 0.01 g; cells with
        # ac <= 0.01 g are assigned the 50 cm cap below
        eq_fail_mask = eq_valid_3D & (ac > 0.01)

        # Use safe (non-NaN) values for the log terms; zero elsewhere
        safe_lna = np.where(eq_fail_mask, ln_ac, 0.0)
        safe_lnS = np.where(eq_fail_mask, ln_Sa, 0.0)
        safe_lnV = np.where(eq_fail_mask, np.log(np.maximum(PGV_3D, eps)), 0.0)

        # Displacement equation (Bray & Macedo 2019, Eq. for median Dn)
        Dn_raw = np.exp(
            a1_D - 2.491*safe_lna - 0.245*safe_lna**2 + 0.344*safe_lna*safe_lnS
            + 2.703*safe_lnS - 0.089*safe_lnS**2
            + a2_D*Ts + a3_D*Ts**2 + 0.607*M
            + a4_D*safe_lnV + a5_D
        )
        # When ac ≤ 0.01, the slope is essentially always failing; assign 50 cm cap
        Dn_raw = np.where(ac <= 0.01, 50.0, Dn_raw)
        Dn_3D  = np.where(eq_fail_mask | ((ac <= 0.01) & eq_valid_3D),
                          np.minimum(Dn_raw, 50.0), np.nan)

    # ── Step 9: Determine which cells fail seismically ────────────────────────
    # A cell is seismically triggered if its Newmark displacement meets the threshold
    eq_fail_3D = eq_valid_3D & (np.where(np.isnan(Dn_3D), 0.0, Dn_3D) >= Seismic_displacement_threshold)

    # ── Step 10: Combine rainfall and seismic failure flags ───────────────────
    # A cell triggers if it fails by rainfall OR by seismic displacement
    combined_3D = rain_fail_3D | eq_fail_3D
    # Collapse depth axis: True if any depth layer triggers at this cell
    fail_any             = np.any(combined_3D, axis=2)
    fail_rows, fail_cols = np.where(fail_any)

    # Return empty arrays if no triggering cell is found at this timestep
    if len(fail_rows) == 0:
        return _EMPTY_ARR.copy(), np.array([], dtype=np.float32), _EMPTY_ARR.copy()

    n_trig = len(fail_rows)

    # Convert (row, col) indices to flat linear indices
    trig_cells  = np.ravel_multi_index((fail_rows, fail_cols), (rows, clms))

    # Extract failure-flag profiles for the triggering cells across all depths
    combined_2D = combined_3D[fail_rows, fail_cols, :]
    rain_2D     = rain_fail_3D[fail_rows, fail_cols, :]
    eq_2D       = eq_fail_3D[fail_rows, fail_cols, :]

    # Identify the reported failure depth (shallowest or deepest) per cell
    if top_bot == 1:
        depth_idx = np.argmax(combined_2D, axis=1)            # shallowest
    else:
        depth_idx = (n_dep - 1) - np.argmax(combined_2D[:, ::-1], axis=1)  # deepest

    # Retrieve actual depth values (m) at the identified depth indices
    trig_depths  = Depths_vector[depth_idx].astype(np.float32)

    # Assign trigger source codes:
    #   1 = rainfall only  |  2 = both rainfall and earthquake  |  3 = earthquake only
    trig_sources = np.zeros(n_trig, dtype=np.int64)
    idx_range    = np.arange(n_trig)
    in_rain = rain_2D[idx_range, depth_idx]
    in_eq   = eq_2D[idx_range, depth_idx]
    trig_sources[in_rain & ~in_eq] = 1
    trig_sources[in_rain &  in_eq] = 2
    trig_sources[~in_rain & in_eq] = 3

    return trig_cells, trig_depths, trig_sources


# =============================================================================
# Function 6 — Export landslide polygons to shapefile
# =============================================================================

def Modeled_landslides_shapefiles(FID_all_times, Cell_size, x, y, t,
                                  Coordinate_Reference_System, output_dir,
                                  file_name=None):
    """
    Write all cumulative landslide polygons up to time step t as a shapefile.

    Each landslide is represented as a single merged polygon (union of all
    its constituent grid-cell bounding boxes), with one row per unique FID.
    Using unary_union instead of one Polygon per cell reduces the number of
    geometric objects by 10–500×, making the shapefile much faster to write.

    Parameters
    ----------
    FID_all_times : ndarray of float (rows, clms)
        Cumulative landslide-ID raster; NaN where no failure has occurred.
        If every cell is NaN, the function returns without writing a file.
    Cell_size : float
        Grid cell size (m); used to build each cell's bounding box.
    x, y : ndarray
        Easting of each column / northing of each row (m) at the cell
        centers, in the same row order as `FID_all_times` (see
        `compute_coordinate_vectors`).
    t : int
        Time step index; used only to name the output file (when
        `file_name` is not given).
    Coordinate_Reference_System : str
        EPSG string for the shapefile's CRS (e.g. ``"EPSG:32610"``).
    output_dir : str or Path
        Directory the shapefile is written into, as
        ``Predicted_landslides_time_{t}.shp``. Must already exist.
    file_name : str, optional
        Name of the shapefile, used instead of
        ``Predicted_landslides_time_{t}.shp``. A steady-state run (constant
        Ru) writes a single shapefile, ``Predicted_landslides.shp``.

    Returns
    -------
    None
        Writes a shapefile (and its `.dbf`/`.shx`/`.prj`/`.cpg` sidecars)
        as a side effect.

    Notes
    -----
    ESRI Shapefile field names are truncated to 10 characters, so the
    ``'landslide_number'`` column is written (and read back) as
    ``'landslide_'``.

    Examples
    --------
    A 2x3 grid holding two 2-cell landslides (IDs 1 and 2), written to a
    temporary directory and read back to confirm the shapefile is correct:

    >>> import numpy as np, geopandas as gpd, tempfile, os
    >>> from landlab.components.Crisis._crisis_core import Modeled_landslides_shapefiles
    >>> FID_all_times = np.array([
    ...     [1.0, 1.0, np.nan],
    ...     [np.nan, 2.0, 2.0],
    ... ])
    >>> x = np.array([0.0, 10.0, 20.0])
    >>> y = np.array([0.0, 10.0])
    >>> outdir = tempfile.mkdtemp()
    >>> Modeled_landslides_shapefiles(FID_all_times, 10.0, x, y, 0, "EPSG:32610", outdir)
    >>> sorted(os.listdir(outdir))
    ['Predicted_landslides_time_0.cpg', 'Predicted_landslides_time_0.dbf', 'Predicted_landslides_time_0.prj', 'Predicted_landslides_time_0.shp', 'Predicted_landslides_time_0.shx']
    >>> gdf = gpd.read_file(os.path.join(outdir, "Predicted_landslides_time_0.shp"))
    >>> list(gdf.columns)
    ['landslide_', 'geometry']
    >>> list(gdf['landslide_'])
    [1.0, 2.0]
    >>> gdf.total_bounds
    array([-5., -5., 25., 15.])

    In a steady-state run (constant Ru, a single time step), the single
    shapefile has no time step in its name:

    >>> outdir = tempfile.mkdtemp()
    >>> Modeled_landslides_shapefiles(FID_all_times, 10.0, x, y, 0, "EPSG:32610", outdir,
    ...                               file_name="Predicted_landslides.shp")
    >>> sorted(f for f in os.listdir(outdir) if f.endswith(".shp"))
    ['Predicted_landslides.shp']
    """
    # Skip writing if no landslides have been delineated yet
    if np.all(np.isnan(FID_all_times)):
        return

    # Get sorted list of unique landslide IDs, excluding NaN (no-landslide cells)
    unique_ids = np.unique(FID_all_times[~np.isnan(FID_all_times)])
    gdfs = []

    for lid in unique_ids:
        # Find the row and column indices of all cells belonging to this landslide
        rows_l, cols_l = np.where(FID_all_times == lid)
        if len(rows_l) == 0:
            continue

        # Half-cell offset: each cell is represented by its bounding box
        h = 0.5 * Cell_size

        # Build a list of axis-aligned bounding boxes, one per landslide cell
        cells = [
            _shapely_box(x[cols_l[ii]] - h, y[rows_l[ii]] - h,
                         x[cols_l[ii]] + h, y[rows_l[ii]] + h)
            for ii in range(len(rows_l))
        ]

        # Merge all cell boxes into a single polygon for this landslide
        merged = _unary_union(cells)

        # Create a one-row GeoDataFrame for this landslide
        gdfs.append(gpd.GeoDataFrame(
            {'landslide_number': [lid], 'geometry': [merged]},
            crs=Coordinate_Reference_System
        ))

    # Combine all landslide GeoDataFrames into one; create empty frame if none
    gdf_all = pd.concat(gdfs, ignore_index=True) if gdfs else gpd.GeoDataFrame(
        columns=['landslide_number', 'geometry'], crs=Coordinate_Reference_System
    )

    # Write to a shapefile named by time step (or by file_name, e.g. for a steady-state run)
    out_path = Path(output_dir) / (file_name or f'Predicted_landslides_time_{t}.shp')
    gdf_all.to_file(str(out_path))


# =============================================================================
# Function 7 — Excel output writers
# =============================================================================

def to_north_to_south_index(linear_index, grid_shape):
    """Convert Landlab linear indices to north-up (DEM file) row order.

    Inside CRISIS, arrays are in Landlab's orientation (row 0 = southern
    edge), so a linear index ``row * ncols + col`` is a Landlab node ID.
    The Excel outputs give cell indices in the row order of the input files
    (north-up, row 0 = northern edge); this converts an index to the same
    cell's index in that row order (row -> nrows - 1 - row, column unchanged).

    Parameters
    ----------
    linear_index : int or array-like of int
        Row-major linear index (index) in Landlab orientation.
    grid_shape : tuple of int
        ``(nrows, ncols)`` of the grid.

    Returns
    -------
    ndarray of int64
        Row-major linear index (indices) of the same cell(s), with rows ordered
        from north to south.

    Examples
    --------
    On a 2 x 3 grid, Landlab node 0 (south-west corner) is cell 3 when the
    grid stores its rows from north to south; node 5 (north-east corner)
    is cell 2.

    >>> from landlab.components.Crisis._crisis_core import to_north_to_south_index
    >>> to_north_to_south_index([0, 1, 2, 3, 4, 5], (2, 3)).tolist()
    [3, 4, 5, 0, 1, 2]
    """
    nrows, ncols = grid_shape
    row, col = np.divmod(np.asarray(linear_index, dtype=np.int64), ncols)
    return (nrows - 1 - row) * ncols + col


def _save_individual_excel(state, output_dir, grid_shape, steady_state=False):
    """Write per-landslide results to an Excel workbook.

    Each sheet corresponds to one output variable; rows = landslides, columns = time steps
    (header row t=0, t=1, ... labels which time step each column is). Sheet names
    carry the units of their values (_m, _m2, _m3); empty cells mean no landslide.

    The Triggering_Cells sheets hold row-major linear indices in the row
    order of the ASCII DEM file (north-up, row 0 = northern edge), converted
    from the Landlab node IDs kept in `state` with `to_north_to_south_index`.

    Parameters
    ----------
    state : dict
        Must contain the ten ``'Individual_{Current,Cumulative}_*_with_time'``
        list-of-lists keys (see `CrisisComponent.state`), each holding one ragged row
        per time step. Rows are padded with `pad_rows` and transposed so
        that time steps become columns, labeled t=0, t=1, ...
    output_dir : str or Path
        Directory the workbook is written into, as
        ``Individual_landslides_with_time.xlsx``. Must already exist.
    grid_shape : tuple of int
        ``(nrows, ncols)`` of the grid, used to convert the triggering-cell
        indices to rows ordered from north to south.
    steady_state : bool, optional
        True for a steady-state run (constant Ru): the workbook is then
        ``Individual_landslides.xlsx``, with one sheet per quantity
        (``Triggering_Cells``, ``Triggering_Depths_m``, ``Triggering_Sources``,
        ``Areas_m2``, ``Volumes_m3``) and a single column holding the
        landslides of the run (no time steps).

    Returns
    -------
    None
        Writes the workbook as a side effect.

    Examples
    --------
    Two time steps on a 2 x 3 grid: 2 landslides triggered at step 0
    (Landlab nodes 0, 1), 1 more at step 1 (node 2) — so the "Current" sheet
    has a ragged 2nd column, padded with an empty cell, while the "Cumulative" sheet
    accumulates over time. Nodes 0, 1, 2 are the southern row, which is
    row 1 (cells 3, 4, 5) when the rows are ordered from north to south.

    >>> import pandas as pd, tempfile, os
    >>> from landlab.components.Crisis._crisis_core import _save_individual_excel
    >>> state = {
    ...     'Individual_Current_triggering_linear_indices_with_time': [[0, 1], [2]],
    ...     'Individual_Current_triggering_depths_with_time':         [[1.0, 2.0], [1.5]],
    ...     'Individual_Current_triggering_sources_with_time':        [[1, 1], [1]],
    ...     'Individual_Current_areas_with_time':                     [[16.0, 8.0], [4.0]],
    ...     'Individual_Current_volumes_with_time':                   [[11.9, 5.0], [2.0]],
    ...     'Individual_Cumulative_triggering_linear_indices_with_time': [[0, 1], [0, 1, 2]],
    ...     'Individual_Cumulative_triggering_depths_with_time':         [[1.0, 2.0], [1.0, 2.0, 1.5]],
    ...     'Individual_Cumulative_triggering_sources_with_time':        [[1, 1], [1, 1, 1]],
    ...     'Individual_Cumulative_areas_with_time':                     [[16.0, 8.0], [16.0, 8.0, 4.0]],
    ...     'Individual_Cumulative_volumes_with_time':                   [[11.9, 5.0], [11.9, 5.0, 2.0]],
    ... }
    >>> outdir = tempfile.mkdtemp()
    >>> _save_individual_excel(state, outdir, grid_shape=(2, 3))
    >>> sorted(os.listdir(outdir))
    ['Individual_landslides_with_time.xlsx']
    >>> xl = pd.ExcelFile(os.path.join(outdir, "Individual_landslides_with_time.xlsx"))
    >>> xl.sheet_names
    ['Triggering_Cells_Current', 'Triggering_Depths_Current_m', 'Triggering_Sources_Current', 'Areas_Current_m2', 'Volumes_Current_m3', 'Triggering_Cells_Cumulative', 'Triggering_Depths_Cumulative_m', 'Triggering_Sources_Cumulative', 'Areas_Cumulative_m2', 'Volumes_Cumulative_m3']
    >>> xl.parse("Triggering_Cells_Current")
       t=0  t=1
    0    3  5.0
    1    4  NaN

    In a steady-state run (constant Ru, a single time step), the workbook
    has no time steps: one sheet per quantity, one column of landslides.

    >>> outdir = tempfile.mkdtemp()
    >>> _save_individual_excel(state, outdir, grid_shape=(2, 3), steady_state=True)
    >>> sorted(os.listdir(outdir))
    ['Individual_landslides.xlsx']
    >>> xl = pd.ExcelFile(os.path.join(outdir, "Individual_landslides.xlsx"))
    >>> xl.sheet_names
    ['Triggering_Cells', 'Triggering_Depths_m', 'Triggering_Sources', 'Areas_m2', 'Volumes_m3']
    >>> xl.parse("Triggering_Cells")["Triggering cell"].tolist()
    [3, 4, 5]
    """
    def to_df(key):
        # pad_rows makes the ragged list rectangular, then transpose so that
        # time steps become columns (few) and landslides become rows (many),
        # which fits Excel's row limit better
        rows = state[key]
        if "triggering_linear_indices" in key:
            # Landlab node IDs -> row-major indices with rows from north to south, like all output files
            rows = [to_north_to_south_index(r, grid_shape).tolist() for r in rows]
        df = pd.DataFrame(pad_rows(rows)).T
        df.columns = [f't={i}' for i in range(df.shape[1])]
        return df

    if steady_state:
        # steady state (constant Ru): no time steps, the landslides of the run
        # (the last entry of the cumulative lists), one sheet per quantity
        def last(key):
            values = list(state[key][-1]) if len(state[key]) else []
            if "triggering_linear_indices" in key:
                values = to_north_to_south_index(values, grid_shape).tolist()
            return values

        path = Path(output_dir) / 'Individual_landslides.xlsx'
        with pd.ExcelWriter(str(path), engine='openpyxl') as writer:
            for key, sheet, column in (
                ('Individual_Cumulative_triggering_linear_indices_with_time', 'Triggering_Cells',    'Triggering cell'),
                ('Individual_Cumulative_triggering_depths_with_time',         'Triggering_Depths_m', 'Triggering depth (m)'),
                ('Individual_Cumulative_triggering_sources_with_time',        'Triggering_Sources',  'Triggering source'),
                ('Individual_Cumulative_areas_with_time',                     'Areas_m2',            'Area (m^2)'),
                ('Individual_Cumulative_volumes_with_time',                   'Volumes_m3',          'Volume (m^3)'),
            ):
                pd.DataFrame({column: last(key)}).to_excel(writer, sheet_name=sheet, index=False)
        return

    path = Path(output_dir) / 'Individual_landslides_with_time.xlsx'
    with pd.ExcelWriter(str(path), engine='openpyxl') as writer:
        kw = dict(index=False, header=True, na_rep='')
        # ── Current (new-at-this-timestep) landslide data ─────────────────────
        to_df('Individual_Current_triggering_linear_indices_with_time'
              ).to_excel(writer, sheet_name='Triggering_Cells_Current',       **kw)
        to_df('Individual_Current_triggering_depths_with_time'
              ).to_excel(writer, sheet_name='Triggering_Depths_Current_m',    **kw)
        to_df('Individual_Current_triggering_sources_with_time'
              ).to_excel(writer, sheet_name='Triggering_Sources_Current',     **kw)
        to_df('Individual_Current_areas_with_time'
              ).to_excel(writer, sheet_name='Areas_Current_m2',               **kw)
        to_df('Individual_Current_volumes_with_time'
              ).to_excel(writer, sheet_name='Volumes_Current_m3',             **kw)
        # ── Cumulative (all-time) landslide data ──────────────────────────────
        to_df('Individual_Cumulative_triggering_linear_indices_with_time'
              ).to_excel(writer, sheet_name='Triggering_Cells_Cumulative',    **kw)
        to_df('Individual_Cumulative_triggering_depths_with_time'
              ).to_excel(writer, sheet_name='Triggering_Depths_Cumulative_m', **kw)
        to_df('Individual_Cumulative_triggering_sources_with_time'
              ).to_excel(writer, sheet_name='Triggering_Sources_Cumulative',  **kw)
        to_df('Individual_Cumulative_areas_with_time'
              ).to_excel(writer, sheet_name='Areas_Cumulative_m2',            **kw)
        to_df('Individual_Cumulative_volumes_with_time'
              ).to_excel(writer, sheet_name='Volumes_Cumulative_m3',          **kw)


def _save_summary_excel(state, output_dir, steady_state=False):
    """
    Write aggregate statistics (total counts, area, volume per time step)
    to a separate Excel workbook for quick overview.

    Parameters
    ----------
    state : dict
        Must contain the six ``'Total_{Current,Cumulative}_*_vector_with_time'``
        1-D arrays (see `CrisisComponent.state`), one entry per time step.
    output_dir : str or Path
        Directory the workbook is written into, as
        ``Total_landslides_with_time.xlsx``. Must already exist.
    steady_state : bool, optional
        True for a steady-state run (constant Ru): the workbook is then
        ``Total_landslides.xlsx``, with a single row (number of landslides,
        total area, total volume) and no time steps.

    Returns
    -------
    None
        Writes the workbook as a side effect: one row per time step
        (Time = t=0, t=1, ...), with the units in the column names.

    Examples
    --------
    >>> import numpy as np, pandas as pd, tempfile, os
    >>> from landlab.components.Crisis._crisis_core import _save_summary_excel
    >>> state = {
    ...     'Total_Current_NL_vector_with_time':        np.array([2.0, 1.0]),
    ...     'Total_Current_Area_vector_with_time':      np.array([24.0, 4.0]),
    ...     'Total_Current_Volume_vector_with_time':    np.array([16.9, 2.0]),
    ...     'Total_Cumulative_NL_vector_with_time':     np.array([2.0, 3.0]),
    ...     'Total_Cumulative_Area_vector_with_time':   np.array([24.0, 28.0]),
    ...     'Total_Cumulative_Volume_vector_with_time': np.array([16.9, 18.9]),
    ... }
    >>> outdir = tempfile.mkdtemp()
    >>> _save_summary_excel(state, outdir)
    >>> sorted(os.listdir(outdir))
    ['Total_landslides_with_time.xlsx']
    >>> summary = pd.read_excel(os.path.join(outdir, "Total_landslides_with_time.xlsx"))
    >>> list(summary.columns)
    ['Time', 'Current Number of Landslides', 'Current Total Area (m^2)', 'Current Total Volume (m^3)', 'Cumulative Number of Landslides', 'Cumulative Total Area (m^2)', 'Cumulative Total Volume (m^3)']
    >>> summary['Time'].tolist()
    ['t=0', 't=1']
    >>> summary['Cumulative Total Volume (m^3)'].tolist()
    [16.9, 18.9]

    In a steady-state run (constant Ru, a single time step), the workbook has
    a single row and no time steps:

    >>> outdir = tempfile.mkdtemp()
    >>> _save_summary_excel(state, outdir, steady_state=True)
    >>> sorted(os.listdir(outdir))
    ['Total_landslides.xlsx']
    >>> total = pd.read_excel(os.path.join(outdir, "Total_landslides.xlsx"))
    >>> list(total.columns), len(total)
    (['Number of Landslides', 'Total Area (m^2)', 'Total Volume (m^3)'], 1)
    >>> total.iloc[0].tolist()
    [3.0, 28.0, 18.9]
    """
    if steady_state:
        # steady state (constant Ru): a single row with the totals of the run
        df = pd.DataFrame({
            'Number of Landslides': [state['Total_Cumulative_NL_vector_with_time'][-1]],
            'Total Area (m^2)':     [state['Total_Cumulative_Area_vector_with_time'][-1]],
            'Total Volume (m^3)':   [state['Total_Cumulative_Volume_vector_with_time'][-1]],
        })
        df.to_excel(str(Path(output_dir) / 'Total_landslides.xlsx'), index=False, engine='openpyxl', na_rep='NaN')
        return

    n_time = len(state['Total_Current_NL_vector_with_time'])
    df = pd.DataFrame({
        'Time':                            [f't={i}' for i in range(n_time)],
        'Current Number of Landslides':    state['Total_Current_NL_vector_with_time'],
        'Current Total Area (m^2)':        state['Total_Current_Area_vector_with_time'],
        'Current Total Volume (m^3)':      state['Total_Current_Volume_vector_with_time'],
        'Cumulative Number of Landslides': state['Total_Cumulative_NL_vector_with_time'],
        'Cumulative Total Area (m^2)':     state['Total_Cumulative_Area_vector_with_time'],
        'Cumulative Total Volume (m^3)':   state['Total_Cumulative_Volume_vector_with_time'],
    })
    path = Path(output_dir) / 'Total_landslides_with_time.xlsx'
    df.to_excel(str(path), index=False, engine='openpyxl', na_rep='NaN')
