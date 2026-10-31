"""""
=============================================================================

  CRISIS — Coupled Regional Rainfall-Induced and Seismic Slope Instability
           Simulations

  Purpose : Landlab Component implementation of the CRISIS forward analysis model.
            Rainfall only, earthquake only, or coupled (rainfall + seismic) triggers

  Authors : Mirna Kassem
            University of California, Berkeley
            kassem_mirna@berkeley.edu

            Prof. Dimitrios Zekkos
            University of California, Berkeley
            zekkos@berkeley.edu

            Code was also written and reviewed by Yuval Shmilovitz (Yuval.Shmilovitz@colorado.edu)

   Documentation: Additional documentation and manual stored at https://github.com/yuvalsmilo/CRISIS-Landlab.git
=============================================================================
"""

# =============================================================================
# Imports
# =============================================================================
# All CRISIS physics functions are imported from the private _crisis_core.py module;
# this file contains interface logic only

from __future__ import annotations
import re
from pathlib import Path
import h5py
import numpy as np
from landlab import Component, RasterModelGrid

Gamma_W = 9.81  # unit weight of water (kN/m³)

# pull in the physics, file-loading and output functions from the core CRISIS module.
# _crisis_core.py sits next to this file: import it relative to the package when
# this file is imported as part of Landlab (landlab.components.crisis), and
# directly when imported from its own folder.
try:
    from . import _crisis_core as _core
except ImportError:
    import _crisis_core as _core

_build_upslope_map                = _core._build_upslope_map
_find_triggering_cells_vectorized = _core._find_triggering_cells_vectorized
Determine_Pseudo_3D_landslides    = _core.Determine_Pseudo_3D_landslides
Modeled_landslides_shapefiles     = _core.Modeled_landslides_shapefiles
_save_individual_excel            = _core._save_individual_excel
_save_summary_excel               = _core._save_summary_excel
load_h5                           = _core.load_h5
read_ascii_header                 = _core.read_ascii_header
lower_left_cell_center            = _core.lower_left_cell_center
compute_coordinate_vectors        = _core.compute_coordinate_vectors
Rupture_Dist                      = _core.Rupture_Dist


# =============================================================================
# Module-level helper
# =============================================================================
# Landlab stores node data as flat 1D arrays; CRISIS physics functions expect
# 2D (nrows, ncols) arrays — this helper reshapes them on the fly.

# small helper: turn a flat Landlab field back into a 2D grid shape
def _node_to_2d(grid, field_name):
    """Return a node field reshaped to (nrows, ncols).

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.crisis.crisis_landlab_component import _node_to_2d
    >>> grid = RasterModelGrid((2, 3))
    >>> _ = grid.add_field("topographic__elevation", np.arange(6.0), at="node")
    >>> _node_to_2d(grid, "topographic__elevation")
    array([[0., 1., 2.],
           [3., 4., 5.]])
    """
    # grab the named field and reshape it to match the grid's row/col layout
    return grid.at_node[field_name].reshape(grid.shape)


# =============================================================================
# Parameter check
# =============================================================================
# No parameter has a default value: every key a run needs must be given in
# params. Which keys are needed depends on the settings (e.g. Ru only with
# Hydrological_model_indicator = 1), so they are checked here, case by case,
# before anything else is done. Every indicator starts at 1.

# keys every run needs, whatever the settings
_ALWAYS_REQUIRED = (
    "dem_file", "output_dir", "Coordinate_Reference_System",
    "Gamma", "z_min", "z_max", "Depth_points",
    "Min_LD_area", "top_or_bottom_indicator",
    "Triggering_Indicator", "Hydrological_model_indicator", "Strength_indicator",
)

# allowed values of the choice keys
_CHOICES = {
    "Triggering_Indicator":          (1, 2),
    "Hydrological_model_indicator":  (1, 2, 3),
    "Strength_variation_with_theta": (1, 2, 3),
    "Strength_indicator":            (1, 2),
    "top_or_bottom_indicator":       (1, 2),
    "Seismic_displacement_model":    (1, 2, 3, 4),
    "Rathje_method":                 (1, 2),
}

# the four corners of the fault plane (Seismic_displacement_model = 3 or 4)
_FAULT_KEYS = tuple(f"fault_{axis}{corner}" for corner in (1, 2, 3, 4) for axis in "XYZ")


def check_params(params):
    """Check that *params* holds every key this run needs, with valid values.

    No parameter has a default value. The keys needed depend on the
    settings:

    - always: ``dem_file``, ``output_dir``, ``Coordinate_Reference_System``,
      ``Gamma``, ``z_min``, ``z_max``, ``Depth_points``, ``Min_LD_area``,
      ``top_or_bottom_indicator``, ``Triggering_Indicator``,
      ``Hydrological_model_indicator`` and ``Strength_indicator``;
    - ``Hydrological_model_indicator = 1``: ``Ru`` (a constant Ru is steady
      state: the run has a single time step, so ``Earthquake_time_point``
      must be 0);
    - ``Hydrological_model_indicator = 2``: ``hydro_files_dir``;
    - ``Hydrological_model_indicator = 2`` or ``3``:
      ``Strength_variation_with_theta``; then ``theta_files_dir`` with
      ``Strength_variation_with_theta = 2``, and ``Theta_s`` and ``Theta_r``
      with ``Strength_variation_with_theta = 2`` or ``3``;
    - ``Strength_indicator = 1``: ``c`` and ``Phi``; ``Strength_indicator = 2``:
      ``cohesion_3d_file`` and ``phi_3d_file``;
    - ``Triggering_Indicator = 2``: ``pga_file``, ``Earthquake_time_point``,
      ``Seismic_displacement_threshold`` and ``Seismic_displacement_model``;
      then ``Rathje_method`` (models 2 and 3), ``M`` (every model except
      Rathje et al. (2014) with ``Rathje_method = 2``), ``pgv_file`` (models
      2 and 3 with ``Rathje_method = 2``, and model 4), and ``Ts`` and the
      fault corners ``fault_X1`` ... ``fault_Z4`` (models 3 and 4).

    A key counts as missing when it is absent, ``None`` or an empty string.
    Every problem is reported at once.

    Parameters
    ----------
    params : dict
        CRISIS parameter dict.

    Raises
    ------
    TypeError
        If *params* is not a dict.
    ValueError
        If a required key is missing, or a choice key (an indicator, a
        model number, ...) has a value that is not allowed. The message
        lists every such key.

    Examples
    --------
    >>> from landlab.components.crisis.crisis_landlab_component import check_params
    >>> params = {
    ...     "dem_file": "dem.asc", "output_dir": "Outputs",
    ...     "Coordinate_Reference_System": "EPSG:32610",
    ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875, "Depth_points": 8,
    ...     "Min_LD_area": 12.0, "top_or_bottom_indicator": 1,
    ...     "Triggering_Indicator": 1, "Hydrological_model_indicator": 1, "Ru": 0.4,
    ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
    ... }
    >>> check_params(params)    # complete: nothing is reported

    ``Ru`` is required only with ``Hydrological_model_indicator = 1``; with
    the HDF5 files (indicator 2), their folder is required instead, and so
    is ``Strength_variation_with_theta``:

    >>> check_params({**params, "Ru": None})
    Traceback (most recent call last):
    ...
    ValueError: Missing or invalid entries in params:
      - Ru is missing (required when Hydrological_model_indicator = 1)
    >>> check_params({**params, "Hydrological_model_indicator": 2})
    Traceback (most recent call last):
    ...
    ValueError: Missing or invalid entries in params:
      - hydro_files_dir is missing (required when Hydrological_model_indicator = 2)
      - Strength_variation_with_theta is missing (required when Hydrological_model_indicator = 2)

    A coupled run with the Bray & Macedo model needs the earthquake keys:

    >>> check_params({**params, "Triggering_Indicator": 2, "Seismic_displacement_model": 4,
    ...               "pga_file": "PGA.h5", "Earthquake_time_point": 0,
    ...               "Seismic_displacement_threshold": 5.0})
    Traceback (most recent call last):
    ...
    ValueError: Missing or invalid entries in params:
      - M is missing (required when Seismic_displacement_model = 4)
      - pgv_file is missing (required when Seismic_displacement_model = 4)
      - Ts is missing (required when Seismic_displacement_model = 4)
      - fault_X1, fault_Y1, fault_Z1, fault_X2, fault_Y2, fault_Z2, fault_X3, fault_Y3, fault_Z3, fault_X4, fault_Y4, fault_Z4 are missing (required when Seismic_displacement_model = 4)

    With a constant Ru (a single time step), the earthquake can only be at
    time step 0:

    >>> check_params({**params, "Triggering_Indicator": 2, "Seismic_displacement_model": 1,
    ...               "pga_file": "PGA.h5", "Earthquake_time_point": 3,
    ...               "Seismic_displacement_threshold": 5.0, "M": 6.5})
    Traceback (most recent call last):
    ...
    ValueError: Missing or invalid entries in params:
      - Earthquake_time_point must be 0 with Hydrological_model_indicator = 1 (a single time step), not 3

    Values that are not allowed are reported too (every indicator starts at
    1):

    >>> check_params({**params, "Strength_indicator": 0})
    Traceback (most recent call last):
    ...
    ValueError: Missing or invalid entries in params:
      - Strength_indicator must be 1 or 2, not 0
    """
    if not isinstance(params, dict):
        raise TypeError(f"params must be a dict, not {type(params).__name__}.")

    problems = []

    def given(key):
        value = params.get(key)
        return value is not None and not (isinstance(value, str) and value.strip() == "")

    def need(keys, reason=None):
        absent = [key for key in keys if not given(key)]
        if absent:
            verb = "is" if len(absent) == 1 else "are"
            when = f" (required when {reason})" if reason else ""
            problems.append(f"{', '.join(absent)} {verb} missing{when}")

    def choice(key):
        """The value of a choice key, or None if it is missing or not allowed."""
        if not given(key):
            return None
        value, allowed = params[key], _CHOICES[key]
        if isinstance(value, bool) or value not in allowed:
            listed = ", ".join(str(a) for a in allowed[:-1]) + f" or {allowed[-1]}"
            problems.append(f"{key} must be {listed}, not {value!r}")
            return None
        return value

    # ── keys of every run ─────────────────────────────────────────────────
    for key in _ALWAYS_REQUIRED:
        need([key])
    choice("top_or_bottom_indicator")

    # ── pore pressure heads and water contents ────────────────────────────
    hydro = choice("Hydrological_model_indicator")
    theta = choice("Strength_variation_with_theta")
    if hydro == 1:
        need(["Ru"], "Hydrological_model_indicator = 1")
    elif hydro in (2, 3):
        if hydro == 2:
            need(["hydro_files_dir"], "Hydrological_model_indicator = 2")
        need(["Strength_variation_with_theta"], f"Hydrological_model_indicator = {hydro}")
        if theta == 2:
            need(["theta_files_dir"], "Strength_variation_with_theta = 2")
        if theta in (2, 3):
            need(["Theta_s", "Theta_r"], f"Strength_variation_with_theta = {theta}")

    # ── shear strength ────────────────────────────────────────────────────
    strength = choice("Strength_indicator")
    if strength == 1:
        need(["c", "Phi"], "Strength_indicator = 1")
    elif strength == 2:
        need(["cohesion_3d_file", "phi_3d_file"], "Strength_indicator = 2")

    # ── earthquake (coupled mode) ─────────────────────────────────────────
    if choice("Triggering_Indicator") == 2:
        need(["pga_file", "Earthquake_time_point", "Seismic_displacement_threshold",
              "Seismic_displacement_model"], "Triggering_Indicator = 2")
        # a constant Ru has a single time step (index 0): the earthquake is at step 0
        if hydro == 1 and given("Earthquake_time_point") and params["Earthquake_time_point"] != 0:
            problems.append(
                "Earthquake_time_point must be 0 with Hydrological_model_indicator = 1 "
                f"(a single time step), not {params['Earthquake_time_point']!r}")
        model = choice("Seismic_displacement_model")
        method = None
        if model in (2, 3):
            need(["Rathje_method"], f"Seismic_displacement_model = {model}")
            method = choice("Rathje_method")
        if model in (1, 3, 4) or (model == 2 and method == 1):
            need(["M"], f"Seismic_displacement_model = {model}"
                        + (f" and Rathje_method = {method}" if model == 2 else ""))
        if (model in (2, 3) and method == 2) or model == 4:
            need(["pgv_file"], f"Seismic_displacement_model = {model}"
                               + (" and Rathje_method = 2" if model in (2, 3) else ""))
        if model in (3, 4):
            need(["Ts"], f"Seismic_displacement_model = {model}")
            need(_FAULT_KEYS, f"Seismic_displacement_model = {model}")

    if problems:
        raise ValueError("Missing or invalid entries in params:\n"
                         + "\n".join(f"  - {p}" for p in problems))


# =============================================================================
# Number of time steps of a set of HDF5 files
# =============================================================================

def _count_time_steps(folder, prefix, Depth_points):
    """Count the time steps of the HDF5 files ``{prefix}_T{t}Z{z}.h5`` in *folder*.

    The time steps must be T0, T1, ... without gaps, and every time step must
    have exactly one file per depth layer, Z1 ... Z{Depth_points} (no layer
    missing, none beyond ``Depth_points``); every time step found is part of
    the run.

    Parameters
    ----------
    folder : str or Path
        Folder holding the files.
    prefix : str
        File name prefix: ``"Pressure_head"`` or ``"Theta"``.
    Depth_points : int
        Number of depth layers.

    Returns
    -------
    int
        Number of time steps.

    Raises
    ------
    FileNotFoundError
        If the folder holds no such file, if a time step is missing (a gap),
        or if a time step is missing one of its depth layers.
    ValueError
        If a time step has depth layers beyond ``Depth_points`` (the files
        and ``Depth_points`` do not match).

    Examples
    --------
    >>> import os, tempfile
    >>> from landlab.components.crisis.crisis_landlab_component import _count_time_steps
    >>> folder = tempfile.mkdtemp()
    >>> for t in range(3):
    ...     for z in (1, 2):
    ...         open(os.path.join(folder, f"Pressure_head_T{t}Z{z}.h5"), "w").close()
    >>> _count_time_steps(folder, "Pressure_head", Depth_points=2)
    3

    Every time step needs all its depth layers:

    >>> _count_time_steps(folder, "Pressure_head", Depth_points=3)
    Traceback (most recent call last):
    ...
    FileNotFoundError: Pressure_head files: time step T0 is missing depth layer(s) Z3 (expected Z1 ... Z3 for every time step).

    and no depth layer beyond ``Depth_points``:

    >>> _count_time_steps(folder, "Pressure_head", Depth_points=1)
    Traceback (most recent call last):
    ...
    ValueError: Pressure_head files: time step T0 has depth layer(s) Z2 beyond Depth_points = 1 (expected exactly Z1 ... Z1 for every time step).

    The time steps must have no gap:

    >>> for z in (1, 2):
    ...     os.remove(os.path.join(folder, f"Pressure_head_T1Z{z}.h5"))
    >>> _count_time_steps(folder, "Pressure_head", Depth_points=2)
    Traceback (most recent call last):
    ...
    FileNotFoundError: Pressure_head files: time step(s) T1 missing (found T0 ... T2; the time steps must be T0, T1, ... without gaps).
    """
    pattern = re.compile(rf"^{re.escape(prefix)}_T(\d+)Z(\d+)\.h5$")
    layers = {}                                    # time step -> set of depth layers found
    for path in Path(folder).iterdir() if Path(folder).is_dir() else ():
        match = pattern.match(path.name)
        if match:
            layers.setdefault(int(match.group(1)), set()).add(int(match.group(2)))
    if not layers:
        raise FileNotFoundError(
            f"No {prefix}_T{{t}}Z{{z}}.h5 files found in '{folder}'."
        )
    last = max(layers)
    gaps = [t for t in range(last + 1) if t not in layers]
    if gaps:
        raise FileNotFoundError(
            f"{prefix} files: time step(s) {', '.join(f'T{t}' for t in gaps)} missing "
            f"(found T0 ... T{last}; the time steps must be T0, T1, ... without gaps)."
        )
    expected = set(range(1, Depth_points + 1))
    for t in range(last + 1):
        # a missing layer is a missing file
        missing = sorted(expected - layers[t])
        if missing:
            raise FileNotFoundError(
                f"{prefix} files: time step T{t} is missing depth layer(s) "
                f"{', '.join(f'Z{z}' for z in missing)} (expected Z1 ... Z{Depth_points} "
                "for every time step)."
            )
        # an extra layer means the files and Depth_points do not match
        extra = sorted(layers[t] - expected)
        if extra:
            raise ValueError(
                f"{prefix} files: time step T{t} has depth layer(s) "
                f"{', '.join(f'Z{z}' for z in extra)} beyond Depth_points = {Depth_points} "
                f"(expected exactly Z1 ... Z{Depth_points} for every time step)."
            )
    return last + 1


# =============================================================================
# CrisisComponent — main class
# =============================================================================
# this is the main class — it implements CRISIS as a native Landlab component

class CrisisComponent(Component):
    """Coupled rainfall-induced and coseismic landslide hazard model (CRISIS).

    Implements a pseudo-3D slope stability analysis to predict landslide polygons
    (with individual location, area, volume) over a Landlab ``RasterModelGrid``.
    Elevation, slope, D8 flow direction and drainage area are read from Landlab
    grid fields. Time-varying pore pressure heads and volumetric water contents
    (both optional) are read from HDF5 files (one per time step and depth
    layer, from an external hydrological model) or from grid fields written by
    a coupled Landlab component.

    Pore pressure heads and water contents have one value per node and depth
    layer. On the grid they are stored as node fields of shape
    ``(number_of_nodes, Depth_points)``: row ``i`` of the field is node ``i``
    (Landlab order, row 0 of the grid = southern edge) and column ``k`` is
    depth layer ``k`` (shallow to deep). Internally they are viewed as
    ``(nrows, ncols, Depth_points)`` arrays by a plain row-major reshape,
    which needs no flip because node ``i`` is grid cell
    ``(i // ncols, i % ncols)``. Only the HDF5 files, which store their rows
    from north to south, are flipped on load.

    The component delegates all physics to ``_crisis_core.py`` and conforms
    natively to the Landlab component API.

    Parameters
    ----------
    grid : RasterModelGrid
        Landlab raster grid.  Must carry elevation, slope, D8 flow direction
        and drainage area fields before the component is instantiated.
    params : dict
        CRISIS model parameters.  All keys and their meanings are documented
        in ``CRISIS_Manual.pdf`` (Section 7).  The most commonly
        used keys are listed below under *Key Parameters*.
    slope, drainage_area, flow_direction : str, optional
        Names of the node fields holding the slope (degrees), the drainage
        area (m²) and the D8 receiver nodes. Default to
        ``"topographic__slope"``, ``"drainage_area"`` and
        ``"flow__receiver_node"``.
    log_callback : callable, optional
        Receives progress strings during the simulation.  Defaults to
        ``print``.
    pressure_head, water_content : str, optional
        Names of the node fields holding the pore pressure heads and the
        volumetric water contents (see *Key Parameters* below). Default to
        ``"soil_water__pressure_head"`` and ``"soil_water__volume_fraction"``.

    See ``__init__`` for the full description of each argument.

    Key Parameters
    --------------
    The most commonly used keys (``CRISIS_Manual.pdf``, Section 7, lists all
    keys, including the earthquake and fault ones). No key has a default
    value: :func:`check_params` lists which keys each setting requires, and
    the component stops with a ``ValueError`` naming every missing one.

    Gamma : float
        Unit weight of soil (kN/m³).
    z_min : float
        Depth of the shallowest depth layer (m). This layer is checked for
        failure like every other layer, unless ``z_min = 0`` (the ground
        surface).
    z_max : float
        Deepest depth layer (m).
    Depth_points : int
        Number of depth discretization layers.
    Min_LD_area : float
        Minimum accepted landslide planar area (m²).
    Triggering_Indicator : {1, 2}
        1 = rainfall only; 2 = rainfall + seismic (then the earthquake keys
        are required too, see :func:`check_params`).
    Hydrological_model_indicator : {1, 2, 3}
        Source of the pore pressure heads:
        1 = constant Ru pore pressure coefficient (steady state);
        2 = time-varying pore pressure head HDF5 files (from an external
        hydrological model);
        3 = time-varying pore pressure heads read from the
        ``soil_water__pressure_head`` grid field (written by a coupled
        Landlab component before each time step).
    Ru : float
        Constant pore pressure coefficient. Required when
        ``Hydrological_model_indicator == 1``.
    hydro_files_dir : str or Path
        Folder with the pore pressure head files ``Pressure_head_T{t}Z{z}.h5``.
        Required when ``Hydrological_model_indicator == 2``.
    Strength_indicator : {1, 2}
        1 = uniform cohesion/friction; 2 = spatially variable.
    c, Phi : float
        Uniform effective cohesion (kPa) and friction angle (degrees).
        Required when ``Strength_indicator == 1``.
    cohesion_3d_file, phi_3d_file : str or Path
        HDF5 files of the effective cohesion (kPa) and friction angle
        (degrees), shape ``(nrows, ncols, Depth_points)``. Required when
        ``Strength_indicator == 2``.
    Strength_variation_with_theta : {1, 2, 3}
        Source of the volumetric water content θ used for Bishop's χ:
        1 = fully saturated (Bishop χ = 1), no water content needed;
        2 = χ = (θ − θr) / (θs − θr), clipped to [0, 1], with θ read from
        the water content HDF5 files;
        3 = the same, with θ read from the ``soil_water__volume_fraction``
        grid field (written by a coupled Landlab component).
        Required when ``Hydrological_model_indicator`` is 2 or 3; not needed,
        and ignored (χ = 1), when ``Hydrological_model_indicator == 1``.
    theta_files_dir : str or Path
        Folder with the water content files ``Theta_T{t}Z{z}.h5``. Required
        when ``Strength_variation_with_theta == 2``.
    top_or_bottom_indicator : {1, 2}
        1 = report shallowest failure depth; 2 = deepest.
    Theta_s, Theta_r : float
        Saturated and residual volumetric water contents (−). Required when
        ``Strength_variation_with_theta`` is 2 or 3.
    Coordinate_Reference_System : str
        EPSG string for shapefile export (e.g. ``"EPSG:32610"``).
    dem_file : str, Path or file-like
        Required. Esri ASCII grid (``.asc``) of the DEM loaded onto the grid:
        a path, or an open text stream holding its content. Its header
        (``ncols``, ``nrows``, ``cellsize``, lower-left corner or center) is
        used to compute the easting / northing of every cell center, and its
        ``cellsize`` is the cell size (m) used for landslide areas, volumes
        and shapefile polygons (read automatically from the header). ``nrows``
        and ``ncols`` must match ``grid.shape``, and ``cellsize`` must equal
        ``grid.dx`` and ``grid.dy``; this is checked.
    output_dir : str or Path
        Folder the output files (shapefiles and Excel workbooks) are written
        into; created if it does not exist.

    Pore pressure heads and water contents given as grid fields (either
    indicator = 3) must be node fields of shape
    ``(number_of_nodes, Depth_points)`` holding floating-point values (any
    float dtype, e.g. float64 or float32, is used as it is).

    The number of time steps is not a parameter: it is 1 with a constant Ru
    (``Hydrological_model_indicator = 1``), the number of time steps of the
    HDF5 files found in ``hydro_files_dir`` (indicator 2), and, with grid
    fields written by a coupled component (indicator 3), the number of times
    :meth:`run_one_step` is called (see :attr:`number_of_time_steps`).

    Attributes
    ----------
    landslide_fid : ndarray, shape (nrows, ncols)
        Cumulative landslide-ID raster.  NaN where no failure has occurred.
        Landlab row order (row 0 = southern edge), like the grid fields.
    failed_cells : ndarray of bool, shape (nrows, ncols)
        True at every cell that has failed in any time step. Landlab row
        order (row 0 = southern edge).
    state : dict
        Full output accumulator (per-step and cumulative landslide results).
    current_time_index : int
        Index of the *next* time step that will be run by ``run_one_step``.
    number_of_time_steps : int or None
        Number of time steps of the run: 1 with a constant Ru, the number of
        time steps of the HDF5 files with indicator 2, ``None`` (no fixed
        number) with indicator 3.

    """

    # ── Component identity ────────────────────────────────────────────────────

    # the name Landlab shows for this component in error messages and docs
    _name = "CrisisComponent"

    # inputs and outputs must be in the SI units declared in _info
    _unit_agnostic = False

    # =============================================================================
    # Field declarations
    # =============================================================================
    # _info declares every grid field this component reads ("intent": "in") or
    # writes ("intent": "out"). Landlab validates required fields at construction
    # time and uses this dict to auto-generate documentation.

    # dictionary listing every input/output field this component needs or makes
    _info = {

        # ── Topographic inputs ────────────────────────────────────────────────
        # All four must exist on the grid before CrisisComponent is instantiated.
        # Slope, D8 flow direction and drainage area are computed beforehand with
        # Landlab (e.g. calc_slope_at_node and a D8 FlowAccumulator). They are
        # marked optional here only because the component accepts custom field
        # names for them; __init__ raises a KeyError if any of them is missing.

        # ground surface height at each grid point, from the DEM
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Surface elevation (DEM).  No-data cells: NaN or values below -900.",
        },

        # slope angle at each grid point, in degrees
        "topographic__slope": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "degrees",
            "mapping": "node",
            "doc": (
                "Slope angle in degrees, computed with Landlab, e.g. "
                "np.degrees(grid.calc_slope_at_node('topographic__elevation'))."
            ),
        },

        # which way water flows from each cell (D8 direction)
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Receiver-node index for each node (D8 routing), e.g. from a Landlab "
                "FlowAccumulator using FlowDirectorD8."
            ),
        },

        "drainage_area": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m^2",
            "mapping": "node",
            "doc": (
                "Upstream contributing area in m², e.g. from a Landlab FlowAccumulator "
                "using FlowDirectorD8."
            ),
        },

        # ── Hydrology (input or output, depending on the indicators) ──────────
        # One value per node and depth layer: shape (number_of_nodes, Depth_points),
        # nodes in Landlab order (row 0 = south), layers ordered shallow to deep.

        # pore pressure head at each node and depth layer
        "soil_water__pressure_head": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": (
                "Pore pressure head at each node and depth layer, shape "
                "(number_of_nodes, Depth_points), layers ordered shallow to deep. "
                "Written by the component when Hydrological_model_indicator is 1 "
                "(from Ru) or 2 (from the HDF5 files); read, and required, when it "
                "is 3."
            ),
        },

        # volumetric water content at each node and depth layer
        "soil_water__volume_fraction": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Volumetric water content at each node and depth layer, shape "
                "(number_of_nodes, Depth_points), layers ordered shallow to deep. "
                "Written by the component when Strength_variation_with_theta is 2 "
                "(from the HDF5 files); read, and required, when it is 3. Not used "
                "when it is 1 or when Hydrological_model_indicator is 1."
            ),
        },

        # ── Outputs ───────────────────────────────────────────────────────────
        # Created on the grid by __init__ and updated after every run_one_step call.

        # ID number telling which landslide (if any) each cell belongs to
        "landslide__failure_id_accumulated": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Cumulative landslide-ID raster. Each failed cell carries the integer ID "
                "of the landslide it belongs to. NaN where no failure has occurred."
            ),
        },

        # simple true/false flag: has this cell ever failed
        "landslide__failed_accumulated": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "True at nodes that have failed in any time step (cumulative).",
        },

        # ID of the landslide each cell belongs to, for landslides new in this step
        "landslide__failure_id": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Landslide ID local to the current time step (1..N) at every cell of "
                "each new landslide; 0 elsewhere."
            ),
        },

        # failure-surface depth at each cell of a new landslide
        "failure__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": (
                "Failure-surface depth at every cell of each landslide that is new in "
                "the current time step; 0 elsewhere."
            ),
        },

        # initiation (triggering) cell of each new landslide
        "triggering__cells": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "True at the triggering cell of each landslide that is new in the current time step.",
        },
    }

    # =============================================================================
    # Constructor  (__init__)
    # =============================================================================
    # this runs once when you create a CrisisComponent — it sets everything up

    def __init__(
            self,
            grid,
            params,
            slope = 'topographic__slope',
            drainage_area = 'drainage_area',
            flow_direction= "flow__receiver_node",
            log_callback=print,
            pressure_head="soil_water__pressure_head",
            water_content="soil_water__volume_fraction"):
        """Initialize the CRISIS component and pre-load all static data.

        Parameters
        ----------
        grid : RasterModelGrid
            Landlab raster grid carrying the required input fields.
        params : dict
            CRISIS model parameters (see class docstring).
        slope : str, optional
            Node field holding slope in degrees. Defaults to
            ``"topographic__slope"``.
        drainage_area : str, optional
            Node field holding drainage area in m². Defaults to
            ``"drainage_area"``.
        flow_direction : str, optional
            Node field holding D8 receiver nodes. Defaults to
            ``"flow__receiver_node"``.
        log_callback : callable, optional
            Progress-message sink.  Defaults to ``print``.
        pressure_head : str, optional
            Node field holding the pore pressure heads (m), shape
            ``(number_of_nodes, Depth_points)``. Written by the component when
            ``Hydrological_model_indicator`` is 1 or 2; must exist on the
            grid when it is 3. Defaults to ``"soil_water__pressure_head"``.
        water_content : str, optional
            Node field holding the volumetric water contents (-), shape
            ``(number_of_nodes, Depth_points)``. Written by the component when
            ``Strength_variation_with_theta`` is 2; must exist on the grid
            when it is 3. Defaults to ``"soil_water__volume_fraction"``.

        Raises
        ------
        TypeError
            If *grid* is not a ``RasterModelGrid``.
        KeyError
            If the slope, drainage-area or flow-direction field does not
            exist on the grid, or if ``Hydrological_model_indicator == 3``
            (``Strength_variation_with_theta == 3``) and the pore pressure head
            (water content) field does not exist on the grid.
        TypeError
            If a pore pressure head / water content field used by the component
            does not hold floating-point values (any float dtype, e.g.
            float32, is accepted and used as it is).
        TypeError
            If *params* is not a dict.
        ValueError
            If a key the settings need is missing from ``params``, or a choice
            key has a value that is not allowed (checked first, by
            :func:`check_params`, which lists every such key); if the header of
            ``dem_file`` does not match the grid (number of rows/columns, or a
            ``cellsize`` different from the grid spacing ``grid.dx`` /
            ``grid.dy``); if a pore pressure head / water content field does
            not have the shape ``(number_of_nodes, Depth_points)``; or, with
            HDF5 files (indicator 2), if a time step has depth layers beyond
            ``Depth_points``, if the water content files do not have the same
            time steps as the pore pressure head files, or if
            ``Earthquake_time_point`` is not one of the time steps of the run.
        FileNotFoundError
            If a file given in ``params`` (e.g. ``dem_file``) cannot be found,
            or, with HDF5 files (indicator 2), if the folder holds no pore
            pressure head (water content) files, a time step is missing (a
            gap), or a time step is missing one of its depth layers.

        Examples
        --------
        Create a 3x100 grid representing a uniformly inclined slope, and
        instantiate CrisisComponent on it (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.current_time_index
        0

        The easting / northing of every cell center is computed from the
        header of ``dem_file``. A header giving the outer corner of the
        lower-left cell (``xllcorner`` / ``yllcorner``) is shifted by half a
        cell; row 0 is the southern-most row, so northing increases with the
        row index.

        >>> params["dem_file"] = io.StringIO("ncols 100\\nnrows 3\\nxllcorner 500000.0\\n"
        ...                                  "yllcorner 4000000.0\\ncellsize 1.0\\n")
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis._x[:3].tolist()
        [500000.5, 500001.5, 500002.5]
        >>> crisis._y.tolist()
        [4000000.5, 4000001.5, 4000002.5]

        A header that does not match the grid is rejected.

        >>> params["dem_file"] = io.StringIO("ncols 50\\nnrows 3\\nxllcorner 500000.0\\n"
        ...                                  "yllcorner 4000000.0\\ncellsize 1.0\\n")
        >>> CrisisComponent(grid, params)
        Traceback (most recent call last):
        ...
        ValueError: dem_file has 3 rows x 50 columns, but the grid has 3 x 100.

        The cell size is read automatically from the ``cellsize`` line of the
        ``dem_file`` header, and must be the grid spacing.

        >>> crisis._Cell_size
        1.0
        >>> params["dem_file"] = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\n"
        ...                                  "yllcenter 0.0\\ncellsize 2.0\\n")
        >>> CrisisComponent(grid, params)
        Traceback (most recent call last):
        ...
        ValueError: dem_file cellsize (2.0) does not match the grid spacing (dx = 1.0, dy = 1.0).

        With a constant Ru (``Hydrological_model_indicator = 1``) there is no
        water content: Bishop χ = 1, and ``Strength_variation_with_theta`` is
        ignored (no water content file is read, no water content field is
        created).

        >>> params["dem_file"] = dem_file
        >>> params["Strength_variation_with_theta"] = 2
        >>> crisis = CrisisComponent(grid, params)
        Hydrological_model_indicator = 1 (constant Ru): Bishop chi = 1, Strength_variation_with_theta = 2 is ignored.
        >>> "soil_water__volume_fraction" in grid.at_node
        False

        A constant Ru is steady state: the run has a single time step (the
        number of time steps is never a parameter).

        >>> params["Strength_variation_with_theta"] = 1
        >>> params["dem_file"] = dem_file
        >>> CrisisComponent(grid, params).number_of_time_steps
        1

        Before anything else, :func:`check_params` checks that every key the
        settings need is given, with an allowed value (no parameter has a
        default value). Every indicator starts at 1; the hydrological model
        indicator must be 1, 2 or 3:

        >>> params["Hydrological_model_indicator"] = 0
        >>> CrisisComponent(grid, params)
        Traceback (most recent call last):
        ...
        ValueError: Missing or invalid entries in params:
          - Hydrological_model_indicator must be 1, 2 or 3, not 0

        and a key needed by the settings, such as the 3D strength files with
        ``Strength_indicator = 2``, must be given:

        >>> params["Hydrological_model_indicator"] = 1
        >>> params["Strength_indicator"] = 2
        >>> CrisisComponent(grid, params)
        Traceback (most recent call last):
        ...
        ValueError: Missing or invalid entries in params:
          - cohesion_3d_file, phi_3d_file are missing (required when Strength_indicator = 2)
        """

        # ── Parameter check ───────────────────────────────────────────────────
        # stop before anything is done if a key this run needs is missing or
        # has a value that is not allowed (no parameter has a default value)
        check_params(params)

        # ── Grid type check ───────────────────────────────────────────────────
        # stop early with a clear error if someone passes the wrong grid type
        if not isinstance(grid, RasterModelGrid):
            raise TypeError("CrisisComponent requires a RasterModelGrid.")

        # ── Layered hydrology fields may hold any floating-point type ─────────
        # Landlab requires an existing field to have the dtype declared in _info
        # (float64). The pore pressure head / water content fields are used in place
        # (read as float32 copies, written by assignment), so any float dtype
        # (e.g. float32 from a coupled component) works. Accept the field's own
        # dtype instead of replacing the field with a float64 copy: a replaced
        # field would leave the coupled component writing into an array the grid no
        # longer holds. Integer fields are rejected by _layer_field_to_3d.
        info = None
        for name in ("soil_water__pressure_head", "soil_water__volume_fraction"):
            if grid.has_field(name, at="node"):
                if info is None:
                    info = {key: dict(entry) for key, entry in self._info.items()}
                info[name]["dtype"] = grid.at_node[name].dtype
        if info is not None:
            self._info = info      # this instance only; the class declaration is unchanged

        # ── Base class initialization ─────────────────────────────────────────
        # super().__init__ registers the grid and validates all _info "in" fields

        # hand the grid to Landlab's base class so it checks our required fields
        super().__init__(grid)

        # ── Store params and logger ───────────────────────────────────────────

        # remember the parameters dict for later use
        self._params = params
        # remember how to print/report progress messages
        self._log    = log_callback

        # ── Simulation clock and grid shape ───────────────────────────────────

        # internal counter: the time step run_one_step will process next
        self._current_time_index = 0

        # store the grid's row and column counts for convenience
        self._rows, self._clms = grid.shape

        # ── Elevation ─────────────────────────────────────────────────────────
        # make our own private copy of the elevation grid
        self._Elev_Mat = _node_to_2d(grid, "topographic__elevation").copy().astype(np.float32)
        # turn "no data" elevation placeholders into proper NaN (missing value)
        self._Elev_Mat[self._Elev_Mat < -900] = np.nan

        # ── Slope, drainage area and D8 flow direction (Landlab grid fields) ──
        # These must already be on the grid, computed beforehand with Landlab
        # (e.g. calc_slope_at_node and a D8 FlowAccumulator). They are never
        # read from external files, and the component does not compute them.
        # Landlab drainage area is already in m², so no unit conversion is needed.

        # remember the field names so run_one_step re-reads the same fields
        self._slope_name          = slope
        self._drainage_area_name  = drainage_area
        self._flow_direction_name = flow_direction

        # stop early with a clear error if any of the three fields is missing
        missing = [name for name in (slope, drainage_area, flow_direction)
                   if not grid.has_field(name, at="node")]
        if missing:
            raise KeyError(
                f"CrisisComponent requires the grid field(s) {missing}. Compute slope "
                "(degrees), D8 flow direction and drainage area with Landlab before "
                "creating the component, e.g. "
                "np.degrees(grid.calc_slope_at_node('topographic__elevation')) and "
                "FlowAccumulator(grid, flow_director='FlowDirectorD8')."
            )

        # point at the grid fields (re-read at every run_one_step)
        self._Slope_Mat = grid.at_node[slope]
        self._Area_Mat  = grid.at_node[drainage_area]

        # convert Landlab receiver nodes to the D8 codes (1-128) used by the core
        self._FD_Mat = self._build_fd_mat_from_landlab(grid, flow_direction)

        # ── Coordinates (easting / northing) and cell size ────────────────────
        # Computed from the header of the ASCII DEM (params["dem_file"]) so the
        # output shapefiles are placed in real-world projected coordinates.
        # The header must describe the same grid as the Landlab grid.

        dem_header = read_ascii_header(params["dem_file"])
        if (dem_header["nrows"], dem_header["ncols"]) != (self._rows, self._clms):
            raise ValueError(
                f"dem_file has {dem_header['nrows']} rows x {dem_header['ncols']} columns, "
                f"but the grid has {self._rows} x {self._clms}."
            )

        # the cell size is read automatically from the DEM header's cellsize;
        # areas and volumes use it, so it must also be the grid spacing
        self._Cell_size = float(dem_header["cellsize"])
        if not (np.isclose(grid.dx, self._Cell_size) and np.isclose(grid.dy, self._Cell_size)):
            raise ValueError(
                f"dem_file cellsize ({self._Cell_size}) does not match the grid spacing "
                f"(dx = {grid.dx}, dy = {grid.dy})."
            )

        # center of the lower-left cell, then one coordinate per column / row
        # (row 0 = southern-most row, matching the Landlab grid)
        xmin, ymin = lower_left_cell_center(dem_header)
        self._x, self._y = compute_coordinate_vectors(
            self._rows, self._clms, xmin, ymin, dem_header["cellsize"]
        )

        # ── Depth discretization vector ───────────────────────────────────────
        # FS and Newmark displacement are evaluated at each depth for every cell

        # build the list of soil depths, shallow to deep (every layer is checked for
        # failure, except a layer at depth 0, which is the ground surface)
        self._Depths_vector = np.linspace(
            params["z_min"], params["z_max"],
            params["Depth_points"], dtype=np.float32,
        )

        # ── Hydrology indicators and their grid fields ────────────────────────
        # Pore pressure heads: 1 = constant Ru, 2 = HDF5 files, 3 = grid field.
        # Water contents (Bishop χ): 1 = χ = 1, 2 = HDF5 files, 3 = grid field;
        # not used with a constant Ru (χ = 1).

        # (both indicators were checked by check_params)
        self._Hydro_indicator = params["Hydrological_model_indicator"]
        if self._Hydro_indicator == 1:
            # a constant Ru has no water content: Bishop χ = 1, so
            # Strength_variation_with_theta is not needed, and ignored if given
            given_theta = params.get("Strength_variation_with_theta")
            if given_theta not in (None, 1):
                self._log(
                    "Hydrological_model_indicator = 1 (constant Ru): Bishop chi = 1, "
                    f"Strength_variation_with_theta = {given_theta} is ignored."
                )
            self._Strength_variation_with_theta = 1
        else:
            self._Strength_variation_with_theta = params["Strength_variation_with_theta"]

        # ── Number of time steps (never a parameter) ──────────────────────────
        # 1: a constant Ru is steady state, a single time step;
        # 2: the time steps of the HDF5 files (all steps found are run);
        # 3: not fixed, as many as run_one_step is called (None)
        if self._Hydro_indicator == 1:
            self._number_of_time_steps = 1
        elif self._Hydro_indicator == 2:
            self._number_of_time_steps = _count_time_steps(
                params["hydro_files_dir"], "Pressure_head", params["Depth_points"])
            if self._Strength_variation_with_theta == 2:
                theta_steps = _count_time_steps(
                    params["theta_files_dir"], "Theta", params["Depth_points"])
                if theta_steps != self._number_of_time_steps:
                    raise ValueError(
                        f"The water content files have {theta_steps} time step(s), but the "
                        f"pore pressure head files have {self._number_of_time_steps}: both must "
                        "cover the same time steps."
                    )
            self._log(
                f"Hydrological_model_indicator = 2: {self._number_of_time_steps} time step(s) "
                f"found in the HDF5 files (T0 ... T{self._number_of_time_steps - 1})."
            )
        else:
            self._number_of_time_steps = None
            # water content files with a coupled component: check them now too
            # (gaps, missing or extra layers); they cover this many time steps
            if self._Strength_variation_with_theta == 2:
                theta_steps = _count_time_steps(
                    params["theta_files_dir"], "Theta", params["Depth_points"])
                self._log(
                    f"Strength_variation_with_theta = 2: {theta_steps} time step(s) found in "
                    f"the water content files (T0 ... T{theta_steps - 1})."
                )

        # the earthquake must happen at one of the time steps of the run
        eq_t = params.get("Earthquake_time_point")
        if (params["Triggering_Indicator"] == 2 and self._number_of_time_steps is not None
                and not 0 <= eq_t < self._number_of_time_steps):
            raise ValueError(
                f"Earthquake_time_point = {eq_t} is not one of the {self._number_of_time_steps} "
                f"time step(s) of the run (0 ... {self._number_of_time_steps - 1})."
            )

        # remember the field names so every time step uses the same fields
        self._pressure_head_name = pressure_head
        self._water_content_name = water_content

        # fields supplied by a coupled component must already be on the grid
        # (_layer_field_to_3d raises if a field is missing or has the wrong shape)
        if self._Hydro_indicator == 3:
            self._layer_field_to_3d(pressure_head)
        if self._Strength_variation_with_theta == 3:
            self._layer_field_to_3d(water_content)

        # fields written by this component: create them (or check an existing
        # one has the right shape); NaN until the first time step fills them
        if self._Hydro_indicator in (1, 2):
            self._ensure_layer_field(pressure_head)
        if self._Strength_variation_with_theta == 2:
            self._ensure_layer_field(water_content)

        # a constant Ru gives the same pressure head at every time step:
        # Ru × (γ / γ_w) × z at every node (the core computes the same values)
        if self._Hydro_indicator == 1:
            Ru_head = (params["Ru"] * params["Gamma"] / Gamma_W) * self._Depths_vector
            self.grid.at_node[pressure_head][:] = Ru_head[np.newaxis, :]

        # ── Static seismic arrays (Triggering_Indicator == 2 only) ───────────
        # PGA and PGV are spatial maps loaded once from HDF5 files (flipped on
        # load by load_h5, like all raster files); rupture distance computed
        # once from the fault plane geometry

        # placeholders for earthquake-related arrays, filled in only if needed
        self._PGA_2D_array     = None
        self._PGV_2D_array     = None
        self._Rupture_2D_array = None

        # only load earthquake data if this run includes seismic triggering
        # (check_params made sure every key this model needs is given)
        if params["Triggering_Indicator"] == 2:

            # load peak ground acceleration map (g) from its HDF5 file
            self._PGA_2D_array = load_h5(params["pga_file"])

            # which displacement model to use for slope failure under shaking,
            # and (Rathje models 2 and 3 only) which Rathje method
            seismic_model = params["Seismic_displacement_model"]
            rathje_method = params["Rathje_method"] if seismic_model in (2, 3) else None

            # PGV is needed by Rathje's PGAV method (models 2/3, method 2) and
            # by Bray & Macedo (model 4)
            if seismic_model == 4 or rathje_method == 2:
                # load peak ground velocity map (cm/s) from its HDF5 file
                self._PGV_2D_array = load_h5(params["pgv_file"])

            # the Rathje flexible and Bray & Macedo models (3, 4) also need the
            # distance from every cell to the fault rupture
            if seismic_model in (3, 4):
                self._Rupture_2D_array = self._compute_rupture_distances(params)

        # ── Spatially variable strength arrays (Strength_indicator == 2 only) ─
        # 3D arrays (rows × cols × Depth_points) loaded once from HDF5

        # placeholders for soil-strength maps that vary by location and depth
        self._Cohesion_3D_array = None
        self._Phi_3D_array      = None

        # only load these if strength is set to vary spatially
        if params["Strength_indicator"] == 2:
            # load cohesion (kPa) and friction angle (degrees) for every cell and depth
            self._Cohesion_3D_array = load_h5(params["cohesion_3d_file"])
            self._Phi_3D_array      = load_h5(params["phi_3d_file"])

        # ── Upslope flow map ──────────────────────────────────────────────────
        # Pre-build {cell: [upslope neighbors]} dict so BFS delineation does not
        # recompute it at every time step (most expensive setup step for large grids)

        # build a lookup of "which cells flow into this cell" once, up front
        self._upslope_map = _build_upslope_map(self._FD_Mat)

        # ── Output directory ──────────────────────────────────────────────────
        # Shapefiles are written to output_dir after each time step; the
        # Excel workbooks are written there by write_outputs() after a run

        # the output folder, created now if it doesn't already exist
        self._output_dir = Path(params["output_dir"])
        self._output_dir.mkdir(parents=True, exist_ok=True)

        # ── State accumulator ─────────────────────────────────────────────────
        # All arrays that change across time steps live here, in the structure
        # expected by the Excel writers in _crisis_core (_save_individual_excel,
        # _save_summary_excel).
        #   "all_times" rasters  — grid-shaped, updated in place each step
        #   "with_time" vectors  — one entry appended per step, written to Excel
        #   (they grow with every run_one_step call, so they always hold exactly
        #   the time steps run so far)

        # this dictionary holds every result that builds up as the simulation runs
        self._state = {

            # ── Cumulative rasters ────────────────────────────────────────────
            # grid of which cells have ever failed (true/false)
            "Fail_Index_all_times":  np.zeros((self._rows, self._clms), dtype=bool),
            # grid of landslide ID numbers (NaN = no landslide there)
            "FID_all_times":         np.full((self._rows, self._clms), np.nan, dtype=np.float32),

            # ── Per-step totals (current = new this step only) ────────────────
            # total new landslide area at each time step
            "Total_Current_Area_vector_with_time":      np.zeros(0, dtype=np.float32),
            # total new landslide volume at each time step
            "Total_Current_Volume_vector_with_time":    np.zeros(0, dtype=np.float32),
            # count of new landslides at each time step
            "Total_Current_NL_vector_with_time":        np.zeros(0, dtype=np.float32),

            # ── Per-step running totals (cumulative = all steps so far) ───────
            # running total area of all landslides so far, per time step
            "Total_Cumulative_Area_vector_with_time":   np.zeros(0, dtype=np.float32),
            # running total volume of all landslides so far, per time step
            "Total_Cumulative_Volume_vector_with_time": np.zeros(0, dtype=np.float32),
            # running total count of landslides so far, per time step
            "Total_Cumulative_NL_vector_with_time":     np.zeros(0, dtype=np.float32),

            # ── Per-step lists — new (current) landslide properties ───────────
            # list of cell locations that triggered new landslides, one entry per step
            "Individual_Current_triggering_linear_indices_with_time": [],
            # list of failure depths for the new landslides, one entry per step
            "Individual_Current_triggering_depths_with_time":         [],
            # list of what caused the new failures (rain/quake), one entry per step
            "Individual_Current_triggering_sources_with_time":        [],
            # list of new landslide areas, one entry per step
            "Individual_Current_areas_with_time":                     [],
            # list of new landslide volumes, one entry per step
            "Individual_Current_volumes_with_time":                   [],

            # ── Per-step snapshots — cumulative landslide registry ────────────
            # snapshot of all triggering cell locations seen so far, per step
            "Individual_Cumulative_triggering_linear_indices_with_time": [],
            # snapshot of all failure depths seen so far, per step
            "Individual_Cumulative_triggering_depths_with_time":         [],
            # snapshot of all trigger sources seen so far, per step
            "Individual_Cumulative_triggering_sources_with_time":        [],
            # snapshot of all landslide areas seen so far, per step
            "Individual_Cumulative_areas_with_time":                     [],
            # snapshot of all landslide volumes seen so far, per step
            "Individual_Cumulative_volumes_with_time":                   [],

            # ── Running master registry — all unique triggering cells ─────────
            # master list of every unique landslide-triggering cell location
            "Individual_Cumulative_triggering_linear_indices": [],
            # master list of failure depths matching the above cells
            "Individual_Cumulative_triggering_depths":         [],
            # master list of trigger sources matching the above cells
            "Individual_Cumulative_triggering_sources":        [],
            # master list of volumes matching the above cells
            "Individual_Cumulative_volumes":                   [],
            # master list of areas matching the above cells
            "Individual_Cumulative_areas":                     [],
        }

        # ── Initialize output fields on the grid ──────────────────────────────

        # create the output fields on the Landlab grid so other components can read them
        self._init_output_fields()

        # ── Extract constant model parameters ─────────────────────────────────
        # These values never change during the simulation, so extract them once
        # at initialization rather than re-extracting them in every run_one_step call

        # ── General model settings ────────────────────────────────────────
        self._Gamma       = params["Gamma"]                    # soil unit weight (kN/m³)
        self._Gamma_W     = Gamma_W                            # unit weight of water (kN/m³)
        # (_Cell_size, the grid cell size in m, was read from the DEM header above)
        self._Min_LD_area = params["Min_LD_area"]              # minimum landslide area (m²)
        self._top_bot     = params["top_or_bottom_indicator"]  # shallowest (1) or deepest (2) failure

        # ── Triggering and strength settings ─────────────────────────────
        self._Strength_indicator            = params["Strength_indicator"]
        self._Triggering_Indicator          = params["Triggering_Indicator"]
        # (_Hydro_indicator and _Strength_variation_with_theta were checked and
        # stored above, with the hydrology fields)
        # a constant Ru is steady state: single-shapefile, no-time-step outputs
        self._steady_state                  = self._Hydro_indicator == 1
        self._CRS                           = params["Coordinate_Reference_System"]

        # ── Settings used only in some cases ──────────────────────────────
        # check_params made sure each one is given when the settings need it;
        # a key the settings do not use is None here and is never used
        self._c       = params.get("c")        # constant cohesion (kPa), Strength_indicator = 1
        self._Phi     = params.get("Phi")      # constant friction angle (degrees), Strength_indicator = 1
        self._Ru      = params.get("Ru")       # constant pore pressure coefficient, indicator 1
        self._Theta_s = params.get("Theta_s")  # saturated water content, theta indicator 2 or 3
        self._Theta_r = params.get("Theta_r")  # residual water content, theta indicator 2 or 3

        # ── Seismic settings (Triggering_Indicator = 2) ───────────────────
        self._Earthquake_time_point          = params.get("Earthquake_time_point")
        self._Seismic_displacement_threshold = params.get("Seismic_displacement_threshold")
        self._Seismic_displacement_model     = params.get("Seismic_displacement_model")
        self._Rathje_method                  = params.get("Rathje_method")    # models 2 and 3
        self._Ts                             = params.get("Ts")               # models 3 and 4
        self._M                              = params.get("M")

    # =============================================================================
    # Public API
    # =============================================================================

    # ─────────────────────────────────────────────────────────────────────────────
    # run_one_step — advance the simulation by one time step
    # ─────────────────────────────────────────────────────────────────────────────

    def run_one_step(self, dt=1.0):
        """Run the CRISIS forward model for a single time step.

        Advances the internal time-step counter by one each call.  The first
        call corresponds to time-step index 0, the second to index 1, and so
        on.  Time-varying HDF5 files (pore pressure heads, water contents) are
        loaded by index; pore pressure head / water content grid fields written
        by a coupled component (indicator = 3) are read as they are at the time
        of the call, so that component must update them before each call.

        Parameters
        ----------
        dt : float, optional
            Duration of this time step (any unit consistent with the model
            clock).  CRISIS does not use *dt* internally; it is accepted for
            compatibility with the standard Landlab ``run_one_step(dt)``
            interface.  Defaults to ``1.0``.

        Raises
        ------
        RuntimeError
            If called more times than the run has time steps (1 with a
            constant Ru, the number of time steps of the HDF5 files with
            indicator 2; no limit with indicator 3).
        FileNotFoundError
            If ``Hydrological_model_indicator == 2`` and a pore pressure head
            file, or ``Strength_variation_with_theta == 2`` and a water
            content file, of the current time step cannot be found.
        KeyError, TypeError
            If a field supplied by a coupled component (indicator = 3) was
            removed from the grid or no longer holds floating-point values.
        ValueError
            If ``Hydrological_model_indicator == 3`` or
            ``Strength_variation_with_theta == 3`` and the field supplied by
            the coupled component no longer has the shape
            ``(number_of_nodes, Depth_points)``.

        Examples
        --------
        Create a 3x100 grid representing a uniformly inclined slope, draining
        to a single outlet on the left edge, and run CRISIS on it.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()

        Set up and run CRISIS.

        >>> Ru = 0.8
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": Ru,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> np.unique(grid.at_node["landslide__failure_id"])
        array([0, 1, 2])
        """

        # ── Re-read topography from the grid fields ───────────────────────────
        # Elevation, slope and drainage area are re-read every step, so updates
        # made to these fields by other components are picked up (the D8 flow
        # direction is read once, in __init__)
        self._Elev_Mat = _node_to_2d(self._grid, "topographic__elevation").copy().astype(np.float32)
        self._Elev_Mat[self._Elev_Mat < -900] = np.nan

        self._Slope_Mat = self.get_topographic_slope()
        self._Slope_Mat = self._Slope_Mat[self._grid.nodes]
        self._Area_Mat = self.get_drainage_area()
        self._Area_Mat = self._Area_Mat[self._grid.nodes]

        # ── Time-step counter and bounds check ────────────────────────────────
        # dt is unused: the internal counter selects the time step (and its files)
        t = self._current_time_index

        # refuse to run past the time steps of the run (checked before the
        # reset, so the fields of the last step are kept)
        # (no limit with a coupled component, indicator 3)
        if self._number_of_time_steps is not None and t >= self._number_of_time_steps:
            raise RuntimeError(
                f"run_one_step called {t + 1} times, but the run has "
                f"{self._number_of_time_steps} time step(s): no time steps left."
            )

        # clear the per-step output fields of the previous step
        self._reset_output_fields()

        # ── Unpack parameters into local variables ────────────────────────────
        # shortcuts to params and grid size, so the code below is easier to read
        params = self._params
        rows   = self._rows
        clms   = self._clms

        # Local shortcuts to the parameters cached in __init__
        Gamma         = self._Gamma
        Gamma_W       = self._Gamma_W
        Cell_size     = self._Cell_size
        Min_LD_area   = self._Min_LD_area
        top_bot       = self._top_bot

        # ── Triggering and strength settings ──────────────────────────────────
        Strength_indicator            = self._Strength_indicator
        Triggering_Indicator          = self._Triggering_Indicator
        Strength_variation_with_theta = self._Strength_variation_with_theta
        Hydro_indicator               = self._Hydro_indicator
        CRS                           = self._CRS

        # Soil strength constants
        c       = self._c
        Phi     = self._Phi
        Ru      = self._Ru
        Theta_s = self._Theta_s
        Theta_r = self._Theta_r

        # Seismic parameters
        Earthquake_time_point          = self._Earthquake_time_point
        Seismic_displacement_threshold = self._Seismic_displacement_threshold
        Seismic_displacement_model     = self._Seismic_displacement_model
        Rathje_method                  = self._Rathje_method
        Ts                             = self._Ts
        M                              = self._M

        # Load hydrology
        Theta_3D, Pressures_3D = self._calc_hydrology(t, params)

        # ── Identify triggering cells ─────────────────────────────────────────
        # Compute FS and/or Newmark displacement for every cell and depth layer.
        # Returns flat arrays of linear indices, failure depths, and trigger sources.

        # run the core math to find which cells fail this step, and at what depth
        Trig_Cells, Trig_Depths, Trig_sources = _find_triggering_cells_vectorized(
            rows, clms,
            self._Elev_Mat, self._Slope_Mat, self._Depths_vector,
            Hydro_indicator, Pressures_3D, Theta_3D, Ru,
            Triggering_Indicator, Strength_indicator, c, Phi,
            Theta_s, Theta_r,
            self._Cohesion_3D_array, self._Phi_3D_array,
            Gamma_W, Gamma, t, Earthquake_time_point,
            self._PGA_2D_array,
            Seismic_displacement_threshold, M,
            Strength_variation_with_theta, top_bot,
            Seismic_displacement_model=Seismic_displacement_model,
            PGV_2D_array=self._PGV_2D_array,
            Ts=Ts,
            Rathje_method=Rathje_method,
            Rupture_2D_array=self._Rupture_2D_array,
        )

        # ── Pull cumulative failure arrays from state ─────────────────────────

        # get the running "failed cells so far" grid out of state
        Fail_Index = self._state["Fail_Index_all_times"]
        # get the running landslide-ID grid out of state
        FID        = self._state["FID_all_times"]

        # clear the per-step cell lists used to fill the per-step output fields
        self.LD_Z = []
        self.LD_Index = []
        self.Trig_Index_Accept = []

        # ── Pseudo-3D BFS landslide delineation ───────────────────────────────
        # Walk upslope from each triggering cell; cells that satisfy the geometric
        # failure criterion form one landslide body each.

        if len(Trig_Cells) > 0:
            # walk uphill from each failing cell to build full landslide shapes
            (Vol, Area, NL,
             Trig_Index_Accept, _,
             Fail_Index, Trig_Depth_Accept,
             Trig_source_Accept, FID,
             LD_Z,
             LD_Index,
             FID_Index) = Determine_Pseudo_3D_landslides(
                self._Elev_Mat, Trig_Cells, self._Area_Mat, self._Slope_Mat,
                self._FD_Mat, Cell_size, Trig_Depths, Trig_sources,
                Fail_Index, FID, Min_LD_area,
                upslope_map=self._upslope_map,
            )

            self.LD_Z = LD_Z
            self.LD_Index = LD_Index
            self.Trig_Index_Accept = Trig_Index_Accept
            self._FID_Index = FID_Index
        else:
            # nothing triggered this step — use empty results
            Vol  = np.array([], dtype=np.float32)
            Area = np.array([], dtype=np.float32)
            NL   = 0
            Trig_Index_Accept  = np.array([], dtype=np.int64)
            Trig_Depth_Accept  = []
            Trig_source_Accept = []

        # ── Write updated failure arrays back into state ───────────────────────
        # save the updated failed-cells grid back into state
        self._state["Fail_Index_all_times"] = Fail_Index
        # save the updated landslide-ID grid back into state
        self._state["FID_all_times"]        = FID

        # ── Update cumulative registry ────────────────────────────────────────
        # Merge newly accepted landslides into the running cumulative lists and
        # append per-timestep snapshots for Excel output

        # fold this step's new landslides into the running master list
        self._update_cumulative_registry(
            t, Trig_Index_Accept, Trig_Depth_Accept,
            Trig_source_Accept, Vol, Area, Cell_size,
        )

        # ── Compute summary scalars ───────────────────────────────────────────
        # Aggregate statistics for this time step; stored in pre-allocated vectors

        # if nothing has ever failed, totals are zero
        if np.all(np.isnan(FID)):
            cum_NL, cum_area_total, cum_vol_total = 0, 0.0, 0.0
        else:
            # otherwise compute the running totals: landslide count, area, volume
            cum_NL         = int(np.nanmax(FID))
            cum_area_total = np.nansum(Fail_Index) * Cell_size ** 2
            cum_vol_total  = np.nansum(self._state["Individual_Cumulative_volumes"])

        # append this step's totals (new landslides, and running cumulative totals)
        for key, value in (
                ("Total_Current_Area_vector_with_time",      np.sum(Area)),
                ("Total_Current_Volume_vector_with_time",    np.sum(Vol)),
                ("Total_Current_NL_vector_with_time",        NL),
                ("Total_Cumulative_Area_vector_with_time",   cum_area_total),
                ("Total_Cumulative_NL_vector_with_time",     cum_NL),
                ("Total_Cumulative_Volume_vector_with_time", cum_vol_total)):
            self._state[key] = np.append(self._state[key], np.float32(value))

        # ── Shapefile export ──────────────────────────────────────────────────
        # One shapefile per time step with all cumulative landslide polygons
        # (nothing is written until at least one landslide has formed); a
        # steady-state run (constant Ru) writes a single Predicted_landslides.shp

        # write out a shapefile of all landslides so far, for this time step
        Modeled_landslides_shapefiles(
            FID, Cell_size, self._x, self._y, t, CRS, str(self._output_dir),
            file_name="Predicted_landslides.shp" if self._steady_state else None,
        )

        # ── Update Landlab node fields ────────────────────────────────────────
        # Flatten 2D state arrays back to 1D for the Landlab grid

        # push the updated 2D results back onto the flat Landlab grid fields
        self._update_output_fields()

        # ── Advance clock ─────────────────────────────────────────────────────

        # move the internal counter to the next time step
        self._current_time_index += 1

    # ─────────────────────────────────────────────────────────────────────────────
    # write_outputs — write the Excel workbooks
    # ─────────────────────────────────────────────────────────────────────────────

    def write_outputs(self):
        """Write the end-of-run Excel workbooks to ``output_dir``.

        Writes ``Individual_landslides_with_time.xlsx`` and
        ``Total_landslides_with_time.xlsx`` from the results accumulated so
        far, with one column (individual) or one row (totals) per time step
        run, labeled t=0, t=1, ... With a coupled component (indicator 3) and
        an earthquake (``Triggering_Indicator = 2``), a notice is printed if
        ``Earthquake_time_point`` was not reached. A steady-state run (constant Ru,
        ``Hydrological_model_indicator = 1``) has no time steps: it writes
        ``Individual_landslides.xlsx`` (one column of landslides per sheet)
        and ``Total_landslides.xlsx`` (a single row) instead. The
        ``Triggering_Cells`` indices are row-major indices in the row order of
        the ASCII DEM file (north-up, row 0 = northern edge):
        ``row * ncols + col``, counting rows from the top. The per-step
        shapefiles are written by :meth:`run_one_step`.

        Examples
        --------
        Write the workbooks of a steady-state run (constant Ru) on a small
        grid that is steep only in its northern half, and check that every
        triggering cell lies in the top (northern) rows of the DEM file.

        >>> import os, tempfile
        >>> import numpy as np, pandas as pd
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> grid = RasterModelGrid((6, 20), xy_spacing=1.0)
        >>> rate = np.where(grid.y_of_node >= 3, 0.9, 0.1)   # steep north half
        >>> _ = grid.add_field("topographic__elevation", rate * grid.x_of_node, at="node")
        >>> _ = grid.add_field("topographic__slope",
        ...     np.degrees(grid.calc_slope_at_node("topographic__elevation")), at="node")
        >>> FlowAccumulator(grid, flow_director="FlowDirectorD8").run_one_step()
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 20\\nnrows 6\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.0, "z_max": 2.0,
        ...     "Depth_points": 9, "Min_LD_area": 4.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.4,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 2.0, "Phi": 30.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> crisis.write_outputs()
        >>> sorted(f for f in os.listdir(folder) if f.endswith((".xlsx", ".shp")))
        ['Individual_landslides.xlsx', 'Predicted_landslides.shp', 'Total_landslides.xlsx']
        >>> cells = pd.read_excel(os.path.join(folder, "Individual_landslides.xlsx"),
        ...     sheet_name="Triggering_Cells")["Triggering cell"]
        >>> file_rows = cells.to_numpy() // 20          # row in the DEM file, 0 = top
        >>> bool(len(file_rows) > 0 and file_rows.max() < 3)
        True
        """
        rows_cols = (self._rows, self._clms)

        # with a coupled component the number of time steps is not fixed, so the
        # earthquake time step may never have been reached: say so
        if (self._Triggering_Indicator == 2 and self._number_of_time_steps is None
                and self._Earthquake_time_point >= self._current_time_index):
            self._log(
                f"Earthquake_time_point = {self._Earthquake_time_point} was not reached "
                f"({self._current_time_index} time step(s) run): no seismic displacement "
                "was computed."
            )

        # per-landslide workbook (triggering-cell indices converted to rows from north to south)
        _save_individual_excel(self._state, self._output_dir, rows_cols,
                               steady_state=self._steady_state)
        # totals workbook (one row per time step, or a single row when steady state)
        _save_summary_excel(self._state, self._output_dir, steady_state=self._steady_state)

    # =============================================================================
    # Properties
    # =============================================================================

    # quick read-only access to the landslide-ID grid
    @property
    def landslide_fid(self):
        """ndarray of float, shape (nrows, ncols) : Cumulative landslide-ID raster.

        Examples
        --------
        Access the cumulative landslide-ID raster after running one time
        step (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> crisis.landslide_fid.shape
        (3, 100)
        >>> np.unique(crisis.landslide_fid[~np.isnan(crisis.landslide_fid)])
        array([1., 2.], dtype=float32)
        """
        return self._state["FID_all_times"]

    # quick read-only access to the failed-cells grid
    @property
    def failed_cells(self):
        """ndarray of bool, shape (nrows, ncols) : True where any failure has occurred.

        Examples
        --------
        Access the cumulative failed-cell mask after running one time step
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> crisis.failed_cells.dtype
        dtype('bool')
        >>> int(crisis.failed_cells.sum())
        97
        """
        return self._state["Fail_Index_all_times"]

    # quick read-only access to the entire state dictionary
    @property
    def state(self):
        """dict : Full output accumulator (per-step and cumulative landslide results).

        Examples
        --------
        Access the state accumulator after running one time step
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> sorted(crisis.state.keys())
        ['FID_all_times', 'Fail_Index_all_times', 'Individual_Cumulative_areas', 'Individual_Cumulative_areas_with_time', 'Individual_Cumulative_triggering_depths', 'Individual_Cumulative_triggering_depths_with_time', 'Individual_Cumulative_triggering_linear_indices', 'Individual_Cumulative_triggering_linear_indices_with_time', 'Individual_Cumulative_triggering_sources', 'Individual_Cumulative_triggering_sources_with_time', 'Individual_Cumulative_volumes', 'Individual_Cumulative_volumes_with_time', 'Individual_Current_areas_with_time', 'Individual_Current_triggering_depths_with_time', 'Individual_Current_triggering_linear_indices_with_time', 'Individual_Current_triggering_sources_with_time', 'Individual_Current_volumes_with_time', 'Total_Cumulative_Area_vector_with_time', 'Total_Cumulative_NL_vector_with_time', 'Total_Cumulative_Volume_vector_with_time', 'Total_Current_Area_vector_with_time', 'Total_Current_NL_vector_with_time', 'Total_Current_Volume_vector_with_time']
        >>> float(crisis.state["Total_Current_NL_vector_with_time"][0])  # new landslides at t = 0
        2.0
        """
        return self._state

    # quick read-only access to the current time-step counter
    @property
    def current_time_index(self):
        """int : Index of the next time step to be run.

        Examples
        --------
        The counter starts at 0 and advances by one after each call to
        :meth:`run_one_step` (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.current_time_index
        0
        >>> crisis.run_one_step()
        >>> crisis.current_time_index
        1
        """
        return self._current_time_index

    # read-only number of time steps of the run (never a parameter)
    @property
    def number_of_time_steps(self):
        """int or None : Number of time steps of the run.

        The number of time steps is never a parameter. It is 1 with a
        constant Ru (``Hydrological_model_indicator = 1``), the number of
        time steps of the HDF5 files found in ``hydro_files_dir`` (indicator
        2), and ``None`` with grid fields written by a coupled component
        (indicator 3): the run then has as many time steps as
        :meth:`run_one_step` is called.

        Examples
        --------
        Same setup as :meth:`run_one_step`: a constant Ru has a single time
        step.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> CrisisComponent(grid, params).number_of_time_steps
        1

        With pore pressure heads written on the grid by a coupled component
        (indicator 3), the number of time steps is not fixed:

        >>> _ = grid.add_field("soil_water__pressure_head",
        ...                    np.zeros((grid.number_of_nodes, 8)), at="node", clobber=True)
        >>> params.update({"Hydrological_model_indicator": 3,
        ...                "Strength_variation_with_theta": 1})
        >>> CrisisComponent(grid, params).number_of_time_steps is None
        True
        """
        return self._number_of_time_steps

    # =============================================================================
    # Private helpers
    # =============================================================================

    # ─────────────────────────────────────────────────────────────────────────────
    # _init_output_fields — create output grid fields with default values
    # ─────────────────────────────────────────────────────────────────────────────

    # sets up the five output fields on the grid before any steps are run
    def _init_output_fields(self):
        """Create output node fields on the grid with default values.

        Examples
        --------
        Called automatically by ``__init__``; the output fields already
        exist on the grid by the time :meth:`run_one_step` is called
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> sorted(f for f in grid.at_node.keys()
        ...        if f.startswith("landslide") or f in ("failure__depth", "triggering__cells"))
        ['failure__depth', 'landslide__failed_accumulated', 'landslide__failure_id', 'landslide__failure_id_accumulated', 'triggering__cells']
        """

        # Each field is only created if it does not already exist on the grid

        # cumulative landslide-ID field, default NaN (no landslide yet)
        if "landslide__failure_id_accumulated" not in self.grid.at_node:
            self.grid.add_full("landslide__failure_id_accumulated", np.nan, at="node", dtype=float)

        # cumulative failed flag, default False
        if "landslide__failed_accumulated" not in self.grid.at_node:
            self.grid.add_zeros("landslide__failed_accumulated", at="node", dtype=bool)

        # per-step failure depth (m), default 0
        if "failure__depth" not in self.grid.at_node:
            self.grid.add_zeros("failure__depth", at="node", dtype=float)

        # per-step triggering-cell flag, default False
        if "triggering__cells" not in self.grid.at_node:
            self.grid.add_zeros("triggering__cells", at="node", dtype=bool)

        # per-step landslide ID, default 0
        if "landslide__failure_id" not in self.grid.at_node:
            self.grid.add_zeros("landslide__failure_id", at="node", dtype=int)

    # ─────────────────────────────────────────────────────────────────────────────
    # _update_output_fields — push 2D results back onto the flat Landlab grid
    # ─────────────────────────────────────────────────────────────────────────────

    # copies our internal 2D results onto the Landlab grid's flat fields
    def _update_output_fields(self):
        """Flatten 2D state arrays to 1D and write them into Landlab node fields.

        Examples
        --------
        Called automatically at the end of :meth:`run_one_step`; the
        accumulated-FID field reflects the state after one time step
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> np.unique(grid.at_node["landslide__failure_id_accumulated"])
        array([ 1.,  2., nan])
        """

        # flatten the landslide-ID grid and write it into the Landlab field
        self.grid.at_node["landslide__failure_id_accumulated"][:] = self._state["FID_all_times"].ravel()
        # flatten the failed-cells grid and write it into the Landlab field
        self.grid.at_node["landslide__failed_accumulated"][:]     = self._state["Fail_Index_all_times"].ravel()

        # per-step fields: only written when new landslides formed in this step
        if len(self.LD_Z) > 0:
            self._grid.at_node['failure__depth'][self.LD_Index] = self.LD_Z
            self._grid.at_node['triggering__cells'][self.Trig_Index_Accept] = True
            self.grid.at_node["landslide__failure_id"][:] = self._FID_Index.flatten()

    # ─────────────────────────────────────────────────────────────────────────────
    # _build_fd_mat_from_landlab — convert Landlab routing to D8 codes
    # ─────────────────────────────────────────────────────────────────────────────

    # converts Landlab's flow-routing info into the D8 direction codes CRISIS expects
    def _build_fd_mat_from_landlab(self, grid, flow_direction):
        """Convert Landlab receiver-node indices to D8 direction codes.

        Uses the same D8 codes (ArcGIS numbering 1-128) as ``_crisis_core``
        (see ``T_UpGrid``), each code labeling the (row, col) offset from a
        cell to its receiver: (0, 1)=1, (1, 1)=2, (1, 0)=4, (1, -1)=8,
        (0, -1)=16, (-1, -1)=32, (-1, 0)=64, (-1, 1)=128. Rows follow
        Landlab's orientation (row 0 = southern edge). Sink / outlet cells
        (receiver == self) are assigned code 0.

        Parameters
        ----------
        grid : RasterModelGrid
            The grid carrying the flow-direction field.
        flow_direction : str
            Name of the node field holding D8 receiver nodes.

        Returns
        -------
        FD_Mat : ndarray of float32, shape (nrows, ncols)
            D8 code of each node (0 = sink / outlet).

        Examples
        --------
        On the uniformly inclined 3x100 grid, every interior cell drains
        one column west (offset (0, -1), D8 code 16) toward the fixed outlet
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> fd_mat = crisis._build_fd_mat_from_landlab(grid, "flow__receiver_node")
        >>> fd_mat.shape
        (3, 100)
        >>> np.unique(fd_mat)
        array([ 0., 16.], dtype=float32)
        >>> fd_mat[1, :5]
        array([ 0., 16., 16., 16., 16.], dtype=float32)
        """

        # ── Setup ─────────────────────────────────────────────────────────────

        # grid dimensions
        nrows, ncols = grid.shape
        # which node each node's water flows into, reshaped to 2D
        receiver     = grid.at_node[flow_direction].reshape(nrows, ncols)
        # output grid of direction codes, starts all zero (no direction / sink)
        FD_Mat       = np.zeros((nrows, ncols), dtype=np.float32)

        # ── Direction code lookup table ───────────────────────────────────────
        # Map (row_offset, col_offset) → D8 code (same table as _crisis_core.T_UpGrid)
        # Offset = receiver_position − source_position

        # lookup table: which direction code matches which row/col shift
        offset_to_code = {
            ( 0,  1):   1,
            ( 1,  1):   2,
            ( 1,  0):   4,
            ( 1, -1):   8,
            ( 0, -1):  16,
            (-1, -1):  32,
            (-1,  0):  64,
            (-1,  1): 128,
        }

        # ── Build flat index arrays ───────────────────────────────────────────
        # Build (row, col) index arrays for every node in the grid

        # build full row/column index grids for every node
        row_idx, col_idx = np.meshgrid(
            np.arange(nrows), np.arange(ncols), indexing="ij"
        )
        # flatten the row indices into a 1D list
        flat_row = row_idx.ravel()
        # flatten the column indices into a 1D list
        flat_col = col_idx.ravel()
        # flatten the receiver-node indices into a 1D list
        flat_rec = receiver.ravel()

        # ── Compute row/col offsets to each receiver ──────────────────────────
        # Convert flat receiver indices to (row, col) and compute offsets

        # turn each receiver's flat index back into a (row, col) pair
        rec_row, rec_col = np.unravel_index(flat_rec, (nrows, ncols))
        # row offset from each cell to its receiver
        dr = rec_row - flat_row
        # column offset from each cell to its receiver
        dc = rec_col - flat_col

        # ── Assign D8 direction codes ─────────────────────────────────────────
        # Assign D8 code; sink cells with no matching offset get code 0

        # go cell by cell and look up its direction code from the offset table
        for i in range(len(flat_row)):
            FD_Mat[flat_row[i], flat_col[i]] = offset_to_code.get(
                (int(dr[i]), int(dc[i])), 0
            )

        # hand back the finished direction-code grid
        return FD_Mat

    # ─────────────────────────────────────────────────────────────────────────────
    # _calc_hydrology — get the water contents and pore pressure heads of one time step
    # ─────────────────────────────────────────────────────────────────────────────

    def _calc_hydrology(self, t, params):
        """Get Theta_3D and Pressures_3D for time step *t*.

        - ``Pressures_3D`` (pore pressure head, m), by
          ``Hydrological_model_indicator``:

          * 1: ``None``; the core derives the pressure head from ``Ru``
            (the ``soil_water__pressure_head`` field was filled with the same
            values by ``__init__``).
          * 2: read from ``Pressure_head_T{t}Z{z}.h5`` in ``hydro_files_dir``,
            each file flipped by ``load_h5``
            (rows from north to south -> Landlab's row order), each layer
            capped at its depth, and written to the pore pressure head field.
          * 3: read from the pore pressure head field (no flip: the field is in
            Landlab's node order), each layer capped at its depth. The field
            itself is not modified.

        - ``Theta_3D`` (volumetric water content), by
          ``Strength_variation_with_theta``:

          * 1, or ``Hydrological_model_indicator == 1``: ``None`` (Bishop χ = 1,
            no water content needed).
          * 2: read from ``Theta_T{t}Z{z}.h5`` in ``theta_files_dir``, each
            file flipped by ``load_h5``, and written to the water content field.
          * 3: read from the water content field (no flip). The field itself
            is not modified.

        Parameters
        ----------
        t : int
            Zero-based time-step index.
        params : dict
            CRISIS parameter dict.

        Returns
        -------
        Theta_3D : ndarray of float32, shape (rows, clms, Depth_points) or None
            None when no water content is used.
        Pressures_3D : ndarray of float32, shape (rows, clms, Depth_points) or None
            None when ``Hydrological_model_indicator == 1``.

        Raises
        ------
        FileNotFoundError
            If a pore pressure head or water content file of time step *t* is
            required but not found.
        ValueError
            If a pore pressure head / water content field supplied by a coupled
            component does not have the shape ``(number_of_nodes, Depth_points)``.

        Examples
        --------
        With ``Strength_variation_with_theta=1`` and
        ``Hydrological_model_indicator=1`` (this example's settings), no water
        content is used (Theta_3D is ``None``) and Pressures_3D is ``None``:
        the core uses ``Ru``. The pore pressure head field holds the Ru pressure
        head ``Ru × (γ / γ_w) × z`` of every layer, at every node
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> theta_3d, pressures_3d = crisis._calc_hydrology(0, params)
        >>> theta_3d is None, pressures_3d is None
        (True, True)
        >>> head = grid.at_node["soil_water__pressure_head"]
        >>> head.shape                                  # (number_of_nodes, Depth_points)
        (300, 8)
        >>> np.allclose(head, 0.8 * 22.0 / 9.81 * crisis._Depths_vector)
        True
        >>> "soil_water__volume_fraction" in grid.at_node
        False

        Pressure heads and water contents from HDF5 files
        (``Hydrological_model_indicator=2``, ``Strength_variation_with_theta=2``).
        The files store their rows from north to south: here the first row
        of every file (the northern edge) has a pressure head of 5 m and a
        water content of 0.4, the other rows -1 m and 0.2. After the flip,
        that row is Landlab's last row (row 2, y = 2 m); each layer's pressure
        head is capped at its depth.

        >>> import h5py
        >>> north_first_head  = np.full((3, 100), -1.0); north_first_head[0] = 5.0
        >>> north_first_theta = np.full((3, 100), 0.2);  north_first_theta[0] = 0.4
        >>> for z in range(1, 9):
        ...     with h5py.File(os.path.join(folder, f"Pressure_head_T0Z{z}.h5"), "w") as f:
        ...         _ = f.create_dataset("data", data=north_first_head)
        ...     with h5py.File(os.path.join(folder, f"Theta_T0Z{z}.h5"), "w") as f:
        ...         _ = f.create_dataset("data", data=north_first_theta)
        >>> params.update({"Hydrological_model_indicator": 2, "hydro_files_dir": folder,
        ...                "Strength_variation_with_theta": 2, "theta_files_dir": folder})
        >>> crisis = CrisisComponent(grid, params)
        Hydrological_model_indicator = 2: 1 time step(s) found in the HDF5 files (T0 ... T0).
        >>> theta_3d, pressures_3d = crisis._calc_hydrology(0, params)
        >>> pressures_3d.shape
        (3, 100, 8)
        >>> pressures_3d[:, 0, 0].tolist(), pressures_3d[:, 0, 7].tolist()  # rows 0 (south) .. 2 (north)
        ([-1.0, -1.0, 0.125], [-1.0, -1.0, 1.875])
        >>> theta_3d[:, 0, 0].tolist()
        [0.20000000298023224, 0.20000000298023224, 0.4000000059604645]

        The capped pressure heads and the water contents are written to the
        grid fields, one row per node. Node 200 is the first node of grid row
        2, the northern edge (y = 2 m):

        >>> float(grid.y_of_node[200])
        2.0
        >>> grid.at_node["soil_water__pressure_head"][200].tolist()
        [0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875]
        >>> grid.at_node["soil_water__pressure_head"][0].tolist() == [-1.0] * 8
        True
        >>> np.round(grid.at_node["soil_water__volume_fraction"][[0, 200], 0], 3).tolist()
        [0.2, 0.4]

        Pressure heads and water contents from grid fields written by a
        coupled Landlab component (``Hydrological_model_indicator=3``,
        ``Strength_variation_with_theta=3``). The fields are in Landlab's node
        order, so they are only reshaped, never flipped: here the pressure
        head of every layer equals the node's y (row 0 = 0 m ... row 2 = 2 m)
        and the water content is 0.1 + 0.1 y.

        >>> layers_of = lambda v: np.repeat(v[:, np.newaxis], 8, axis=1)  # same value in all 8 layers
        >>> grid.at_node["soil_water__pressure_head"][:] = layers_of(grid.y_of_node)
        >>> grid.at_node["soil_water__volume_fraction"][:] = layers_of(0.1 + 0.1 * grid.y_of_node)
        >>> params.update({"Hydrological_model_indicator": 3, "Strength_variation_with_theta": 3})
        >>> crisis = CrisisComponent(grid, params)
        >>> theta_3d, pressures_3d = crisis._calc_hydrology(0, params)
        >>> pressures_3d[:, 0, 0].tolist(), pressures_3d[:, 0, 7].tolist()  # capped at 0.125 / 1.875 m
        ([0.0, 0.125, 0.125], [0.0, 1.0, 1.875])
        >>> np.round(theta_3d[:, 50, 3].astype(float), 3).tolist()
        [0.1, 0.2, 0.3]

        The fields themselves are left as the other component wrote them
        (the pressure heads are capped in a copy):

        >>> float(grid.at_node["soil_water__pressure_head"].max())
        2.0
        """

        # ── Shortcuts ─────────────────────────────────────────────────────────
        # number of depth layers
        Depth_points = params["Depth_points"]
        # where the pore pressure heads come from: 1 = Ru, 2 = HDF5 files, 3 = grid field
        Hydro_ind    = self._Hydro_indicator
        # where the water contents come from: 1 = not used, 2 = HDF5 files, 3 = grid field
        # (already set to 1 by __init__ for a constant Ru)
        Theta_var    = self._Strength_variation_with_theta

        # ── File paths (only used with indicator = 2, the HDF5 files) ─────────

        # folder to look in for pore pressure head files
        def hydro_dir():
            return Path(params["hydro_files_dir"])

        # folder to look in for water content files
        def theta_dir():
            return Path(params["theta_files_dir"])

        # ── Volumetric water content — Theta_3D ───────────────────────────────
        # Theta_var == 1: Bishop χ = 1, no water content needed
        Theta_3D = None

        if Theta_var == 2:
            # load water content HDF5 files from an external hydrological model
            folder  = theta_dir()
            missing = [d + 1 for d in range(Depth_points)
                       if not (folder / f"Theta_T{t}Z{d + 1}.h5").exists()]
            if missing:
                raise FileNotFoundError(
                    f"Strength_variation_with_theta = 2, but the water content file "
                    f"Theta_T{t}Z{missing[0]}.h5 was not found in '{folder}'. Expected "
                    f"Theta_T{t}Z1.h5 ... Theta_T{t}Z{Depth_points}.h5 for time step {t}."
                )
            # load one file per depth layer into a 3D array
            # (load_h5 flips each file, rows from north to south, to Landlab's row order)
            Theta_3D = np.zeros((self._rows, self._clms, Depth_points), dtype=np.float32)
            for d in range(Depth_points):
                Theta_3D[:, :, d] = load_h5(folder / f"Theta_T{t}Z{d + 1}.h5")
            # store this step's water contents on the grid (no flip: same row order)
            self._write_layer_field(self._water_content_name, Theta_3D)

        elif Theta_var == 3:
            # water contents written on the grid by a coupled component: Landlab
            # node order already, so reshape only (no flip); astype makes a copy,
            # so the field is not modified
            Theta_3D = self._layer_field_to_3d(self._water_content_name).astype(np.float32)

        # ── Pore pressure head — Pressures_3D ────────────────────────────────
        # Hydro_ind == 1: None, the core uses Ru (field filled once by __init__).
        # Pressure is capped at the layer depth — water table cannot exceed surface.

        # start with no pressure data; only fill it in if time-varying hydrology is on
        Pressures_3D = None

        if Hydro_ind == 2:
            # load pore pressure head HDF5 files from an external hydrological model
            folder  = hydro_dir()
            missing = [d + 1 for d in range(Depth_points)
                       if not (folder / f"Pressure_head_T{t}Z{d + 1}.h5").exists()]
            if missing:
                raise FileNotFoundError(
                    f"Hydrological_model_indicator = 2, but the pore pressure head file "
                    f"Pressure_head_T{t}Z{missing[0]}.h5 was not found in '{folder}'. Expected "
                    f"Pressure_head_T{t}Z1.h5 ... Pressure_head_T{t}Z{Depth_points}.h5 for "
                    f"time step {t}."
                )
            # load one file per depth layer, capping pressure at that layer's depth
            # (load_h5 flips each file, rows from north to south, to Landlab's row order)
            Pressures_3D = np.zeros((self._rows, self._clms, Depth_points), dtype=np.float32)
            for d in range(Depth_points):
                P2D = load_h5(folder / f"Pressure_head_T{t}Z{d + 1}.h5")
                Pressures_3D[:, :, d] = np.minimum(P2D, self._Depths_vector[d])
            # store this step's capped pressure heads on the grid (no flip: same row order)
            self._write_layer_field(self._pressure_head_name, Pressures_3D)

        elif Hydro_ind == 3:
            # pore pressure heads written on the grid by a coupled component: Landlab
            # node order already, so reshape only (no flip). Each layer is capped
            # at its depth, like indicator = 2 (HDF5 files); np.minimum returns a new array,
            # so the field is not modified.
            P3D = self._layer_field_to_3d(self._pressure_head_name).astype(np.float32)
            Pressures_3D = np.minimum(P3D, self._Depths_vector[np.newaxis, np.newaxis, :])

        # ── Return both arrays ────────────────────────────────────────────────

        # hand back both the water content and pore pressure head arrays
        return Theta_3D, Pressures_3D

    # ─────────────────────────────────────────────────────────────────────────────
    # Layered node fields — (number_of_nodes, Depth_points) <-> (nrows, ncols, Depth_points)
    # ─────────────────────────────────────────────────────────────────────────────
    # Landlab numbers nodes row by row from the south-west corner: node i is grid
    # cell (row, col) = (i // ncols, i % ncols), with row 0 = southern edge. A field
    # with one row per node therefore becomes a (nrows, ncols, Depth_points) array in
    # CRISIS's row order by a plain row-major (C-order) reshape — no flip.

    # turns a layered node field into a (nrows, ncols, Depth_points) array
    def _layer_field_to_3d(self, name):
        """Return the layered node field *name* as a (nrows, ncols, Depth_points) view.

        The field must have one row per node (Landlab order) and one column
        per depth layer (shallow to deep): shape ``(number_of_nodes,
        Depth_points)``. A 1D field is accepted when ``Depth_points == 1``.
        Any floating-point dtype is accepted (e.g. float32 from another
        component); the field is used as it is, never replaced.
        The array returned is a view of the field (no copy, no flip):
        element ``[row, col, k]`` is the field value of node
        ``row * ncols + col`` in layer ``k``.

        Parameters
        ----------
        name : str
            Name of the node field.

        Returns
        -------
        ndarray, shape (nrows, ncols, Depth_points)

        Raises
        ------
        KeyError
            If the field does not exist on the grid.
        TypeError
            If the field does not hold floating-point values.
        ValueError
            If the field does not have the shape ``(number_of_nodes, Depth_points)``.

        Examples
        --------
        Same setup as :meth:`run_one_step`, with pore pressure heads given on
        the grid by a coupled component (``Hydrological_model_indicator = 3``):
        here, each node's pressure head is its node number in layer 0 and
        minus its node number in layer 1 (``Depth_points = 2``).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.5, "z_max": 1.0,
        ...     "Depth_points": 2, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 3,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> node = np.arange(grid.number_of_nodes, dtype=float)
        >>> _ = grid.add_field("soil_water__pressure_head",
        ...                    np.stack([node, -node], axis=1), at="node")
        >>> crisis = CrisisComponent(grid, params)
        >>> head_3d = crisis._layer_field_to_3d("soil_water__pressure_head")
        >>> head_3d.shape
        (3, 100, 2)

        Element ``[row, col, k]`` is node ``row * ncols + col``, layer ``k``:
        row 0 is the southern edge (y = 0), row 2 the northern edge (y = 2 m).

        >>> head_3d[0, 0].tolist(), head_3d[2, 5].tolist()
        ([0.0, -0.0], [205.0, -205.0])
        >>> float(grid.y_of_node[0]), float(grid.y_of_node[205])
        (0.0, 2.0)

        A field that is missing, or that does not have one column per depth
        layer, is rejected.

        >>> crisis._layer_field_to_3d("no_such_field")
        Traceback (most recent call last):
        ...
        KeyError: "Grid field 'no_such_field' not found. With Hydrological_model_indicator = 3 (pore pressure heads) or Strength_variation_with_theta = 3 (water contents), add it to the grid before creating CrisisComponent, with shape (number_of_nodes, Depth_points) = (300, 2)."
        >>> _ = grid.add_field("three_layers", np.zeros((300, 3)), at="node")
        >>> crisis._layer_field_to_3d("three_layers")
        Traceback (most recent call last):
        ...
        ValueError: Grid field 'three_layers' has shape (300, 3), but must have shape (number_of_nodes, Depth_points) = (300, 2): one row per node, one column per depth layer (shallow to deep).
        >>> _ = grid.add_field("int_layers", np.zeros((300, 2), dtype=int), at="node")
        >>> crisis._layer_field_to_3d("int_layers")
        Traceback (most recent call last):
        ...
        TypeError: Grid field 'int_layers' has dtype int64, but must hold floating-point values (e.g. float64 or float32).

        A float32 field (e.g. written by a coupled component) is accepted as
        it is: it is not converted or replaced, so the coupled component keeps
        writing into the same array that CRISIS reads.

        >>> head32 = np.stack([node, -node], axis=1).astype(np.float32)
        >>> _ = grid.add_field("soil_water__pressure_head", head32, at="node", clobber=True)
        >>> crisis = CrisisComponent(grid, params)
        >>> grid.at_node["soil_water__pressure_head"].dtype
        dtype('float32')
        >>> head32[205] = [9.0, -9.0]           # the other component updates its array
        >>> crisis._layer_field_to_3d("soil_water__pressure_head")[2, 5].tolist()
        [9.0, -9.0]
        """
        n_nodes = self.grid.number_of_nodes
        n_dep   = len(self._Depths_vector)

        # the field must exist on the grid
        if name not in self.grid.at_node:
            raise KeyError(
                f"Grid field {name!r} not found. With Hydrological_model_indicator = 3 "
                "(pore pressure heads) or Strength_variation_with_theta = 3 (water contents), "
                "add it to the grid before creating CrisisComponent, with shape "
                f"(number_of_nodes, Depth_points) = ({n_nodes}, {n_dep})."
            )
        values = self.grid.at_node[name]

        # any floating-point dtype is fine (float64, float32, ...); integers are not
        if not np.issubdtype(values.dtype, np.floating):
            raise TypeError(
                f"Grid field {name!r} has dtype {values.dtype}, but must hold "
                "floating-point values (e.g. float64 or float32)."
            )

        # a single depth layer may be stored as a plain 1D field
        if n_dep == 1 and values.shape == (n_nodes,):
            values = values[:, np.newaxis]

        # one row per node, one column per depth layer
        if values.shape != (n_nodes, n_dep):
            raise ValueError(
                f"Grid field {name!r} has shape {values.shape}, but must have shape "
                f"(number_of_nodes, Depth_points) = ({n_nodes}, {n_dep}): one row per "
                "node, one column per depth layer (shallow to deep)."
            )

        # node i -> (i // ncols, i % ncols): row-major reshape, no flip, no copy
        return values.reshape(self._rows, self._clms, n_dep)

    # writes a (nrows, ncols, Depth_points) array into a layered node field
    def _write_layer_field(self, name, array_3d):
        """Write a (nrows, ncols, Depth_points) array into the layered node field *name*.

        The reverse of :meth:`_layer_field_to_3d`: cell ``(row, col)`` goes to
        node ``row * ncols + col`` (row-major reshape, no flip). The field is
        updated in place.

        Examples
        --------
        Called by :meth:`_calc_hydrology` for indicator = 2 (HDF5 files); see its
        examples. On its own (same setup as :meth:`_layer_field_to_3d`):

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> import h5py
        >>> for z in (1, 2):   # pore pressure head files of one time step, 2 depth layers
        ...     with h5py.File(os.path.join(folder, f"Pressure_head_T0Z{z}.h5"), "w") as f:
        ...         _ = f.create_dataset("data", data=np.zeros((nrows, ncols)))
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.5, "z_max": 1.0,
        ...     "Depth_points": 2, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 2, "hydro_files_dir": folder,
        ...     "Strength_variation_with_theta": 1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        Hydrological_model_indicator = 2: 1 time step(s) found in the HDF5 files (T0 ... T0).
        >>> a3d = np.zeros((3, 100, 2)); a3d[2, 5] = [7.0, 8.0]   # row 2 (north), column 5
        >>> crisis._write_layer_field("soil_water__pressure_head", a3d)
        >>> head = grid.at_node["soil_water__pressure_head"]
        >>> np.flatnonzero(head[:, 0]).tolist(), head[205].tolist()   # node 2 * 100 + 5
        ([205], [7.0, 8.0])
        """
        field = self.grid.at_node[name]
        field[:] = np.asarray(array_3d).reshape(field.shape)

    # creates a layered node field written by the component (or checks an existing one)
    def _ensure_layer_field(self, name):
        """Create the layered node field *name* if missing, or check its shape.

        A missing field is created, filled with NaN until the first time step
        writes it; an existing field is kept and its shape is checked (see
        :meth:`_layer_field_to_3d`).

        Examples
        --------
        Called by ``__init__`` for the fields the component writes. With
        ``Hydrological_model_indicator = 2``, the pore pressure head field is
        created, empty (NaN), until the first time step reads the HDF5 files
        (same setup as :meth:`_write_layer_field`):

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> import h5py
        >>> for z in (1, 2):   # pore pressure head files of one time step, 2 depth layers
        ...     with h5py.File(os.path.join(folder, f"Pressure_head_T0Z{z}.h5"), "w") as f:
        ...         _ = f.create_dataset("data", data=np.zeros((nrows, ncols)))
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.5, "z_max": 1.0,
        ...     "Depth_points": 2, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 2, "hydro_files_dir": folder,
        ...     "Strength_variation_with_theta": 1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        Hydrological_model_indicator = 2: 1 time step(s) found in the HDF5 files (T0 ... T0).
        >>> head = grid.at_node["soil_water__pressure_head"]
        >>> head.shape, bool(np.isnan(head).all())
        ((300, 2), True)
        """
        if name in self.grid.at_node:
            self._layer_field_to_3d(name)          # raises if the shape is wrong
        else:
            shape = (self.grid.number_of_nodes, len(self._Depths_vector))
            self.grid.add_field(name, np.full(shape, np.nan), at="node")

    # ─────────────────────────────────────────────────────────────────────────────
    # _compute_rupture_distances — distance from every cell to the fault plane
    # ─────────────────────────────────────────────────────────────────────────────

    # works out how far each cell is from the fault rupture plane
    def _compute_rupture_distances(self, params):
        """Compute closest-distance-to-rupture for Rathje flexible / Bray models.

        Parameters
        ----------
        params : dict
            CRISIS parameter dict, with the fault corners ``fault_X1`` ...
            ``fault_Z4``.

        Returns
        -------
        ndarray of float32, shape (rows, clms)
            Closest distance (km) from every cell to the fault plane.

        Examples
        --------
        On the same 3x100 grid used by :meth:`run_one_step` (cell centers at
        x = 0..99 m, y = 0..2 m), place a vertical fault plane along x=0
        running from y=0 to y=50 m and from the surface to 1000 m depth, and
        compute each cell's distance (km) to it. Column 0 lies on the plane;
        columns 1 and 2 lie above its top edge (elevations 0.9 m and 1.8 m),
        so their distance is to that edge: sqrt(1² + 0.9²) and sqrt(2² + 1.8²) m.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> import h5py
        >>> pga_file = os.path.join(folder, "PGA.h5")  # coupled mode needs a PGA raster
        >>> with h5py.File(pga_file, "w") as f:
        ...     _ = f.create_dataset("data", data=np.full((nrows, ncols), 0.3))
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.35,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ...     "Triggering_Indicator": 2, "Seismic_displacement_model": 3,
        ...     "pga_file": pga_file, "Earthquake_time_point": 0,
        ...     "Seismic_displacement_threshold": 5.0, "M": 6.5,
        ...     "Rathje_method": 1, "Ts": 0.3,
        ...     "fault_X1": 0.0, "fault_Y1": 0.0,  "fault_Z1": 0.0,
        ...     "fault_X2": 0.0, "fault_Y2": 50.0, "fault_Z2": 0.0,
        ...     "fault_X3": 0.0, "fault_Y3": 50.0, "fault_Z3": 1000.0,
        ...     "fault_X4": 0.0, "fault_Y4": 0.0,  "fault_Z4": 1000.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        Computing fault rupture distances...
        >>> rupture = crisis._compute_rupture_distances(params)
        Computing fault rupture distances...
        >>> rupture.shape
        (3, 100)
        >>> np.round(rupture[:, :3], 3)
        array([[0.   , 0.001, 0.003],
               [0.   , 0.001, 0.003],
               [0.   , 0.001, 0.003]], dtype=float32)
        """

        # ── Define fault plane corners ────────────────────────────────────────
        # 4-corner fault plane; Z values are depths (positive downward)

        # the four corners of the fault plane, as x/y/z points
        corners = np.array([
            [params["fault_X1"], params["fault_Y1"], params["fault_Z1"]],
            [params["fault_X2"], params["fault_Y2"], params["fault_Z2"]],
            [params["fault_X3"], params["fault_Y3"], params["fault_Z3"]],
            [params["fault_X4"], params["fault_Y4"], params["fault_Z4"]],
        ], dtype=np.float64)

        # ── Prepare 2D coordinate and elevation grids ─────────────────────────
        # Cell-center coordinates computed in __init__ (float64, to avoid
        # precision loss with large projected coordinate values)

        # elevation grid at high precision
        _elev        = self._Elev_Mat.astype(np.float64)
        # grid size shortcuts
        _rows, _clms = _elev.shape

        # stretch the 1D x (per column) and y (per row) vectors into 2D grids
        _x2d = np.broadcast_to(self._x[np.newaxis, :], (_rows, _clms))
        _y2d = np.broadcast_to(self._y[:, np.newaxis], (_rows, _clms))

        # replace NaN (no-data) elevations with 0 (sea level) for the distance calculation
        _z2d = np.where(np.isnan(_elev), 0.0, _elev)

        # ── Compute and return rupture distances ──────────────────────────────

        # let the user know this (possibly slow) calculation is starting
        self._log("Computing fault rupture distances...")
        # compute and return the distance from every cell to the fault plane
        return Rupture_Dist(_x2d, _y2d, _z2d, corners).astype(np.float32)

    # ─────────────────────────────────────────────────────────────────────────────
    # _update_cumulative_registry — merge new landslides into the master records
    # ─────────────────────────────────────────────────────────────────────────────

    # adds this time step's newly accepted landslides into the running master records
    def _update_cumulative_registry(self, t, Trig_Index_Accept, Trig_Depth_Accept,
                                    Trig_source_Accept, Vol, Area, Cell_size):
        """Merge newly accepted landslides into the cumulative tracking registry.

        Parameters
        ----------
        t : int
            Current time-step index.
        Trig_Index_Accept : ndarray of int64
            Linear indices of the triggering cells of the new landslides.
        Trig_Depth_Accept : list of float
            Their triggering depths (m).
        Trig_source_Accept : list of int
            Their triggering sources (1 = rainfall, 2 = rainfall + earthquake,
            3 = earthquake).
        Vol : ndarray of float32
            Their volumes (m³).
        Area : ndarray of float32
            Their areas (m²).
        Cell_size : float
            Cell size (m), from the DEM header.

        Examples
        --------
        Re-applying the same triggering cells from time step 0 at time
        step 1 leaves the cumulative registry unchanged, since duplicate
        cell indices are deduplicated (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> sorted(crisis.state["Individual_Cumulative_triggering_linear_indices"].tolist())
        [100, 176]
        >>> Trig_Index_Accept  = crisis.Trig_Index_Accept
        >>> Trig_Depth_Accept  = crisis.state["Individual_Current_triggering_depths_with_time"][-1]
        >>> Trig_source_Accept = crisis.state["Individual_Current_triggering_sources_with_time"][-1]
        >>> Vol  = crisis.state["Individual_Current_volumes_with_time"][-1]
        >>> Area = crisis.state["Individual_Current_areas_with_time"][-1]
        >>> crisis._update_cumulative_registry(1, Trig_Index_Accept, Trig_Depth_Accept, Trig_source_Accept, Vol, Area, dx)
        >>> sorted(crisis.state["Individual_Cumulative_triggering_linear_indices"].tolist())
        [100, 176]
        """

        # ── Shortcuts ─────────────────────────────────────────────────────────

        # current landslide-ID grid
        FID   = self._state["FID_all_times"]
        # shortcut to the state dictionary
        state = self._state

        # the master list of triggering cell locations seen so far
        cum_idx = state["Individual_Cumulative_triggering_linear_indices"]

        # ── Merge new landslides into the running registry ────────────────────
        # While the registry is empty: initialize it directly.
        # Afterwards: concatenate and deduplicate — np.unique keeps the first
        # occurrence of each linear index so re-triggered cells are not double-counted.

        # while the registry is still empty, start it with this step's results
        if len(cum_idx) == 0:
            state["Individual_Cumulative_triggering_linear_indices"] = Trig_Index_Accept.copy()
            state["Individual_Cumulative_triggering_depths"]         = list(Trig_Depth_Accept)
            state["Individual_Cumulative_triggering_sources"]        = list(Trig_source_Accept)
            state["Individual_Cumulative_volumes"]                   = Vol.copy()

        elif len(Trig_Index_Accept) > 0:
            # later steps: merge old and new data together
            merged_idx = np.concatenate((cum_idx, Trig_Index_Accept))
            merged_dep = np.concatenate((state["Individual_Cumulative_triggering_depths"],
                                         Trig_Depth_Accept))
            merged_src = np.concatenate((state["Individual_Cumulative_triggering_sources"],
                                         Trig_source_Accept))
            merged_vol = np.concatenate((state["Individual_Cumulative_volumes"], Vol))

            # drop duplicate cells, keeping only the first time each one appeared
            # (np.unique returns them sorted, so the registry is ordered by linear index)
            _, unique_idx = np.unique(merged_idx, return_index=True)
            state["Individual_Cumulative_triggering_linear_indices"] = merged_idx[unique_idx]
            state["Individual_Cumulative_triggering_depths"]         = merged_dep[unique_idx]
            state["Individual_Cumulative_triggering_sources"]        = merged_src[unique_idx]
            state["Individual_Cumulative_volumes"]                   = merged_vol[unique_idx]

        # ── Recompute cumulative areas from the FID raster ────────────────────
        # Area of each landslide = number of cells carrying its ID × cell area

        # grab the up-to-date master list of triggering cells
        cum_trig = state["Individual_Cumulative_triggering_linear_indices"]

        if len(cum_trig) > 0:
            # only look at cells that actually belong to a landslide
            valid_fid = FID[~np.isnan(FID)]
            # count how many cells belong to each landslide ID
            fids_uniq, counts = np.unique(valid_fid, return_counts=True)
            # turn that cell count into a real area in square meters, per landslide ID
            fid_to_area = {int(f): cnt * Cell_size ** 2 for f, cnt in zip(fids_uniq, counts)}
            # look up which landslide ID each triggering cell belongs to
            fid_of_trig = FID.flat[cum_trig]
            # look up the area for each triggering cell's landslide
            cum_areas   = np.array([fid_to_area.get(int(f), 0.0) for f in fid_of_trig],
                                   dtype=np.float32)
        else:
            # no triggering cells yet, so no areas to report
            cum_areas = np.array([], dtype=np.float32)

        # save the freshly computed cumulative areas
        state["Individual_Cumulative_areas"] = cum_areas

        # ── Append per-timestep snapshots for Excel output ────────────────────
        # one entry per column in the final Excel workbook

        # record this step's new triggering cell locations
        state["Individual_Current_triggering_linear_indices_with_time"].append(Trig_Index_Accept)
        # record this step's new failure depths
        state["Individual_Current_triggering_depths_with_time"].append(Trig_Depth_Accept)
        # record this step's new trigger sources
        state["Individual_Current_triggering_sources_with_time"].append(Trig_source_Accept)
        # record this step's new landslide areas
        state["Individual_Current_areas_with_time"].append(Area)
        # record this step's new landslide volumes
        state["Individual_Current_volumes_with_time"].append(Vol)

        # save a snapshot of the full master list of triggering cells, for this step
        state["Individual_Cumulative_triggering_linear_indices_with_time"].append(
            state["Individual_Cumulative_triggering_linear_indices"].copy())
        # save a snapshot of the full master list of failure depths, for this step
        state["Individual_Cumulative_triggering_depths_with_time"].append(
            state["Individual_Cumulative_triggering_depths"].copy())
        # save a snapshot of the full master list of trigger sources, for this step
        state["Individual_Cumulative_triggering_sources_with_time"].append(
            state["Individual_Cumulative_triggering_sources"].copy())
        # save a snapshot of the cumulative areas, for this step
        state["Individual_Cumulative_areas_with_time"].append(cum_areas)
        # save a snapshot of the full master list of volumes, for this step
        state["Individual_Cumulative_volumes_with_time"].append(
            state["Individual_Cumulative_volumes"].copy())

    # ─────────────────────────────────────────────────────────────────────────────
    # get_topographic_slope / get_drainage_area — current topography from the grid
    # ─────────────────────────────────────────────────────────────────────────────

    def get_topographic_slope(self, slope=None):
        """Return the current slope array (degrees) from the grid field.

        Parameters
        ----------
        slope : str, optional
            Name of the Landlab node field to read the slope from.
            Defaults to the ``slope`` field name given to the constructor.

        Returns
        -------
        ndarray of float, shape (number_of_nodes,)

        Examples
        --------
        The grid already carries a ``topographic__slope`` field, so it is
        returned directly (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.get_topographic_slope()[:3]
        array([41.9872125, 41.9872125, 41.9872125])
        """

        return self._grid.at_node[slope or self._slope_name]

    def get_drainage_area(self, drainage_area=None):
        """Return the current drainage area array (m^2) from the grid field.

        Parameters
        ----------
        drainage_area : str, optional
            Name of the Landlab node field to read the drainage area from.
            Defaults to the ``drainage_area`` field name given to the
            constructor.

        Returns
        -------
        ndarray of float, shape (number_of_nodes,)

        Examples
        --------
        The grid already carries a ``drainage_area`` field (created by
        ``FlowAccumulator``), so it is returned directly
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.get_drainage_area()[:3]
        array([0., 0., 0.])
        """

        return self._grid.at_node[drainage_area or self._drainage_area_name]

    # ─────────────────────────────────────────────────────────────────────────────
    # _reset_output_fields — clear the per-step output fields
    # ─────────────────────────────────────────────────────────────────────────────

    def _reset_output_fields(self):
        """Clear the per-time-step output fields back to their default values.

        Called at the start of every :meth:`run_one_step` call so that
        ``triggering__cells``, ``failure__depth``, and
        ``landslide__failure_id`` reflect only the current time step, not
        cells triggered in a previous one.

        Examples
        --------
        After one time step, ``triggering__cells`` has some ``True``
        entries; resetting clears them back to their defaults
        (same setup as :meth:`run_one_step`).

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import CrisisComponent
        >>> nrows, ncols = 3, 100
        >>> dx = 1.0  # grid spacing (m)
        >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
        >>> slope = 0.9  # rise per unit run
        >>> elevation = slope * grid.x_of_node
        >>> _ = grid.add_field("topographic__elevation", elevation, at="node")
        >>> s = grid.add_zeros("topographic__slope", at="node")
        >>> s[:] = np.degrees(grid.calc_slope_at_node("topographic__elevation"))
        >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> outlet_id = grid.nodes_at_left_edge[1]
        >>> grid.status_at_node[outlet_id] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid, runoff_rate=1)
        >>> fa.run_one_step()
        >>> import os, tempfile
        >>> folder = tempfile.mkdtemp()
        >>> import io
        >>> dem_file = io.StringIO("ncols 100\\nnrows 3\\nxllcenter 0.0\\nyllcenter 0.0\\ncellsize 1.0\\n")  # ASCII DEM header of this grid, in memory
        >>> params = {
        ...     "output_dir": folder,
        ...     "Coordinate_Reference_System": "EPSG:32610", "dem_file": dem_file,
        ...     "Gamma": 22.0, "z_min": 0.125, "z_max": 1.875,
        ...     "Depth_points": 8, "Min_LD_area": 12.0,
        ...     "top_or_bottom_indicator": 1, "Triggering_Indicator": 1,
        ...     "Hydrological_model_indicator": 1, "Ru": 0.8,
        ...     "Strength_variation_with_theta": 1,
        ...     "Theta_s": 0.45, "Theta_r": 0.1,
        ...     "Strength_indicator": 1, "c": 20.0, "Phi": 37.0,
        ... }
        >>> crisis = CrisisComponent(grid, params)
        >>> crisis.run_one_step()
        >>> crisis._reset_output_fields()
        >>> bool(grid.at_node["triggering__cells"].any())
        False
        >>> float(grid.at_node["failure__depth"].sum())
        0.0
        >>> int(grid.at_node["landslide__failure_id"].sum())
        0
        """

        # reset the three per-step output fields to their defaults
        self._grid.at_node['triggering__cells'][:] = False
        self._grid.at_node['failure__depth'].fill(0.0)
        self._grid.at_node["landslide__failure_id"].fill(0)