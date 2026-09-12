#! /usr/bin/env python

"""Base class and runner for a grid-based Landlab model.

Model authors subclass :class:`Model` and implement :meth:`~Model.update`
to advance their components by a supplied time step. The base class constructs the
grid and clock, while :class:`ModelRunner` advances time and runs scheduled output
events.

The following landscape-evolution model combines uniform uplift, stream-power
erosion, and linear hillslope diffusion:

*(Greg Tucker, University of Colorado Boulder)*

Examples
--------
>>> import numpy as np
>>> from landlab.components import FlowAccumulator
>>> from landlab.components import LinearDiffuser
>>> from landlab.components import StreamPowerEroder
>>> from landlab.core.model import Model

>>> class LandscapeEvolutionModel(Model):
...     DEFAULT_PARAMS = {
...         "grid": {
...             "source": "create",
...             "create_grid": {
...                 "RasterModelGrid": {
...                     "shape": (5, 5),
...                     "xy_spacing": 1.0,
...                 },
...             },
...         },
...     }
...
...     def __init__(self, grid, *, clock, params):
...         super().__init__(grid, clock=clock, params=params)
...         rng = np.random.default_rng()
...         elevation = grid.add_zeros("topographic__elevation", at="node")
...         elevation[grid.core_nodes] = rng.uniform(size=len(grid.core_nodes))
...
...         self._uplift_rate = params["model"]["parameters"]["uplift_rate"]
...         self._flow_accumulator = FlowAccumulator(
...             grid, **params["model"]["components"]["flow_accumulator"]
...         )
...         self._flow_accumulator.run_one_step()
...         self._eroder = StreamPowerEroder(
...             grid, **params["model"]["components"]["eroder"]
...         )
...         self._diffuser = LinearDiffuser(
...             grid, **params["model"]["components"]["diffuser"]
...         )
...
...     def update(self, dt):
...         elevation = self.grid.at_node["topographic__elevation"]
...         elevation[self.grid.core_nodes] += self._uplift_rate * dt
...         self._diffuser.run_one_step(dt)
...         self._flow_accumulator.run_one_step()
...         self._eroder.run_one_step(dt)
...
...     def report(self, current_time):
...         print(f"model time: {current_time:g}")
...

>>> model = LandscapeEvolutionModel.from_params(
...     {
...         "clock": {"start": 0.0, "stop": 2.0, "step": 1.0},
...         "model": {
...             "parameters": {"uplift_rate": 0.001},
...             "components": {
...                 "flow_accumulator": {"flow_director": "D8"},
...                 "eroder": {"K_sp": 0.01},
...                 "diffuser": {"linear_diffusivity": 0.1},
...             },
...         },
...         "events": {
...             "report": {"times": [0.0, 1.0, 2.0]},
...         },
...     }
... )

>>> model.run()
model time: 0
model time: 1
model time: 2
>>> model.current_time
2.0
>>> np.all(np.isfinite(model.grid.at_node["topographic__elevation"]))
True
"""

from __future__ import annotations

import os
import tomllib
from typing import Any
from typing import ClassVar
from typing import Self

import numpy as np
from requireit import require_contains
from requireit import require_nonnegative
from requireit import require_one_of

from landlab.core.model_parameter_loader import load_params
from landlab.core.model_runner import Clock
from landlab.core.model_runner import ModelRunner
from landlab.core.model_runner import _build_events
from landlab.grid.base import ModelGrid
from landlab.io.legacy_vtk import write_legacy_vtk
from landlab.io.native_landlab import save_grid
from landlab.io.netcdf import write_netcdf

__all__ = ["Model"]


def _merge_params(
    user: dict[str, Any],
    *,
    defaults: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Merge parameters with defaults, returning a new nested dictionary.

    Merge default parameters into the user-parameter dictionary, adding
    defaults where user values are absent. Nested dictionaries are merged
    recursively, except for ``grid``, which is treated as a single value.

    Parameters
    ----------
    user : dict
        dict containing names and values of user-defined parameters
    defaults : dict, optional
        dict containing all parameter names and their default values

    Returns
    -------
    merged : dict
        The merged parameters.

    Examples
    --------
    >>> user = {"a": 1, "d": {"da": 4}, "e": 5, "grid": {"RasterModelGrid": []}}
    >>> defaults = {"a": 2, "b": 3, "d": {"db": 6}, "grid": {"HexModelGrid": []}}
    >>> merged = _merge_params(user, defaults=defaults)
    >>> merged["a"] == user["a"]
    True
    >>> merged["b"] == defaults["b"]
    True
    >>> sorted(merged["d"].items())
    [('da', 4), ('db', 6)]

    >>> merged["grid"]
    {'RasterModelGrid': []}
    """
    defaults = {} if defaults is None else defaults

    merged = {**defaults, **user}
    for k, v in merged.items():
        if isinstance(v, dict):
            default_value = defaults.get(k)
            if k == "grid" or not isinstance(default_value, dict):
                default_value = None

            merged[k] = _merge_params(v, defaults=default_value)

    return merged


def _resolve_array_filepaths(params: dict[str, Any]) -> dict[str, Any]:
    """Return new parameters with array filepath references resolved.

    Dictionary values containing an ``"_filepath"`` key are replaced by
    arrays loaded from the referenced files. Nested parameter dictionaries
    are processed recursively.

    Parameters
    ----------
    params : dict
        Parameter dictionary that may contain array filepath references.

    Returns
    -------
    resolved : dict
        A new parameter dictionary containing the resolved arrays.
    """
    resolved = {}
    for key, value in params.items():
        if isinstance(value, dict):
            if "_filepath" in value:
                resolved[key] = np.load(value["_filepath"])
            else:
                resolved[key] = _resolve_array_filepaths(value)
        else:
            resolved[key] = value
    return resolved


class _FilenameSequence:
    def __init__(self, base_name: str, *, ndigits: int = 0, ext: str = "") -> None:
        self._base_name = base_name
        self._ndigits = require_nonnegative(ndigits, name="ndigits")
        self._ext = ext
        self._frame = 0

    def __next__(self) -> str:
        self._frame += 1
        return self._build_filename()

    def __iter__(self) -> Self:
        return self

    def _build_filename(self) -> str:
        return f"{self._base_name}" f"{self._frame:0{self._ndigits}d}" f"{self._ext}"


class _GridSaver:
    EXTENSIONS = {
        "grid": ".grid",
        "netcdf": ".nc",
        "vtk": ".vtk",
    }

    def __init__(
        self,
        grid: ModelGrid,
        base_name: str,
        *,
        fmt="grid",
        ndigits: int = 4,
    ) -> None:
        fmt = require_one_of(fmt, allowed=_GridSaver.EXTENSIONS, name="fmt")
        self._filenames = _FilenameSequence(
            base_name, ndigits=ndigits, ext=self.EXTENSIONS[fmt]
        )

        self._grid = grid

        self._write = getattr(self, f"_write_{fmt}")

    def __call__(self, time: float) -> None:
        self.save()

    def save(self) -> str:
        filename = next(self._filenames)
        self._write(filename)
        return filename

    def _write_grid(self, filename: str) -> None:
        save_grid(self._grid, filename, clobber=True)

    def _write_netcdf(self, filename: str) -> None:
        write_netcdf(filename, self._grid)

    def _write_vtk(self, filename: str) -> None:
        write_legacy_vtk(filename, self._grid, clobber=True)


class Model:
    """Base class for a time-dependent, grid-based Landlab model.

    ``Model`` provides configuration constructors, scheduled reporting and
    output, and model time management. Subclasses define the model physics by
    constructing their components and implementing :meth:`update`. They may
    override :meth:`plot`, :meth:`report`, and :meth:`save` to customize the
    corresponding scheduled events.

    Parameters
    ----------
    grid : ModelGrid
        Grid shared by the model's components.
    clock : Clock
        Start time, stop time, and default time-step duration.
    params : dict
        Model parameters. The ``events`` section configures the scheduled
        events.

    See Also
    --------
    Clock
        Definition of the model time domain.
    ModelRunner
        Time-stepping and event orchestration.
    """

    # Default parameters, to be overridden in derived classes
    DEFAULT_PARAMS: ClassVar[dict[str, Any]] = {}

    def __init__(
        self,
        grid: ModelGrid,
        *,
        clock: Clock,
        params: dict[str, Any],
    ) -> None:
        """Initialize the model.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab `ModelGrid`.
        clock : Clock
            Start time, stop time, and default time-step duration.
        params : dict
            Dictionary containing names and values of model parameters
        """
        self.grid = grid
        self.params = params

        event_params = params.get("events", {})
        events = _build_events(
            event_params,
            clock=clock,
            actions={
                "plot": self.plot,
                "report": self.report,
                "save": self.save,
            },
        )

        if "save" in events:
            if events["save"].schedule.is_due(clock.start):
                events["save"].schedule.advance()

        save_params = event_params.get("save", {})
        self._saver = _GridSaver(
            grid,
            save_params.get("base_name", "model-output"),
            fmt=save_params.get("format", "grid"),
            ndigits=save_params.get("ndigits", 4),
        )
        self._runner = ModelRunner(self, clock=clock, events=events)

    @property
    def run_duration(self) -> float:
        return self._runner.run_duration

    @property
    def dt(self) -> float:
        return self._runner.dt

    @classmethod
    def from_file(cls, input_file: str) -> Self:
        """Create a model from parameters stored in a YAML or TOML file.

        The file contents are loaded into a parameter dictionary and passed to
        :meth:`from_params`. Files with a ``.toml`` extension are read as TOML;
        all other files are read as YAML.

        Parameters
        ----------
        input_file : str
            Name of the parameter file.

        Returns
        -------
        Model
            Model constructed from the parameters in ``input_file``.
        """
        if os.path.splitext(input_file)[1].lower() == ".toml":
            with open(input_file, "rb") as fp:
                params = tomllib.load(fp)
        else:
            with open(input_file) as fp:
                params = load_params(fp)
        return cls.from_params(params=params)

    @classmethod
    def from_params(cls, params: dict[str, Any] | None = None) -> Self:
        """Create a model from a parameter dictionary.

        User parameters are merged with :attr:`DEFAULT_PARAMS`, references to
        arrays stored in files are resolved, and the model grid and clock are
        constructed before the class is initialized.

        Parameters
        ----------
        params : dict, optional
            Model parameters that override :attr:`DEFAULT_PARAMS`.

        Returns
        -------
        Model
            Model constructed from the merged parameters.
        """
        params = {} if params is None else params

        params = _merge_params(params, defaults=cls.DEFAULT_PARAMS)
        params = _resolve_array_filepaths(params)

        params = require_contains(params, required=("clock", "grid"), name="params")

        grid = _setup_grid(params["grid"])
        clock = Clock(**params["clock"])
        return cls(grid, clock=clock, params=params)

    @property
    def current_time(self) -> float:
        return self._runner.current_time

    def report(self, current_time: float) -> None:
        """Issue a text update on status."""
        print(f"time = {current_time}")

    def plot(self, current_time: float = 0.0) -> None:
        """Virtual function for plotting; to be overridden."""
        raise NotImplementedError("plot")

    def save(self, current_time: float) -> None:
        """Save a grid."""
        self._saver(current_time)

    def update(self, dt: float) -> None:
        """
        Advance the model by one time step of duration dt.

        The derived class should override this function.
        """
        pass

    def update_until(self, update_to_time: float, dt: float) -> None:
        """Advance the model to an absolute model time.

        This method advances the model without running scheduled events.

        Parameters
        ----------
        update_to_time : float
            Model time to which the model should advance. If this is not later
            than the current time, the model is unchanged.
        dt : float
            Maximum time-step duration.
        """
        self._runner.update_until(update_to_time, dt=dt)

    def run(self, run_duration: float | None = None, dt: float | None = None) -> None:
        """Advance the model while running scheduled events.

        Parameters
        ----------
        run_duration : float, optional
            Duration of the run. By default, advance from the current time to the
            stop time of the model clock.
        dt : float, optional
            Maximum time-step duration. By default, use the step specified by the
            model clock.
        """
        self._runner.run(run_duration, dt=dt)


def _setup_grid(params: dict) -> ModelGrid:
    """Load or create the grid.

    Parameters
    ----------
    params : dict
        Dictionary containing parameters related grid setup.

    Notes
    -----
    Must include an item "source" for which the valid values are
    "create" (create a new grid), "file" (read a grid from file), or
    "grid_object" (indicating that a grid object is included
    directly in the parameter dictionary).

    If "create", then there must be an item "create_grid" that
    contains a dict in which the key is the name of the grid type
    ("RasterModelGrid", "HexModelGrid") and the value is a dict
    containing the names and values for the grid object's
    parameters (such as a tuple for "shape", etc.)

    If "file", then there must be an item "grid_file_name" that
    contains the file name as a string.

    If "grid_object", then there must be an item called "grid_object"
    containing the grid object!

    Examples
    --------
    >>> p = {"source": "create"}
    >>> p["create_grid"] = {"RasterModelGrid": {"shape": (4, 5), "xy_spacing": 2.0}}
    >>> grid = _setup_grid(params=p)
    >>> grid.shape
    (4, 5)

    >>> from landlab import RasterModelGrid
    >>> p = {"source": "grid_object"}
    >>> p["grid_object"] = RasterModelGrid((3, 3))
    >>> grid = _setup_grid(params=p)
    >>> grid.shape
    (3, 3)
    """
    from requireit import require_one_of

    from landlab import create_grid
    from landlab.io.native_landlab import load_grid

    source = require_one_of(
        params["source"], allowed=("create", "file", "grid_object"), name="source"
    )

    if source == "create":
        return create_grid(params, section="create_grid")

    if source == "file":
        return load_grid(params["grid_file_name"])

    if source == "grid_object" and isinstance(params["grid_object"], ModelGrid):
        return params["grid_object"]

    raise ValueError("grid source must be one of 'create', 'file', or a grid instance")
