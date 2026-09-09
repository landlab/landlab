#! /usr/bin/env python

"""Base class and runner for a grid-based Landlab model.

Model authors subclass :class:`LandlabModel` and implement :meth:`~LandlabModel.update`
to advance their components by a supplied time step. The base class constructs the
grid and clock, while :class:`ModelRunner` advances time and runs scheduled output
events.

The following landscape-evolution model combines uniform uplift, stream-power
erosion, and linear hillslope diffusion:

*(Greg Tucker, University of Colorado Boulder)*

Examples
--------
>>> import numpy as np
>>> from landlab.components import FastscapeEroder
>>> from landlab.components import FlowAccumulator
>>> from landlab.components import LinearDiffuser
>>> class LandscapeEvolutionModel(LandlabModel):
...     def __init__(self, grid, *, clock, params):
...         super().__init__(grid, clock=clock, params=params)
...         self._uplift_rate = params["uplift_rate"]
...         elevation = grid.add_zeros("topographic__elevation", at="node")
...         elevation[:] = 0.01 * grid.node_y
...         self._flow_router = FlowAccumulator(grid, flow_director="D8")
...         self._flow_router.run_one_step()
...         self._eroder = FastscapeEroder(grid, K_sp=params["erodibility"])
...         self._diffuser = LinearDiffuser(
...             grid, linear_diffusivity=params["diffusivity"]
...         )
...
...     def update(self, dt):
...         elevation = self.grid.at_node["topographic__elevation"]
...         elevation[self.grid.core_nodes] += self._uplift_rate * dt
...         self._flow_router.run_one_step()
...         self._eroder.run_one_step(dt)
...         self._diffuser.run_one_step(dt)
...
...     def report(self, current_time):
...         print(f"model time: {current_time:g}")
...
>>> model = LandscapeEvolutionModel.from_params(
...     {
...         "grid": {
...             "source": "create",
...             "create_grid": {
...                 "RasterModelGrid": {
...                     "shape": (5, 5),
...                     "xy_spacing": 1.0,
...                 },
...             },
...         },
...         "clock": {"start": 0.0, "stop": 2.0, "step": 1.0},
...         "uplift_rate": 0.001,
...         "erodibility": 0.01,
...         "diffusivity": 0.1,
...         "output": {
...             "report_times": [0.0, 1.0, 2.0],
...             "plot_times": [],
...             "save_times": [],
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
from collections.abc import Callable
from collections.abc import Iterator
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import count
from typing import Any
from typing import Self

import numpy as np
from requireit import require_less_than
from requireit import require_nonnegative
from requireit import require_one_of
from requireit import require_positive
from requireit import require_sorted

from landlab.core.component_utils import iter_adaptive_time_steps
from landlab.core.component_utils import iter_time_steps
from landlab.core.model_parameter_loader import load_params
from landlab.grid.base import ModelGrid
from landlab.io.legacy_vtk import write_legacy_vtk
from landlab.io.native_landlab import save_grid
from landlab.io.netcdf import write_netcdf


def merge_params(
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
    >>> merged = merge_params(user, defaults=defaults)
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

            merged[k] = merge_params(v, defaults=default_value)

    return merged


def resolve_array_filepaths(params: dict[str, Any]) -> dict[str, Any]:
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
                resolved[key] = resolve_array_filepaths(value)
        else:
            resolved[key] = value
    return resolved


@dataclass(frozen=True, slots=True)
class Clock:
    """Define the time domain and default time step for a model run.

    Parameters
    ----------
    start : float, optional
        Initial model time.
    stop : float, optional
        Final model time. It must be greater than ``start``.
    step : float, optional
        Positive, finite default time-step duration.

    Examples
    --------
    >>> clock = Clock(start=2.0, stop=8.0, step=0.5)
    >>> clock.duration
    6.0
    """

    start: float = 0.0
    stop: float = np.inf
    step: float = 1.0

    def __post_init__(self) -> None:
        require_less_than(self.start, self.stop, name="start")
        require_positive(self.step, name="step")
        if np.isinf(self.step):
            raise ValueError("step must be finite")

    @property
    def duration(self) -> float:
        return self.stop - self.start


class _PauseSchedule:
    """Track the next pause in a sequence of times.

    Parameters
    ----------
    schedule : float or sequence of float
        Constant interval between pauses, or a sequence of absolute times
        at which to pause.
    start : float, optional
        Earliest time in the schedule. For a constant interval, this is
        also the first pause. Explicit times before ``start`` are skipped.
    stop : float, optional
        Latest time in the schedule. The stop time is included. With the
        default of infinity, a constant-interval schedule is unbounded.

    Examples
    --------
    >>> schedule = _PauseSchedule(1.0, start=0.0, stop=4.0)
    >>> schedule.next_pause
    0.0
    >>> schedule.advance()
    1.0
    >>> schedule = _PauseSchedule([0.0, 0.5, 2.0, 4.0], start=0.5, stop=4.0)
    >>> schedule.next_pause
    0.5
    """

    def __init__(
        self,
        schedule: float | Sequence[float],
        *,
        start: float = 0.0,
        stop: float = np.inf,
    ) -> None:
        self._times = _iter_pause_times(schedule=schedule, start=start, stop=stop)
        self._next_pause = next(self._times, np.inf)

    @property
    def next_pause(self) -> float:
        return self._next_pause

    def is_due(self, time: float) -> bool:
        return time >= self._next_pause

    def advance(self) -> float:
        self._next_pause = next(self._times, np.inf)
        return self._next_pause


@dataclass(slots=True)
class _Event:
    schedule: _PauseSchedule
    action: Callable[[float], None]

    @property
    def next_time(self) -> float:
        return self.schedule.next_pause

    def run_if_due(self, time: float) -> None:
        if self.schedule.is_due(time):
            self.action(time)
            self.schedule.advance()


def _iter_pause_times(
    schedule: float | Sequence[float],
    *,
    start: float = 0.0,
    stop: float = np.inf,
) -> Iterator[float]:
    if isinstance(schedule, (float, int)):
        require_positive(schedule, name="pause interval")
        if not np.isfinite(schedule):
            raise ValueError("pause interval must be finite")
        for step in count():
            next_pause = start + step * schedule
            if next_pause > stop:
                break
            yield next_pause
    else:
        require_sorted(schedule, strict=True, name="schedule")
        for next_pause in schedule:
            if next_pause < start:
                continue
            if next_pause > stop:
                break
            yield next_pause


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


class ModelRunner:
    def __init__(
        self,
        model: LandlabModel,
        *,
        clock: Clock,
        events: Mapping[str, _Event] | None = None,
    ) -> None:
        self._model = model
        self._clock = clock
        self._current_time = clock.start
        self._events = {} if events is None else dict(events)

    @property
    def current_time(self) -> float:
        return self._current_time

    @property
    def run_duration(self) -> float:
        return self._clock.duration

    @property
    def dt(self) -> float:
        return self._clock.step

    def update_until(
        self,
        update_to_time: float,
        dt: float,
    ) -> None:
        duration = update_to_time - self.current_time
        if duration <= 0.0:
            return

        for this_dt in iter_time_steps(duration, dt=dt):
            self._model.update(this_dt)
            self._current_time += this_dt
        self._current_time = update_to_time

    def run(
        self,
        run_duration: float | None = None,
        dt: float | None = None,
    ) -> None:
        if run_duration is None:
            run_duration = self._clock.stop - self.current_time
        if dt is None:
            dt = self._clock.step

        self._run_scheduled_actions()
        for time_until_pause in iter_adaptive_time_steps(
            run_duration, calc_dt=self._time_to_next_pause
        ):
            self.update_until(self.current_time + time_until_pause, dt)
            self._run_scheduled_actions()

    def _time_to_next_pause(self) -> float:
        return (
            min((event.next_time for event in self._events.values()), default=np.inf)
            - self.current_time
        )

    def _run_scheduled_actions(self) -> None:
        for event in self._events.values():
            event.run_if_due(self.current_time)


class LandlabModel:
    """
    Base class for a generic Landlab grid-based model.

    Examples
    --------
    >>> from landlab.utils.model_base import LandlabModel
    >>> class MyModel(LandlabModel):
    ...     pass
    ...
    >>> p = {"grid": {"source": "create"}}
    >>> p["grid"]["create_grid"] = {
    ...     "RasterModelGrid": {"shape": (4, 5), "xy_spacing": 2.0}
    ... }
    >>> model = MyModel.from_params(p)
    >>> model.grid.shape
    (4, 5)
    """

    # Default parameters, to be overridden in derived classes
    DEFAULT_PARAMS = {
        "grid": {
            "source": "create",
            "create_grid": {
                "RasterModelGrid": {
                    "shape": (5, 5),
                    "xy_spacing": 1.0,
                },
            },
        },
        "clock": {"start": 0.0, "stop": 2.0, "step": 1.0},
        "output": {
            "plot_times": 10.0,  # float or list
            "save_times": 10.0,  # float or list
            "report_times": 1.0,  # float or list
            "save_path": "model_output",
            "clobber": True,
            "fields": None,
            "plot_to_file": True,
        },
    }

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
        params : dict
            Dictionary containing names and values of model parameters
        """
        self.grid = grid
        self.params = params

        output_params = params["output"]
        self._saver = _GridSaver(
            grid,
            output_params["save_path"],
            fmt=output_params.get("format", "grid"),
            ndigits=4,
        )

        events = _build_events(
            output_params,
            clock=clock,
            actions={
                "plot": self.plot,
                "report": self.report,
                "save": self.save,
            },
        )

        if events["save"].schedule.is_due(clock.start):
            events["save"].schedule.advance()

        self._runner = ModelRunner(self, clock=clock, events=events)

    @property
    def run_duration(self) -> float:
        return self._runner.run_duration

    @property
    def dt(self) -> float:
        return self._runner.dt

    @classmethod
    def from_file(cls, input_file: str) -> Self:
        """Create a model from an input file.

        Parameters
        ----------
        input_file : str
            Name of a YAML or TOML file containing model parameters. TOML files
            are identified by a ``.toml`` extension; other files are read as YAML.
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
        params = {} if params is None else params

        params = merge_params(params, defaults=cls.DEFAULT_PARAMS)
        params = resolve_array_filepaths(params)

        grid = setup_grid(params["grid"])
        clock = Clock(**params["clock"])
        return cls(grid, clock=clock, params=params)

    @property
    def current_time(self) -> float:
        return self._runner.current_time

    def report(self, current_time: float) -> None:
        """Issue a text update on status."""
        print(self.__class__.__name__, "time =", current_time)

    def plot(self, current_time: float = 0.0) -> None:
        """Virtual function for plotting; to be overridden."""
        print("Base class placeholder for plot() at time", current_time)

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
        """Iterate up to given time, using time-step duration dt."""
        self._runner.update_until(update_to_time, dt=dt)

    def run(self, run_duration: float | None = None, dt: float | None = None) -> None:
        """Run the model for given duration, or self.run_duration if none
        given.

        Includes file output of images and model state at user-specified
        intervals.
        """
        self._runner.run(run_duration, dt=dt)


def _build_events(
    params: dict[str, Any],
    *,
    clock: Clock,
    actions: Mapping[str, Callable[[float], None]],
) -> dict[str, _Event]:
    start, stop = clock.start, clock.stop

    events = {
        name: _Event(
            _PauseSchedule(params[f"{name}_times"], start=start, stop=stop),
            action=action,
        )
        for name, action in actions.items()
    }

    return events


def setup_grid(params: dict) -> ModelGrid:
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
    >>> grid = setup_grid(params=p)
    >>> grid.shape
    (4, 5)

    >>> from landlab import RasterModelGrid
    >>> p = {"source": "grid_object"}
    >>> p["grid_object"] = RasterModelGrid((3, 3))
    >>> grid = setup_grid(params=p)
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
