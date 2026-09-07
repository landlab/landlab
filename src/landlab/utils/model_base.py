#! /usr/bin/env python

# # Base class for a grid-based Landlab model
#
# This code defines LandlabModel, a Python class that is designed to make it easier
# to create a standalone model code using Landlab. The model developer writes
# a class that inherits from LandlabModel and adds the functionality needed to
# implement their model. LandlabModel provides code to handle formatted user input,
# in the form of either a Python dictionary or the name of a yaml-format input
# file (given as a string). The LandlabModel __init__() method will combine the
# user inputs with a set of default parameter values defined in the derived
# class header (for parameters whose value has not been specified by the user).
# For model execution, the user simply needs to override the built-in update()
# method. LandlabModel calls this via a built-in run() method (which runs the model
# from start to finish) and a built-in update_until() method (which calls
# update() until the either the run is complete or it is time to pause and
# generate output).
#
# *(Greg Tucker, University of Colorado Boulder)*
#
from collections.abc import Iterator
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import count
from typing import Any
from typing import Self

import numpy as np
from requireit import require_less_than
from requireit import require_positive
from requireit import require_sorted

from landlab.core.component_utils import iter_adaptive_time_steps
from landlab.core.component_utils import iter_time_steps
from landlab.core.model_parameter_loader import load_params
from landlab.grid.base import ModelGrid


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
        self._clock = clock
        self._current_time = self._clock.start

        self.setup_for_output(self.params, self._clock)

    @property
    def run_duration(self) -> float:
        return self._clock.duration

    @property
    def dt(self) -> float:
        return self._clock.step

    @classmethod
    def from_file(cls, input_file: str) -> Self:
        """Create a model from an input file.

        Parameters
        ----------
        input_file : str
            Name of yaml-format file containing names and values
        """
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
        return self._current_time

    def setup_for_output(self, params: dict, clock: Clock) -> None:
        """
        Setup variables for control of plotting and saving.

        Parameters
        ----------
        params : dict
            Parameter dictionary. Must include a key ``output`` with a dictionary
        containing values for ``plot_times``, ``save_times``, and ``report_times``.
        Each of these should be either a ``float`` or a ``list``. If a list, the value
        is interpreted as a list of model times for plotting, saving, or reporting.
        If a single float, the value is interpreted as the (regular) time
        interval (in model time) for plotting, saving, or reporting.
            Should also contain a key ``clock`` as a dictionary that has values
        for ``start`` and ``stop``.

        Notes
        -----
        The "format" parameter can be "grid" (native Landlab grid format), "netcdf",
        or "vtk". The default is "grid".
        """
        op_params = params["output"]

        self._plot_schedule = _PauseSchedule(
            op_params["plot_times"], start=clock.start, stop=clock.stop
        )
        self._save_schedule = _PauseSchedule(
            op_params["save_times"], start=clock.start, stop=clock.stop
        )
        self._report_schedule = _PauseSchedule(
            op_params["report_times"], start=clock.start, stop=clock.stop
        )

        if self._save_schedule.is_due(clock.start):
            self._save_schedule.advance()

        self.ndigits_for_save_files = 4
        self.save_num = 0  # current save file frame number
        self.save_path = op_params["save_path"]
        if op_params["plot_to_file"]:
            self.ndigits_for_plot_files = 4
            self.plot_num = 0  # current plot image frame number
        self.display_params = params

        if "format" in op_params:
            save_fmt = op_params["format"]
        else:
            save_fmt = "grid"

        if save_fmt == "grid":
            self.save_state = self.save_state_grid_format
        elif save_fmt == "vtk":
            self.save_state = self.save_state_vtk_format
        elif save_fmt == "netcdf":
            self.save_state = self.save_state_netcdf_format
        else:
            print("Unrecognized save format '" + save_fmt + "'.")
            print("Valid formats are: grid, vtk, netcdf")
            raise ValueError

    def report(self, current_time: float) -> None:
        """Issue a text update on status."""
        print(self.__class__.__name__, "time =", current_time)

    def plot(self, current_time: float = 0.0) -> None:
        """Virtual function for plotting; to be overridden."""
        print("Base class placeholder for plot() at time", current_time)

    def save_state_grid_format(
        self, save_path: str, save_num: int, ndigits: int
    ) -> None:
        """
        Save the grid and its fields in native Landlab format.

        Override this function to add to or modify what gets saved.
        """
        from landlab.io.native_landlab import save_grid

        save_grid(
            self.grid, save_path + str(save_num).zfill(ndigits) + ".grid", clobber=True
        )

    def save_state_vtk_format(
        self, save_path: str, save_num: int, ndigits: int
    ) -> None:
        """
        Save grid fields in legacy VTK format.

        Override this function to add to or modify what gets saved.
        """
        from landlab.io.legacy_vtk import write_legacy_vtk

        write_legacy_vtk(
            save_path + str(save_num).zfill(ndigits) + ".vtk", self.grid, clobber=True
        )

    def save_state_netcdf_format(
        self, save_path: str, save_num: int, ndigits: int
    ) -> None:
        """
        Save grid fields in NetCDF format.

        Override this function to add to or modify what gets saved.
        """
        from landlab.io.netcdf import write_netcdf

        write_netcdf(save_path + str(save_num).zfill(ndigits) + ".nc", self.grid)

    def update(self, dt: float) -> None:
        """
        Advance the model by one time step of duration dt.

        The derived class should override this function.
        """
        pass

    def update_until(self, update_to_time: float, dt: float) -> None:
        """Iterate up to given time, using time-step duration dt."""
        duration = update_to_time - self.current_time
        if duration <= 0.0:
            return

        for this_dt in iter_time_steps(duration, dt=dt):
            self.update(this_dt)
            self._current_time += this_dt

    def run(self, run_duration: float | None = None, dt: float | None = None) -> None:
        """Run the model for given duration, or self.run_duration if none
        given.

        Includes file output of images and model state at user-specified
        intervals.
        """
        if run_duration is None:
            run_duration = self.run_duration
        if dt is None:
            dt = self.dt

        self._run_scheduled_actions()
        for time_until_pause in iter_adaptive_time_steps(
            run_duration, calc_dt=self._time_to_next_pause
        ):
            self.update_until(self.current_time + time_until_pause, dt)
            self._run_scheduled_actions()

    def _time_to_next_pause(self) -> float:
        return (
            min(
                self._plot_schedule.next_pause,
                self._save_schedule.next_pause,
                self._report_schedule.next_pause,
            )
            - self.current_time
        )

    def _run_scheduled_actions(self) -> None:
        if self._report_schedule.is_due(self.current_time):
            self.report(self.current_time)
            self._report_schedule.advance()

        if self._plot_schedule.is_due(self.current_time):
            self.plot(self.current_time)
            self._plot_schedule.advance()

        if self._save_schedule.is_due(self.current_time):
            self.save_num += 1
            self.save_state(self.save_path, self.save_num, self.ndigits_for_save_files)
            self._save_schedule.advance()


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
