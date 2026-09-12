from __future__ import annotations

from collections.abc import Callable
from collections.abc import Iterator
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import count
from typing import Any
from typing import Protocol

import numpy as np
from requireit import require_contains
from requireit import require_less_than
from requireit import require_positive
from requireit import require_sorted

from landlab.core.component_utils import iter_adaptive_time_steps
from landlab.core.component_utils import iter_time_steps

__all__ = ["Clock", "ModelRunner"]


class _TimeSteppable(Protocol):
    def update(self, dt: float) -> None: ...


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


class ModelRunner:
    """Advance a model through time and run scheduled events.

    A model runner owns the current model time. It advances the model by calling
    its ``update`` method with time-step durations no greater than the requested
    step. Scheduled events are run at their specified absolute model times.

    Parameters
    ----------
    model : _TimeSteppable
        Object with an ``update(dt)`` method that advances its state.
    clock : Clock
        Time domain and default time step for the run.
    events : mapping of str to _Event, optional
        Named events to run according to their schedules.

    Examples
    --------
    >>> from landlab.core.model_runner import ModelRunner

    >>> class MyModel:
    ...     def __init__(self):
    ...         self.elapsed = 0.0
    ...
    ...     def update(self, dt):
    ...         self.elapsed += dt
    ...
    >>> model = MyModel()
    >>> runner = ModelRunner(model, clock=Clock(start=1.0, stop=3.5, step=1.0))
    >>> runner.run()
    >>> runner.current_time
    3.5
    >>> model.elapsed
    2.5
    """

    def __init__(
        self,
        model: _TimeSteppable,
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


def _build_events(
    params: dict[str, Any],
    *,
    clock: Clock,
    actions: Mapping[str, Callable[[float], None]],
) -> dict[str, _Event]:
    require_contains(actions, required=params, name="actions")

    start, stop = clock.start, clock.stop
    events = {
        name: _Event(
            _PauseSchedule(event_config["times"], start=start, stop=stop),
            action=actions[name],
        )
        for name, event_config in params.items()
    }

    return events
