from dataclasses import FrozenInstanceError
from itertools import islice
from unittest.mock import Mock

import numpy as np
import pytest
from requireit import ValidationError

from landlab.core.model_runner import Clock
from landlab.core.model_runner import ModelRunner
from landlab.core.model_runner import _build_events
from landlab.core.model_runner import _Event
from landlab.core.model_runner import _iter_pause_times
from landlab.core.model_runner import _PauseSchedule


def test_clock_defaults():
    clock = Clock()
    assert clock.start == 0.0
    assert np.isinf(clock.stop)
    assert clock.step == 1.0
    assert np.isinf(clock.duration)


def test_clock():
    clock = Clock(start=1.0, stop=10, step=3)
    assert clock.start == 1.0
    assert clock.stop == 10
    assert clock.step == 3
    assert clock.duration == 9.0


@pytest.mark.parametrize("attr", ("start", "stop", "step"))
def test_clock_is_unchanging(attr):
    clock = Clock(start=1.0, stop=10, step=3)
    with pytest.raises(FrozenInstanceError):
        setattr(clock, attr, 999)


@pytest.mark.parametrize("step", (0.0, -1.0, np.inf, np.nan))
def test_clock_bad_step(step):
    with pytest.raises(ValueError, match="^step must"):
        Clock(step=step)


@pytest.mark.parametrize("start, stop", ((1.1, 1.0), (2.0, 2.0)))
def test_clock_start_greater_than_stop(start, stop):
    with pytest.raises(ValidationError, match="^start must be"):
        Clock(start=start, stop=stop)


def test_build_events_pairs_actions_with_schedules():
    actions = {"report": Mock(), "save": Mock()}
    events = _build_events(
        {
            "report": {"times": [1.0, 2.0]},
            "save": {"times": 0.5},
        },
        clock=Clock(start=1.0, stop=2.0),
        actions=actions,
    )

    assert set(events) == {"report", "save"}
    assert events["report"].action is actions["report"]
    assert events["report"].next_time == 1.0
    assert events["save"].action is actions["save"]
    assert events["save"].next_time == 1.0


def test_build_events_missing_action_raises():
    params = {
        "report": {"times": [1.0, 2.0]},
        "save": {"times": 0.5},
    }
    actions = {"report": Mock()}
    with pytest.raises(ValidationError, match="^actions must contain save"):
        _build_events(params, clock=Clock(start=1.0, stop=2.0), actions=actions)


def test_iter_pause_times_with_constant_interval():
    times = _iter_pause_times(2.0, start=1.0, stop=7.0)

    assert list(times) == [1.0, 3.0, 5.0, 7.0]


def test_iter_pause_times_with_explicit_times():
    times = _iter_pause_times([0.0, 1.0, 2.5, 4.0, 6.0], start=1.0, stop=4.0)

    assert list(times) == [1.0, 2.5, 4.0]


def test_iter_pause_times_can_be_unbounded():
    times = _iter_pause_times(0.5, start=1.0)
    actual = list(islice(times, 4))

    assert actual == [1.0, 1.5, 2.0, 2.5]


@pytest.mark.parametrize("interval", [0.0, -1.0])
def test_iter_pause_times_rejects_nonpositive_interval(interval):
    with pytest.raises(ValidationError, match="^pause interval must be"):
        next(_iter_pause_times(interval))


@pytest.mark.parametrize("interval", [np.inf, np.nan])
def test_iter_pause_times_rejects_nonfinite_interval(interval):
    with pytest.raises((ValueError, ValidationError), match="^pause interval must"):
        next(_iter_pause_times(interval))


@pytest.mark.parametrize("schedule", [[0.0, 2.0, 1.0], [0.0, 1.0, 1.0]])
def test_iter_pause_times_requires_strictly_increasing_times(schedule):
    with pytest.raises(ValidationError, match="^schedule must be"):
        next(_iter_pause_times(schedule))


def test_pause_schedule_starts_at_first_pause():
    schedule = _PauseSchedule([0.0, 1.0, 2.5])

    assert schedule.next_pause == 0.0


def test_pause_schedule_reports_when_pause_is_due():
    schedule = _PauseSchedule([1.0, 2.0])

    assert not schedule.is_due(0.5)
    assert schedule.is_due(1.0)
    assert schedule.is_due(1.5)


def test_pause_schedule_advance_returns_next_pause():
    schedule = _PauseSchedule([1.0, 2.0])

    assert schedule.advance() == 2.0
    assert schedule.next_pause == 2.0


def test_pause_schedule_is_infinite_when_exhausted():
    schedule = _PauseSchedule([1.0])

    assert np.isinf(schedule.advance())
    assert np.isinf(schedule.next_pause)
    assert not schedule.is_due(1e9)


def test_pause_schedule_honors_start_and_stop():
    schedule = _PauseSchedule([0.0, 1.0, 2.0, 3.0], start=1.0, stop=2.0)

    assert schedule.next_pause == 1.0
    assert schedule.advance() == 2.0
    assert np.isinf(schedule.advance())


def test_pause_schedule_with_constant_interval():
    schedule = _PauseSchedule(0.5, start=1.0, stop=2.0)

    assert schedule.next_pause == 1.0
    assert schedule.advance() == 1.5
    assert schedule.advance() == 2.0
    assert np.isinf(schedule.advance())


def test_event_reports_next_scheduled_time():
    event = _Event(_PauseSchedule([1.0, 2.0]), action=lambda time: None)

    assert event.next_time == 1.0


def test_event_does_not_run_before_it_is_due():
    action = Mock()
    event = _Event(_PauseSchedule([1.0, 2.0]), action=action)

    event.run_if_due(0.5)

    action.assert_not_called()
    assert event.next_time == 1.0


def test_event_runs_action_when_due():
    action = Mock()
    event = _Event(_PauseSchedule([1.0, 2.0]), action=action)

    event.run_if_due(1.0)

    action.assert_called_once_with(1.0)


def test_event_advances_after_running_action():
    event = _Event(_PauseSchedule([1.0, 2.0]), action=lambda time: None)

    event.run_if_due(1.0)

    assert event.next_time == 2.0


def test_event_is_exhausted_after_its_last_action():
    action = Mock()
    event = _Event(_PauseSchedule([1.0]), action=action)

    event.run_if_due(1.0)
    event.run_if_due(2.0)

    action.assert_called_once_with(1.0)
    assert np.isinf(event.next_time)


def test_model_runner_uses_clock():
    runner = ModelRunner(Mock(), clock=Clock(start=1.0, stop=5.0, step=0.25))

    assert runner.current_time == 1.0
    assert runner.run_duration == 4.0
    assert runner.dt == 0.25


def test_model_runner_update_until_advances_model_and_time():
    model = Mock()
    runner = ModelRunner(model, clock=Clock(start=1.0, stop=5.0, step=1.0))

    runner.update_until(3.5, dt=1.0)

    actual_steps = [call.args[0] for call in model.update.call_args_list]

    assert len(actual_steps)
    assert actual_steps == pytest.approx([2.5 / 3.0] * 3)
    assert runner.current_time == 3.5


def test_model_runner_update_until_ignores_past_time():
    model = Mock()
    runner = ModelRunner(model, clock=Clock(start=1.0, stop=5.0))

    runner.update_until(0.5, dt=0.25)

    model.update.assert_not_called()
    assert runner.current_time == 1.0


def test_model_runner_run_without_events():
    model = Mock()
    runner = ModelRunner(model, clock=Clock(start=1.0, stop=3.5, step=1.0))

    runner.run()

    actual_steps = [call.args[0] for call in model.update.call_args_list]

    assert len(actual_steps)
    assert actual_steps == pytest.approx([2.5 / 3.0] * 3)
    assert runner.current_time == 3.5


def test_model_runner_run_stops_at_clock_stop_after_partial_update():
    model = Mock()
    runner = ModelRunner(model, clock=Clock(start=1.0, stop=5.0, step=1.0))
    runner.update_until(2.0, dt=1.0)
    assert runner.current_time == 2.0

    model.reset_mock()

    runner.run()

    actual_steps = [call.args[0] for call in model.update.call_args_list]
    assert actual_steps == [1.0, 1.0, 1.0]
    # assert model.update.call_args_list == [call(1.0), call(1.0), call(1.0)]
    assert runner.current_time == 5.0


def test_model_runner_runs_scheduled_events():
    model = Mock()
    action = Mock()
    event = _Event(_PauseSchedule([1.0, 2.0, 3.0]), action=action)
    runner = ModelRunner(
        model,
        clock=Clock(start=1.0, stop=3.0, step=0.75),
        events={"report": event},
    )

    runner.run()

    actual_times = [call.args[0] for call in action.call_args_list]
    assert actual_times == [1.0, 2.0, 3.0]
    assert runner.current_time == 3.0


def test_model_runner_copies_events():
    action = Mock()
    events = {"report": _Event(_PauseSchedule([1.0]), action=action)}
    runner = ModelRunner(Mock(), clock=Clock(start=1.0, stop=2.0), events=events)
    events.clear()

    runner.run()

    action.assert_called_once_with(1.0)
