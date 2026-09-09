from contextlib import chdir
from dataclasses import FrozenInstanceError
from itertools import islice
from unittest.mock import Mock
from unittest.mock import patch

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from requireit import ValidationError

from landlab import RasterModelGrid
from landlab.io.native_landlab import save_grid
from landlab.utils.model_base import Clock
from landlab.utils.model_base import LandlabModel
from landlab.utils.model_base import _build_events
from landlab.utils.model_base import _Event
from landlab.utils.model_base import _FilenameSequence
from landlab.utils.model_base import _GridSaver
from landlab.utils.model_base import _iter_pause_times
from landlab.utils.model_base import _PauseSchedule
from landlab.utils.model_base import merge_params
from landlab.utils.model_base import resolve_array_filepaths
from landlab.utils.model_base import setup_grid


@pytest.fixture
def model_params():
    return {
        "grid": {
            "source": "create",
            "create_grid": {
                "RasterModelGrid": {
                    "shape": (4, 5),
                    "xy_spacing": (2.0, 4.0),
                }
            },
        },
        "clock": {"start": 1.0, "stop": 5.0, "step": 0.5},
        "output": {
            "plot_times": 10.0,
            "save_times": 10.0,
            "report_times": 10.0,
            "save_path": "model-output",
            "clobber": True,
            "fields": None,
            "plot_to_file": False,
        },
    }


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


def test_build_events_pairs_actions_with_schedules():
    actions = {"report": Mock(), "save": Mock()}
    events = _build_events(
        {"report_times": [1.0, 2.0], "save_times": 0.5},
        clock=Clock(start=1.0, stop=2.0),
        actions=actions,
    )

    assert set(events) == {"report", "save"}
    assert events["report"].action is actions["report"]
    assert events["report"].next_time == 1.0
    assert events["save"].action is actions["save"]
    assert events["save"].next_time == 1.0


@pytest.mark.parametrize("base", ("foobar", "foo.bar", ""))
@pytest.mark.parametrize("ext", (".nc", "vtk", ""))
def test_filename_sequence_numbers_files(base, ext):
    filenames = _FilenameSequence(base, ndigits=4, ext=ext)

    assert list(islice(filenames, 4)) == [
        f"{base}0001{ext}",
        f"{base}0002{ext}",
        f"{base}0003{ext}",
        f"{base}0004{ext}",
    ]


def test_filename_sequence_defaults_to_unpadded_names_without_extension():
    filenames = _FilenameSequence("foo")

    assert next(filenames) == "foo1"
    assert next(filenames) == "foo2"


def test_filename_sequence_is_an_iterator():
    filenames = _FilenameSequence("frog")

    assert iter(filenames) is filenames


def test_filename_sequence_padding_is_a_minimum_width():
    filenames = _FilenameSequence("f", ndigits=1)

    assert list(islice(filenames, 10)) == [
        "f1",
        "f2",
        "f3",
        "f4",
        "f5",
        "f6",
        "f7",
        "f8",
        "f9",
        "f10",
    ]


def test_filename_sequence_rejects_negative_ndigits():
    with pytest.raises(ValidationError, match="^ndigits must be"):
        _FilenameSequence("frame", ndigits=-1)


@pytest.mark.parametrize(
    ("fmt", "writer_name", "ext"),
    [
        ("grid", "save_grid", ".grid"),
        ("netcdf", "write_netcdf", ".nc"),
        ("vtk", "write_legacy_vtk", ".vtk"),
    ],
)
def test_grid_saver_uses_writer_for_format(fmt, writer_name, ext):
    grid = RasterModelGrid((3, 4))
    saver = _GridSaver(grid, "foo-output", fmt=fmt, ndigits=3)

    with patch(f"landlab.utils.model_base.{writer_name}") as writer:
        filename = saver.save()

    assert filename == f"foo-output001{ext}"
    writer.assert_called_once()


def test_grid_saver_advances_filename():
    saver = _GridSaver(RasterModelGrid((3, 4)), "frame", ndigits=2)

    with patch("landlab.utils.model_base.save_grid"):
        assert saver.save() == "frame01.grid"
        assert saver.save() == "frame02.grid"


def test_grid_saver_is_callable():
    saver = _GridSaver(RasterModelGrid((3, 4)), "frame")

    with patch("landlab.utils.model_base.save_grid") as writer:
        assert saver(10.0) is None
        assert writer.call_count == 1


def test_grid_saver_rejects_unknown_format():
    with pytest.raises(ValidationError, match="^fmt must be one of"):
        _GridSaver(RasterModelGrid((3, 4)), "frame", fmt="foobar")


def test_merge_params():
    user = {"a": 1, "dict": {"user": 2}}
    defaults = {"a": 0, "b": 3, "dict": {"default": 4}}

    actual = merge_params(user, defaults=defaults)

    assert actual == {
        "a": 1,
        "b": 3,
        "dict": {"user": 2, "default": 4},
    }


def test_merge_params_copies_nested_dicts():
    user = {
        "user_only": {"nested": {"value": 1}},
        "merged": {"user": {"value": 2}},
    }
    defaults = {
        "default_only": {"nested": {"value": 3}},
        "merged": {"default": {"value": 4}},
    }

    actual = merge_params(user, defaults=defaults)

    assert actual == {
        "user_only": {"nested": {"value": 1}},
        "default_only": {"nested": {"value": 3}},
        "merged": {"user": {"value": 2}, "default": {"value": 4}},
    }
    assert actual["user_only"] is not user["user_only"]
    assert actual["user_only"]["nested"] is not user["user_only"]["nested"]
    assert actual["default_only"] is not defaults["default_only"]
    assert actual["default_only"]["nested"] is not defaults["default_only"]["nested"]
    assert actual["merged"] is not user["merged"]
    assert actual["merged"]["user"] is not user["merged"]["user"]
    assert actual["merged"]["default"] is not defaults["merged"]["default"]


def test_merge_params_does_not_merge_grid_dict():
    user = {"grid": {"RasterModelGrid": {"shape": (3, 4)}}}
    defaults = {"grid": {"HexModelGrid": {"shape": (5, 6)}}}

    actual = merge_params(user, defaults=defaults)

    assert actual["grid"] == user["grid"]
    assert actual["grid"] is not user["grid"]
    assert actual["grid"]["RasterModelGrid"] is not user["grid"]["RasterModelGrid"]


def test_merge_params_preserves_grid_instance():
    grid = RasterModelGrid((3, 4))

    params = {"grid": {"source": "grid_object", "grid_object": grid}}
    actual = merge_params(params)

    assert actual["grid"]["grid_object"] is grid
    assert actual["grid"] is not params["grid"]


def test_merge_params_dict_overrides_non_dict_default():
    actual = merge_params({"value": {"dict": 1}}, defaults={"value": 0})

    assert actual == {"value": {"dict": 1}}


def test_model_init_uses_in_memory_grid_and_params(model_params):
    clock = Clock(start=0.0, stop=100.0, step=0.25)
    model_params.pop("grid")

    grid = RasterModelGrid((3, 4))

    model = LandlabModel(grid, clock=clock, params=model_params)

    assert model.grid is grid
    assert model.params is model_params
    assert model.current_time == 0.0
    assert model.run_duration == 100.0
    assert model.dt == 0.25


def test_model_from_params(model_params):
    model_params["grid"]["create_grid"]["RasterModelGrid"] = {
        "shape": (40, 50),
        "xy_spacing": (0.5, 8.0),
        "xy_of_lower_left": (-16.0, 32.0),
    }
    model = LandlabModel.from_params(model_params)

    assert isinstance(model.grid, RasterModelGrid)
    assert model.grid.shape == (40, 50)
    assert model.grid.dx == 0.5
    assert model.grid.dy == 8.0
    assert model.grid.xy_of_lower_left == (-16.0, 32.0)
    assert model.params == model_params
    assert model.params is not model_params


def test_model_from_params_uses_defaults():
    model = LandlabModel.from_params()

    assert isinstance(model.grid, RasterModelGrid)
    assert model.grid.shape == (5, 5)
    assert model.current_time == 0.0
    assert model.run_duration == 2.0
    assert model.dt == 1.0


def test_model_from_params_returns_subclass(model_params):
    class FrogModel(LandlabModel):
        pass

    assert isinstance(FrogModel.from_params(model_params), FrogModel)


def test_model_from_file(tmp_path):
    input_file = tmp_path / "model.yaml"
    input_file.write_text("""
grid:
  source: create
  create_grid:
    RasterModelGrid:
      shape: [3, 4]
      xy_spacing: [2.0, 4.0]
clock:
  start: 2.0
  stop: 8.0
  step: 0.25
""")

    model = LandlabModel.from_file(input_file)

    assert isinstance(model.grid, RasterModelGrid)
    assert model.grid.shape == (3, 4)
    assert model.grid.dx == 2.0
    assert model.grid.dy == 4.0
    assert model.current_time == 2.0
    assert model.run_duration == 6.0
    assert model.dt == 0.25


def test_setup_grid_creates_grid():
    params = {
        "source": "create",
        "create_grid": {
            "RasterModelGrid": {
                "shape": (4, 5),
                "xy_spacing": (2.0, 4.0),
            }
        },
    }

    grid = setup_grid(params)

    assert isinstance(grid, RasterModelGrid)
    assert grid.shape == (4, 5)
    assert grid.dx == 2.0
    assert grid.dy == 4.0


def test_setup_grid_loads_grid(tmp_path):
    original = RasterModelGrid((3, 4), xy_spacing=(2.0, 4.0))
    original.add_zeros("topographic__elevation", at="node")
    path = tmp_path / "model.grid"
    save_grid(original, path)

    actual = setup_grid({"source": "file", "grid_file_name": path})

    assert isinstance(actual, RasterModelGrid)
    assert actual.shape == original.shape
    assert actual.spacing == original.spacing
    assert "topographic__elevation" in actual.at_node


def test_setup_grid_uses_existing_grid():
    expected = RasterModelGrid((3, 4))

    actual = setup_grid({"source": "grid_object", "grid_object": expected})

    assert actual is expected


def test_setup_grid_rejects_unknown_source():
    with pytest.raises(ValidationError, match="^source must be one of"):
        setup_grid({"source": "unknown"})


def test_setup_grid_rejects_non_grid_object():
    with pytest.raises(ValueError, match="^grid source must be"):
        setup_grid({"source": "grid_object", "grid_object": object()})


def test_resolve_array_filepaths(tmp_path):
    expected_1d = np.arange(3) / 2
    expected_2d = np.arange(4).reshape((2, 2)) / 2
    expected_col = np.arange(10).reshape((-1, 1)) / 4

    np.save(tmp_path / "test1", expected_1d)
    np.save(tmp_path / "test2", expected_2d)
    np.save(tmp_path / "test3", expected_col)

    p = {
        "a": 123,
        "b": {"c": 456, "d": {"_filepath": "test1.npy"}},
        "e": {"_filepath": "test2.npy"},
        "f": {"_filepath": "test3.npy"},
    }
    with chdir(tmp_path):
        actual = resolve_array_filepaths(p)

    assert actual is not p
    assert actual["b"] is not p["b"]
    assert p["b"]["d"] == {"_filepath": "test1.npy"}
    assert p["e"] == {"_filepath": "test2.npy"}
    assert p["f"] == {"_filepath": "test3.npy"}
    assert_array_equal(actual["b"]["d"], expected_1d)
    assert_array_equal(actual["e"], expected_2d)
    assert_array_equal(actual["f"], expected_col)


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
