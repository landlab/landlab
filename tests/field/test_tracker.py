import pytest
from numpy.testing import assert_array_equal

from landlab.field.errors import BadOperationError
from landlab.field.errors import MissingTrackedFieldError
from landlab.field.errors import NothingToTrackError
from landlab.field.tracker import FieldTracker
from landlab.field.tracker import open_tracker
from landlab.grid.raster import RasterModelGrid


def test_tracker_open_close():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("z", at="node")

    tracker = FieldTracker(grid)

    names = tracker.open()
    assert tracker.tracking == names
    grid.at_node["z"] += 4.0
    tracker.close()

    assert tracker.tracking is None

    assert names == ("at_node:z",)

    actual = dict(tracker.fields)
    expected = {
        "at_node:z": (
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
        )
    }
    assert actual.keys() == expected.keys()
    assert len(actual["at_node:z"]) == len(expected["at_node:z"])
    assert_array_equal(actual["at_node:z"][0], expected["at_node:z"][0])
    assert_array_equal(actual["at_node:z"][1], expected["at_node:z"][1])


def test_tracker_open_context():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("z", at="node")

    with open_tracker(grid) as tracker:
        assert tracker.tracking == ("at_node:z",)
        grid.at_node["z"] += 4.0
    assert tracker.tracking is None

    actual = dict(tracker.fields)
    expected = {
        "at_node:z": (
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
        )
    }
    assert actual.keys() == expected.keys()
    assert len(actual["at_node:z"]) == len(expected["at_node:z"])
    assert_array_equal(actual["at_node:z"][0], expected["at_node:z"][0])
    assert_array_equal(actual["at_node:z"][1], expected["at_node:z"][1])


def test_tracker_checkpoint():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("z", at="node")

    with open_tracker(grid) as tracker:
        grid.at_node["z"] += 4.0
        tracker.checkpoint()
        grid.at_node["z"] += 2.0

    actual = dict(tracker.fields)
    expected = {
        "at_node:z": (
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
            [7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0],
        )
    }
    assert actual.keys() == expected.keys()
    assert len(actual["at_node:z"]) == len(expected["at_node:z"])
    assert_array_equal(actual["at_node:z"][0], expected["at_node:z"][0])
    assert_array_equal(actual["at_node:z"][1], expected["at_node:z"][1])
    assert_array_equal(actual["at_node:z"][2], expected["at_node:z"][2])


def test_nothing_to_track_error():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("z", at="node")

    with pytest.raises(NothingToTrackError):
        with open_tracker(grid, include="at_link*"):
            pass


def test_bad_operation_error():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("z", at="node")

    tracker = FieldTracker(grid)
    with pytest.raises(BadOperationError):
        tracker.checkpoint()

    with open_tracker(grid) as tracker:
        pass
    with pytest.raises(BadOperationError):
        tracker.checkpoint()


def test_ok_to_close_unopened():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("z", at="node")

    tracker = FieldTracker(grid)
    assert tracker.tracking is None
    tracker.close()

    with open_tracker(grid) as tracker:
        pass
    assert tracker.tracking is None
    tracker.close()


def test_missing_field_error():
    grid = RasterModelGrid((3, 4))
    tracker = FieldTracker(grid)

    grid.add_ones("z", at="node")
    with pytest.raises(MissingTrackedFieldError):
        with open_tracker(grid) as tracker:
            grid.at_node.pop("z")
    assert tracker.tracking is None

    grid.add_ones("z", at="node")
    tracker.open()
    grid.at_node.pop("z")
    with pytest.raises(MissingTrackedFieldError):
        tracker.close()
    assert tracker.tracking is None


def test_reduce():
    grid = RasterModelGrid((3, 4))
    grid.add_ones("z", at="node")

    with open_tracker(grid) as tracker:
        grid.at_node["z"] += 5.0

    def reducer(x, y):
        return x - y

    expected = {
        "at_node:z": [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    }

    actual = dict(tracker.reduce(reducer))
    assert actual.keys() == expected.keys()
    assert_array_equal(actual["at_node:z"], expected["at_node:z"])

    with open_tracker(grid) as tracker:
        grid.at_node["z"] += 2.0
        tracker.checkpoint()
        grid.at_node["z"] += 3.0

    actual = dict(tracker.reduce(reducer))
    assert actual.keys() == expected.keys()
    assert_array_equal(actual["at_node:z"], expected["at_node:z"])
