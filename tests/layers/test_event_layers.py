import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.layers import EventLayers
from landlab.layers.eventlayers import _BlockSlice
from landlab.layers.eventlayers import _valid_keywords_or_raise


def test_EventLayersMixIn():
    grid = RasterModelGrid((4, 4))
    assert hasattr(grid, "event_layers")
    assert grid.event_layers.number_of_layers == 0
    assert grid.event_layers.number_of_stacks == 4


def test_setitem_with_scalar():
    layers = EventLayers(5)
    layers.add(1.0, age=3.0)
    layers.add(2.0, age=4.0)

    truth = np.array([[3.0, 3.0, 3.0, 3.0, 3.0], [4.0, 4.0, 4.0, 4.0, 4.0]])
    assert_array_equal(layers["age"], truth)

    layers["age"] = 2.0
    truth = np.array([[2.0, 2.0, 2.0, 2.0, 2.0], [2.0, 2.0, 2.0, 2.0, 2.0]])
    assert_array_equal(layers["age"], truth)


def test_set_item_with_1d():
    layers = EventLayers(5)
    layers.add(1.0, age=3.0)
    layers.add(2.0, age=4.0)

    truth = np.array([[3.0, 3.0, 3.0, 3.0, 3.0], [4.0, 4.0, 4.0, 4.0, 4.0]])
    assert_array_equal(layers["age"], truth)

    layers["age"] = [4.0, 7.0]

    truth = np.array([[4.0, 4.0, 4.0, 4.0, 4.0], [7.0, 7.0, 7.0, 7.0, 7.0]])
    assert_array_equal(layers["age"], truth)


def test_set_item_with_2d():
    layers = EventLayers(5)
    layers.add(1.0, age=3.0)
    layers.add(2.0, age=4.0)

    truth = np.array([[3.0, 3.0, 3.0, 3.0, 3.0], [4.0, 4.0, 4.0, 4.0, 4.0]])
    assert_array_equal(layers["age"], truth)

    layers["age"] = [[4.0, 4.0, 4.0, 4.0, 4.0], [7.0, 7.0, 7.0, 7.0, 7.0]]

    truth = np.array([[4.0, 4.0, 4.0, 4.0, 4.0], [7.0, 7.0, 7.0, 7.0, 7.0]])
    assert_array_equal(layers["age"], truth)


def test__str__():
    layers = EventLayers(5)
    layers.add(1.0, age=3.0)
    vals = str(layers)
    assert vals.splitlines() == [
        "number_of_layers: 1",
        "number_of_stacks: 5",
        "tracking: age",
    ]


def test__repr__():
    layers = EventLayers(5)
    layers.add(1.0, age=3.0)
    vals = repr(layers)
    assert vals == "EventLayers(5)"


def test_adding_untracked_layer():
    layers = EventLayers(3)
    layers.add(1.0, type=3.0, size="sand")
    layers.add([0.0, 0.0, 1.0], type=3.0, size="sand")
    with pytest.raises(ValueError):
        layers.add([1.0], type=3.0, size="sand", spam="eggs")


def test_reduce_with_no_args():
    layers = EventLayers(3)

    layers.add(1.5)
    layers.reduce()
    assert_array_equal(layers.dz, [[1.5, 1.5, 1.5]])

    layers.add(2.5)
    layers.reduce()
    assert_array_equal(layers.dz, [[4.0, 4.0, 4.0]])


def test_reduce_with_start():
    layers = EventLayers(3)

    layers.add([3.0, 4.0, 5.0])
    layers.add([1.0, 2.0, 0.5])
    layers.add([1.0, 2.0, 0.5])
    layers.reduce(1, 3)

    assert_array_equal(layers.dz, [[3.0, 4.0, 5.0], [2.0, 4.0, 1.0]])


def test_reduce_with_start_and_stop():
    layers = EventLayers(3)

    layers.add([3.0, 4.0, 5.0])
    layers.add([1.0, 2.0, 0.5])
    layers.add([1.0, 2.0, 0.5])
    layers.add([2.0, 5.0, 6.0])

    layers.reduce(1, 3)
    assert_array_equal(layers.dz, [[3.0, 4.0, 5.0], [2.0, 4.0, 1.0], [2.0, 5.0, 6.0]])


def test_reduce_with_start_and_stop_and_step():
    layers = EventLayers(3)

    layers.add([3.0, 4.0, 5.0])
    layers.add([1.0, 2.0, 0.5])
    layers.add([1.0, 2.0, 0.5])
    layers.add([2.0, 5.0, 6.0])

    layers.reduce(1, 3, 2)
    assert_array_equal(layers.dz, [[3.0, 4.0, 5.0], [2.0, 4.0, 1.0], [2.0, 5.0, 6.0]])


def test_reduce_with_stop_equals_start():
    layers = EventLayers(3)

    layers.add([3.0, 4.0, 5.0])
    layers.add([1.0, 2.0, 0.5])
    layers.add([1.0, 2.0, 0.5])
    layers.add([2.0, 5.0, 6.0])

    layers.reduce(1, 1)
    assert_array_equal(
        layers.dz, [[3.0, 4.0, 5.0], [1.0, 2.0, 0.5], [1.0, 2.0, 0.5], [2.0, 5.0, 6.0]]
    )


def test_reduce_with_stop_less_than_start():
    layers = EventLayers(3)

    layers.add([3.0, 4.0, 5.0])
    layers.add([1.0, 2.0, 0.5])
    layers.add([1.0, 2.0, 0.5])
    layers.add([2.0, 5.0, 6.0])

    with pytest.raises(ValueError):
        layers.reduce(3, 1)


def test_reduce_with_one_layer():
    layers = EventLayers(3)
    layers.add(1.0)
    layers.reduce()
    assert_array_equal(layers.dz, [[1, 1, 1]])


def test_reduce_with_no_layers():
    layers = EventLayers(3)
    layers.reduce()
    assert_array_equal(layers.dz, np.empty((0, 3)))


def test_reduce_with_attrs():
    layers = EventLayers(3)
    layers.add([1, 1, 1], age=0.0)
    layers.add([1, 2, 5], age=1.0)
    layers.add([2, 2, 2], age=2.0)
    layers.reduce(age=np.sum)
    assert_array_equal(layers["age"], [[3.0, 3.0, 3.0]])


def test_reduce_with_reducer():
    layers = EventLayers(3)
    layers.add([2, 2, 2], age=1.0)
    layers.add([2, 2, 2], age=3.0)
    layers.add([2, 2, 2], age=4.0)
    layers.reduce(1, 3, age=np.mean)
    assert_array_equal(layers["age"], [[1.0, 1.0, 1.0], [3.5, 3.5, 3.5]])


def test_reduce_with_start_too_big():
    layers = EventLayers(3)
    layers.add(1)
    layers.reduce(2, 4)
    assert_array_equal(layers.dz, [[1, 1, 1]])


def test_reduce_with_unknown_property():
    layers = EventLayers(3)
    for layer in range(4):
        layers.add(layer)
    with pytest.raises(TypeError):
        layers.reduce(1, 3, age=np.mean)


def test_reduce_with_missing_property():
    layers = EventLayers(3)
    for layer in range(4):
        layers.add(layer, age=layer)
    with pytest.raises(TypeError):
        layers.reduce(1, 3)


@pytest.mark.parametrize("args", [(), (4,), (0, 4), (0, 4, 4), (-4, 4, 4), (4, 0, -4)])
def test_reduce_args(args):
    layers = EventLayers(3)
    for layer in range(4):
        layers.add(layer)
    layers.reduce(4)
    assert_array_equal(layers.dz, [[6, 6, 6]])


def test_reduce_partial_block():
    layers = EventLayers(3)
    for layer in range(6):
        layers.add(layer)
    layers.reduce(1, 6, 4)
    assert_array_equal(layers.dz, [[0, 0, 0], [10, 10, 10], [5, 5, 5]])


def test_block_slice_no_args():
    block = _BlockSlice()
    assert block.start == 0
    assert block.stop is None
    assert block.step is None
    assert block.indices(4) == (0, 4, 4)


def test_block_slice_one_arg():
    block = _BlockSlice(4)
    assert block.start == 0
    assert block.stop == 4
    assert block.step is None
    assert block.indices(4) == (0, 4, 4)
    assert block.indices(3) == (0, 3, 3)


def test_block_slice_two_args():
    block = _BlockSlice(1, 4)
    assert block.start == 1
    assert block.stop == 4
    assert block.step is None
    assert block.indices(4) == (1, 4, 3)
    assert block.indices(5) == (1, 4, 3)
    assert block.indices(3) == (1, 3, 2)


def test_block_slice_three_args():
    block = _BlockSlice(1, 7, 2)
    assert block.start == 1
    assert block.stop == 7
    assert block.step == 2
    assert block.indices(8) == (1, 7, 2)

    block = _BlockSlice(1, 7, 7)
    assert block.start == 1
    assert block.stop == 7
    assert block.step == 6
    assert block.indices(8) == (1, 7, 6)


@pytest.mark.parametrize("args", [(0, 7, 2), (0, 8, 2), (0, 8, 21)])
@pytest.mark.parametrize("n_rows", [7, 8])
def test_block_slice_full_blocks(args, n_rows):
    start, stop, step = _BlockSlice(*args).indices(n_rows)
    assert (stop - start) % step == 0


def test_block_slice_stop_less_than_start():
    with pytest.raises(ValueError):
        _BlockSlice(5, 4)


def test_valid_keywords():
    _valid_keywords_or_raise({"foo": 0}, optional=["foo", "bar"])
    with pytest.raises(TypeError):
        _valid_keywords_or_raise({"foo": 0}, required=["foo", "bar"])
    with pytest.raises(TypeError):
        _valid_keywords_or_raise(["baz"], optional=["foo", "bar"])


@pytest.mark.parametrize("as_type", [list, tuple, set])
def test_valid_keywords_as_types(as_type):
    _valid_keywords_or_raise(as_type(["foo"]), optional=as_type(["foo", "bar"]))
    with pytest.raises(TypeError):
        _valid_keywords_or_raise(as_type(["foo"]), required=as_type(["foo", "bar"]))
    with pytest.raises(TypeError):
        _valid_keywords_or_raise(as_type(["baz"]), optional=as_type(["foo", "bar"]))


def test_valid_keywords_empty():
    _valid_keywords_or_raise([])
    _valid_keywords_or_raise([], optional=["foo"])
    with pytest.raises(TypeError):
        _valid_keywords_or_raise([], required=["foo"])
