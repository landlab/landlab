from __future__ import print_function

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.layers import EventLayers


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
    assert vals == "number_of_layers: 1\nnumber_of_stacks: 5\ntracking: age"


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
