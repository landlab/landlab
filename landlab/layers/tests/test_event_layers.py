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
    layers.add(1., age=3.)
    layers.add(2., age=4.)

    truth = np.array([[3., 3., 3., 3., 3.], [4., 4., 4., 4., 4.]])
    assert_array_equal(layers["age"], truth)

    layers["age"] = 2.
    truth = np.array([[2., 2., 2., 2., 2.], [2., 2., 2., 2., 2.]])
    assert_array_equal(layers["age"], truth)


def test_set_item_with_1d():
    layers = EventLayers(5)
    layers.add(1., age=3.)
    layers.add(2., age=4.)

    truth = np.array([[3., 3., 3., 3., 3.], [4., 4., 4., 4., 4.]])
    assert_array_equal(layers["age"], truth)

    layers["age"] = [4., 7.]

    truth = np.array([[4., 4., 4., 4., 4.], [7., 7., 7., 7., 7.]])
    assert_array_equal(layers["age"], truth)


def test_set_item_with_2d():
    layers = EventLayers(5)
    layers.add(1., age=3.)
    layers.add(2., age=4.)

    truth = np.array([[3., 3., 3., 3., 3.], [4., 4., 4., 4., 4.]])
    assert_array_equal(layers["age"], truth)

    layers["age"] = [[4., 4., 4., 4., 4.], [7., 7., 7., 7., 7.]]

    truth = np.array([[4., 4., 4., 4., 4.], [7., 7., 7., 7., 7.]])
    assert_array_equal(layers["age"], truth)


def test__str__():
    layers = EventLayers(5)
    layers.add(1., age=3.)
    vals = str(layers)
    assert vals == "number_of_layers: 1\nnumber_of_stacks: 5\ntracking: age"


def test__repr__():
    layers = EventLayers(5)
    layers.add(1., age=3.)
    vals = repr(layers)
    assert vals == "EventLayers(5)"


def test_adding_untracked_layer():
    layers = EventLayers(3)
    layers.add(1., type=3., size="sand")
    layers.add([0., 0., 1.], type=3., size="sand")
    with pytest.raises(ValueError):
        layers.add([1.], type=3., size="sand", spam="eggs")
