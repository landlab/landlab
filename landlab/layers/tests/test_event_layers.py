from __future__ import print_function

import numpy as np
from numpy.testing import assert_array_equal

from landlab.layers import EventLayersMixIn, EventLayers


def test_EventLayersMixIn():
    pass


def test_setitem_with_scalar():
    layers = EventLayers(5)
    layers.add(1., age=3.)
    layers.add(2., age=4.)

    truth = np.array([[ 3.,  3.,  3.,  3.,  3.],
                      [ 4.,  4.,  4.,  4.,  4.]])
    assert_array_equal(layers['age'], truth)

    layers['age'] = 2.
    truth = np.array([[ 2.,  2.,  2.,  2.,  2.],
                      [ 2.,  2.,  2.,  2.,  2.]])
    assert_array_equal(layers['age'], truth)


def test_set_item_with_1d():
    layers = EventLayers(5)
    layers.add(1., age=3.)
    layers.add(2., age=4.)

    truth = np.array([[ 3.,  3.,  3.,  3.,  3.],
                      [ 4.,  4.,  4.,  4.,  4.]])
    assert_array_equal(layers['age'], truth)

    layers['age'] = [4., 7.]

    truth = np.array([[ 4.,  4.,  4.,  4.,  4.],
                      [ 7.,  7.,  7.,  7.,  7.]])
    assert_array_equal(layers['age'], truth)


def test_set_item_with_2d():
    layers = EventLayers(5)
    layers.add(1., age=3.)
    layers.add(2., age=4.)

    truth = np.array([[ 3.,  3.,  3.,  3.,  3.],
                      [ 4.,  4.,  4.,  4.,  4.]])
    assert_array_equal(layers['age'], truth)

    layers['age'] = [[ 4.,  4.,  4.,  4.,  4.],
                     [ 7.,  7.,  7.,  7.,  7.]]

    truth = np.array([[ 4.,  4.,  4.,  4.,  4.],
                      [ 7.,  7.,  7.,  7.,  7.]])
    assert_array_equal(layers['age'], truth)
