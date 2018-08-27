from __future__ import print_function
import pytest

from landlab.layers import MaterialLayers, MaterialLayersMixIn

def test_MaterialLayersMixIn():
    pass


def test_adding_ignored_layer():
    layers = MaterialLayers(3)
    layers.add(1., type=3., size='sand')
    layers.add([0., 0., 1.], type=3., size='sand')
    with pytest.raises(ValueError):
        layers.add([1.], type=3., size='sand', spam='eggs')
