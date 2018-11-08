from __future__ import print_function
import pytest

from landlab.layers import MaterialLayers
from landlab import RasterModelGrid

def test_MaterialLayersMixIn():
    mg = RasterModelGrid(4,4)
    ml = mg.material_layers
    assert hasattr(mg, 'material_layers') == True
    assert ml.number_of_layers == 0
    assert ml.number_of_stacks == 4

def test_adding_untracked_layer():
    layers = MaterialLayers(3)
    layers.add(1., type=3., size='sand')
    layers.add([0., 0., 1.], type=3., size='sand')
    with pytest.raises(ValueError):
        layers.add([1.], type=3., size='sand', spam='eggs')
