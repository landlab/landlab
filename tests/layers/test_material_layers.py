import pytest

from landlab import RasterModelGrid
from landlab.layers import MaterialLayers


def test_MaterialLayersMixIn():
    grid = RasterModelGrid((4, 4))
    assert hasattr(grid, "material_layers")
    assert grid.material_layers.number_of_layers == 0
    assert grid.material_layers.number_of_stacks == 4


def test_adding_untracked_layer():
    layers = MaterialLayers(3)
    layers.add(1.0, type=3.0, size="sand")
    layers.add([0.0, 0.0, 1.0], type=3.0, size="sand")
    with pytest.raises(ValueError):
        layers.add([1.0], type=3.0, size="sand", spam="eggs")
