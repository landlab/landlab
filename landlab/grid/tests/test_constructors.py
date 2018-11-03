import numpy as np
from landlab import RasterModelGrid, HexModelGrid, CLOSED_BOUNDARY

def test_RMG_from_dict():
    params = {'shape': (10,20),
              'spacing': 25,
              'bc': {"right": "closed",
                     "top": "closed",
                     "left": "closed",
                     "bottom": "closed"},
               'origin': (35, 55)}

    mg = RasterModelGrid.from_dict(params)

    # assert things.
    assert mg.shape == mg.shape
    assert mg.dx == 25
    assert (mg.x_of_node.min(), mg.y_of_node) == (35, 55)
    assert np.all(mg.status_at_node[mg.boundary_nodes]==CLOSED_BOUNDARY) == True
