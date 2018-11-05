import numpy as np
from landlab import RasterModelGrid, HexModelGrid, CLOSED_BOUNDARY

def test_RMG_from_dict():
    params = {'shape': (10,20),
              'spacing': 25,
              'bc': {"right": "closed",
                     "top": "closed",
                     "left": "closed",
                     "bottom": "closed"},
              'origin': (35, 55),
              "axis_name": ('spam', 'eggs'),
              "axis_units": ('smoot', 'parsec')}

    mg = RasterModelGrid.from_dict(params)

    # assert things.
    assert mg.shape == mg.shape
    assert mg.dx == 25
    assert (mg.y_of_node.min(), mg.x_of_node.min()) == (35, 55)
    assert np.all(mg.status_at_node[mg.boundary_nodes]==CLOSED_BOUNDARY) == True
    assert mg.axis_units == ('smoot', 'parsec')
    assert mg.axis_name == ('spam', 'eggs')
