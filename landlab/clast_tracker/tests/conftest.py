import pytest
import numpy as np
from landlab import RasterModelGrid
from landlab.clast_tracker import ClastCollection

S=0.3

grid_east = RasterModelGrid((5,5))
_ = grid_east.add_field('node',
                          'topographic__elevation',
                          grid_east.node_x*-S+1.2)
grid_east.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                              left_is_closed=True,
                                              right_is_closed=False,
                                              top_is_closed=True)

grid_north = RasterModelGrid((5,5))
_ = grid_north.add_field('node',
                          'topographic__elevation',
                          grid_north.node_y*-S+1.2)
grid_north.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                               left_is_closed=True,
                                               right_is_closed=True,
                                               top_is_closed=False)

grid_west = RasterModelGrid((5,5))
_ = grid_west.add_field('node',
                          'topographic__elevation',
                          grid_west.node_x*S)
grid_west.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                              left_is_closed=False,
                                              right_is_closed=True,
                                              top_is_closed=True)

grid_south = RasterModelGrid((5,5))
_ = grid_south.add_field('node',
                          'topographic__elevation',
                          grid_south.node_y*S)
grid_south.set_closed_boundaries_at_grid_edges(bottom_is_closed=False,
                                               left_is_closed=True,
                                               right_is_closed=True,
                                               top_is_closed=True)


@pytest.fixture
def cc_east():
    return ClastCollection(grid_east,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)

@pytest.fixture
def cc_north():
    return ClastCollection(grid_north,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)

@pytest.fixture
def cc_west():
    return ClastCollection(grid_west,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)

@pytest.fixture
def cc_south():
    return ClastCollection(grid_south,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)



