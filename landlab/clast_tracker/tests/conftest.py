import pytest
import numpy as np
from landlab import RasterModelGrid
from landlab.clast_tracker import ClastCollection

S=0.3

@pytest.fixture
def grid_east():
    grid_east = RasterModelGrid((5,5))
    grid_east.add_field('node',
                        'topographic__elevation',
                        grid_east.node_x*-S+1.2);
    grid_east.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                                  left_is_closed=True,
                                                  right_is_closed=False,
                                                  top_is_closed=True)
    return grid_east

@pytest.fixture
def grid_north():
    grid_north = RasterModelGrid((5,5))
    grid_north.add_field('node',
                         'topographic__elevation',
                         grid_north.node_y*-S+1.2);
    grid_north.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                                   left_is_closed=True,
                                                   right_is_closed=True,
                                                   top_is_closed=False)
    return grid_north

@pytest.fixture
def grid_west():
    grid_west = RasterModelGrid((5,5))
    grid_west.add_field('node',
                        'topographic__elevation',
                        grid_west.node_x*S);
    grid_west.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                                  left_is_closed=False,
                                                  right_is_closed=True,
                                                  top_is_closed=True)
    return grid_west

@pytest.fixture
def grid_south():
    grid_south = RasterModelGrid((5,5))
    grid_south.add_field('node',
                         'topographic__elevation',
                         grid_south.node_y*S);
    grid_south.set_closed_boundaries_at_grid_edges(bottom_is_closed=False,
                                                   left_is_closed=True,
                                                   right_is_closed=True,
                                                   top_is_closed=True)
    return grid_south



#from landlab import imshow_grid
#%matplotlib auto 
#imshow_grid(grid_east, 'topographic__elevation')





@pytest.fixture
def cc_east(grid_east):
    return ClastCollection(grid_east,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)

@pytest.fixture
def cc_north(grid_north):
    return ClastCollection(grid_north,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)

@pytest.fixture
def cc_west(grid_west):
    return ClastCollection(grid_west,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)

@pytest.fixture
def cc_south(grid_south):
    return ClastCollection(grid_south,
                           clast_x=np.ones(1)*2,
                           clast_y=np.ones(1)*2,
                           clast_elev=np.ones(1)*6,
                           clast_radius=np.ones(1)*0.5)



