import pytest
import numpy as np
from landlab import RasterModelGrid
from landlab.clast_tracker import ClastCollection

S = 0.3

# NO SLOPE:
@pytest.fixture
def grid_flat():
    grid_flat = RasterModelGrid((5, 5))
    grid_flat.add_field(
        "node",
        "topographic__elevation",
        np.ones(grid_flat.number_of_nodes),
    )
    grid_flat.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_flat


# SLOPES TO ROW/COL:
@pytest.fixture
def grid_east():
    grid_east = RasterModelGrid((5, 5))
    grid_east.add_field(
        "node",
        "topographic__elevation",
        grid_east.node_x * -S + 1.2,
    )
    grid_east.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=False,
        top_is_closed=True,
    )
    return grid_east


@pytest.fixture
def grid_north():
    grid_north = RasterModelGrid((5, 5))
    grid_north.add_field(
        "node",
        "topographic__elevation",
        grid_north.node_y * -S + 1.2,
    )
    grid_north.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=False,
    )
    return grid_north


@pytest.fixture
def grid_west():
    grid_west = RasterModelGrid((5, 5))
    grid_west.add_field(
        "node",
        "topographic__elevation",
        grid_west.node_x * S,
    )
    grid_west.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=False,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_west


@pytest.fixture
def grid_south():
    grid_south = RasterModelGrid((5, 5))
    grid_south.add_field(
        "node",
        "topographic__elevation",
        grid_south.node_y * S,
    )
    grid_south.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=False,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_south


# SLOPES TO DIAGONALS
@pytest.fixture
def grid_ne():
    grid_ne = RasterModelGrid((5, 5))
    grid_ne.add_field(
        "node",
        "topographic__elevation",
        (grid_ne.node_y + grid_ne.node_x) * -S,
    )
    grid_ne.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_ne


@pytest.fixture
def grid_nw():
    grid_nw = RasterModelGrid((5, 5))
    grid_nw.add_field(
        "node",
        "topographic__elevation",
        (grid_nw.node_x - grid_nw.node_y) * S,
    )
    grid_nw.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_nw


@pytest.fixture
def grid_sw():
    grid_sw = RasterModelGrid((5, 5))
    grid_sw.add_field(
        "node",
        "topographic__elevation",
        (grid_sw.node_x + grid_sw.node_y) * S,
    )
    grid_sw.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_sw


@pytest.fixture
def grid_se():
    grid_se = RasterModelGrid((5, 5))
    grid_se.add_field(
        "node",
        "topographic__elevation",
        (grid_se.node_y - grid_se.node_x) * S,
    )
    grid_se.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_se


# INTERMEDIATE SLOPES:
@pytest.fixture
def grid_ene():
    grid_ene = RasterModelGrid((5, 5))
    grid_ene.add_field(
        "node",
        "topographic__elevation",
        (grid_ene.node_y + 2 * grid_ene.node_x) * -S,
    )
    grid_ene.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_ene


@pytest.fixture
def grid_nne():
    grid_nne = RasterModelGrid((5, 5))
    grid_nne.add_field(
        "node",
        "topographic__elevation",
        (2 * grid_nne.node_y + grid_nne.node_x) * -S,
    )
    grid_nne.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_nne


@pytest.fixture
def grid_nnw():
    grid_nnw = RasterModelGrid((5, 5))
    grid_nnw.add_field(
        "node",
        "topographic__elevation",
        (grid_nnw.node_x - 2 * grid_nnw.node_y) * S,
    )
    grid_nnw.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_nnw


@pytest.fixture
def grid_wnw():
    grid_wnw = RasterModelGrid((5, 5))
    grid_wnw.add_field(
        "node",
        "topographic__elevation",
        (2 * grid_wnw.node_x - grid_wnw.node_y) * S,
    )
    grid_wnw.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_wnw


@pytest.fixture
def grid_wsw():
    grid_wsw = RasterModelGrid((5, 5))
    grid_wsw.add_field(
        "node",
        "topographic__elevation",
        (2 * grid_wsw.node_x + grid_wsw.node_y) * S,
    )
    grid_wsw.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_wsw


@pytest.fixture
def grid_ssw():
    grid_ssw = RasterModelGrid((5, 5))
    grid_ssw.add_field(
        "node",
        "topographic__elevation",
        (grid_ssw.node_x + 2 * grid_ssw.node_y) * S,
    )
    grid_ssw.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_ssw


@pytest.fixture
def grid_sse():
    grid_sse = RasterModelGrid((5, 5))
    grid_sse.add_field(
        "node",
        "topographic__elevation",
        (2 * grid_sse.node_y - grid_sse.node_x) * S,
    )
    grid_sse.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_sse


@pytest.fixture
def grid_ese():
    grid_ese = RasterModelGrid((5, 5))
    grid_ese.add_field(
        "node",
        "topographic__elevation",
        (grid_ese.node_y - 2 * grid_ese.node_x) * S,
    )
    grid_ese.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return grid_ese


@pytest.fixture
def cc_flat(grid_flat):
    return ClastCollection(
        grid_flat,
        clast_x=np.ones(2) * [2, 3],
        clast_y=np.ones(2) * 2,
        clast_elev=np.ones(2) * 6,
        clast_radius=np.ones(2) * 0.5,
    )


@pytest.fixture
def cc_east(grid_east):
    return ClastCollection(
        grid_east,
        clast_x=np.ones(2) * [2, 3],
        clast_y=np.ones(2) * 2,
        clast_elev=np.ones(2) * 6,
        clast_radius=np.ones(2) * 0.5,
    )


@pytest.fixture
def cc_north(grid_north):
    return ClastCollection(
        grid_north,
        clast_x=np.ones(2) * [2, 3],
        clast_y=np.ones(2) * 2,
        clast_elev=np.ones(2) * 6,
        clast_radius=np.ones(2) * 0.5,
    )


@pytest.fixture
def cc_west(grid_west):
    return ClastCollection(
        grid_west,
        clast_x=np.ones(2) * [2, 3],
        clast_y=np.ones(2) * 2,
        clast_elev=np.ones(2) * 6,
        clast_radius=np.ones(2) * 0.5,
    )


@pytest.fixture
def cc_south(grid_south):
    return ClastCollection(
        grid_south,
        clast_x=np.ones(2) * [2, 3],
        clast_y=np.ones(2) * 2,
        clast_elev=np.ones(2) * 6,
        clast_radius=np.ones(2) * 0.5,
    )


@pytest.fixture
def cc_ne(grid_ne):
    return ClastCollection(
        grid_ne,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_nw(grid_nw):
    return ClastCollection(
        grid_nw,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_sw(grid_sw):
    return ClastCollection(
        grid_sw,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_se(grid_se):
    return ClastCollection(
        grid_se,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_ene(grid_ene):
    return ClastCollection(
        grid_ene,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_nne(grid_nne):
    return ClastCollection(
        grid_nne,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_nnw(grid_nnw):
    return ClastCollection(
        grid_nnw,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_wnw(grid_wnw):
    return ClastCollection(
        grid_wnw,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_wsw(grid_wsw):
    return ClastCollection(
        grid_wsw,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_ssw(grid_ssw):
    return ClastCollection(
        grid_ssw,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_sse(grid_sse):
    return ClastCollection(
        grid_sse,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )


@pytest.fixture
def cc_ese(grid_ese):
    return ClastCollection(
        grid_ese,
        clast_x=np.ones(3) * [1, 2, 3],
        clast_y=np.ones(3) * 2,
        clast_elev=np.ones(3) * 6,
        clast_radius=np.ones(3) * 0.5,
    )
