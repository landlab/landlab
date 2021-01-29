import numpy as np
import pytest

from landlab import HexModelGrid
from landlab.grid.mappers import map_link_vector_components_to_node


@pytest.mark.parametrize("angle", (0.0, 60.0, 120.0, 180.0, 240.0, 300.0,))
@pytest.mark.parametrize("node_layout", ("rect", "hex"))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
def test_along_links(orientation, node_layout, angle):

    if orientation == "horizontal":
        east = np.deg2rad(0.0)
        nne = np.deg2rad(60.0)
        nnw = np.deg2rad(120.0)
        angle = np.deg2rad(angle)
    else:
        east = np.deg2rad(330.0)
        nne = np.deg2rad(30.0)
        nnw = np.deg2rad(90.0)
        angle = np.deg2rad(angle - 30.0)

    grid = HexModelGrid((3, 3), node_layout=node_layout, orientation=orientation)
    value_at_link = np.empty(grid.number_of_links)
    is_core_node = grid.status_at_node == grid.BC_NODE_IS_CORE

    east_links = np.isclose(grid.angle_of_link, east)
    nne_links = np.isclose(grid.angle_of_link, nne)
    nnw_links = np.isclose(grid.angle_of_link, nnw)

    value_at_link[east_links] = np.cos(angle - east)
    value_at_link[nne_links] = np.cos(angle - nne)
    value_at_link[nnw_links] = np.cos(angle - nnw)

    vx, vy = map_link_vector_components_to_node(grid, value_at_link)

    assert np.sqrt(vx[is_core_node] ** 2 + vy[is_core_node] ** 2) == pytest.approx(1.0)

    assert vx[~is_core_node] == pytest.approx(0.0)
    assert vy[~is_core_node] == pytest.approx(0.0)

    assert vx[is_core_node] == pytest.approx(np.cos(angle))
    assert vy[is_core_node] == pytest.approx(np.sin(angle))
