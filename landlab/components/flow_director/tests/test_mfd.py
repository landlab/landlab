
import pytest
import numpy as np

from landlab import RasterModelGrid
from landlab.components.flow_director import flow_direction_mfd


def test_bad_argument_mfd():
    mg = RasterModelGrid((5, 5), spacing=(1, 1))
    z = mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")

    neighbors_at_node = mg.adjacent_nodes_at_node
    links_at_node = mg.links_at_node
    active_link_dir_at_node = mg.active_link_dirs_at_node
    link_slope = np.arctan(mg.calc_grad_at_link(z))
    slopes_to_neighbors_at_node = link_slope[links_at_node] * active_link_dir_at_node

    with pytest.raises(ValueError):
        flow_direction_mfd.flow_directions_mfd(
            z,
            neighbors_at_node,
            links_at_node,
            active_link_dir_at_node,
            link_slope,
            partition_method="foo",
        )
