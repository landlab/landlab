import pytest

from landlab import RadialModelGrid


@pytest.mark.parametrize("n_rings", (2, 3, 4))
@pytest.mark.parametrize("nodes_in_first_ring", (8, 16, 32))
def test_as_dataset(n_rings, nodes_in_first_ring):
    xy_of_center = (1.0, -1.0)
    grid = RadialModelGrid(n_rings, nodes_in_first_ring, xy_of_center=xy_of_center)
    ds = grid.as_dataset()
    assert ds.n_rings == n_rings
    assert ds.nodes_in_first_ring == nodes_in_first_ring
    assert tuple(ds.xy_of_center) == xy_of_center
