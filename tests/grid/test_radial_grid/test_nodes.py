from numpy.testing import assert_array_equal
from pytest import approx

from landlab import RadialModelGrid


def test_radius_at_node():
    grid = RadialModelGrid(2, 8)
    assert grid.radius_at_node == approx(
        [
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            1.0,
            2.0,
            2.0,
            1.0,
            1.0,
            2.0,
            1.0,
            0.0,
            1.0,
            2.0,
            1.0,
            1.0,
            2.0,
            2.0,
            1.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
        ]
    )


def test_spacing_of_rings():
    grid = RadialModelGrid(2, spacing=1.0)
    assert grid.spacing_of_rings == approx(1.0)
    assert isinstance(grid.spacing_of_rings, float)

    grid = RadialModelGrid(2, spacing=11.0)
    assert grid.spacing_of_rings == approx(11.0)
    assert isinstance(grid.spacing_of_rings, float)


def test_number_of_rings():
    grid = RadialModelGrid(2)
    assert grid.number_of_rings == 2
    assert isinstance(grid.number_of_rings, int)

    grid = RadialModelGrid(4)
    assert grid.number_of_rings == 4


def test_nodes_per_ring():
    grid = RadialModelGrid(3, 8)
    assert_array_equal(grid.number_of_nodes_in_ring, [8, 16, 32])
