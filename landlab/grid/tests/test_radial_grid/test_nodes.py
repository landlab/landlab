from numpy.testing import assert_array_equal

from landlab import RadialModelGrid


def test_radius_at_node():
    grid = RadialModelGrid(2)
    assert_array_equal(grid.radius_at_node,
                       [ 2.,  2.,  2.,  2.,  1.,  2.,  1.,  1.,  2.,  2.,
                         0.,  2.,  1.,  1.,  2.,  1.,  2.,  2.,  2.,  2.])


def test_spacing_of_shells():
    grid = RadialModelGrid(2)
    assert grid.spacing_of_shells == 1.
    assert isinstance(grid.spacing_of_shells, float)

    grid = RadialModelGrid(2, dr=11)
    assert grid.spacing_of_shells == 11.
    assert isinstance(grid.spacing_of_shells, float)


def test_number_of_shells():
    grid = RadialModelGrid(2)
    assert grid.number_of_shells == 3
    assert isinstance(grid.number_of_shells, int)

    grid = RadialModelGrid(4)
    assert grid.number_of_shells == 5


def test_nodes_per_shell():
    grid = RadialModelGrid(2)
    assert_array_equal(grid.nodes_per_shell, [1, 6, 13])
