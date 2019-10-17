from numpy.testing import assert_array_equal

from landlab import HexModelGrid


def test_number_of_patches():
    grid = HexModelGrid((4, 3))
    assert grid.number_of_patches == 19

    grid = HexModelGrid((3, 4))
    assert grid.number_of_patches == 14

    grid = HexModelGrid((4, 3), node_layout="rect")
    assert grid.number_of_patches == 12

    grid = HexModelGrid((3, 4), node_layout="rect")
    assert grid.number_of_patches == 12


def test_nodes_at_path():
    grid = HexModelGrid((3, 2))
    assert_array_equal(
        grid.nodes_at_patch,
        [[3, 0, 1], [3, 2, 0], [4, 3, 1], [5, 2, 3], [6, 3, 4], [6, 5, 3]],
    )


def test_links_at_patch():
    grid = HexModelGrid((3, 2))
    assert (
        grid.links_at_patch
        == [[3, 2, 0], [5, 1, 2], [6, 3, 4], [8, 7, 5], [10, 9, 6], [11, 8, 9]]
    ).all()
    # assert_array_equal(grid.links_at_patch, [[ 3,  2,  0],
    #                                          [ 5,  1,  2],
    #                                          [ 6,  3,  4],
    #                                          [ 8,  7,  5],
    #                                          [10,  9,  6],
    #                                          [11,  8,  9]])
