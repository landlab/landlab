import numpy as np
from numpy.testing import assert_array_equal

from landlab import FramedVoronoiGrid


def test_rect_number_of_patches():
    grid = FramedVoronoiGrid((4, 3))
    assert grid.number_of_patches == 12

    grid = FramedVoronoiGrid((3, 4))
    assert grid.number_of_patches == 12


def test_rect_nodes_at_patch():
    grid = FramedVoronoiGrid((3, 3))
    assert_array_equal(
        grid.nodes_at_patch,
        [
            [3, 0, 1],
            [4, 1, 2],
            [4, 3, 1],
            [5, 4, 2],
            [6, 3, 4],
            [7, 4, 5],
            [7, 6, 4],
            [8, 7, 5],
        ],
    )


def test_rect_links_at_patch():
    grid = FramedVoronoiGrid((3, 3))
    assert_array_equal(
        grid.links_at_patch,
        np.array(
            [
                [3, 2, 0],
                [5, 4, 1],
                [7, 3, 4],
                [8, 5, 6],
                [10, 9, 7],
                [12, 11, 8],
                [14, 10, 11],
                [15, 12, 13],
            ]
        ),
    )
