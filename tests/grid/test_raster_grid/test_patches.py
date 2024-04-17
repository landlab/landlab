import numpy as np
from pytest import approx
from pytest import raises

from landlab import RasterModelGrid


class TestPatchesAtNode:
    patch_values = np.array(
        [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [2, 1, -1, -1],
            [3, 2, -1, -1],
            [-1, 3, -1, -1],
            [4, -1, -1, 0],
            [5, 4, 0, 1],
            [6, 5, 1, 2],
            [7, 6, 2, 3],
            [-1, 7, 3, -1],
            [8, -1, -1, 4],
            [9, 8, 4, 5],
            [10, 9, 5, 6],
            [11, 10, 6, 7],
            [-1, 11, 7, -1],
            [-1, -1, -1, 8],
            [-1, -1, 8, 9],
            [-1, -1, 9, 10],
            [-1, -1, 10, 11],
            [-1, -1, 11, -1],
        ]
    )

    patch_mask = np.array(
        [
            [False, True, True, True],
            [False, False, True, True],
            [False, False, True, True],
            [False, False, True, True],
            [True, False, True, True],
            [False, True, True, False],
            [False, False, False, False],
            [False, False, False, False],
            [False, False, False, False],
            [True, False, False, True],
            [False, True, True, False],
            [False, False, False, False],
            [False, False, False, False],
            [False, False, False, False],
            [True, False, False, True],
            [True, True, True, False],
            [True, True, False, False],
            [True, True, False, False],
            [True, True, False, False],
            [True, True, False, True],
        ]
    )

    def test_patches_at_node(self):
        rmg = RasterModelGrid((4, 5))
        patches_out = rmg.patches_at_node
        assert np.all(patches_out == self.patch_values)

    def test_nodes_at_patch(self):
        rmg = RasterModelGrid((4, 5))
        nodes_out = rmg.nodes_at_patch
        assert np.all(
            nodes_out
            == np.array(
                [
                    [6, 5, 0, 1],
                    [7, 6, 1, 2],
                    [8, 7, 2, 3],
                    [9, 8, 3, 4],
                    [11, 10, 5, 6],
                    [12, 11, 6, 7],
                    [13, 12, 7, 8],
                    [14, 13, 8, 9],
                    [16, 15, 10, 11],
                    [17, 16, 11, 12],
                    [18, 17, 12, 13],
                    [19, 18, 13, 14],
                ]
            )
        )


class TestSlopesAtPatches:
    def test_slopes_at_patches(self):
        rmg = RasterModelGrid((4, 5))
        rmg.at_node["topographic__elevation"] = rmg.node_x.copy()
        slopes_out = rmg.calc_slope_at_node()
        assert np.all(slopes_out == approx(np.full(20, np.pi / 4.0, dtype=float)))

    def test_slopes_at_patches_comps(self):
        rmg = RasterModelGrid((4, 5))
        rmg.at_node["topographic__elevation"] = -rmg.node_y
        slopes_out = rmg.calc_slope_at_node(
            rmg.at_node["topographic__elevation"], return_components=True
        )
        assert np.all(slopes_out[0] == approx(np.full(20, np.pi / 4.0, dtype=float)))
        assert np.all(
            slopes_out[1][1] == approx(np.full(20, -np.pi / 4.0, dtype=float))
        )
        assert np.all(slopes_out[1][0] == approx(np.zeros(20, dtype=float)))


def test_number_of_patches():
    grid = RasterModelGrid((4, 5))
    assert grid.number_of_patches == 12


def test_patches_at_node():
    grid = RasterModelGrid((3, 3))
    assert np.all(
        grid.patches_at_node
        == np.array(
            [
                [0, -1, -1, -1],
                [1, 0, -1, -1],
                [-1, 1, -1, -1],
                [2, -1, -1, 0],
                [3, 2, 0, 1],
                [-1, 3, 1, -1],
                [-1, -1, -1, 2],
                [-1, -1, 2, 3],
                [-1, -1, 3, -1],
            ],
            dtype=int,
        )
    )


def test_links_at_patch():
    grid = RasterModelGrid((3, 4))
    assert np.all(
        grid.links_at_patch
        == np.array(
            [
                [4, 7, 3, 0],
                [5, 8, 4, 1],
                [6, 9, 5, 2],
                [11, 14, 10, 7],
                [12, 15, 11, 8],
                [13, 16, 12, 9],
            ],
            dtype=int,
        )
    )

    with raises(ValueError):
        grid.links_at_patch[0] = -1


def test_patches_at_link():
    grid = RasterModelGrid((3, 4))
    assert np.all(
        grid.patches_at_link
        == np.array(
            [
                [-1, 0],
                [-1, 1],
                [-1, 2],
                [0, -1],
                [1, 0],
                [2, 1],
                [-1, 2],
                [0, 3],
                [1, 4],
                [2, 5],
                [3, -1],
                [4, 3],
                [5, 4],
                [-1, 5],
                [3, -1],
                [4, -1],
                [5, -1],
            ],
            dtype=int,
        )
    )

    with raises(ValueError):
        grid.patches_at_link[0] = -1
