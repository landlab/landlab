import numpy as np
from numpy.testing import assert_array_equal

from landlab.graph.sort.ext.remap_element import offset_to_sorted_block
from landlab.graph.sort.intpair import map_pairs_to_values
from landlab.graph.sort.intpair import map_rolling_pairs_to_values
from landlab.graph.sort.intpair import pair_isin


def test_pair_isin_one_pair():
    src = np.asarray(
        [[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]], dtype=int
    )
    pairs = np.asarray([[4, 8]], dtype=int)
    out = pair_isin(src, pairs)

    assert tuple(out) == (False,)


def test_pair_isin_1():
    src = np.asarray(
        [[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]], dtype=int
    )
    pairs = np.asarray(
        [
            [5, 2],
            [5, 4],
            [5, 8],
            [1, 2],
            [1, 4],
            [1, 0],
            [1, 3],
            [2, 4],
            [0, 6],
            [0, 3],
            [7, 3],
            [7, 4],
            [7, 6],
            [7, 8],
            [3, 6],
            [3, 4],
            [4, 8],
        ],
        dtype=int,
    )
    out = pair_isin(src, pairs, sorter=np.argsort(src[:, 0]))

    assert tuple(out) == (
        True,
        False,
        True,
        True,
        False,
        True,
        False,
        False,
        False,
        True,
        False,
        False,
        True,
        True,
        True,
        False,
        False,
    )


def test_pair_isin():
    src = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)

    pairs = np.asarray([[1, 1], [3, 1], [1, 2], [5, 1]], dtype=int)
    out = pair_isin(src, pairs)

    assert np.all(out == [True, True, True, False])


def test_pair_isin_with_out_keyword():
    src = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)

    pairs = np.asarray([[1, 1], [3, 1], [1, 2], [5, 1]], dtype=int)
    out = np.empty(len(pairs), dtype=bool)
    rtn = pair_isin(src, pairs, out=out)

    assert rtn is out
    assert np.all(out == [True, True, True, False])


def test_map_pairs():
    src = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    data = np.arange(len(src), dtype=int)

    pairs = np.asarray([[1, 1], [3, 1]], dtype=int)
    out = map_pairs_to_values((src, data), pairs)

    assert np.all(out == [1, 3])


def test_map_pairs_with_out_keyword():
    src = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    data = np.arange(len(src), dtype=int)

    pairs = np.asarray([[1, 1], [3, 1]], dtype=int)
    out = np.empty(len(pairs), dtype=int)
    rtn = map_pairs_to_values((src, data), pairs, out=out)

    assert rtn is out
    assert np.all(out == [1, 3])


def test_map_pairs_with_sorter_keyword():
    src = np.asarray([[0, 1], [1, 1], [3, 1], [2, 1], [4, 1]], dtype=int)
    data = np.array([0, 1, 3, 2, 4], dtype=int)

    pairs = np.asarray([[1, 1], [3, 1]], dtype=int)
    out = map_pairs_to_values((src, data), pairs, sorter=np.argsort(src[:, 0]))

    assert np.all(out == [1, 3])


def test_map_pairs_at_end():
    src = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    data = np.arange(len(src), dtype=int)

    pairs = np.asarray([[1, 1], [3, 1], [4, 1], [1, 4], [4, 2]], dtype=int)
    out = map_pairs_to_values((src, data), pairs)

    assert np.all(out == [1, 3, 4, 4, -1])


def test_map_pairs_transposed():
    src = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    data = np.arange(len(src), dtype=int)

    pairs = np.asarray([[1, 0], [1, 3]], dtype=int)
    out = map_pairs_to_values((src, data), pairs)

    assert np.all(out == [0, 3])


def test_map_pairs_all_missing():
    src = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    data = np.arange(len(src), dtype=int)

    pairs = np.asarray([[5, 1], [1, 42]], dtype=int)
    out = map_pairs_to_values((src, data), pairs)

    assert np.all(out == [-1, -1])


def test_map_pairs_big_data():
    n_src_pairs = 2**20
    src = np.random.randint(n_src_pairs * 2, size=n_src_pairs * 2, dtype=int).reshape(
        (-1, 2)
    )
    src = src[np.lexsort((src[:, 1], src[:, 0]))]
    data = np.arange(n_src_pairs, dtype=int)

    ids = np.random.randint(n_src_pairs, size=n_src_pairs // 10)
    pairs = src[ids]
    out = map_pairs_to_values((src, data), pairs)

    assert_array_equal(out, ids)


def test_map_rolling_pairs():
    src = np.asarray([[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]], dtype=int)
    data = np.asarray([0, 1, 2, 3, 4], dtype=int)

    pairs = np.asarray([[0, 1, 2, 3], [0, 2, 3, 4]], dtype=int)
    out = map_rolling_pairs_to_values((src, data), pairs)

    assert np.all(out == [[0, 1, 2, -1], [-1, 2, 3, 4]])


def test_map_rolling_pairs_with_out_keyword():
    src = np.asarray([[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]], dtype=int)
    data = np.asarray([0, 1, 2, 3, 4], dtype=int)

    pairs = np.asarray([[0, 1, 2, 3], [0, 2, 3, 4]], dtype=int)
    out = np.empty_like(pairs)
    rtn = map_rolling_pairs_to_values((src, data), pairs, out=out)

    assert rtn is out
    assert np.all(out == [[0, 1, 2, -1], [-1, 2, 3, 4]])


def test_map_rolling_pairs_with_size_of_row():
    src = np.asarray([[0, 1], [1, 2], [2, 3], [2, 0], [3, 4], [4, 0]], dtype=int)
    data = np.asarray([0, 1, 2, 3, 4, 5], dtype=int)

    pairs = np.asarray([[0, 1, 2, -1], [0, 2, -1, -1]], dtype=int)
    out = map_rolling_pairs_to_values((src, data), pairs, size_of_row=[3, 2])

    assert np.all(out == [[0, 1, 3, -1], [3, 3, -1, -1]])


def test_offset_2d():
    pairs = np.asarray([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    out = np.empty(len(pairs), dtype=int)
    offset_to_sorted_block(pairs, out)

    assert np.all(out == [0, 1, 2, 3, 4])


def test_offset_2d_with_missing_at_start():
    pairs = np.asarray([[1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    out = np.empty(5, dtype=int)
    offset_to_sorted_block(pairs, out)

    assert np.all(out == [0, 0, 1, 2, 3])


def test_offset_2d_with_missing_in_middle():
    pairs = np.asarray([[0, 1], [1, 1], [3, 1], [4, 1]], dtype=int)
    out = np.empty(5, dtype=int)
    offset_to_sorted_block(pairs, out)

    assert np.all(out == [0, 1, 2, 2, 3])


def test_offset_with_different_strides():
    pairs = np.arange(5, dtype=int).reshape((-1, 1))
    out = np.empty(len(pairs), dtype=int)
    for _ in range(5):
        offset_to_sorted_block(pairs, out)
        assert np.all(out == [0, 1, 2, 3, 4])
        pairs = np.hstack((pairs, pairs))


def test_offset_with_sort_out():
    pairs = np.asarray([[1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    out = np.empty(4, dtype=int)
    offset_to_sorted_block(pairs, out)

    assert np.all(out == [0, 0, 1, 2])


def test_offset_with_long_out():
    pairs = np.asarray([[1, 1], [2, 1], [3, 1], [4, 1]], dtype=int)
    out = np.empty(7, dtype=int)
    offset_to_sorted_block(pairs, out)

    assert np.all(out == [0, 0, 1, 2, 3, 4, 4])


def test_offset_to_sorted_block_with_negatives():
    pairs = np.array(
        [
            [-1, 0],
            [-1, 8],
            [-1, 7],
            [0, 3],
            [0, 2],
            [1, 4],
            [1, 7],
            [2, 8],
            [3, 5],
            [4, 5],
            [5, 6],
            [6, 7],
            [6, 8],
        ]
    )
    out = np.empty(len(pairs), dtype=int)
    offset_to_sorted_block(pairs, out)
    assert tuple(out) == (3, 5, 7, 8, 9, 10, 11, 13, 13, 13, 13, 13, 13)


def test_offset_to_sorted_block_all_negatives():
    pairs = np.array([[-1, 0], [-1, 8], [-1, 7]])
    out = np.empty(len(pairs) + 1, dtype=int)
    offset_to_sorted_block(pairs, out)
    assert tuple(out) == (3, 3, 3, 3)
