import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.graph.sort.ext._deprecated_sparse import offset_to_sorted_block
from landlab.graph.sort.ext.intpair import fill_offsets_to_sorted_blocks
from landlab.graph.sort.intpair import IntPairCollection
from landlab.graph.sort.intpair import IntPairMapping
from landlab.graph.sort.intpair import map_pairs_to_values
from landlab.graph.sort.intpair import map_rolling_pairs_to_values
from landlab.graph.sort.intpair import pair_isin


def test_fill_offsets():

    array = [0, 2, 6, 7]
    offsets = np.full(10, -2)

    fill_offsets_to_sorted_blocks(np.asarray(array), offsets)

    assert_array_equal(offsets, [0, 1, 1, 2, 2, 2, 2, 3, 4, 4])

    array = [0, 2, 2, 2, 6, 7, 7]
    offsets = np.full(10, -2)
    fill_offsets_to_sorted_blocks(np.asarray(array), offsets)

    assert_array_equal(offsets, [0, 1, 1, 4, 4, 4, 4, 5, 7, 7])

    actual = []
    for i in range(len(offsets) - 1):
        actual += [i] * (offsets[i + 1] - offsets[i])

    assert actual == array

    array = [0, 2, 6, 7]
    offsets = np.full(5, -2)
    fill_offsets_to_sorted_blocks(np.asarray(array), offsets)

    assert_array_equal(offsets, [0, 1, 1, 2, 2])

    actual = []
    for i in range(len(offsets) - 1):
        actual += [i] * (offsets[i + 1] - offsets[i])

    assert actual == [0, 2]


def test_collection_contains():
    pairs = IntPairCollection(
        [[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]]
    )
    assert (4, 8) not in pairs
    assert (3, 0) in pairs
    assert (0, 3) in pairs

    assert pairs.contains_pairs([(4, 8)]) == (False,)
    assert pairs.contains_pairs((4, 8)) == (False,)
    assert_array_equal(
        pairs.contains_pairs([(4, 8), (7, 8), (0, 1)]), (False, True, True)
    )

    pairs = IntPairCollection([[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]])
    assert_array_equal(
        pairs.contains_pairs([[1, 1], [3, 1], [1, 2], [5, 1]]),
        (True, True, True, False),
    )


def test_collection_with_sorter():
    src = np.asarray([[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]])
    pairs = IntPairCollection(src, sorter=np.argsort(src[:, 0]))

    assert_array_equal(
        pairs.contains_pairs(
            [
                *([5, 2], [5, 4], [5, 8], [1, 2], [1, 4]),
                *([1, 0], [1, 3], [2, 4], [0, 6], [0, 3]),
                *([7, 3], [7, 4], [7, 6], [7, 8], [3, 6]),
                *([3, 4], [4, 8]),
            ]
        ),
        [
            *(True, False, True, True, False),
            *(True, False, False, False, True),
            *(False, False, True, True, True),
            *(False, False),
        ],
    )


def test_mapping():
    pairs = IntPairMapping(
        [[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], values=[0, 10, 20, 30, 40]
    )

    assert_array_equal(pairs.get_items([(1, 1), (3, 1)]), (10, 30))


def test_mapping_with_sorter_keyword():
    pairs = IntPairMapping(
        [[0, 1], [1, 1], [3, 1], [2, 1], [4, 1]],
        values=[0, 10, 30, 20, 40],
        sorter=[0, 1, 3, 2, 4],
    )

    assert_array_equal(pairs.get_items([(1, 1), (3, 1)]), (10, 30))


def test_mapping_pairs_at_end():
    pairs = IntPairMapping(
        [[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], values=[0, 10, 20, 30, 40]
    )
    assert_array_equal(
        pairs.get_items([[1, 1], [3, 1], [4, 1], [1, 4], [4, 2]]), [10, 30, 40, 40, -1]
    )


def test_mapping_pairs_transposed():
    pairs = IntPairMapping(
        [[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], values=[0, 10, 20, 30, 40]
    )

    assert_array_equal(pairs.get_items([(1, 0), (1, 3)]), [0, 30])


def test_mapping_pairs_all_missing():
    pairs = IntPairMapping(
        [[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]], values=[0, 10, 20, 30, 40]
    )
    assert_array_equal(pairs.get_items([(5, 1), (1, 42)]), (-1, -1))


def test_mapping_pairs_big_data():
    grid = RasterModelGrid((2**10, 2**10))

    n_pairs = len(grid.nodes_at_link)

    pairs = IntPairMapping(grid.nodes_at_link, values=np.arange(n_pairs) * 10)

    ids = np.random.randint(n_pairs, size=n_pairs // 10)

    assert_array_equal(pairs.get_items(grid.nodes_at_link[ids]), ids * 10)


def test_mapping_rolling_pairs():
    pairs = IntPairMapping(
        [[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]], values=[0, 10, 20, 30, 40]
    )

    assert_array_equal(
        pairs.get_items([[0, 1, 2, 3], [0, 2, 3, 4]], wraparound=True),
        [[0, 10, 20, -1], [-1, 20, 30, 40]],
    )

    assert_array_equal(
        pairs.get_items([[0, 1, 2, 3], [0, 2, 3, 4]], wraparound=False),
        [[0, 10, 20], [-1, 20, 30]],
    )


def test_mapping_rolling_pairs_with_jagged_rows():
    pairs = IntPairMapping(
        [[0, 1], [1, 2], [2, 3], [2, 0], [3, 4], [4, 0]], values=[0, 10, 20, 30, 40, 50]
    )

    assert_array_equal(
        pairs.get_items([[0, 1, 2, -1], [0, 2, -1, -1]], wraparound=True),
        [[0, 10, 30, -1], [30, 30, -1, -1]],
    )
    assert_array_equal(
        pairs.get_items([[0, 1, 2, -1], [0, 2, -1, -1]], wraparound=False),
        [[0, 10, -1], [30, -1, -1]],
    )


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
