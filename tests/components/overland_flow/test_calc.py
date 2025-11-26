import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components.overland_flow._calc import calc_grad_at_link
from landlab.components.overland_flow._calc import zero_out_dry_links


@pytest.mark.parametrize("dx", (1.0, 2.0, -1.0))
@pytest.mark.parametrize("dy", (1.0, 2.0, -1.0))
def test_calc_grad_non_uniform_spacing(dx, dy):
    grid = RasterModelGrid((3, 4), xy_spacing=(dx, dy))

    length_of_link = np.empty(grid.number_of_links, dtype=float)
    length_of_link[grid.horizontal_links] = dx
    length_of_link[grid.vertical_links] = dy

    z = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [8, 9, 10, 11],
    ]
    expected = np.asarray(
        [
            *[1.0, 1.0, 1.0],
            *[4.0, 4.0, 4.0, 4.0],
            *[1.0, 1.0, 1.0],
            *[4.0, 4.0, 4.0, 4.0],
            *[1.0, 1.0, 1.0],
        ]
    )
    expected /= length_of_link

    where = np.arange(grid.number_of_links, dtype=int)
    actual = np.empty(grid.number_of_links, dtype=float)

    ret = calc_grad_at_link(
        np.asarray(z, dtype=float).ravel(),
        length_of_link=length_of_link,
        nodes_at_link=grid.nodes_at_link.astype(int),
        where=where,
        out=actual,
    )

    assert_allclose(actual, expected)
    assert ret is actual


@pytest.mark.parametrize(
    "here",
    (
        [
            *[1, 1, 1],
            *[1, 1, 1, 1],
            *[1, 1, 1],
            *[1, 1, 1, 1],
            *[1, 1, 1],
        ],
        [
            *[0, 0, 0],
            *[0, 0, 0, 0],
            *[0, 0, 0],
            *[0, 0, 0, 0],
            *[0, 0, 0],
        ],
        [
            *[1, 1, 1],
            *[1, 0, 0, 0],
            *[0, 1, 1],
            *[0, 0, 0, 0],
            *[0, 1, 0],
        ],
    ),
)
@pytest.mark.parametrize("float_dtype", [np.float32, np.float64])
@pytest.mark.parametrize("id_dtype", [np.int32, np.int64])
def test_calc_grad_with_mask(here, float_dtype, id_dtype):
    grid = RasterModelGrid((3, 4), xy_spacing=(1.0, 1.0))
    z = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [8, 9, 10, 11],
    ]
    expected = np.asarray(
        [
            *[1.0, 1.0, 1.0],
            *[4.0, 4.0, 4.0, 4.0],
            *[1.0, 1.0, 1.0],
            *[4.0, 4.0, 4.0, 4.0],
            *[1.0, 1.0, 1.0],
        ]
    )

    here = np.asarray(here, dtype=bool)
    where = np.flatnonzero(here)
    actual = np.full(grid.number_of_links, -999, dtype=float_dtype)

    ret = calc_grad_at_link(
        np.asarray(z, dtype=float_dtype).ravel(),
        length_of_link=grid.length_of_link.astype(float_dtype),
        nodes_at_link=grid.nodes_at_link.astype(id_dtype),
        where=where.astype(id_dtype),
        out=actual,
    )

    assert np.all(here) or np.allclose(actual[~here], -999)
    assert np.all(~here) or np.allclose(actual[here], expected[here])

    assert ret is actual


@pytest.mark.parametrize(
    "h, where, initial_out, expected",
    [
        (
            [0.1, 0.0, -0.3, 0.2],
            [0, 1, 2, 3],
            [5.0, 4.0, 3.0, 2.0],
            [5.0, 0.0, 0.0, 2.0],
        ),
        (
            [0.1, 0.0, -0.3, 0.2],
            [1, 2],
            [5.0, 4.0, 3.0, 2.0],
            [5.0, 0.0, 0.0, 2.0],
        ),
        (
            [0.1, 0.5, 0.2, 0.3],
            [0, 2],
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
        ),
        (
            [0.0, -1.0, 0.0],
            [0, 1, 2],
            [1.0, 2.0, 3.0],
            [0.0, 0.0, 0.0],
        ),
    ],
)
def test_zero_out_dry_links(h, where, initial_out, expected):
    actual = np.array(initial_out)

    rtn = zero_out_dry_links(
        np.asarray(h, dtype=float), where=np.asarray(where, dtype=int), out=actual
    )

    assert rtn is actual
    assert_array_equal(actual, expected)


def test_zero_out_dry_links_does_not_touch_indices_not_in_where():
    h = [0.0, -1.0, 0.5, -0.1]
    where = [0, 3]

    actual = np.full_like(h, -999)
    expected = [0.0, -999, -999, 0.0]

    rtn = zero_out_dry_links(
        np.asarray(h, dtype=float),
        where=np.asarray(where, dtype=int),
        out=actual,
    )

    assert rtn is actual
    assert_array_equal(actual, expected)


def test_zero_out_dry_links_nowhere():
    h = [0.0, -1.0, 0.5]
    where = []

    actual = np.full_like(h, -999)
    expected = actual.copy()

    rtn = zero_out_dry_links(
        np.asarray(h, dtype=float),
        where=np.asarray(where, dtype=int),
        out=actual,
    )

    assert rtn is actual
    assert_array_equal(actual, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_zero_out_dry_links_dtype(dtype):
    h = [0.0, 1.0]
    where = [0, 1]

    actual = np.asarray([5.0, 6.0], dtype=dtype)
    expected = [0.0, 6.0]
    rtn = zero_out_dry_links(
        np.asarray(h, dtype=dtype),
        where=np.asarray(where, dtype=int),
        out=actual,
    )

    assert rtn.dtype == dtype
    assert_array_equal(actual, expected)
