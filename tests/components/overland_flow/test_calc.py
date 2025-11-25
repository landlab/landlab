import numpy as np
import pytest
from numpy.testing import assert_allclose

from landlab import RasterModelGrid
from landlab.components.overland_flow._calc import calc_grad_at_link


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
