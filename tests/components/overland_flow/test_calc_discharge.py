import numpy as np
from hypothesis import assume
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal

from landlab.components.overland_flow._calc import calc_discharge_at_links

valid_q = st.floats(min_value=-50.0, max_value=50.0)
valid_mannings = st.floats(min_value=0.0, max_value=5.0)
valid_h = st.floats(min_value=1e-6, max_value=20.0, exclude_min=True)
valid_slope = st.floats(min_value=-2.0, max_value=2.0)
valid_g = st.floats(min_value=0.0, max_value=20.0)
valid_dt = st.floats(min_value=0.0, max_value=20.0)

invalid_h = st.one_of(
    st.floats(max_value=0.0, allow_nan=False, allow_infinity=False), st.just(np.nan)
)


def array_like(array, *, fill_with=None):
    return hnp.arrays(
        dtype=array.dtype,
        shape=array.shape,
        elements=hnp.from_dtype(array.dtype) if fill_with is None else fill_with,
    )


def _reference(
    *,
    q_at_link,
    q_mean_at_link,
    water_slope_at_link,
    h_at_link,
    mannings_at_link,
    g,
    dt,
):
    numerator = q_mean_at_link - g * dt * h_at_link * water_slope_at_link

    denominator = 1.0 + g * dt * mannings_at_link**2 * abs(q_at_link) / h_at_link ** (
        7 / 3
    )

    return numerator / denominator


@given(
    size=st.integers(min_value=0, max_value=64),
    data=st.data(),
)
def test_calc_discharge_at_links_teeny_h(size, data):
    g, dt = 9.81, 10.0

    out = np.full(size, np.nan)
    where = np.arange(size, dtype=int)

    q = data.draw(array_like(out, fill_with=valid_q), label="discharge")
    q_mean = data.draw(array_like(out, fill_with=valid_q), label="mean discharge")
    slope = data.draw(array_like(out, fill_with=valid_slope), label="slope")
    mannings = data.draw(array_like(out, fill_with=valid_mannings), label="mannings")

    tiny_h = st.floats(min_value=0.0, max_value=1e-6, exclude_min=True)
    h = data.draw(array_like(out, fill_with=tiny_h), label="water depth")

    calc_discharge_at_links(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=g,
        dt=dt,
        where=where,
        out=out,
    )

    expected = np.zeros_like(out)
    assert_almost_equal(out, expected)


@given(
    g=valid_g,
    dt=valid_dt,
    size=st.integers(min_value=0, max_value=64),
    data=st.data(),
)
def test_calc_discharge_at_links_matches_ref(g, dt, size, data):
    out = np.full(size, np.nan)
    where = np.arange(size, dtype=int)

    q = data.draw(array_like(out, fill_with=valid_q), label="discharge")
    q_mean = data.draw(array_like(out, fill_with=valid_q), label="mean discharge")
    slope = data.draw(array_like(out, fill_with=valid_slope), label="slope")
    h = data.draw(array_like(out, fill_with=valid_h), label="water depth")
    mannings = data.draw(array_like(out, fill_with=valid_mannings), label="mannings")

    expected = _reference(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=g,
        dt=dt,
    )
    actual = calc_discharge_at_links(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=g,
        dt=dt,
        where=where,
        out=out,
    )
    assert actual is out
    assert_almost_equal(out, expected)


@given(
    g=valid_g,
    dt=valid_dt,
    size=st.integers(min_value=0, max_value=64),
    data=st.data(),
)
def test_calc_discharge_at_links_symmetrical(g, dt, size, data):
    out = np.full(size, np.nan)
    where = np.arange(size, dtype=int)

    q = data.draw(array_like(out, fill_with=valid_q), label="discharge")
    q_mean = data.draw(array_like(out, fill_with=valid_q), label="mean discharge")
    slope = data.draw(array_like(out, fill_with=valid_slope), label="slope")
    h = data.draw(array_like(out, fill_with=valid_h), label="water depth")
    mannings = data.draw(array_like(out, fill_with=valid_mannings), label="mannings")

    rtn = calc_discharge_at_links(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=g,
        dt=dt,
        where=where,
        out=out,
    )
    assert rtn is out
    assert not np.any(np.isnan(out))

    expected = -out

    q *= -1
    q_mean *= -1
    slope *= -1

    out = np.full(size, np.nan)
    calc_discharge_at_links(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=g,
        dt=dt,
        where=where,
        out=out,
    )
    assert not np.any(np.isnan(out))

    assert_array_equal(out, expected)


@given(
    size=st.integers(min_value=0, max_value=64),
    bad_h=invalid_h,
    data=st.data(),
)
def test_calc_discharge_at_links_writes_only_where_h_positive(size, bad_h, data):
    g, dt = 9.81, 1.0

    h = np.full(size, bad_h, dtype=float)
    mask = data.draw(hnp.arrays(np.bool_, size, elements=st.booleans()), label="mask")
    h[mask] = 1.0

    good = h > 0
    bad = (h <= 0) | np.isnan(h)
    assume(np.any(good) and np.any(bad))

    q = data.draw(array_like(h, fill_with=valid_q), label="discharge")
    q_mean = data.draw(array_like(h, fill_with=valid_q), label="mean discharge")
    slope = data.draw(array_like(h, fill_with=valid_slope), label="slope")
    mannings = data.draw(array_like(h, fill_with=valid_mannings), label="mannings")

    before = data.draw(array_like(h, fill_with=valid_q), label="initial water depth")
    out = before.copy()

    calc_discharge_at_links(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=g,
        dt=dt,
        where=np.arange(size, dtype=int),
        out=out,
    )

    assert_array_equal(out[bad], before[bad])


@given(
    size=st.integers(min_value=0, max_value=64),
    data=st.data(),
)
def test_calc_discharge_at_links_write_only_where(size, data):
    mask = data.draw(hnp.arrays(np.bool_, size, elements=st.booleans()), label="mask")
    where = np.flatnonzero(mask)
    before = np.full(size, np.nan)

    q = data.draw(array_like(before, fill_with=valid_q), label="discharge")
    q_mean = data.draw(array_like(before, fill_with=valid_q), label="mean discharge")
    slope = data.draw(array_like(before, fill_with=valid_slope), label="slope")
    h = data.draw(array_like(before, fill_with=valid_h), label="water depth")
    mannings = data.draw(array_like(before, fill_with=valid_mannings), label="mannings")

    out = calc_discharge_at_links(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=9.81,
        dt=1.0,
        where=where,
        out=before.copy(),
    )

    assert np.all(np.isnan(out[~mask]))
    assert np.all(out[mask] != before[mask])


@given(
    size=st.integers(min_value=0, max_value=64),
    data=st.data(),
)
def test_calc_discharge_at_links_frictionless(size, data):
    g, dt = 9.81, 10.0
    mannings = np.zeros(size, dtype=float)

    q = data.draw(array_like(mannings, fill_with=valid_q), label="discharge")
    q_mean = data.draw(array_like(mannings, fill_with=valid_q), label="mean discharge")
    slope = data.draw(array_like(mannings, fill_with=valid_slope), label="slope")
    h = data.draw(array_like(mannings, fill_with=valid_h), label="water depth")

    actual = np.full(size, np.nan)
    where = np.arange(size, dtype=int)

    calc_discharge_at_links(
        q_at_link=q,
        q_mean_at_link=q_mean,
        water_slope_at_link=slope,
        h_at_link=h,
        mannings_at_link=mannings,
        g=g,
        dt=dt,
        where=where,
        out=actual,
    )

    assert_almost_equal(actual, q_mean - g * dt * h * slope)
