"""Tests for StreamPowerEroder's sp_type/m_sp/n_sp/a_sp/b_sp/c_sp validation.

The exponent derivation itself is delegated to four module-level, public
functions (``stream_power_exponents``, ``total_stream_power_exponents``,
``unit_stream_power_exponents``, ``shear_stress_stream_power_exponents``),
so a user can precompute ``m``/``n`` and pass them directly via
``sp_type="set_mn"`` instead of specifying ``sp_type``/``a_sp``/``b_sp``/
``c_sp`` on the component itself. Those are tested directly, independent
of any grid.

Validation raises either ``requireit.ValidationError`` (any check that
delegates to a ``requireit`` validator: ``require_one_of``,
``require_nonnegative``, `require_between``, and the local ``require_none``
/``require_not_none`` helpers) or a plain ``ValueError`` (the two checks
that are hand-written rather than delegated: ``m_sp``/``n_sp`` must stay
at their defaults when deriving from ``a_sp``/``b_sp``/``c_sp``, and
``channel_width_field`` must be ``None`` when ``sp_type == "Total"``).
"""

import numpy as np
import pytest
from requireit import ValidationError

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import StreamPowerEroder
from landlab.components.stream_power.stream_power import (
    shear_stress_stream_power_exponents,
)
from landlab.components.stream_power.stream_power import stream_power_exponents
from landlab.components.stream_power.stream_power import total_stream_power_exponents
from landlab.components.stream_power.stream_power import unit_stream_power_exponents


@pytest.fixture
def grid():
    grid = RasterModelGrid((5, 5))
    grid.add_zeros("topographic__elevation", at="node")
    grid.at_node["topographic__elevation"][:] = grid.node_x + grid.node_y
    flow_accumulator = FlowAccumulator(grid, flow_director="D8")
    flow_accumulator.run_one_step()
    return grid


class TestStreamPowerExponents:
    """The set_mn identity case: m_sp/n_sp are used, and validated, as-is."""

    def test_pass_through(self):
        assert stream_power_exponents(0.3, 0.7) == (0.3, 0.7)

    def test_negative_m_raises(self):
        with pytest.raises(ValidationError):
            stream_power_exponents(-1.0, 1.0)

    def test_negative_n_raises(self):
        with pytest.raises(ValidationError):
            stream_power_exponents(0.5, -1.0)


class TestTotalStreamPowerExponents:
    def test_formula(self):
        m, n = total_stream_power_exponents(a_sp=2.0, c_sp=3.0)
        assert m == 2.0 * 3.0
        assert n == 2.0

    def test_negative_a_sp_raises(self):
        with pytest.raises(ValidationError):
            total_stream_power_exponents(a_sp=-1.0, c_sp=1.0)

    def test_negative_c_sp_raises(self):
        with pytest.raises(ValidationError):
            total_stream_power_exponents(a_sp=1.0, c_sp=-1.0)


class TestUnitStreamPowerExponents:
    def test_formula(self):
        m, n = unit_stream_power_exponents(a_sp=2.0, b_sp=0.25, c_sp=3.0)
        assert m == 2.0 * 3.0 * (1.0 - 0.25)
        assert n == 2.0

    def test_negative_a_sp_raises(self):
        with pytest.raises(ValidationError):
            unit_stream_power_exponents(a_sp=-1.0, b_sp=0.5, c_sp=1.0)

    def test_negative_c_sp_raises(self):
        with pytest.raises(ValidationError):
            unit_stream_power_exponents(a_sp=1.0, b_sp=0.5, c_sp=-1.0)

    def test_negative_b_sp_raises(self):
        with pytest.raises(ValidationError):
            unit_stream_power_exponents(a_sp=1.0, b_sp=-0.1, c_sp=1.0)

    def test_b_sp_greater_than_one_raises(self):
        with pytest.raises(ValidationError):
            unit_stream_power_exponents(a_sp=1.0, b_sp=1.5, c_sp=1.0)

    def test_b_sp_of_one_is_the_inclusive_boundary(self):
        m, n = unit_stream_power_exponents(a_sp=1.0, b_sp=1.0, c_sp=1.0)
        assert m == 0.0
        assert n == 1.0


class TestShearStressStreamPowerExponents:
    def test_formula(self):
        m, n = shear_stress_stream_power_exponents(a_sp=2.0, b_sp=0.25, c_sp=3.0)
        assert m == 2.0 * 2.0 * 3.0 * (1.0 - 0.25) / 3.0
        assert n == 2.0 * 2.0 / 3.0

    def test_negative_a_sp_raises(self):
        with pytest.raises(ValidationError):
            shear_stress_stream_power_exponents(a_sp=-1.0, b_sp=0.5, c_sp=1.0)

    def test_b_sp_out_of_range_raises(self):
        with pytest.raises(ValidationError):
            shear_stress_stream_power_exponents(a_sp=1.0, b_sp=1.5, c_sp=1.0)


class TestSetMN:
    """sp_type='set_mn' (the default) uses m_sp/n_sp directly."""

    def test_defaults(self, grid):
        sp = StreamPowerEroder(grid)
        assert sp._m == 0.5
        assert sp._n == 1.0

    def test_custom_m_and_n(self, grid):
        sp = StreamPowerEroder(grid, m_sp=0.3, n_sp=0.7)
        assert sp._m == 0.3
        assert sp._n == 0.7

    def test_negative_m_raises(self, grid):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, m_sp=-1.0)

    def test_negative_n_raises(self, grid):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, n_sp=-1.0)

    @pytest.mark.parametrize("param", ("a_sp", "b_sp", "c_sp"))
    def test_abc_forbidden_with_set_mn(self, grid, param):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type="set_mn", **{param: 1.0})

    def test_precomputed_exponents_can_be_passed_through_set_mn(self, grid):
        """The whole point of the public *_exponents functions: a user can
        compute m/n themselves and bypass sp_type/a_sp/b_sp/c_sp entirely.
        """
        m, n = unit_stream_power_exponents(a_sp=1.0, b_sp=0.5, c_sp=1.0)
        sp = StreamPowerEroder(grid, sp_type="set_mn", m_sp=m, n_sp=n)
        assert sp._m == m
        assert sp._n == n


class TestSpTypeValidation:
    def test_unrecognized_sp_type_raises(self, grid):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type="not_a_type", a_sp=1.0, c_sp=1.0)

    @pytest.mark.parametrize("sp_type", ("Total", "Unit", "Shear_stress"))
    @pytest.mark.parametrize("mn_param", ("m_sp", "n_sp"))
    def test_m_n_must_stay_default_for_non_set_mn(self, grid, sp_type, mn_param):
        kwargs = {"a_sp": 1.0, "b_sp": 0.5, "c_sp": 1.0, mn_param: 0.3}
        with pytest.raises(ValueError):
            StreamPowerEroder(grid, sp_type=sp_type, **kwargs)

    @pytest.mark.parametrize("sp_type", ("Total", "Unit", "Shear_stress"))
    def test_a_sp_is_required(self, grid, sp_type):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type=sp_type, b_sp=0.5, c_sp=1.0)

    @pytest.mark.parametrize("sp_type", ("Total", "Unit", "Shear_stress"))
    def test_negative_a_sp_raises(self, grid, sp_type):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type=sp_type, a_sp=-1.0, b_sp=0.5, c_sp=1.0)

    @pytest.mark.parametrize("sp_type", ("Unit", "Shear_stress"))
    def test_negative_b_sp_raises(self, grid, sp_type):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type=sp_type, a_sp=1.0, b_sp=-1.0, c_sp=1.0)

    @pytest.mark.parametrize("sp_type", ("Total", "Unit", "Shear_stress"))
    def test_negative_c_sp_raises(self, grid, sp_type):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type=sp_type, a_sp=1.0, b_sp=0.5, c_sp=-1.0)


class TestExponentDerivation:
    """StreamPowerEroder dispatches to the matching *_exponents function."""

    def test_total(self, grid):
        sp = StreamPowerEroder(grid, sp_type="Total", a_sp=2.0, c_sp=3.0)
        assert sp._m == 2.0 * 3.0
        assert sp._n == 2.0

    def test_unit(self, grid):
        sp = StreamPowerEroder(grid, sp_type="Unit", a_sp=2.0, b_sp=0.25, c_sp=3.0)
        assert sp._m == 2.0 * 3.0 * (1.0 - 0.25)
        assert sp._n == 2.0

    def test_shear_stress(self, grid):
        sp = StreamPowerEroder(
            grid, sp_type="Shear_stress", a_sp=2.0, b_sp=0.25, c_sp=3.0
        )
        assert sp._m == 2.0 * 2.0 * 3.0 * (1.0 - 0.25) / 3.0
        assert sp._n == 2.0 * 2.0 / 3.0


class TestChannelWidthFieldInteraction:
    """channel_width_field's interaction with b_sp, gated on sp_type."""

    def test_forbidden_for_total(self, grid):
        with pytest.raises(ValueError):
            StreamPowerEroder(
                grid,
                sp_type="Total",
                a_sp=1.0,
                c_sp=1.0,
                channel_width_field=np.ones(grid.number_of_nodes),
            )

    @pytest.mark.parametrize("sp_type", ("Unit", "Shear_stress"))
    def test_b_sp_required_without_a_width_field(self, grid, sp_type):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type=sp_type, a_sp=1.0, c_sp=1.0)

    @pytest.mark.parametrize("sp_type", ("Unit", "Shear_stress"))
    def test_b_sp_defaults_to_zero_when_width_field_given(self, grid, sp_type):
        sp = StreamPowerEroder(
            grid,
            sp_type=sp_type,
            a_sp=1.0,
            c_sp=1.0,
            channel_width_field=np.ones(grid.number_of_nodes),
        )
        # b_sp defaulted to 0.0, so m matches the b_sp=0 formula.
        if sp_type == "Unit":
            assert sp._m == 1.0 * 1.0 * (1.0 - 0.0)
        else:
            assert sp._m == 2.0 * 1.0 * 1.0 * (1.0 - 0.0) / 3.0

    @pytest.mark.parametrize("sp_type", ("Unit", "Shear_stress"))
    def test_b_sp_and_width_field_can_both_be_given(self, grid, sp_type):
        sp = StreamPowerEroder(
            grid,
            sp_type=sp_type,
            a_sp=1.0,
            b_sp=0.5,
            c_sp=1.0,
            channel_width_field=np.ones(grid.number_of_nodes),
        )
        if sp_type == "Unit":
            assert sp._m == 1.0 * 1.0 * (1.0 - 0.5)
        else:
            assert sp._m == 2.0 * 1.0 * 1.0 * (1.0 - 0.5) / 3.0


def _b_sp_kwargs(sp_type):
    """b_sp is only accepted for Unit/Shear_stress -- must stay unset for Total."""
    return {} if sp_type == "Total" else {"b_sp": 0.5}


class TestDischargeFieldInteraction:
    """discharge_field's interaction with c_sp, gated on sp_type."""

    @pytest.mark.parametrize("sp_type", ("Total", "Unit", "Shear_stress"))
    def test_c_sp_required_with_default_discharge_field(self, grid, sp_type):
        with pytest.raises(ValidationError):
            StreamPowerEroder(grid, sp_type=sp_type, a_sp=1.0, **_b_sp_kwargs(sp_type))

    @pytest.mark.parametrize("sp_type", ("Total", "Unit", "Shear_stress"))
    def test_c_sp_defaults_to_one_with_custom_discharge_field(self, grid, sp_type):
        # FlowAccumulator's D8 director already populates this field.
        sp = StreamPowerEroder(
            grid,
            sp_type=sp_type,
            a_sp=1.0,
            discharge_field="surface_water__discharge",
            **_b_sp_kwargs(sp_type),
        )
        expected_n = 2.0 / 3.0 if sp_type == "Shear_stress" else 1.0
        assert sp._n == expected_n

    @pytest.mark.parametrize("sp_type", ("Total", "Unit", "Shear_stress"))
    def test_c_sp_and_custom_discharge_field_can_both_be_given(self, grid, sp_type):
        sp = StreamPowerEroder(
            grid,
            sp_type=sp_type,
            a_sp=1.0,
            c_sp=1.0,
            discharge_field="surface_water__discharge",
            **_b_sp_kwargs(sp_type),
        )
        if sp_type == "Total":
            assert sp._m == 1.0 * 1.0
        elif sp_type == "Unit":
            assert sp._m == 1.0 * 1.0 * (1.0 - 0.5)
        else:
            assert sp._m == 2.0 * 1.0 * 1.0 * (1.0 - 0.5) / 3.0
