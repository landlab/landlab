#! /usr/bin/env python
"""Unit tests for landlab.components.flexure.Flexure1D."""
import numpy as np
import pytest
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
    assert_array_less,
)

from landlab import RasterModelGrid
from landlab.components.flexure import Flexure1D

(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e3, 10e3), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_name(flex1d):
    """Test component name exists and is a string."""
    assert isinstance(flex1d.name, str)


def test_input_var_names(flex1d):
    """Test input_var_names is a tuple of strings."""
    assert isinstance(flex1d.input_var_names, tuple)
    for name in flex1d.input_var_names:
        assert isinstance(name, str)


def test_output_var_names(flex1d):
    """Test output_var_names is a tuple of strings."""
    assert isinstance(flex1d.output_var_names, tuple)
    for name in flex1d.output_var_names:
        assert isinstance(name, str)


def test_var_units(flex1d):
    """Test input/output var units."""
    assert isinstance(flex1d.units, tuple)
    for name, units in flex1d.units:
        assert name in flex1d.input_var_names + flex1d.output_var_names
        assert isinstance(units, str)


def test_var_mapping(flex1d):
    """Test input/output var mappings."""
    assert isinstance(flex1d._var_mapping, dict)
    for name in flex1d.input_var_names + flex1d.output_var_names:
        assert name in flex1d._var_mapping
        assert isinstance(flex1d._var_mapping[name], str)
        assert flex1d._var_mapping[name] in (
            "node",
            "link",
            "patch",
            "corner",
            "face",
            "cell",
        )


def test_var_doc(flex1d):
    """Test input/output var docs."""
    assert isinstance(flex1d._var_doc, dict)
    for name in flex1d.input_var_names + flex1d.output_var_names:
        assert name in flex1d._var_doc
        assert isinstance(flex1d._var_doc[name], str)


def test_calc_airy():
    """Test airy isostasy."""
    flex = Flexure1D(RasterModelGrid((3, 5)), method="airy")
    flex.load_at_node[:] = flex.gamma_mantle

    assert_array_equal(flex.dz_at_node, 0.0)
    flex.update()
    assert_array_equal(flex.dz_at_node, 1.0)


def test_with_method_flexure():
    n = 101
    i_mid = (n - 1) // 2
    flex = Flexure1D(RasterModelGrid((3, n)), method="flexure")
    flex.load_at_node[1, i_mid] = 1.0

    flex.update()

    assert np.argmax(flex.dz_at_node[1]) == i_mid
    assert np.all(
        flex.dz_at_node[:, i_mid::-1] == pytest.approx(flex.dz_at_node[:, i_mid:])
    )


def test_run_one_step():
    """Test the run_one_step method."""
    flex = Flexure1D(RasterModelGrid((3, 5)), method="airy")
    flex.load_at_node[:] = flex.gamma_mantle

    assert_array_equal(flex.dz_at_node, 0.0)
    flex.run_one_step()
    assert_array_equal(flex.dz_at_node, 1.0)


def test_with_one_row():
    """Test calculating on one row."""
    flex = Flexure1D(RasterModelGrid((3, 5)), method="airy", rows=1)
    flex.load_at_node[:] = flex.gamma_mantle

    assert_array_equal(flex.dz_at_node, 0.0)
    flex.update()
    assert_array_equal(flex.dz_at_node[0], 0.0)
    assert_array_equal(flex.dz_at_node[1], 1.0)
    assert_array_equal(flex.dz_at_node[2], 0.0)


def test_with_two_row():
    """Test calculating on one row."""
    flex = Flexure1D(RasterModelGrid((3, 5)), method="airy", rows=(0, 2))
    flex.load_at_node[:] = -flex.gamma_mantle

    assert_array_equal(flex.dz_at_node, 0.0)
    flex.update()
    assert_array_equal(flex.dz_at_node[0], -1.0)
    assert_array_equal(flex.dz_at_node[1], 0.0)
    assert_array_equal(flex.dz_at_node[2], -1.0)


def test_field_is_updated():
    """Test the output field is updated."""
    flex = Flexure1D(RasterModelGrid((3, 5)), method="airy", rows=(0, 2))
    flex.load_at_node[:] = -flex.gamma_mantle

    assert_array_equal(flex.dz_at_node, 0.0)
    flex.update()

    dz = flex.grid.at_node["lithosphere_surface__increment_of_elevation"]
    assert_array_equal(flex.dz_at_node.flatten(), dz)


def test_calc_flexure():
    """Test calc_flexure function."""
    x = np.arange(100.0)
    loads = np.ones(100)
    dz = Flexure1D.calc_flexure(x, loads, 1.0, 1.0)

    assert_array_less(0.0, dz)
    assert isinstance(dz, np.ndarray)
    assert dz.shape == loads.shape
    assert dz.dtype == loads.dtype


def test_calc_flexure_with_out_keyword():
    """Test calc_flexure out keyword."""
    x = np.arange(100.0)
    loads = np.ones(100)
    buffer = np.empty_like(x)
    dz = Flexure1D.calc_flexure(x, loads, 1.0, 1.0, out=buffer)
    assert np.may_share_memory(dz, buffer)


def test_calc_flexure_with_multiple_rows():
    """Test calc_flexure with multiple rows of loads."""
    x = np.arange(100.0) * 1e3
    loads = np.ones(500).reshape((5, 100))
    dz = Flexure1D.calc_flexure(x, loads, 1e4, 1.0)

    assert_array_less(0.0, dz)
    assert isinstance(dz, np.ndarray)
    assert dz.shape == loads.shape
    assert dz.dtype == loads.dtype

    for row in range(5):
        assert_array_almost_equal(dz[0], dz[row])


DEPENDS_ON = {
    "eet": ("rigidity", "alpha"),
    "youngs": ("rigidity", "alpha"),
    "rho_water": ("gamma_mantle", "alpha"),
    "rho_mantle": ("gamma_mantle", "alpha"),
    "gravity": ("gamma_mantle", "alpha"),
}


def test_setters(flexure_keyword):
    EPS = 1e-6
    flex = Flexure1D(RasterModelGrid((3, 5)))
    val_before = {}
    for name in DEPENDS_ON[flexure_keyword]:
        val_before[name] = 1.0 * getattr(flex, name)
    for name in DEPENDS_ON[flexure_keyword]:
        setattr(
            flex, flexure_keyword, getattr(flex, flexure_keyword) * (1.0 + EPS) + EPS
        )
        assert val_before[name] != getattr(flex, name)


def test_method_keyword():
    """Test using the method keyword."""
    flex = Flexure1D(RasterModelGrid((3, 5)), method="airy")
    assert flex.method == "airy"
    flex = Flexure1D(RasterModelGrid((3, 5)), method="flexure")
    assert flex.method == "flexure"
    with pytest.raises(ValueError):
        Flexure1D(RasterModelGrid((3, 5)), method="Flexure")


def test_flexure_keywords(flexure_keyword):
    flex = Flexure1D(RasterModelGrid((3, 5)), **{flexure_keyword: 1.0})
    assert getattr(flex, flexure_keyword) == 1.0
    assert isinstance(getattr(flex, flexure_keyword), float)
    with pytest.raises(ValueError):
        flex = Flexure1D(RasterModelGrid((3, 5)), **{flexure_keyword: -1})


def test_x_at_node():
    """Test x_at_node is reshaped and shares memory with the grid."""
    flex = Flexure1D(RasterModelGrid((3, 5)))

    assert_array_equal(
        flex.x_at_node,
        [
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [0.0, 1.0, 2.0, 3.0, 4.0],
        ],
    )


def test_dz_at_node():
    """Test dz_at_node is reshaped and shares memory with its field."""
    flex = Flexure1D(RasterModelGrid((3, 5)))

    vals = flex.grid.at_node["lithosphere_surface__increment_of_elevation"]
    assert_array_equal(vals, 0.0)

    assert np.may_share_memory(vals, flex.dz_at_node)
    assert flex.dz_at_node.shape == (3, 5)


def test_load_at_node():
    """Test load_at_node is reshaped and shares memory with its field."""
    flex = Flexure1D(RasterModelGrid((3, 5)))

    vals = flex.grid.at_node["lithosphere__increment_of_overlying_pressure"]
    assert_array_equal(vals, 0.0)

    assert np.may_share_memory(vals, flex.load_at_node)
    assert flex.load_at_node.shape == (3, 5)


def test_x_is_contiguous():
    """Test that x_at_node is contiguous."""
    flex = Flexure1D(RasterModelGrid((3, 5)))
    assert flex.x_at_node.flags["C_CONTIGUOUS"]


def test_dz_is_contiguous():
    """Test that dz_at_node is contiguous."""
    flex = Flexure1D(RasterModelGrid((3, 5)))
    assert flex.dz_at_node.flags["C_CONTIGUOUS"]


def test_load_is_contiguous():
    """Test that load_at_node is contiguous."""
    flex = Flexure1D(RasterModelGrid((3, 5)))
    assert flex.load_at_node.flags["C_CONTIGUOUS"]


def test_subside_loads():
    flex = Flexure1D(RasterModelGrid((3, 5)), method="airy")
    dz_airy = flex.subside_loads([0.0, 0.0, flex.gamma_mantle, 0.0, 0])
    assert dz_airy.shape == flex.grid.shape
    assert np.all(dz_airy == [0.0, 0.0, 1.0, 0.0, 0])

    flex = Flexure1D(RasterModelGrid((3, 5)), method="flexure")
    dz_flexure = flex.subside_loads([0.0, 0.0, flex.gamma_mantle, 0.0, 0])

    assert dz_flexure.shape == flex.grid.shape
    assert np.argmax(dz_flexure) == 2
    assert np.all(dz_flexure[:, 2] < dz_airy[:, 2])
    assert np.all(dz_flexure[:, 2::-1] == pytest.approx(dz_flexure[:, 2:]))
