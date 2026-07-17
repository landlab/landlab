# tests/helper_functions/test_stability_displacement.py

import numpy as np

from landlab import RasterModelGrid
from landlab.components import ShallowLandslider


def make_grid(n=5, spacing=10.0):
    mg = RasterModelGrid((n, n), xy_spacing=spacing)
    mg.add_ones("topographic__elevation", at="node")
    mg.add_ones("soil__depth", at="node")
    return mg


def test_factor_of_safety_matches_formula(monkeypatch):
    """Validate the FoS implementation against the analytic formula."""
    mg = make_grid()
    slope_rad = np.deg2rad(30.0)

    def constant_slope(elevs=None, **kwargs):
        slope = np.ones(mg.number_of_nodes) * slope_rad
        if kwargs.get("return_components", False):
            return slope, (slope, slope)
        return slope

    monkeypatch.setattr(mg, "calc_slope_at_node", constant_slope)
    monkeypatch.setattr(
        mg,
        "calc_aspect_at_node",
        lambda **kwargs: np.zeros(mg.number_of_nodes),
    )

    coh = 1000.0
    phi_rad = np.deg2rad(30.0)
    gamma_s = 15e3
    gamma_w = 9.8e3
    sub = 0.0

    comp = ShallowLandslider(
        mg,
        cohesion_eff=coh,
        angle_int_frict=30.0,
        submerged_soil_proportion=sub,
        update_soil=False,
    )

    fos = comp._factor_of_safety(
        mg,
        coh,
        phi_rad,
        submerged_soil_proportion=sub,
        soil_unit_weight=gamma_s,
        water_unit_weight=gamma_w,
    )

    soil_depth = mg.at_node["soil__depth"]
    psi = sub * gamma_w * soil_depth
    slope = np.ones(mg.number_of_nodes) * slope_rad

    expected = (
        (coh - psi * np.tan(phi_rad)) / (gamma_s * soil_depth * np.sin(slope))
    ) + (np.tan(phi_rad) / np.tan(slope))

    assert np.allclose(fos, expected, rtol=1e-6, atol=1e-9)


def test_compute_stability_uses_configured_submerged_proportion(monkeypatch):
    """The public pipeline forwards saturation to the factor-of-safety equation."""
    mg = make_grid()
    slope = np.full(mg.number_of_nodes, np.deg2rad(30.0))
    monkeypatch.setattr(mg, "calc_slope_at_node", lambda **kwargs: slope)
    monkeypatch.setattr(
        mg, "calc_aspect_at_node", lambda **kwargs: np.zeros(mg.number_of_nodes)
    )
    comp = ShallowLandslider(
        mg,
        cohesion_eff=1000.0,
        angle_int_frict=30.0,
        submerged_soil_proportion=0.8,
    )

    comp._compute_stability()
    expected = comp._factor_of_safety(
        mg,
        1000.0,
        np.deg2rad(30.0),
        submerged_soil_proportion=0.8,
    )

    assert np.allclose(comp.results["factor_of_safety"], expected, equal_nan=True)


def test_critical_transient_acceleration(monkeypatch):
    mg = make_grid()
    slope_rad = np.deg2rad(20)

    def constant_slope(elevs=None, **kwargs):
        slope = np.ones(mg.number_of_nodes) * slope_rad
        if kwargs.get("return_components", False):
            return slope, (slope, slope)
        return slope

    monkeypatch.setattr(mg, "calc_slope_at_node", constant_slope)
    monkeypatch.setattr(
        mg,
        "calc_aspect_at_node",
        lambda **kwargs: np.zeros(mg.number_of_nodes),
    )

    g = 9.81
    a_h = np.ones(mg.number_of_nodes) * 0.3 * g
    a_v = np.ones(mg.number_of_nodes) * 0.1 * g

    phi = np.deg2rad(30)
    coh = 500.0
    gamma_s = 15e3
    gamma_w = 9.8e3
    sub = 0.0

    comp = ShallowLandslider(
        mg,
        cohesion_eff=coh,
        angle_int_frict=30,
        submerged_soil_proportion=sub,
    )
    ac, aslip, adiff = comp._critical_transient_acceleration(
        mg,
        coh,
        phi,
        sub,
        a_h=a_h,
        a_v=a_v,
        soil_unit_weight=gamma_s,
        water_unit_weight=gamma_w,
        g=g,
    )

    soil_depth = mg.at_node["soil__depth"]
    psi = sub * gamma_w * soil_depth

    a_c_simple = (
        np.tan(phi)
        * (g * np.cos(slope_rad) - a_v * np.cos(slope_rad) - a_h * np.sin(slope_rad))
        + ((g * coh) - (psi * g * np.tan(phi))) / (gamma_s * soil_depth)
        - g * np.sin(slope_rad)
    )
    a_c_simple[mg.boundary_nodes] = 0
    a_s = a_h * np.cos(slope_rad) - a_v * np.sin(slope_rad)

    assert np.allclose(ac, a_c_simple, rtol=1e-6)
    assert np.allclose(aslip, a_s)
    assert np.allclose(adiff, a_s - a_c_simple)


def test_newmark_displacement_and_masking():
    mg = make_grid()
    comp = ShallowLandslider(
        mg,
        cohesion_eff=10,
        angle_int_frict=30,
        compute_displacement=True,
    )

    # construct simple labels
    labels = np.zeros(mg.shape, dtype=int)
    labels[2:4, 2:4] = 1
    diff = np.zeros(mg.number_of_nodes)
    active_idx = np.where(labels.ravel() == 1)[0]
    diff[active_idx] = 2.0  # m/s2

    disp = comp._calculate_newmark_displacement(
        a_difference_1d=diff,
        selected_labels_2d=labels,
        time_shaking_2d=np.ones(mg.shape) * 3.0,
    )

    unlabeled = np.where(labels.ravel() == 0)[0]
    assert np.all(np.isnan(disp[unlabeled]))
    assert np.allclose(disp[active_idx], 0.5 * 2.0 * 9.0)


def test_newmark_displacement_threshold_behavior():
    mg = make_grid()
    comp = ShallowLandslider(
        mg,
        cohesion_eff=10,
        angle_int_frict=30,
        compute_displacement=True,
        displacement_threshold=5.0,
    )

    comp.run_one_step()
    disp = mg.at_node["landslide__newmark_displacement"]

    # All high displacement nodes must exceed threshold
    for idx in comp._high_disp_nodes:
        assert disp[idx] > 5.0


def test_critical_acceleration_sets_boundary_to_zero(monkeypatch):
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)

    monkeypatch.setattr(
        mg,
        "calc_slope_at_node",
        lambda elevs=None, **kwargs: np.ones(mg.number_of_nodes) * np.deg2rad(10),
    )

    ac, *_ = comp._critical_transient_acceleration(
        mg,
        10,
        np.deg2rad(30),
        0.0,
        a_h=np.zeros(mg.number_of_nodes),
        a_v=np.zeros(mg.number_of_nodes),
    )

    assert np.all(ac[mg.boundary_nodes] == 0)
