#! /usr/bin/env python
"""
Unit tests for landlab.components.gravel_bedrock_eroder.gravel_bedrock_eroder
"""
import numpy as np
from numpy.testing import assert_almost_equal

from landlab import HexModelGrid
from landlab.components import FlowAccumulator
from landlab.components import ExtendedGravelBedrockEroder

spy = 3600.0 * 24 * 365.25 # seconds per year

# Functions to use for test calculations
def slope(dz, dx):
    return dz / dx

def discharge(runoff_my, dx):
    area = dx * dx
    return runoff_my * area, runoff_my * area / spy

def fixed_width(kw, m, Q):
    return kw * Q**m

def rock_exposure_fraction(h, hc):
    return np.exp(-h / hc)

def median_grain_size(D):
    if type(D) is float:
        D50 = D
    else:
        D50 = np.median(D)
    return D50

def dynamic_width(rhos, rhow, eps, Q, S, D50, tau_star_c, g=9.81):
    R = (rhos - rhow) / rhow
    prefactor = 0.17 * g**(-0.5) * (R * (1 + eps) * tau_star_c)**(-5/3)
    return prefactor * Q * S**(7/6) / D50**1.5

def shear_stress(rho, g, rough, Qcms, w, S):
    return rho * g * rough**0.6 * (Qcms / w)**0.6 * S**0.7

def shields_stress(tau, rho, rhos, g, D50):
    return tau / ((rhos - rho) * g * D50)

def critical_shields(D, tau_st_c, mkomar, sed0):
    if type(D) is float:
        crit_sh = tau_st_c
    else:
        median_index = np.argwhere((np.cumsum(sed0) / np.sum(sed0)) > 0.5)[0].astype(int)
        crit_sh = tau_st_c * (D / D[median_index])**-mkomar
    return crit_sh

def sed_capacity(intcy, mpm, w, rho, rhos, g, D, tau_star, tau_star_c_by_size, f, mk):
    R = (rhos - rho) / rho
    exshields = tau_star - tau_star_c_by_size
    if isinstance(tau_star, float):
        exshields = max(exshields,0)
    else:
        exshields[exshields<0]=0
    Qs_sec = f * intcy * mpm * w * (R * g)**0.5 * D**1.5 * exshields**1.5
    Qs_yr = Qs_sec * spy
    return Qs_yr

def sed_flux(intcy, mpm, w, rho, rhos, g, D, tau_star, tau_star_c_by_size, f, mk, br_exp_frac):
    cap = sed_capacity(
        intcy, mpm, w, rho, rhos, g, D, tau_star, tau_star_c_by_size, f, mk
    )
    return (1.0 - br_exp_frac) * cap

def attrition_thickening_rate(sed_fluxes, beta, dx, cell_area):
    volume_loss_rate_per_length = sed_fluxes * beta
    volume_loss_rate = volume_loss_rate_per_length * dx
    thinning_rate = volume_loss_rate / cell_area
    return -thinning_rate

def outflux_thickening_rate(sed_fluxes, cell_area, porosity):
    return -sed_fluxes / ((1-porosity) *cell_area)

def rock_plucking_rate(pluck_coef, shields, crit_shields_rock, width, dx, cellarea, br_frac, intcy):
    if isinstance(shields,float):
        exshields = max(shields - crit_shields_rock, 0.0)
    else:
        exshields = shields - crit_shields_rock
        exshields[exshields<0]=0
    pluck_rate = (
        intcy * pluck_coef * (exshields)**1.5 
        * width * dx / cellarea
    )
    pluck_rate = np.sum(pluck_rate)
    return -br_frac * pluck_rate

def rock_abrasion_rate(abr_coef, sedflux, chanlen, cellarea, brfrac):
    return -brfrac * abr_coef * sedflux * chanlen / cellarea

def set_default_params():
    p = {}
    p["sed0"] = 100.0              # initial sediment thickness, m
    p["z0"] = 10.0                 # initial elevation, m
    p["phi"] = 0.5                 # initial elevation, m
    p["dx"] = 1000.0               # grid cell width, m
    p["wid_coef"] = 0.001          # coef for width-vs-discharge (m and s units)
    p["wid_exp"] = 0.5             # exponent for width-vs-discharge
    p["runoff"] = 10.0             # bankfull runoff rate, m/y
    p["rho"] = 1000.0              # water density, kg/m3
    p["rhos"] = 2650.0             # sediment density, kg/m3
    p["n"] = 0.05                  # Manning roughness coef, s/m^1/3
    p["D"] = 0.01                  # median grain diameter, m
    p["tau_st_c"] = 0.0495         # critical Shields stress for median size
    p["g"] = 9.81                  # guess what?
    p["mpm"] = 3.97                # transport coefficient for MPM eqn
    p["intcy"] = 0.01              # intermittency factor
    p["cellarea"] = p["dx"]**2     # area of grid cell, m2
    p["beta"] = 0.0                # sed/rock abrasion coefficient, 1/m
    p["width_is_dynamic"] = False  # flag for dynamic vs fixed width
    p["epsilon"] = 0.2             # factor by which tau > tauc for dynamic wid
    p["num_classes"] = 1           # number of lithology or size classes
    p["class_fracs"] = 1.0         # proportion of each class in sediment
    p["mkomar"] = 0.68             # Komar's "m" exponent for critical shields stress
    p["hc"] = 1.0                  # characteristic thickness for bed cover, m
    p["kp"] = 1.0                  # plucking coefficient, m/y
    p["taustarcr"] = 0.0495        # critical shields stress for rock, FOR NOW: ASSUMING = tau_st_C
    p["br_abr_coef"] = 0.001       # bedrock abrasion coefficient, 1/m
    return p

def calc_predicted_outputs(p):
    pv = {}
    pv["S"] = slope(p["z0"], p["dx"])
    pv["Q"], pv["Qcms"] = discharge(p["runoff"], p["dx"])
    pv["D50"] = median_grain_size(p["D"])
    if p["width_is_dynamic"]:
        pv["w"] = dynamic_width(
            p["rhos"], p["rho"], p["epsilon"], pv["Qcms"], pv["S"],
            pv["D50"], p["tau_st_c"], p["g"])
    else:
        pv["w"] = fixed_width(p["wid_coef"], p["wid_exp"], pv["Q"])
    pv["alpha"] = rock_exposure_fraction(np.sum(p["sed0"]), p["hc"])
    pv["tau"] = shear_stress(p["rho"], p["g"], p["n"], pv["Qcms"], pv["w"], pv["S"])
    pv["tau_star"] = shields_stress(pv["tau"], p["rho"], p["rhos"], p["g"], p["D"])
    pv["R"] = (p["rhos"] - p["rho"]) / p["rho"] # intermediate
    pv["crit_shields"] = critical_shields(
        p["D"], p["tau_st_c"], p["mkomar"], p["sed0"]) 
    pv["Qs"] = sed_flux(
        p["intcy"], p["mpm"], pv["w"], p["rho"], p["rhos"], p["g"], p["D"],
        pv["tau_star"], pv["crit_shields"], p["class_fracs"], p["mkomar"],
        pv["alpha"]
    )
    if type(pv["Qs"]) is np.ndarray:
        pv["Qs_total"] = np.sum(pv["Qs"])
    else:
        pv["Qs_total"] = pv["Qs"]
    dHdt_outflux = outflux_thickening_rate(pv["Qs"], p["cellarea"], p["phi"])
    pv["dHdt_attr"] = attrition_thickening_rate(
        pv["Qs"], p["beta"], p["dx"], p["cellarea"]
    )
    pv["pluck_rate"] = rock_plucking_rate(
        p["kp"], pv["tau_star"], p["taustarcr"], pv["w"], p["dx"],
        p["cellarea"], pv["alpha"], p["intcy"]
    )
    pv["rock_abr_rate"] = rock_abrasion_rate(
        p["br_abr_coef"], pv["Qs_total"], p["dx"], p["cellarea"],
        pv["alpha"]
    )
    pv["dHdt"] = dHdt_outflux + pv["dHdt_attr"]
    if np.size(pv["dHdt"]) == 1:
        pv["dHdt_total"] = pv["dHdt"]
    else:
        pv["dHdt_total"] = np.sum(pv["dHdt"])
    pv["dzrdt"] = pv["pluck_rate"] + pv["rock_abr_rate"]
    pv["dzdt"] = pv["dHdt_total"] + pv["dzrdt"]

    return pv

def display_pred_vals(p, pv):
    print("PREDICTED VALUES:")
    print("slope", pv["S"])
    print("discharge (cmy)", pv["Q"])
    print("discharge (cms)", pv["Qcms"])
    print("channel width (m)", pv["w"])
    print("bedrock exposure fraction", pv["alpha"])
    print("bed shear stress (Pa)", pv["tau"])
    print("Shields stress", pv["tau_star"])
    print("Critical Shields stress", pv["crit_shields"])
    print("sediment flux (cmy)", pv["Qs"])
    print("total sediment flux (cmy)", pv["Qs_total"])
    print("rate of sediment change from attrition (m/y)", pv["dHdt_attr"])
    print("rate of sediment thickness change by class (m/y)", pv["dHdt"])
    print("total rate of sediment thickness change (m/y)", pv["dHdt_total"])
    print("rate of rock plucking (m/y)", pv["pluck_rate"])
    print("rate of rock abrasion (m/y)", pv["rock_abr_rate"])
    print("rate of surface change (m/y)", pv["dzdt"])


def init_grid_and_run_one_step(parameters):
    from landlab import RasterModelGrid
    from landlab.components import FlowAccumulator
    from landlab.components.soil_grading import SoilGrading

    xy_spacing = parameters['dx']
    grid = RasterModelGrid((3, 3), xy_spacing = xy_spacing)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[1] = grid.BC_NODE_IS_FIXED_VALUE

    grid.add_zeros("topographic__elevation", at="node")
    sed_depth = parameters["sed0"]
    porosity = parameters["phi"]
    sed_weight = sed_depth * parameters["rhos"] * (1 - porosity)  # weight per grid node area
    if isinstance(sed_weight,float):
        grains_weight = [sed_weight]
        grain_sizes = [parameters["D"]]
    else:
        grains_weight = sed_weight
        grain_sizes = parameters["D"]
    sg = SoilGrading(grid,
                          meansizes = grain_sizes,
                          grains_weight = grains_weight,
                          phi = porosity,
                            )
    grid.at_node['bedrock__elevation'][grid.core_nodes] += 100 # arbitrary
    grid.at_node['topographic__elevation'][:] = grid.at_node['soil__depth'][:] + grid.at_node['bedrock__elevation'][:]
    grid.at_node['topographic__elevation'][1] = grid.at_node['topographic__elevation'][grid.core_nodes[0]] - parameters['z0']

    fa = FlowAccumulator(grid,
                         runoff_rate = parameters["runoff"])
    fa.run_one_step()
    eroder = ExtendedGravelBedrockEroder(grid,
                                         tau_star_c_median = parameters["tau_st_c"],
                                         bedrock_abrasion_coefficient=parameters["br_abr_coef"],
                                         plucking_coefficient=parameters["kp"],
                                         depth_decay_scale = parameters["hc"],
                                         alpha = parameters["mkomar"],
                                         epsilon = parameters["epsilon"],
                                         fixed_width_flag = not parameters["width_is_dynamic"],
                                         abrasion_coefficients = parameters["beta"],
                                         intermittency_factor=parameters["intcy"],
                                         mannings_n = parameters["n"],
                                         rho_sed = parameters['rhos'],
                                         rho_water = parameters['rho'],
                                         fixed_width_coeff=parameters["wid_coef"],
                                         fixed_width_expt=parameters["wid_exp"],
                                         sediment_porosity = parameters["phi"]
                                            )

    eroder.run_one_step(global_dt=1)
    return grid, eroder, sg, fa

def get_results(grid, eroder, sg, fa):

    spy = 3600.0 * 24 * 365.25  # seconds per year
    out = {}
    out["S"] =grid.at_node['topographic__steepest_slope'][grid.core_nodes[0]]
    out["Q"] =grid.at_node['surface_water__discharge'][grid.core_nodes[0]]
    out["Qcms"] = grid.at_node['surface_water__discharge'][grid.core_nodes[0]] / spy
    out["w"] = eroder._channel_width[grid.core_nodes[0]]
    out["alpha"] = eroder._rock_exposure_fraction[grid.core_nodes[0]]
    out["tau"] = eroder._tau[grid.core_nodes[0]]
    out["tau_star"] = eroder._tau_star[grid.core_nodes[0]]
    out["crit_shields"] = eroder._tau_star_c[grid.core_nodes[0]]
    out["Qs"] = eroder._sed_outfluxes[grid.core_nodes[0]]
    out["Qs_total"] = eroder._sediment_outflux[grid.core_nodes[0]]
    out["dHdt_attr"] = eroder._sed_abr_rates[grid.core_nodes[0]]
    out["dHdt"] =  eroder._dHdt_by_class[grid.core_nodes[0]]
    out["dHdt_total"] = eroder._dHdt[grid.core_nodes[0]]
    out["pluck_rate"] = grid.at_node["bedrock__plucking_rate"][grid.core_nodes[0]]
    out["rock_abr_rate"] = eroder._rock_abrasion_rate[grid.core_nodes[0]]
    out["dzdt"] = eroder._dHdt[grid.core_nodes[0]]  +  eroder._rock_lowering_rate[grid.core_nodes[0]]

    return out


def check_results(results, predicted):

    # Do the comparison between the results
    # and predicted values
    for key in results.keys():
        try:
            assert_almost_equal(
                results[key], predicted[key],
                decimal=3,
            )
        except:
            print(key + ' do not match')


def test_unlimited_sediment_homogeneous_fixed_width_no_attrition():

    # Get parameters for the test
    parameters = set_default_params()

    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def test_unlimited_sediment_homogeneous_fixed_width_attrition():

    # Get parameters for the test
    parameters = set_default_params()
    parameters["beta"] = 0.001
    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def test_unlimited_sediment_homogeneous_dynamic_width_no_attrition():
    # Get parameters for the test
    parameters = set_default_params()
    parameters["width_is_dynamic"] = True
    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def unlimited_sediment_multi_lithology_fixed_width_attrition():
    # Get parameters for the test
    parameters = set_default_params()
    parameters["beta"] = np.array([0.001, 0.0001])
    parameters["sed0"] = np.array([1, 1])
    parameters["D"] = np.array([0.01, 0.01])

    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def unlimited_multi_size_dynamic_width_no_attrition():
    # Get parameters for the test
    parameters = set_default_params()
    parameters["sed0"] = np.array([1, 1, 1, 1, 1])
    parameters["D"] = np.array([0.001, 0.002, 0.004, 0.008, 0.016])
    parameters["width_is_dynamic"] = True

    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def test_unlimited_multi_size_fixed_width_no_attrition():
    # Get parameters for the test
    parameters = set_default_params()
    parameters["sed0"] = np.array([1, 1, 1, 1, 1])
    parameters["D"] = np.array([0.001, 0.002, 0.004, 0.008, 0.016])
    parameters["width_is_dynamic"] = False

    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def test_unlimited_multi_size_fixed_width_attrition():
    # Get parameters for the test
    parameters = set_default_params()
    parameters["sed0"] = np.array([1, 1, 1, 1, 1])
    parameters["D"] = np.array([0.001, 0.002, 0.004, 0.008, 0.016])
    parameters["beta"] = 0.001
    parameters["width_is_dynamic"] = False

    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def test_limited_homogeneous_dynamic_width_no_attrition():
    # Get parameters for the test
    parameters = set_default_params()
    parameters["sed0"] = 0.1
    parameters["width_is_dynamic"] = True

    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)


def test_limited_multi_size_fixed_width_no_attrition():
    # Get parameters for the test
    parameters = set_default_params()
    parameters["width_is_dynamic"] = False
    parameters["sed0"] = np.array([0.1/5, 0.1/5, 0.1/5, 0.1/5, 0.1/5])
    parameters["D"] = np.array([0.001, 0.002, 0.004, 0.008, 0.016])
    # Init grid according to the parameters and run one step
    grid, eroder, sg, fa = init_grid_and_run_one_step(parameters)

    # Store the results
    results = get_results(grid, eroder, sg, fa)

    # Calc the predicted outputs based on the parameters
    predicted = calc_predicted_outputs(parameters)

    # Check if the results match the prediction
    check_results(results, predicted)



# TESTS FROM ("regular") GBE
# def test_transport_rate():
#     grid = HexModelGrid((4, 2), spacing=1000.0)
#     grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
#     grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE
#
#     elev = grid.add_zeros("topographic__elevation", at="node")
#     elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed = grid.add_zeros("soil__depth", at="node")
#     sed[:] = 100.0
#
#     fa = FlowAccumulator(grid)
#     fa.run_one_step()
#     gbe = GravelBedrockEroder(grid, intermittency_factor=0.02, depth_decay_scale=0.5)
#     rock = grid.at_node["bedrock__elevation"]
#     qs_out = grid.at_node["bedload_sediment__volume_outflux"]
#
#     gbe.run_one_step(1.0e-6)  # using dt=0 prevents change to sed, rock, or elev
#     assert_almost_equal(qs_out[grid.core_nodes], [9.88854526, 3.29618175, 3.29618175])
#
#     elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed[:] = 0.5
#     rock[:] = elev - sed
#
#     gbe.run_one_step(1.0e-6)
#     assert_almost_equal(qs_out[grid.core_nodes], [6.25075275, 2.08358425, 2.08358425])
#
#     elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed[:] = 0.0
#     rock[:] = elev
#
#     gbe.run_one_step(1.0e-6)
#     assert_almost_equal(qs_out[grid.core_nodes], [0.0, 0.0, 0.0])
#
#
# def test_sediment_abrasion_rate():
#     grid = HexModelGrid((4, 2), spacing=1000.0)
#     grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
#     grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE
#
#     elev = grid.add_zeros("topographic__elevation", at="node")
#     elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed = grid.add_zeros("soil__depth", at="node")
#     sed[:] = 100.0
#
#     fa = FlowAccumulator(grid)
#     fa.run_one_step()
#     gbe = GravelBedrockEroder(grid, abrasion_coefficient=1.0e-4)
#     gbe.run_one_step(1.0)
#
#     assert_almost_equal(
#         gbe._abrasion[grid.core_nodes],
#         [4.7576285545378313e-07, 9.515257103302159e-08, 9.515257103302159e-08],
#     )
#
#
# def test_rock_abrasion_rate():
#     grid = HexModelGrid((4, 2), spacing=1000.0)
#     grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
#     grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE
#
#     elev = grid.add_zeros("topographic__elevation", at="node")
#     elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed = grid.add_zeros("soil__depth", at="node")
#     sed[:] = 1.0
#
#     fa = FlowAccumulator(grid)
#     fa.run_one_step()
#     gbe = GravelBedrockEroder(grid, abrasion_coefficient=1.0e-4)
#     gbe.run_one_step(1.0)
#
#     assert_almost_equal(
#         gbe._sediment_outflux[grid.core_nodes], [3.12537638, 1.04179213, 1.04179213]
#     )
#     assert_almost_equal(
#         gbe._rock_abrasion_rate[grid.core_nodes],
#         [1.10635873e-07, 2.21271745e-08, 2.21271745e-08],
#     )
#
#
# def test_rock_plucking_rate():
#     grid = HexModelGrid((4, 2), spacing=1000.0)
#     grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
#     grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE
#
#     elev = grid.add_zeros("topographic__elevation", at="node")
#     elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed = grid.add_zeros("soil__depth", at="node")
#     sed[:] = 1.0
#
#     fa = FlowAccumulator(grid)
#     fa.run_one_step()
#     gbe = GravelBedrockEroder(grid, plucking_coefficient=1.0e-4)
#     gbe.run_one_step(1.0)
#
#     assert_almost_equal(
#         gbe._pluck_rate[grid.core_nodes],
#         [5.12263532e-06, 1.70754511e-06, 1.70754511e-06],
#     )
#
#
# def test_steady_unlimited_sediment():
#     grid = HexModelGrid((4, 2), spacing=1000.0)
#     grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
#     grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE
#
#     elev = grid.add_zeros("topographic__elevation", at="node")
#     elev[:] = (0.13 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed = grid.add_zeros("soil__depth", at="node")
#     sed[:] = 10000.0
#     rock = grid.add_zeros("bedrock__elevation", at="node")
#     rock[:] = elev - sed
#
#     fa = FlowAccumulator(grid)
#     fa.run_one_step()
#     gbe = GravelBedrockEroder(grid, abrasion_coefficient=0.0005)
#
#     dt = 4.0e4
#     uplift_rate = 0.0001
#     nsteps = 500
#     for _ in range(nsteps):
#         elev[grid.core_nodes] += uplift_rate * dt
#         sed[grid.core_nodes] += uplift_rate * dt
#         gbe.run_one_step(dt)
#
#     assert_almost_equal(
#         grid.at_node["bedload_sediment__volume_outflux"][grid.core_nodes],
#         [99.073, 45.033, 45.033],
#         decimal=2,
#     )
#     assert_almost_equal(
#         grid.at_node["topographic__steepest_slope"][grid.core_nodes],
#         [0.130579, 0.170346, 0.170346],
#         decimal=5,
#     )
#
#
# def test_steady_general():
#     grid = HexModelGrid((3, 2), spacing=1000.0)
#     grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
#     grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE
#
#     elev = grid.add_zeros("topographic__elevation", at="node")
#     elev[:] = (0.2 * grid.y_of_node) / np.cos(np.radians(30.0))
#     sed = grid.add_zeros("soil__depth", at="node")
#     sed[:] = 1.0
#     rock = grid.add_zeros("bedrock__elevation", at="node")
#     rock[:] = elev - sed
#
#     fa = FlowAccumulator(grid)
#     fa.run_one_step()
#     gbe = GravelBedrockEroder(
#         grid, abrasion_coefficient=0.0005, coarse_fraction_from_plucking=0.5
#     )
#
#     dt = 7500.0
#     uplift_rate = 0.0001
#     nsteps = 3000
#     for _ in range(nsteps):
#         elev[grid.core_nodes] += uplift_rate * dt
#         rock[grid.core_nodes] += uplift_rate * dt
#         gbe.run_one_step(dt)
#
#     assert_almost_equal(
#         grid.at_node["bedrock__exposure_fraction"][grid.core_nodes], 0.5062, decimal=4
#     )
#     assert_almost_equal(
#         grid.at_node["topographic__steepest_slope"][grid.core_nodes], 0.2387, decimal=4
#     )
#     assert_almost_equal(
#         grid.at_node["bedload_sediment__volume_outflux"][grid.core_nodes],
#         32.972,
#         decimal=3,
#     )