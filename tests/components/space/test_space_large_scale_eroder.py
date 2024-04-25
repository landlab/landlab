import copy as cp

import numpy as np
import pytest
from numpy import testing

from landlab import FieldError
from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import PriorityFloodFlowRouter
from landlab.components import SpaceLargeScaleEroder

try:
    PriorityFloodFlowRouter.load_richdem()
except ModuleNotFoundError:
    with_richdem = False
else:
    with_richdem = True


def test_inputFields_flowRouter():
    """
    SpaceLargeScaleEroder should throw an error when topograhy is not equal to the sum of
    bedrock and soil thickness
    """
    # %%
    # Make a raster model grid and create a plateau
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    br[:] = mg.x_of_node + mg.y_of_node
    soil = mg.add_zeros("soil__depth", at="node")
    z[:] = br + soil
    fa = FlowAccumulator(mg, flow_director="D8")
    fa.run_one_step()

    # make plateau at 10m
    br += 10

    # Instanciate the slider
    with pytest.raises(AssertionError):
        _ = SpaceLargeScaleEroder(mg)


# %%
def test_inputFields_soil():
    """
    SpaceLargeScaleEroder should throw an error when the soil__depth field is not provided
    """
    # %%
    mg = RasterModelGrid((5, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("bedrock__elevation", at="node")
    fa = FlowAccumulator(mg, flow_director="D8")
    fa.run_one_step()

    # Instanciate the slider
    with pytest.raises(FieldError):
        _ = SpaceLargeScaleEroder(mg)


# %%
def test_inputFields_bedrock():
    """
    SpaceLargeScaleEroder should instanciate the bedrock__elevation field
    when it is not provided
    """
    # %%
    mg = RasterModelGrid((5, 5))
    _ = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    fa = FlowAccumulator(mg, flow_director="D8")
    fa.run_one_step()
    _ = SpaceLargeScaleEroder(mg)

    assert "bedrock__elevation" in mg.at_node


# %%
def test_properties_phi_fraction_fines_LS():
    """
    SpaceLargeScaleEroder should throw an error when phi/fraction_fines_LS < 0 or phi > 0
    """
    # %%
    mg = RasterModelGrid((5, 5))
    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    br[:] = mg.x_of_node + mg.y_of_node
    soil = mg.add_zeros("soil__depth", at="node")
    z[:] = br + soil
    fa = FlowAccumulator(mg, flow_director="D8")
    fa.run_one_step()

    # Instanciate the slider
    with pytest.raises(ValueError):
        _ = SpaceLargeScaleEroder(mg, phi=-0.2)
    # Instanciate the slider
    with pytest.raises(ValueError):
        _ = SpaceLargeScaleEroder(mg, phi=1.2)
    # Instanciate the slider
    with pytest.raises(ValueError):
        _ = SpaceLargeScaleEroder(mg, F_f=-0.2)
    # Instanciate the slider
    with pytest.raises(ValueError):
        _ = SpaceLargeScaleEroder(mg, F_f=1.2)


# %%


def test_route_to_multiple_error_raised():
    # %%
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    br[:] = mg.x_of_node + mg.y_of_node
    soil = mg.add_zeros("soil__depth", at="node")
    z[:] = br + soil
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        SpaceLargeScaleEroder(
            mg,
            K_sed=0.1,
            K_br=0.1,
            F_f=0.5,
            phi=0.1,
            H_star=1.0,
            v_s=0.001,
            m_sp=1.0,
            n_sp=0.5,
            sp_crit_sed=0,
            sp_crit_br=0,
        )


# %%
def test_soil_field_already_on_grid():
    # %%
    """
    Test that an existing soil grid field is not changed by instantiating
    SpaceLargeScaleEroder.
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")
    soil += 1.0  # add 1m of soil everywehre

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 10000 + mg.node_x / 10000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    br[:] = z - soil

    # Create a D8 flow handler
    FlowAccumulator(mg, flow_director="D8")

    # Instantiate SpaceLargeScaleEroder
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=0.01,
        K_br=0.01,
        F_f=0.0,
        phi=0.0,
        v_s=0.001,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ensure that 'soil__depth' field is everywhere equal to 1.0 m.
    testing.assert_array_equal(
        np.ones(mg.number_of_nodes),
        sp._soil__depth,
        err_msg="SpaceLargeScaleEroder soil depth field test failed",
        verbose=True,
    )

    # %% Check getters
    assert sp.K_br == pytest.approx(0.01)
    assert sp.K_sed == pytest.approx(0.01)
    assert sp.fraction_fines == pytest.approx(0.0)
    assert sp.sediment_porosity == pytest.approx(0.0)
    assert sp.settling_velocity == pytest.approx(0.001)
    assert sp.drainage_area_exp == pytest.approx(0.5)
    assert sp.slope_exp == pytest.approx(1.0)

    # sediment erosion is zero before running the component
    testing.assert_array_equal(np.zeros(mg.number_of_nodes), sp.Es)
    # rock erosion is zero before running the component
    testing.assert_array_equal(np.zeros(mg.number_of_nodes), sp.Er)

    # %% Check setters
    sp.K_br = 0.02
    assert sp.K_br == pytest.approx(0.02)
    sp.K_sed = 0.02
    assert sp.K_sed == pytest.approx(0.02)

    with pytest.raises(AttributeError):
        sp.Es = np.zeros(mg.number_of_nodes)

    with pytest.raises(AttributeError):
        sp.Er = np.zeros(mg.number_of_nodes)


# %%


def test_br_field_already_on_grid():
    # %%
    """
    Test that an existing bedrock elevation grid field is not changed by
    instantiating SpaceLargeScaleEroder.
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    br += 1.0  # make bedrock elevation 5m below surface
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 10000 + mg.node_x / 10000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    z[:] = br + soil

    # Create a D8 flow handler
    FlowAccumulator(mg, flow_director="D8")

    # Instantiate SpaceLargeScaleEroder
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=0.01,
        K_br=0.01,
        F_f=0.0,
        phi=0.0,
        v_s=0.001,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ensure that 'bedrock__elevation' field is everywhere equal to 1.0 m.
    testing.assert_array_equal(
        np.ones(mg.number_of_nodes),
        sp._bedrock__elevation,
        err_msg="SpaceLargeScaleEroder bedrock field test failed",
        verbose=True,
    )


@pytest.mark.parametrize("k_sed", (0.00001, [0.00001] * 25))
def test_matches_detachment_solution(k_sed):
    # %%
    """
    Test that model matches the detachment-limited analytical solution
    for slope/area relationship at steady state: S=(U/K_br)^(1/n)*A^(-m/n).
    """

    # %% set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 10000 + mg.node_x / 10000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    br[:] = z[:] - soil[:]

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director="D8")

    # Parameter values for detachment-limited test
    K_br = 0.01
    U = 0.0001
    dt = 1.0
    F_f = 1.0  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=k_sed,
        K_br=K_br,
        F_f=F_f,
        phi=0.1,
        H_star=1.0,
        v_s=0.001,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (2000x1-year timesteps).
    for _ in range(2000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt  # m
        br[mg.core_nodes] = z[mg.core_nodes] - soil[mg.core_nodes]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(U / K_br, 1.0 / n_sp) * np.power(
        mg.at_node["drainage_area"][mg.core_nodes], -m_sp / n_sp
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder detachment-limited test failed",
        verbose=True,
    )


# %%
def test_matches_detachment_solution_n_gr_1():
    # %%
    """
    Test that model matches the detachment-limited analytical solution
    for slope/area relationship at steady state: S=(U/K_br)^(1/n)*A^(-m/n).
    """

    # %% set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 10000 + mg.node_x / 10000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    br[:] = z[:] - soil[:]

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director="D8")

    # Parameter values for detachment-limited test
    K_br = 0.01
    U = 0.0001
    dt = 1.0
    F_f = 1.0  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.1

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=0.00001,
        K_br=K_br,
        F_f=F_f,
        phi=0.1,
        H_star=1.0,
        v_s=0.001,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (2000x1-year timesteps).
    for _ in range(4000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt  # m
        br[mg.core_nodes] = z[mg.core_nodes] - soil[mg.core_nodes]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(U / K_br, 1.0 / n_sp) * np.power(
        mg.at_node["drainage_area"][mg.core_nodes], -m_sp / n_sp
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder detachment-limited test failed",
        verbose=True,
    )


@pytest.mark.slow
def test_matches_transport_solution():
    # %%
    """
    Test that model matches the transport-limited analytical solution
    for slope/area relationship at steady state: S=((U * v_s) / (K_sed * A^m)
    + U / (K_sed * A^m))^(1/n).

    Also test that model matches the analytical solution for steady-state
    sediment flux: Qs = U * A * (1 - phi).
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    soil[:] += 100.0  # initial soil depth of 100 m
    br[:] = z[:]
    z[:] += soil[:]

    # Create a D8 flow handler
    fa = FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    # Parameter values for detachment-limited test
    K_sed = 0.01
    U = 0.0001
    dt = 1.0
    F_f = 1.0  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0
    v_s = 0.5
    phi = 0.5

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=K_sed,
        K_br=0.01,
        F_f=F_f,
        phi=phi,
        H_star=1.0,
        v_s=v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (5000x1-year timesteps).
    for _ in range(5000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        br[mg.core_nodes] += U * dt  # m
        soil[0] = (
            100.0  # enforce constant soil depth at boundary to keep lowering steady
        )
        z[:] = br[:] + soil[:]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(
        (
            (U * v_s * (1 - phi))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        )
        + (
            (U * (1 - phi))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        ),
        1.0 / n_sp,
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder transport-limited slope-area test failed",
        verbose=True,
    )

    # compare numerical and analytical sediment flux solutions
    num_sedflux = mg.at_node["sediment__outflux"][mg.core_nodes]
    analytical_sedflux = U * mg.at_node["drainage_area"][mg.core_nodes] * (1 - phi)

    # test for match with anakytical sediment flux
    testing.assert_array_almost_equal(
        num_sedflux,
        analytical_sedflux,
        decimal=8,
        err_msg="SpaceLargeScaleEroder transport-limited sediment flux test failed",
        verbose=True,
    )
    # %%


@pytest.mark.slow
def test_matches_bedrock_alluvial_solution():
    # %%
    """
    Test that model matches the bedrock-alluvial analytical solution
    for slope/area relationship at steady state:
    S=((U * v_s * (1 - F_f)) / (K_sed * A^m) + U / (K_br * A^m))^(1/n).

    Also test that the soil depth everywhere matches the bedrock-alluvial
    analytical solution at steady state:
    H = -H_star * ln(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s))).
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    soil[:] += 0.0  # initial condition of no soil depth.
    br[:] = z[:]
    z[:] += soil[:]

    # Create a D8 flow handler
    fa = FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    # Parameter values for detachment-limited test
    K_br = 0.002
    K_sed = 0.002
    U = 0.0001
    dt = 10.0
    F_f = 0.2  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0
    v_s = 0.25
    H_star = 0.1

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=K_sed,
        K_br=K_br,
        F_f=F_f,
        phi=0.0,
        H_star=H_star,
        v_s=v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (10000x1-year timesteps).
    for _ in range(10000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        br[mg.core_nodes] += U * dt  # m
        soil[0] = 0.0  # enforce 0 soil depth at boundary to keep lowering steady
        z[:] = br[:] + soil[:]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(
        (
            (U * v_s * (1 - F_f))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        )
        + (U / (K_br * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))),
        1.0 / n_sp,
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial slope-area test failed",
        verbose=True,
    )

    # compare numerical and analytical soil depth solutions
    num_h = mg.at_node["soil__depth"][mg.core_nodes]
    analytical_h = -H_star * np.log(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s)))

    # test for match with analytical sediment depth
    testing.assert_array_almost_equal(
        num_h,
        analytical_h,
        decimal=5,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial soil thickness test failed",
        verbose=True,
    )
    # %%


def test_can_run_with_hex():
    """Test that model can run with hex model grid."""
    # %%
    # Set up a 5x5 grid with open boundaries and low initial elevations.
    mg = HexModelGrid((7, 7))
    z = mg.add_zeros("topographic__elevation", at="node")
    _ = mg.add_zeros("soil__depth", at="node")
    z[:] = 0.01 * mg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director="FlowDirectorSteepest")

    # Parameter values for test 1
    U = 0.001
    dt = 10.0

    # Create the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=0.00001,
        K_br=0.00000000001,
        F_f=0.5,
        phi=0.1,
        H_star=1.0,
        v_s=0.001,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state.
    for _ in range(2000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt


@pytest.mark.skipif(not with_richdem, reason="richdem is not installed")
def test_matches_detachment_solution_PF():
    # %%
    """
    Test that model matches the detachment-limited analytical solution
    for slope/area relationship at steady state: S=(U/K_br)^(1/n)*A^(-m/n).
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 10000 + mg.node_x / 10000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    br[:] = z[:] - soil[:]

    fa = PriorityFloodFlowRouter(
        mg, surface="topographic__elevation", flow_metric="D8", suppress_out=True
    )
    fa.run_one_step()

    # Parameter values for detachment-limited test
    K_br = 0.01
    U = 0.0001
    dt = 1.0
    F_f = 1.0  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=0.00001,
        K_br=K_br,
        F_f=F_f,
        phi=0.1,
        H_star=1.0,
        v_s=0.001,
        v_s_lake=0.001,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (2000x1-year timesteps).
    for _ in range(2000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt  # m
        br[mg.core_nodes] = z[mg.core_nodes] - soil[mg.core_nodes]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(U / K_br, 1.0 / n_sp) * np.power(
        mg.at_node["drainage_area"][mg.core_nodes], -m_sp / n_sp
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder detachment-limited test failed",
        verbose=True,
    )


# %%
@pytest.mark.slow
@pytest.mark.skipif(not with_richdem, reason="richdem is not installed")
def test_matches_transport_solution_PF():
    """
    Test that model matches the transport-limited analytical solution
    for slope/area relationship at steady state: S=((U * v_s) / (K_sed * A^m)
    + U / (K_sed * A^m))^(1/n).

    Also test that model matches the analytical solution for steady-state
    sediment flux: Qs = U * A * (1 - phi).
    """
    # %%
    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    soil[:] += 100.0  # initial soil depth of 100 m
    br[:] = z[:]
    z[:] += soil[:]

    # Create a D8 flow handler
    fa = PriorityFloodFlowRouter(
        mg, surface="topographic__elevation", flow_metric="D8", suppress_out=True
    )
    fa.run_one_step()

    # Parameter values for detachment-limited test
    K_sed = 0.01
    U = 0.0001
    dt = 1.0
    F_f = 1.0  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0
    v_s = 0.5
    phi = 0.5

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=K_sed,
        K_br=0.01,
        F_f=F_f,
        phi=phi,
        H_star=1.0,
        v_s=v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (5000x1-year timesteps).
    for _ in range(5000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        br[mg.core_nodes] += U * dt  # m
        soil[0] = (
            100.0  # enforce constant soil depth at boundary to keep lowering steady
        )
        z[:] = br[:] + soil[:]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(
        (
            (U * v_s * (1 - phi))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        )
        + (
            (U * (1 - phi))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        ),
        1.0 / n_sp,
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder transport-limited slope-area test failed",
        verbose=True,
    )

    # compare numerical and analytical sediment flux solutions
    num_sedflux = mg.at_node["sediment__outflux"][mg.core_nodes]
    analytical_sedflux = U * mg.at_node["drainage_area"][mg.core_nodes] * (1 - phi)

    # test for match with anakytical sediment flux
    testing.assert_array_almost_equal(
        num_sedflux,
        analytical_sedflux,
        decimal=8,
        err_msg="SpaceLargeScaleEroder transport-limited sediment flux test failed",
        verbose=True,
    )


# %%
@pytest.mark.slow
@pytest.mark.skipif(not with_richdem, reason="richdem is not installed")
def test_matches_bedrock_alluvial_solution_PF():
    """
    Test that model matches the bedrock-alluvial analytical solution
    for slope/area relationship at steady state:
    S=((U * v_s * (1 - F_f)) / (K_sed * A^m) + U / (K_br * A^m))^(1/n).

    Also test that the soil depth everywhere matches the bedrock-alluvial
    analytical solution at steady state:
    H = -H_star * ln(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s))).
    """
    # %%
    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    soil[:] += 0.0  # initial condition of no soil depth.
    br[:] = z[:]
    z[:] += soil[:]

    # Create a D8 flow handler
    fa = PriorityFloodFlowRouter(
        mg, surface="topographic__elevation", flow_metric="D8", suppress_out=True
    )
    fa.run_one_step()

    # Parameter values for detachment-limited test
    K_br = 0.002
    K_sed = 0.002
    U = 0.0001
    dt = 10.0
    F_f = 0.2  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0
    v_s = 0.25
    H_star = 0.1

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=K_sed,
        K_br=K_br,
        F_f=F_f,
        phi=0.0,
        H_star=H_star,
        v_s=v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (10000x1-year timesteps).
    for _ in range(10000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        br[mg.core_nodes] += U * dt  # m
        soil[0] = 0.0  # enforce 0 soil depth at boundary to keep lowering steady
        z[:] = br[:] + soil[:]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(
        (
            (U * v_s * (1 - F_f))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        )
        + (U / (K_br * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))),
        1.0 / n_sp,
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial slope-area test failed",
        verbose=True,
    )

    # compare numerical and analytical soil depth solutions
    num_h = mg.at_node["soil__depth"][mg.core_nodes]
    analytical_h = -H_star * np.log(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s)))

    # test for match with analytical sediment depth
    testing.assert_array_almost_equal(
        num_h,
        analytical_h,
        decimal=5,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial soil thickness test failed",
        verbose=True,
    )
    # %%


# %%
@pytest.mark.slow
@pytest.mark.skipif(not with_richdem, reason="richdem is not installed")
def test_matches_bedrock_alluvial_solution_PF_extended_range():
    """
    Test that model matches the bedrock-alluvial analytical solution
    for slope/area relationship at steady state:
    S=((U * v_s * (1 - F_f)) / (K_sed * A^m) + U / (K_br * A^m))^(1/n).

    Also test that the soil depth everywhere matches the bedrock-alluvial
    analytical solution at steady state:
    H = -H_star * ln(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s))).
    """
    # %%
    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    np.random.seed(seed=5000)
    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )

    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    soil[:] += 0.0  # initial condition of no soil depth.
    br[:] = z[:]
    z[:] += soil[:]

    # Create a D8 flow handler
    fa = PriorityFloodFlowRouter(
        mg, surface="topographic__elevation", flow_metric="D8", suppress_out=True
    )
    fa.run_one_step()

    # Parameter values for detachment-limited test
    K_br = 0.005
    K_sed = 0.01
    U = 0.0001
    dt = 10.0
    F_f = 0.0
    m_sp = 0.5
    n_sp = 1.0
    v_s = 1.0
    H_star = 1.0

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=K_sed,
        K_br=K_br,
        F_f=F_f,
        phi=0.0,
        H_star=H_star,
        v_s=v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (10000x1-year timesteps).
    for _ in range(10000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        br[mg.core_nodes] += U * dt  # m
        soil[0] = 0.0  # enforce 0 soil depth at boundary to keep lowering steady
        z[:] = br[:] + soil[:]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(
        (
            (U * v_s * (1 - F_f))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        )
        + (U / (K_br * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))),
        1.0 / n_sp,
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial slope-area test failed",
        verbose=True,
    )

    # compare numerical and analytical soil depth solutions
    num_h = mg.at_node["soil__depth"][mg.core_nodes]
    analytical_h = -H_star * np.log(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s)))

    # test for match with analytical sediment depth
    testing.assert_array_almost_equal(
        num_h,
        analytical_h,
        decimal=5,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial soil thickness test failed",
        verbose=True,
    )


# %%
@pytest.mark.slow
@pytest.mark.skipif(not with_richdem, reason="richdem is not installed")
def test_matches_bedrock_alluvial_solution_PF_high_v_high_hstar():
    """
    Test that model matches the bedrock-alluvial analytical solution
    for slope/area relationship at steady state:
    S=((U * v_s * (1 - F_f)) / (K_sed * A^m) + U / (K_br * A^m))^(1/n).
    Also test that the soil depth everywhere matches the bedrock-alluvial
    analytical solution at steady state:
    H = -H_star * ln(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s))).
    Also test that the sediment flux at the outlet is:
    Qs = U * A.
    """
    # %%
    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    np.random.seed(seed=5000)
    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )

    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    soil[:] += 0.0  # initial condition of no soil depth.
    br[:] = z[:]
    z[:] += soil[:]

    # Create a D8 flow handler
    fa = PriorityFloodFlowRouter(
        mg, surface="topographic__elevation", flow_metric="D8", suppress_out=True
    )
    fa.run_one_step()

    # Parameter values for detachment-limited test
    K_br = 0.005
    K_sed = 0.01
    U = 0.0001
    dt = 10.0
    F_f = 0.0
    m_sp = 0.5
    n_sp = 1.0
    v_s = 10.0
    H_star = 2.0
    phi = 0.0

    # Instantiate the SpaceLargeScaleEroder component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=K_sed,
        K_br=K_br,
        F_f=F_f,
        phi=phi,
        H_star=H_star,
        v_s=v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )

    # ... and run it to steady state (25000x10-year timesteps).
    for _ in range(25000):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        br[mg.core_nodes] += U * dt  # m
        soil[0] = 0.0  # enforce 0 soil depth at boundary to keep lowering steady
        z[:] = br[:] + soil[:]

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(
        (
            (U * v_s * (1 - F_f))
            / (K_sed * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))
        )
        + (U / (K_br * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))),
        1.0 / n_sp,
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial slope-area test failed",
        verbose=True,
    )

    # compare numerical and analytical soil depth solutions
    num_h = mg.at_node["soil__depth"][mg.core_nodes]
    analytical_h = -H_star * np.log(1 - (v_s / (K_sed / (K_br * (1 - F_f)) + v_s)))

    # test for match with analytical sediment depth
    testing.assert_array_almost_equal(
        num_h,
        analytical_h,
        decimal=5,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial soil thickness test failed",
        verbose=True,
    )

    # compare numerical and analytical sediment flux solutions
    num_sedflux = mg.at_node["sediment__outflux"][mg.core_nodes]
    analytical_sedflux = U * mg.at_node["drainage_area"][mg.core_nodes] * (1 - phi)

    # test for match with analytical sediment flux
    testing.assert_array_almost_equal(
        num_sedflux,
        analytical_sedflux,
        decimal=5,
        err_msg="SpaceLargeScaleEroder bedrock-alluvial sediment flux test failed",
        verbose=True,
    )


# %%
@pytest.mark.slow
def test_MassBalance():
    # %%
    # set up a 15x15 grid with one open outlet node and low initial elevations.
    nr = 15
    nc = 15
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")
    br = mg.add_zeros("bedrock__elevation", at="node")
    soil = mg.add_zeros("soil__depth", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )
    soil[:] += 0.0  # initial condition of no soil depth.
    br[:] = z[:]
    z[:] += soil[:]

    # Create a D8 flow handler
    fa = FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    # Parameter values for detachment-limited test
    K_br = 0.002
    K_sed = 0.002
    U = 0.0001
    dt = 10.0
    F_f = 0.2  # all detached rock disappears; detachment-ltd end-member
    m_sp = 0.5
    n_sp = 1.0
    v_s = 0.25
    H_star = 0.1

    # Instantiate the Space component...
    sp = SpaceLargeScaleEroder(
        mg,
        K_sed=K_sed,
        K_br=K_br,
        F_f=F_f,
        phi=0.0,
        H_star=H_star,
        v_s=v_s,
        m_sp=m_sp,
        n_sp=n_sp,
        sp_crit_sed=0,
        sp_crit_br=0,
    )
    # Get values before run
    z = mg.at_node["topographic__elevation"]
    br = mg.at_node["bedrock__elevation"]
    H = mg.at_node["soil__depth"]
    cores = mg.core_nodes
    area = mg.cell_area_at_node

    # ... and run it to steady state (10000x1-year timesteps).
    for _ in range(10000):
        fa.run_one_step()
        soil_B = cp.deepcopy(H)
        bed_B = cp.deepcopy(br)
        vol_SSY_riv, V_leaving_riv = sp.run_one_step(dt=dt)
        diff_MB = (
            np.sum((bed_B[cores] - br[cores]) * area[cores])
            + np.sum((soil_B[cores] - H[cores]) * area[cores]) * (1 - sp._phi)
            - vol_SSY_riv * dt
            - V_leaving_riv
        )

        br[mg.core_nodes] += U * dt  # m
        soil[0] = 0.0  # enforce 0 soil depth at boundary to keep lowering steady
        z[:] = br[:] + soil[:]

        # Test Every iteration
        testing.assert_array_almost_equal(
            z[cores],
            br[cores] + H[cores],
            decimal=5,
            err_msg="Topography does not equal sum of bedrock and soil! Decrease timestep",
            verbose=True,
        )
        testing.assert_array_less(
            abs(diff_MB),
            1e-8 * mg.number_of_nodes,
            err_msg=(
                "Mass balance error SpaceLargeScaleEroder! Try to resolve by "
                "becreasing timestep"
            ),
            verbose=True,
        )
