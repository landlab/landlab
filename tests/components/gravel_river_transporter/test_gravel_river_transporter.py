#! /usr/bin/env python
"""
Unit tests for landlab.components.gravel_river_transporter.gravel_river_transporter
"""
from numpy.testing import assert_allclose
from numpy.testing import assert_equal
from numpy.testing import assert_raises

from landlab import HexModelGrid
from landlab import RadialModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import GravelRiverTransporter


def test_analytical_solution_one_cell_basic_solver():
    # parameters
    dx = 1000.0  # node spacing, m
    uplift_rate = 0.0001  # rate of relative uplift, m/y
    bankfull_runoff_rate = 10.0  # bankfull runoff rate, m/y
    abrasion_coef = 1.0 / 2000.0  # abrasion coefficient, 1/m
    trans_coef = 0.041  # transport coefficient, -
    intermittency_fac = 0.01  # intermittency factor, -
    porosity = 0.35
    dt = 10000.0  # time step duration, y

    grid = RasterModelGrid((3, 3), xy_spacing=dx)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")

    fa = FlowAccumulator(grid, runoff_rate=bankfull_runoff_rate)
    fa.run_one_step()
    transporter = GravelRiverTransporter(
        grid,
        intermittency_factor=intermittency_fac,
        transport_coefficient=trans_coef,
        abrasion_coefficient=abrasion_coef,
        sediment_porosity=porosity,
        solver="explicit",
    )

    for _ in range(200):
        fa.run_one_step()
        elev[grid.core_nodes] += uplift_rate * dt
        transporter.run_one_step(dt)

    # Analytical solutions for sediment flux, slope, and abrasion rate
    Qs_pred = (uplift_rate * dx * dx * (1.0 - porosity)) / (abrasion_coef * dx / 2 + 1)
    S_pred = (
        Qs_pred / (trans_coef * intermittency_fac * bankfull_runoff_rate * dx * dx)
    ) ** (6.0 / 7.0)
    Ea_pred = abrasion_coef * (Qs_pred / 2.0) / dx  # factor of 2 is cell midpoint

    assert_allclose(transporter._sediment_outflux[4], Qs_pred, rtol=1.0e-4)
    assert_allclose(transporter._slope[4], S_pred, rtol=1.0e-4)
    assert_allclose(transporter._abrasion[4], Ea_pred, rtol=1.0e-4)


def test_analytical_solution_one_cell_matrix_solver():
    # parameters
    dx = 1000.0  # node spacing, m
    uplift_rate = 0.0001  # rate of relative uplift, m/y
    bankfull_runoff_rate = 10.0  # bankfull runoff rate, m/y
    abrasion_coef = 1.0 / 2000.0  # abrasion coefficient, 1/m
    trans_coef = 0.041  # transport coefficient, -
    intermittency_fac = 0.01  # intermittency factor, -
    porosity = 0.35
    dt = 10000.0  # time step duration, y

    grid = RasterModelGrid((3, 3), xy_spacing=dx)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")

    fa = FlowAccumulator(grid, runoff_rate=bankfull_runoff_rate)
    fa.run_one_step()
    transporter = GravelRiverTransporter(
        grid,
        intermittency_factor=intermittency_fac,
        transport_coefficient=trans_coef,
        abrasion_coefficient=abrasion_coef,
        sediment_porosity=porosity,
        solver="matrix",
    )

    for _ in range(200):
        fa.run_one_step()
        elev[grid.core_nodes] += uplift_rate * dt
        transporter.run_one_step(dt)

    # Analytical solution for slope (matrix solver does not compute sediment
    # flux or abrasion rate)
    Qs_pred = (uplift_rate * dx * dx * (1.0 - porosity)) / (abrasion_coef * dx / 2 + 1)
    S_pred = (
        Qs_pred / (trans_coef * intermittency_fac * bankfull_runoff_rate * dx * dx)
    ) ** (6.0 / 7.0)
    assert_allclose(transporter._slope[4], S_pred, rtol=1.0e-4)


def test_analytical_solution_four_cells_basic_solver():
    # parameters
    dx = 1000.0  # node spacing, m
    initial_slope = 0.0001  # starting topographic slope, m/m
    bankfull_runoff_rate = 10.0  # r, m/y
    trans_coef = 0.041  # this is kQ
    intermittency_fac = 0.01  # this is I
    abrasion_length_scale = 2000.0  # this is the inverse of beta, in m
    dt = 10000.0  # time step, years
    uplift_rate = 0.0001  # U, m/y
    porosity = 0.35  # lambda_p
    nsteps = 400  # number of time steps

    # Derived parameters
    abrasion_coef = 1.0 / abrasion_length_scale  # this is beta, 1/m

    grid = RasterModelGrid((3, 6), xy_spacing=dx)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[11] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = initial_slope * (grid.x_of_node.max() - grid.x_of_node)

    fa = FlowAccumulator(grid, runoff_rate=bankfull_runoff_rate)
    fa.run_one_step()
    transporter = GravelRiverTransporter(
        grid, abrasion_coefficient=abrasion_coef, solver="explicit"
    )

    for _ in range(nsteps):
        fa.run_one_step()
        elev[grid.core_nodes] += uplift_rate * dt
        transporter.run_one_step(dt)

    # Compare with analytical solutions in upper cell
    Qs_pred = (uplift_rate * dx * dx * (1.0 - porosity)) / (abrasion_coef * dx / 2 + 1)
    assert_allclose(transporter._sediment_outflux[7], Qs_pred, rtol=1.0e-3)
    S_pred = (
        Qs_pred / (trans_coef * intermittency_fac * bankfull_runoff_rate * dx * dx)
    ) ** (6.0 / 7.0)
    assert_allclose(transporter._slope[7], S_pred, rtol=1.0e-3)
    Ea_pred = abrasion_coef * (Qs_pred / 2) / dx
    assert_allclose(transporter._abrasion[7], Ea_pred, rtol=1.0e-3)


def test_analytical_solution_four_cells_matrix_solver():
    # parameters
    dx = 1000.0  # node spacing, m
    initial_slope = 0.0001  # starting topographic slope, m/m
    bankfull_runoff_rate = 10.0  # r, m/y
    trans_coef = 0.041  # this is kQ
    intermittency_fac = 0.01  # this is I
    abrasion_length_scale = 2000.0  # this is the inverse of beta, in m
    dt = 10000.0  # time step, years
    uplift_rate = 0.0001  # U, m/y
    porosity = 0.35  # lambda_p
    nsteps = 400  # number of time steps

    # Derived parameters
    abrasion_coef = 1.0 / abrasion_length_scale  # this is beta, 1/m

    grid = RasterModelGrid((3, 6), xy_spacing=dx)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[11] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = initial_slope * (grid.x_of_node.max() - grid.x_of_node)

    fa = FlowAccumulator(grid, runoff_rate=bankfull_runoff_rate)
    fa.run_one_step()
    transporter = GravelRiverTransporter(
        grid, abrasion_coefficient=abrasion_coef, solver="matrix"
    )

    for _ in range(nsteps):
        fa.run_one_step()
        elev[grid.core_nodes] += uplift_rate * dt
        transporter.run_one_step(dt)

    # Compare with analytical solution for slope in upper cell
    Qs_pred = (uplift_rate * dx * dx * (1.0 - porosity)) / (abrasion_coef * dx / 2 + 1)
    S_pred = (
        Qs_pred / (trans_coef * intermittency_fac * bankfull_runoff_rate * dx * dx)
    ) ** (6.0 / 7.0)
    assert_allclose(transporter._slope[7], S_pred, rtol=1.0e-3)


def test_exception_handling():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    FlowAccumulator(grid)
    assert_raises(
        ValueError, GravelRiverTransporter, grid, solver="not one we actually have"
    )


def test_link_length_handling():
    # Hex grid with fixed and uniform flow-link length
    grid = HexModelGrid((5, 3), spacing=2.0)
    grid.add_zeros("topographic__elevation", at="node")
    fa = FlowAccumulator(grid)
    transporter = GravelRiverTransporter(grid)
    assert_allclose(transporter._flow_link_length_over_cell_area, 0.57735, rtol=1.0e-6)

    # Raster grid with possibility of flow on diagonals
    grid = RasterModelGrid((3, 4), xy_spacing=2.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE
    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[5] = 1.0
    elev[6] = 2.0
    fa = FlowAccumulator(grid, flow_director="FlowDirectorD8")
    fa.run_one_step()
    assert_equal(grid.at_node["flow__link_to_receiver_node"][5:7], [17, 8])
    transporter = GravelRiverTransporter(grid, abrasion_coefficient=1.0)
    transporter.run_one_step(1.0)
    assert_equal(transporter._flow_length_is_variable, True)
    assert_equal(transporter._grid_has_diagonals, True)
    assert_allclose(
        transporter._flow_link_length_over_cell_area, [0.707107, 0.5], rtol=1.0e-6
    )

    # Radial grid with variable link length but no diagonals
    grid = RadialModelGrid(n_rings=1, nodes_in_first_ring=5, spacing=2.0)
    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[2] = 1.0
    elev[0] = -1.0
    fa = FlowAccumulator(grid)
    fa.run_one_step()
    transporter = GravelRiverTransporter(grid, abrasion_coefficient=1.0)
    transporter.run_one_step(1.0)
    assert_equal(transporter._flow_length_is_variable, True)
    assert_equal(transporter._grid_has_diagonals, False)
    assert_allclose(transporter._flow_link_length_over_cell_area, [0.5505524])
