#! /usr/bin/env python
"""
Unit tests for landlab.components.gravel_river_transporter.gravel_river_transporter
"""
from numpy.testing import assert_allclose, assert_raises

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, GravelRiverTransporter


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

    # Analytical solutions for sediment flux, slope, and abrasion rate
    Qs_pred = (uplift_rate * dx * dx * (1.0 - porosity)) / (abrasion_coef * dx / 2 + 1)
    S_pred = (
        Qs_pred / (trans_coef * intermittency_fac * bankfull_runoff_rate * dx * dx)
    ) ** (6.0 / 7.0)
    Ea_pred = abrasion_coef * (Qs_pred / 2.0) / dx  # factor of 2 is cell midpoint

    # assert_allclose(transporter._sediment_outflux[4], Qs_pred, rtol=1.0e-4)
    assert_allclose(transporter._slope[4], S_pred, rtol=1.0e-4)
    # assert_allclose(transporter._abrasion[4], Ea_pred, rtol=1.0e-4)


def test_exception_handling():
    grid = RasterModelGrid((3, 3))
    elev = grid.add_zeros("topographic__elevation", at="node")
    FlowAccumulator(grid)
    assert_raises(
        ValueError, GravelRiverTransporter, grid, solver="not one we actually have"
    )
