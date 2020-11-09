#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 17:25:27 2020

@author: dylanward
"""

import numpy as np
import pytest

from landlab import HexModelGrid, RasterModelGrid
from landlab.components import ExponentialWeathererIntegrated


def test_no_soil_depth_field():
    # should pass if the code fails correctly on instantiation
    # given no established soil__depth field on the model grid
    mg = RasterModelGrid((5, 5))
    # omitting: mg.add_zeros("soil__depth", at="node")
    with pytest.raises(KeyError):
        ExponentialWeathererIntegrated(mg)


def test_create_weatherer_and_change_rate():
    # Unit test from original exp weatherer
    grid = RasterModelGrid((3, 3), 1.0)
    grid.add_zeros("soil__depth", at="node")

    ew = ExponentialWeathererIntegrated(grid, soil_production__maximum_rate=0.0001)
    ew.maximum_weathering_rate = 0.0004
    assert ew.maximum_weathering_rate == 0.0004


def test_basic_production_rate_func():
    # This tests drop-in compatiblity with the existing ExponentialWeatherer component
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    expw2 = ExponentialWeathererIntegrated(mg)
    expw2.calc_soil_prod_rate()
    np.allclose(mg.at_node["soil_production__rate"][mg.core_nodes], 1.0)


def test_run_step_no_dt():
    # This tests drop-in compatiblity with the existing ExponentialWeatherer component
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("soil_production__rate", at="node")

    expw2 = ExponentialWeathererIntegrated(mg)
    expw2.run_one_step()

    errors = []
    if not np.allclose(
        mg.at_node["soil_production__dt_weathered_depth"][mg.core_nodes], 0.0
    ):
        errors.append("error in weathered depth")
    if not np.allclose(mg.at_node["soil_production__rate"][mg.core_nodes], 1.0):
        errors.append("error in basic production rate")
    # assert no error message has been registered, else print messages
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


def test_with_big_dt():
    # With a timestep
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("soil__depth", at="node")

    expw2 = ExponentialWeathererIntegrated(mg)
    dt = 1000.0
    expw2.run_one_step(dt)

    errors = []
    if not np.allclose(
        mg.at_node["soil_production__dt_weathered_depth"][mg.core_nodes], 6.9088
    ):
        errors.append("error in weathered depth")
    if not np.allclose(
        mg.at_node["soil_production__dt_produced_depth"][mg.core_nodes], 6.9088
    ):
        errors.append("error in produced depth")
    if not np.allclose(dt * mg.at_node["soil_production__rate"][mg.core_nodes], 1000.0):
        errors.append("error in basic euler integration")
    # assert no error message has been registered, else print messages
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


def test_density_contrast():
    # Different densities
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("soil__depth", at="node")
    dt = 1000.0
    expw3 = ExponentialWeathererIntegrated(
        mg,
        soil_production__maximum_rate=1.0,
        soil_production__decay_depth=1.0,
        soil_production__expansion_factor=1.3,
    )

    expw3.run_one_step(dt)
    errors = []
    # replace assertions by conditions
    if not np.allclose(
        mg.at_node["soil_production__dt_weathered_depth"][mg.core_nodes], 5.5161
    ):
        errors.append("error in weathered depth")
    if not np.allclose(
        mg.at_node["soil_production__dt_produced_depth"][mg.core_nodes], 7.1709
    ):
        errors.append("error in produced depth")
    if not np.allclose(dt * mg.at_node["soil_production__rate"][mg.core_nodes], 1000.0):
        errors.append("error in basic euler integration")
    # assert no error message has been registered, else print messages
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


def test_hex_grid():
    mg = HexModelGrid((5, 5))
    mg.add_zeros("soil__depth", at="node")
    dt = 1000.0
    expw3 = ExponentialWeathererIntegrated(
        mg,
        soil_production__maximum_rate=1.0,
        soil_production__decay_depth=1.0,
        soil_production__expansion_factor=1.3,
    )

    expw3.run_one_step(dt)
    errors = []
    # replace assertions by conditions
    if not np.allclose(
        mg.at_node["soil_production__dt_weathered_depth"][mg.core_nodes], 5.5161
    ):
        errors.append("error in weathered depth")
    if not np.allclose(
        mg.at_node["soil_production__dt_produced_depth"][mg.core_nodes], 7.1709
    ):
        errors.append("error in produced depth")
    if not np.allclose(dt * mg.at_node["soil_production__rate"][mg.core_nodes], 1000.0):
        errors.append("error in basic euler integration")
    # assert no error message has been registered, else print messages
    assert not errors, "errors occured:\n{}".format("\n".join(errors))
