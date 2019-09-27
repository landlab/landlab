#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:26:31 2019

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_equal

from landlab import RasterModelGrid
from landlab.components import GroundwaterDupuitPercolator


def test_simple_water_table():
    """Test a one-node steady simulation.

    Notes
    -----
    The analytical solution for one interior cell with one open boundary is as
    follows. The incoming recharge must equal the outgoing discharge. The
    incoming recharge is R, and the surface area is 1 m2, so the incoming
    volume per unit time is R m/s x 1 m x 1 m. The outgoing discharge is equal
    to the conductivity, K (m/s), times the thickness at the boundary, Hb,
    times the hydraulic gradient, which in this case is (H - 0) / dx = H.
    The thickness at the boundary in this test is H/2. Therefore:

        K H^2 / 2 = R, or

        H = sqrt( 2 R / K ).

    With R = 10^-8 m/s and K = 10^-2 m/s, we should have

        H ~ 0.00141 m.
    """
    boundaries = {"top": "closed", "left": "closed", "bottom": "closed"}
    rg = RasterModelGrid((3, 3), bc=boundaries)
    gdp = GroundwaterDupuitPercolator(
        rg, recharge_rate=1.0e-8, hydraulic_conductivity=0.01
    )
    for i in range(12):
        gdp.run_one_step(5.0e4)

    assert_equal(np.round(gdp._thickness[4], 5), 0.00141)

    # Re-instantiate to test the case when the necessary fields already exist
    gdp = GroundwaterDupuitPercolator(rg)
