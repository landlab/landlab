#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:39:32 2017

@author: gtucker
"""

from landlab import RasterModelGrid
from landlab.components import (DepthDependentCubicDiffuser,
                                ExponentialWeatherer)


def test_4x7_grid_vs_analytical_solution():
    """Test against known analytical solution."""
    
    # Create a 4-row by 7-column grid with 10 m spacing
    mg = RasterModelGrid((4, 7), 10.0)

    # Close off top and bottom (N and S) boundaries so it becomes a 1D problem
    mg.set_closed_boundaries_at_grid_edges(False, True, False, True)

    # Create an elevation field, initially zero
    z = mg.add_zeros('node', 'topographic__elevation')

    # Instantiate components, and set their parameters. Note that traditional
    # diffusivity, D, is D = SCE x H*, where SCE is soil-creep efficiency.
    # Here we want D = 0.01 m2/yr and H* = 0,.5 m, so we set SCE = 0.02.
    diffuser = DepthDependentCubicDiffuser(mg, soil_creep_efficiency=0.02,
                                           slope_crit=0.8,
                                           soil_transport_decay_depth=0.5)
    weatherer = ExponentialWeatherer(mg, max_soil_production_rate=0.0002,
                                     soil_production_depth_decay=0.5)
    
    
    
if __name__ == '__main__':
    test_4x7_grid_vs_analytical_solution()

    