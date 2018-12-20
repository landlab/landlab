#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:39:32 2017

@author: gtucker
"""
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import DepthDependentTaylorDiffuser, ExponentialWeatherer


def test_4x7_grid_vs_analytical_solution():
    """Test against known analytical solution."""

    # Create a 4-row by 7-column grid with 10 m spacing
    mg = RasterModelGrid((4, 7), xy_spacing=10.0)

    # Close off top and bottom (N and S) boundaries so it becomes a 1D problem
    mg.set_closed_boundaries_at_grid_edges(False, True, False, True)

    # Create an elevation field, initially zero
    z = mg.add_zeros("node", "topographic__elevation")

    # Instantiate components, and set their parameters. Note that traditional
    # diffusivity, D, is D = SCE x H*, where SCE is soil-creep efficiency.
    # Here we want D = 0.01 m2/yr and H* = 0,.5 m, so cwe set SCE = 0.02.
    diffuser = DepthDependentTaylorDiffuser(
        mg, linear_diffusivity=0.01, slope_crit=0.8, soil_transport_decay_depth=0.5
    )
    weatherer = ExponentialWeatherer(
        mg, soil_production__maximum_rate=0.0002, soil_production__decay_depth=0.5
    )

    # Get a reference to bedrock elevation field
    z_bedrock = mg.at_node["bedrock__elevation"]

    # Estimate a reasonable time step. Here we take advantage of the fact that
    # we know the final slope at the outer links will be about 1.33. Stability
    # for the cubic term involves an "effective D" parameter, Deff, that should
    # be Deff = D (S / Sc)^2. (see notebook calcs)
    baselevel_rate = 0.0001
    dt = 250.0

    # Run for 750 ky
    for i in range(3000):

        z[mg.core_nodes] += baselevel_rate * dt
        z_bedrock[mg.core_nodes] += baselevel_rate * dt

        weatherer.calc_soil_prod_rate()
        diffuser.run_one_step(dt)

    # Test: these numbers represent equilibrium. See Jupyter notebook for
    # calculations.
    my_nodes = mg.nodes[2, :]
    assert_array_equal(
        np.round(z[my_nodes], 1), np.array([0.0, 4.0, 6.7, 7.7, 6.7, 4.0, 0.0])
    )
    assert_array_equal(
        np.round(mg.at_node["soil__depth"][8:13], 2),
        np.array([0.35, 0.35, 0.35, 0.35, 0.35]),
    )


def test_raise_stability_error():
    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros("node", "soil__depth")
    z = mg.add_zeros("node", "topographic__elevation")
    BRz = mg.add_zeros("node", "bedrock__elevation")
    z += mg.node_x.copy() ** 2
    BRz = z.copy() - 1.0
    soilTh[:] = z - BRz
    expweath = ExponentialWeatherer(mg)
    DDdiff = DepthDependentTaylorDiffuser(mg)
    expweath.calc_soil_prod_rate()
    with pytest.raises(RuntimeError):
        DDdiff.soilflux(10, if_unstable="raise")


def test_raise_kwargs_error():
    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros("node", "soil__depth")
    z = mg.add_zeros("node", "topographic__elevation")
    BRz = mg.add_zeros("node", "bedrock__elevation")
    z += mg.node_x.copy() ** 2
    BRz = z.copy() - 1.0
    soilTh[:] = z - BRz
    with pytest.raises(TypeError):
        DepthDependentTaylorDiffuser(mg, diffusivity=1)


def test_infinite_taylor_error():

    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros("node", "soil__depth")
    z = mg.add_zeros("node", "topographic__elevation")
    BRz = mg.add_zeros("node", "bedrock__elevation")
    z += mg.node_x.copy() ** 4
    BRz = z.copy() - 1.0
    soilTh[:] = z - BRz
    expweath = ExponentialWeatherer(mg)
    DDdiff = DepthDependentTaylorDiffuser(mg, nterms=400)
    expweath.calc_soil_prod_rate()
    with pytest.raises(RuntimeError):
        DDdiff.soilflux(10)


# def test_warn():
#    mg = RasterModelGrid((5, 5))
#    soilTh = mg.add_zeros('node', 'soil__depth')
#    z = mg.add_zeros('node', 'topographic__elevation')
#    BRz = mg.add_zeros('node', 'bedrock__elevation')
#    z += mg.node_x.copy()**2
#    BRz = z.copy() - 1.0
#    soilTh[:] = z - BRz
#    expweath = ExponentialWeatherer(mg)
#    DDdiff = DepthDependentTaylorDiffuser(mg)
#    expweath.calc_soil_prod_rate()
#
#    with warnings.catch_warnings(record=True) as w:
#    # Cause all warnings to always be triggered.
#        warnings.simplefilter("always")
#        # Trigger a warning.
#        DDdiff.soilflux(dt=10, if_unstable='warn')
#        # Verify some things
#        assert len(w) == 1
#        assert issubclass(w[-1].category, RuntimeWarning)

if __name__ == "__main__":
    test_4x7_grid_vs_analytical_solution()
