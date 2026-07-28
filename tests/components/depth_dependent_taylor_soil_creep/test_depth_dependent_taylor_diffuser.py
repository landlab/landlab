#!/usr/bin/env python2
"""
Created on Fri Mar  3 10:39:32 2017

@author: gtucker
"""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import DepthDependentTaylorDiffuser
from landlab.components import ExponentialWeatherer


@pytest.mark.slow
def test_4x7_grid_vs_analytical_solution():
    """Test against known analytical solution."""

    # Create a 4-row by 7-column grid with 10 m spacing
    mg = RasterModelGrid((4, 7), xy_spacing=10.0)

    # Close off top and bottom (N and S) boundaries so it becomes a 1D problem
    mg.set_closed_boundaries_at_grid_edges(False, True, False, True)

    # Create an elevation field, initially zero
    z = mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")

    # Instantiate components, and set their parameters. Note that traditional
    # diffusivity, D, is D = SCE x H*, where SCE is soil-creep efficiency.
    # Here we want D = 0.01 m2/yr and H* = 0,.5 m, so cwe set SCE = 0.02.
    weatherer = ExponentialWeatherer(
        mg, soil_production_maximum_rate=0.0002, soil_production_decay_depth=0.5
    )

    diffuser = DepthDependentTaylorDiffuser(
        mg, soil_transport_velocity=0.01, slope_crit=0.8, soil_transport_decay_depth=0.5
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
    for _ in range(3000):
        z[mg.core_nodes] += baselevel_rate * dt
        z_bedrock[mg.core_nodes] += baselevel_rate * dt

        weatherer.calc_soil_prod_rate()
        diffuser.run_one_step(dt)

    # Test: these numbers represent equilibrium. See Jupyter notebook for
    # calculations.
    my_nodes = mg.nodes[2, :]
    assert_array_equal(
        np.round(z[my_nodes], 1), np.array([0.0, 6.2, 10.7, 12.6, 10.7, 6.2, 0.0])
    )
    assert_array_equal(
        np.round(mg.at_node["soil__depth"][8:13], 2),
        np.array([0.35, 0.35, 0.35, 0.35, 0.35]),
    )


def test_4x7_grid_vs_analytical_solution_array_K():
    """Test against known analytical solution using K as an array of float."""

    # Create a 4-row by 7-column grid with 10 m spacing
    mg = RasterModelGrid((4, 7), xy_spacing=10.0)

    # Close off top and bottom (N and S) boundaries so it becomes a 1D problem
    mg.set_closed_boundaries_at_grid_edges(False, True, False, True)

    # Create an elevation field, initially zero
    z = mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")

    # Create an array of K
    k_array = np.full(mg.number_of_links, 0.01)
    # Instantiate components, and set their parameters. Note that traditional
    # diffusivity, D, is D = SCE x H*, where SCE is soil-creep efficiency.
    # Here we want D = 0.01 m2/yr and H* = 0,.5 m, so cwe set SCE = 0.02.
    weatherer = ExponentialWeatherer(
        mg, soil_production_maximum_rate=0.0002, soil_production_decay_depth=0.5
    )

    diffuser = DepthDependentTaylorDiffuser(
        mg,
        soil_transport_velocity=k_array,
        slope_crit=0.8,
        soil_transport_decay_depth=0.5,
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
    for _ in range(3000):
        z[mg.core_nodes] += baselevel_rate * dt
        z_bedrock[mg.core_nodes] += baselevel_rate * dt

        weatherer.calc_soil_prod_rate()
        diffuser.run_one_step(dt)

    # Test: these numbers represent equilibrium. See Jupyter notebook for
    # calculations.
    my_nodes = mg.nodes[2, :]
    assert_array_equal(
        np.round(z[my_nodes], 1), np.array([0.0, 6.2, 10.7, 12.6, 10.7, 6.2, 0.0])
    )
    assert_array_equal(
        np.round(mg.at_node["soil__depth"][8:13], 2),
        np.array([0.35, 0.35, 0.35, 0.35, 0.35]),
    )


def test_raise_stability_error():
    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros("soil__depth", at="node")
    z = mg.add_zeros("topographic__elevation", at="node")
    BRz = mg.add_zeros("bedrock__elevation", at="node")
    z += mg.node_x.copy() ** 2
    BRz = z.copy() - 1.0
    soilTh[:] = z - BRz
    expweath = ExponentialWeatherer(mg)
    DDdiff = DepthDependentTaylorDiffuser(mg, if_unstable="raise")
    expweath.calc_soil_prod_rate()
    with pytest.raises(RuntimeError):
        DDdiff.soilflux(10)


def test_raise_kwargs_error():
    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros("soil__depth", at="node")
    z = mg.add_zeros("topographic__elevation", at="node")
    BRz = mg.add_zeros("bedrock__elevation", at="node")
    z += mg.node_x.copy() ** 2
    BRz = z.copy() - 1.0
    soilTh[:] = z - BRz
    with pytest.raises(TypeError):
        DepthDependentTaylorDiffuser(mg, diffusivity=1)


def test_infinite_taylor_error():
    mg = RasterModelGrid((5, 5))
    soilTh = mg.add_zeros("soil__depth", at="node")
    z = mg.add_zeros("topographic__elevation", at="node")
    BRz = mg.add_zeros("bedrock__elevation", at="node")
    z += mg.node_x.copy() ** 4
    BRz = z.copy() - 1.0
    soilTh[:] = z - BRz
    expweath = ExponentialWeatherer(mg)
    DDdiff = DepthDependentTaylorDiffuser(mg, nterms=400)
    expweath.calc_soil_prod_rate()
    with pytest.raises(RuntimeError):
        DDdiff.soilflux(10)


def test_variable_K_matches_uniform_scalar():
    mg_scalar = RasterModelGrid((5, 5))

    dt = 1
    soilTh_scalar = mg_scalar.add_zeros("soil__depth", at="node")
    z_scalar = mg_scalar.add_zeros("topographic__elevation", at="node")
    BRz_scalar = mg_scalar.add_zeros("bedrock__elevation", at="node")
    z_scalar += mg_scalar.node_x.copy() ** 2
    BRz_scalar = z_scalar.copy() - 1.0
    soilTh_scalar[:] = z_scalar - BRz_scalar
    expweath_scalar = ExponentialWeatherer(mg_scalar)
    expweath_scalar.calc_soil_prod_rate()
    DDdiff_scalar = DepthDependentTaylorDiffuser(mg_scalar, soil_transport_velocity=2.0)
    DDdiff_scalar.run_one_step(dt)

    mg_array = RasterModelGrid((5, 5))
    soilTh = mg_array.add_zeros("soil__depth", at="node")
    z = mg_array.add_zeros("topographic__elevation", at="node")
    BRz = mg_array.add_zeros("bedrock__elevation", at="node")
    z += mg_array.node_x.copy() ** 2
    BRz = z.copy() - 1.0
    soilTh[:] = z - BRz
    expweath = ExponentialWeatherer(mg_array)
    expweath.calc_soil_prod_rate()
    k_array = np.full(mg_array.number_of_links, 2.0)
    # Initializing component with a float

    # Initializing component with an array of float
    DDdiff_array = DepthDependentTaylorDiffuser(
        mg_array, soil_transport_velocity=k_array
    )
    DDdiff_array.run_one_step(dt)

    np.testing.assert_allclose(
        mg_scalar.at_node["topographic__elevation"],
        mg_array.at_node["topographic__elevation"],
    )


# def test_warn():
#    mg = RasterModelGrid((5, 5))
#    soilTh = mg.add_zeros("soil__depth", at="node")
#    z = mg.add_zeros("topographic__elevation", at="node")
#    BRz = mg.add_zeros("bedrock__elevation", at="node")
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
