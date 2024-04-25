"""
Created on Mon Nov 08 13:50:43 2021

@author: BenjaminCampforts


Doc tests and unit tests for bedrock landslides using the BedrockLandslider component.
"""

import numpy as np
import pytest
from numpy import testing

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import BedrockLandslider
from landlab.components import PriorityFloodFlowRouter

try:
    PriorityFloodFlowRouter.load_richdem()
except ModuleNotFoundError:
    pytestmark = pytest.mark.skip(reason="richdem is not installed")


def test_inputFields_flowRouter():
    """
    BedrockLandslider should throw an error when topograhy is not equal to the sum of
    bedrock and soil thickness
    """
    # Make a raster model grid and create a plateau
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make plateau at 10m
    b += 10

    PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # instantiate the slider
    with pytest.raises(RuntimeError):
        BedrockLandslider(mg)


def test_inputFields_soil():
    """
    BedrockLandslider should throw an error when the soil__depth field is not provided
    """
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    # mg.add_zeros("soil__depth", at='node')
    mg.add_zeros("bedrock__elevation", at="node")

    PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # instantiate the slider
    with pytest.raises(FieldError):
        BedrockLandslider(mg)


def test_inputFields_bedrock():
    """
    BedrockLandslider should create the bedrock__elevation field when not provided
    """
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    # mg.add_zeros("bedrock__elevation", at='node')

    PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)
    BedrockLandslider(mg)

    # instantiate the slider
    assert "bedrock__elevation" in mg.at_node


def test_properties_phi_fraction_fines_LS():
    """
    BedrockLandslider should throw an error when phi/fraction_fines_LS < 0 or phi > 0
    """
    # Make a raster model grid and create a plateau
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    mg.add_zeros("bedrock__elevation", at="node")

    PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # instantiate the slider
    with pytest.raises(ValueError):
        BedrockLandslider(mg, phi=-0.2)
    # instantiate the slider
    with pytest.raises(ValueError):
        BedrockLandslider(mg, phi=1.2)
    # instantiate the slider
    with pytest.raises(ValueError):
        BedrockLandslider(mg, fraction_fines_LS=-0.2)
    # instantiate the slider
    with pytest.raises(ValueError):
        BedrockLandslider(mg, fraction_fines_LS=1.2)


def test_sliding_plain():
    """
    Test if BedrockLandslider maths are sound and follow the Culmann theory
    """
    # Make a raster model grid and create a plateau
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make plateau at 10m
    b += 10
    # Lower boundary cell to 0
    b[2] = 0
    z[:] = b + s
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)
    hy = BedrockLandslider(mg, landslides_return_time=1)

    # run flow director and BedrockLandslider for one timestep
    fd.run_one_step()
    hy.run_one_step(dt=1)
    # After one timestep, we can predict eactly where the landslide will occur.
    # The return time is set to 1 year so that probability for sliding is 100%
    # The angle of internal driction is 1m/m, the topographical gradient is 10 m/m
    # At cardinal cells, the sliding plain will be at (1+10)/2 = 5.5 m/m.
    # With a dx of 1, the cardial cell next to the critical sliding node must
    # be 5.5 m and the diagonal one at 5.5 * sqrt(2) = 7.8 m
    err_msg = "Error in the calculation of the sliding plain"
    testing.assert_almost_equal(
        [5.5 * np.sqrt(2), 5.5, 5.5 * np.sqrt(2)], z[6:9], decimal=5, err_msg=err_msg
    )


def test_sliding_evolution():
    """
    Test that if BedrockLandslider is run for a long enough time, slopes evolve to the
    angle of internal friction
    """
    # Make a raster model grid and create a plateau
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make plateau at 10m
    b += 10
    # Lower boundary cell to 0
    b[2] = 0
    z[:] = b + s
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # Instanciate the slider, but this time provide the node at which landslides
    # should occur (node 2)
    hy = BedrockLandslider(
        mg,
        angle_int_frict=1,
        threshold_slope=1,
        landslides_return_time=0.01,
        landslides_on_boundary_nodes=True,
        critical_sliding_nodes=[2],
    )

    # run flow director and BedrockLandslider for long enough so that surface evolves
    # towards the angle of internal friction timestep

    for _ in range(50):
        fd.run_one_step()
        hy.run_one_step(dt=1)

    # After 50 iterations, landslide erosion starting fromnode 2 evolves to a
    # known slope, equal to the angle of internal friction, multiplied by the
    # distance from node 2 giving:
    topo_cal = np.array(
        [
            [2.0, 1.0, 0.0, 1.0, 2.0],
            [2.23606798, 1.41421356, 1.0, 1.41421356, 2.23606798],
            [2.82842712, 2.23606798, 2.0, 2.23606798, 2.82842712],
            [3.60555128, 3.16227766, 3.0, 3.16227766, 3.60555128],
            [4.47213595, 4.12310563, 4.0, 4.12310563, 4.47213595],
        ]
    ).flatten()

    err_msg = "Slope is not evolving to theoretical value"
    testing.assert_almost_equal(topo_cal, z, decimal=5, err_msg=err_msg)


def test_boundary_nodes():
    # %%
    """
    Test that if BedrockLandslider cannot make or initate landslides at boundary nodes,
    it doesn't
    """
    # Make a raster model grid and create a plateau
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make plateau at 10m
    b += 10
    # Lower boundary cell to 0
    b[2] = 0
    z[:] = b + s
    ini_topo_bd = z[mg.boundary_nodes]
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # Instanciate the slider, but this time provide the node at which landslides
    # should occur (node 2)
    hy = BedrockLandslider(
        mg,
        angle_int_frict=1,
        landslides_return_time=0.01,
        landslides_on_boundary_nodes=False,
    )

    # run flow director and BedrockLandslider for long enough so that surface evolves
    # towards the angle of internal friction timestep

    for _ in range(5):
        fd.run_one_step()
        hy.run_one_step(dt=1)
    # Final topo at boundary nodes should be inital topo
    err_msg = "No landslide erosion allowed at boundary nodes"

    testing.assert_almost_equal(
        ini_topo_bd, z[mg.boundary_nodes], decimal=5, err_msg=err_msg
    )


# %%
def test_mass_balance_noporosity():
    """
    Test is mass balance is conserved during sliding events.
    First test for soils with zero porosity and where none of the sediment is
    evacuated as suspended sediment.
    """
    # Make a raster model grid and create a plateau
    n_rows = 9
    n_columns = 9
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make plateau at 10m
    b += 10
    # Elevate central node
    b[40] = 10
    # Create slope at margin so that sediment will be evacuated
    b[mg.boundary_nodes] = 0
    z[:] = b + s
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # Instanciate the slider, but this time provide the node at which landslides
    # should occur (node 2)
    hy = BedrockLandslider(
        mg,
        angle_int_frict=1,
        landslides_return_time=1,
        landslides_on_boundary_nodes=True,
        phi=0,
        fraction_fines_LS=0,
        verbose_landslides=True,
    )

    # Assume equal bulk density for rock and soil
    total_vol_before = np.sum(b) + np.sum(s)

    # Keep track of total volumes evacuated
    vol_SSY_tot = 0
    V_leaving_tot = 0
    # Run for some time
    for _ in range(5):
        fd.run_one_step()
        vol_SSY, V_leaving = hy.run_one_step(dt=1)
        vol_SSY_tot += vol_SSY
        V_leaving_tot += V_leaving

    # Check if list with landslide properties is generated
    assert isinstance(hy.landslides_size, list)
    assert isinstance(hy.landslides_volume, list)
    assert isinstance(hy.landslides_volume_bed, list)
    assert isinstance(hy.landslides_volume_sed, list)

    total_vol_after = np.sum(b) + np.sum(s) + vol_SSY_tot + V_leaving_tot
    # Check mass balance
    err_msg = "Mass balance not okey"
    testing.assert_almost_equal(
        total_vol_before, total_vol_after, decimal=5, err_msg=err_msg
    )


def test_mass_balance_porosity():
    """
    Test is mass balance is conserved during sliding events.
    Test for soils with given porosity phi and where none of the sediment is
    evacuated as suspended sediment.
    """
    # Make a raster model grid and create a plateau
    n_rows = 9
    n_columns = 9
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make plateau at 10m
    b += 10
    # Elevate central node
    b[40] = 10
    # Create slope at margin so that sediment will be evacuated
    b[mg.boundary_nodes] = 0
    z[:] = b + s
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # Instanciate the slider, but this time provide the node at which landslides
    # should occur (node 2)
    phi = 0.3
    hy = BedrockLandslider(
        mg,
        angle_int_frict=1,
        landslides_return_time=1,
        landslides_on_boundary_nodes=True,
        phi=phi,
        fraction_fines_LS=0,
    )

    # Assume equal bulk density for rock and soil
    total_vol_before = np.sum(b) + np.sum(s) * (1 - phi)

    # Keep track of total volumes evacuated
    vol_SSY_tot = 0
    V_leaving_tot = 0
    # Run for some time
    for _ in range(5):
        fd.run_one_step()
        vol_SSY, V_leaving = hy.run_one_step(dt=1)
        vol_SSY_tot += vol_SSY
        V_leaving_tot += V_leaving

    total_vol_after = np.sum(b) + np.sum(s) * (1 - phi) + vol_SSY_tot + V_leaving_tot
    # Check mass balance
    err_msg = "Mass balance not okey"
    testing.assert_almost_equal(
        total_vol_before, total_vol_after, decimal=5, err_msg=err_msg
    )


def test_mass_balance_porosity_suspension():
    """
    Test is mass balance is conserved during sliding events.
    Test for soils with given porosity phi and where none of the sediment is
    evacuated as suspended sediment.
    """
    # %%
    # Make a raster model grid and create a plateau
    n_rows = 9
    n_columns = 9
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make plateau at 10m
    b += 10
    # Elevate central node
    b[40] = 10
    # Create slope at margin so that sediment will be evacuated
    b[mg.boundary_nodes] = 0
    z[:] = b + s
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)

    # Instanciate the slider, but this time provide the node at which landslides
    # should occur (node 2)
    phi = 0.3
    fraction_fines_LS = 0.5
    hy = BedrockLandslider(
        mg,
        angle_int_frict=1,
        landslides_return_time=1,
        landslides_on_boundary_nodes=True,
        phi=phi,
        fraction_fines_LS=fraction_fines_LS,
    )

    # Assume equal bulk density for rock and soil
    total_vol_before = np.sum(b) + np.sum(s) * (1 - phi)

    # Keep track of total volumes evacuated
    vol_SSY_tot = 0
    V_leaving_tot = 0
    # Run for some time
    for _ in range(5):
        fd.run_one_step()
        vol_SSY, V_leaving = hy.run_one_step(dt=1)
        vol_SSY_tot += vol_SSY
        V_leaving_tot += V_leaving

    total_vol_after = np.sum(b) + np.sum(s) * (1 - phi) + vol_SSY_tot + V_leaving_tot
    # Check mass balance
    err_msg = "Mass balance not okey"
    testing.assert_almost_equal(
        total_vol_before, total_vol_after, decimal=5, err_msg=err_msg
    )
