"""
Created on Mon Tuesday Oct 27 2026

@author: PaulMorgan
Built on the architecture of BedrockLandslider by Benjamin Campforts
and repeats many of the same tests


Doc tests and unit tests for soil landslides using the SoilLandslider component.
"""

import numpy as np
import pytest
from numpy import testing

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import SoilLandsliderGeo
from landlab.components import PriorityFloodFlowRouter

try:
    PriorityFloodFlowRouter.load_richdem()
except ModuleNotFoundError:
    pytestmark = pytest.mark.skip(reason="richdem is not installed")


def test_inputFields_flowRouter():
    """
    SoilLandslider should throw an error when topograhy is not equal to the sum of
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
        SoilLandsliderGeo(mg)


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
        SoilLandsliderGeo(mg)


def test_inputFields_bedrock():
    """
    SoilLandslider should create the bedrock__elevation field when not provided
    """
    n_rows = 5
    n_columns = 5
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("soil__depth", at="node")
    # mg.add_zeros("bedrock__elevation", at='node')

    PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)
    SoilLandsliderGeo(mg)

    # instantiate the slider
    assert "bedrock__elevation" in mg.at_node


def test_properties_phi_fraction_fines_LS():
    """
    SoilLandlsider should throw an error when phi/fraction_fines_LS < 0 or phi > 0
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
        SoilLandsliderGeo(mg, phi=-0.2)
    # instantiate the slider
    with pytest.raises(ValueError):
        SoilLandsliderGeo(mg, phi=1.2)
    # instantiate the slider
    with pytest.raises(ValueError):
        SoilLandsliderGeo(mg, fraction_fines_LS=-0.2)
    # instantiate the slider
    with pytest.raises(ValueError):
        SoilLandsliderGeo(mg, fraction_fines_LS=1.2)


def test_mass_balance_noporosity():
    """
    Test is mass balance is conserved during sliding events.
    First test for soils with zero porosity and where none of the sediment is
    evacuated as suspended sediment.
    """
    # Make a raster model grid and create a slope
    n_rows = 9
    n_columns = 9
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")

    # make the top half of the grid a slope, sloping down from the north
    mid_row = n_rows // 2          
    y_mid = mid_row * spacing          # y-coordinate of the middle row
    y = mg.y_of_node
    slope = 1   # 45 degrees = gradient of 1 m/m
    b[:] = np.where(y < y_mid, 1.0, 1.0 + (y - y_mid) * slope)
    #add soil depth
    s += 1
    #get the topography to match
    z[:] = b + s
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)
    fd.run_one_step()

    # Instanciate the slider, keep fraction fines LS at 0 suspended sediment should be 0 
    # Allow several landslides to occur
    sls = SoilLandsliderGeo(
        mg,
        angle_int_frict=0.8,
        landslides_return_time=10,
        landslides_on_boundary_nodes=True,
        phi=0,
        fraction_fines_LS=0,
        verbose_landslides=True,
    )

    # all the soil before 
    total_vol_before = np.sum(b) + np.sum(s)

    # Keep track of total volumes evacuated SSY-suspended sediment yield
    vol_SSY_tot = 0
    V_leaving_tot = 0
    
    # Run for some time
    for _ in range(5):
        fd.run_one_step()
        vol_SSY, V_leaving = sls.run_one_step(dt=1)
        vol_SSY_tot += vol_SSY
        V_leaving_tot += V_leaving

    # Check if list with landslide properties is generated
    assert isinstance(sls.landslides_size, list)
    assert isinstance(sls.landslides_volume, list)
    assert isinstance(sls.landslides_volume_sed, list)

    total_vol_after = np.sum(b) + np.sum(s) + vol_SSY_tot + V_leaving_tot
    # Check mass balance
    err_msg = "Mass balance not okey"
    testing.assert_almost_equal(
        total_vol_before, total_vol_after, decimal=5, err_msg=err_msg
    )



def test_factor_of_safety():
    """
    Test if the factor of safety function in SoilLandslider from data gathered from 
    the grid is equivalent to FS values calculated using the equation. 

    """
    
    #first set the parameters that we will use
    #marked as Theoretical
    theo_slope_degrees = 45.0 # degrees
    theo_soildepth = 1.0  # [m]
    theo_relativewetness = 1 # Rw fraction of 1
    theo_cohesion = 4000 # cohesion_eff [Pa]
    theo_frictionangle = 0.8  # angle_int_frict, tan(phi) [m/m]
    theo_rho_s = 2000.0  # [kg/m^3]
    theo_rho_h2o = 1000.0  # [kg/m^3]
    theo_g = 9.81  # [m/s^2]
    
    # then calculate a theoretical factor of safety for a cell with these values
    #follow Eq. 1 in the SoilLandslider paper
    theoretical_FS = (
        theo_cohesion + theo_g * theo_soildepth * np.cos(np.deg2rad(theo_slope_degrees))**2
        * (theo_rho_s - theo_relativewetness * theo_rho_h2o)
        * np.tan(np.arctan(theo_frictionangle))
    ) / (
        theo_rho_s * theo_g * theo_soildepth
        * np.cos(np.deg2rad(theo_slope_degrees)) * np.sin(np.deg2rad(theo_slope_degrees))
    )
    
    # now create a specific grid to match the above scenario
    # Make a raster model grid and create a slope on the left half of the grid
    # need 5 rows to avoid some edge effects
    n_rows = 5
    n_columns = 13
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")
    
    
    # make the left half of the grid a slope
    mid_column = n_columns // 2 
    x_mid = mid_column * spacing  # x-coordinate of the middle column
    x = mg.x_of_node
     
    grid_slope = np.tan(np.deg2rad(45))  # 45 degrees or a gradient of 1.0 m/m
    b[:] = np.where(x < x_mid, 1.0 + (x_mid - x) * grid_slope, 1.0) # set the gridslope values
    
    # add soil depth
    s += 1
    #set to topography 
    z[:] = b + s
    
    #run the flow router
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)
    fd.run_one_step()
    
    #instantiate the soillandslier, make sure the input parameters match up
    sls = SoilLandsliderGeo(
        mg,
        angle_int_frict=theo_frictionangle,
        cohesion_eff=theo_cohesion,
        rho_s=theo_rho_s,
        rho_h2o=theo_rho_h2o,
        Rw_mode="constant",
        Rw=theo_relativewetness,
        landslides_return_time=100000000000000,
        landslides_on_boundary_nodes=False,
        verbose_landslides=True,
    )
    sls.run_one_step(dt=1)
    
    # now read in the calculated FS value from a specific cell
    # at cell 29 slope should be 1 and H  should be 1
    calculated_FS=mg.at_node["factor_of_safety"][29]
    
    err_msg = "FS values don't match"
    testing.assert_almost_equal(
        theoretical_FS, calculated_FS, decimal=5, err_msg=err_msg
    )
    
def test_stable_stopper_factor():
    """
    Test if the stable stopper factor can halt landslide propagation as intended. 
    
    """    

    # start by building a grid the same as the test factor of safety grid
    # Make a raster model grid and create a slope on the left half of the grid
    # need 5 rows to avoid some edge effects
    n_rows = 5
    n_columns = 13
    spacing = 1
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
    z = mg.add_zeros("topographic__elevation", at="node")
    s = mg.add_zeros("soil__depth", at="node")
    b = mg.add_zeros("bedrock__elevation", at="node")
    
    
    # make the left half of the grid a slope
    mid_column = n_columns // 2 
    x_mid = mid_column * spacing  # x-coordinate of the middle column
    x = mg.x_of_node
     
    grid_slope = np.tan(np.deg2rad(45))  # 45 degrees or a gradient of 1.0 m/m
    b[:] = np.where(x < x_mid, 1.0 + (x_mid - x) * grid_slope, 1.0) # set the gridslope values
    
    # add soil depth
    s += 1
    #set to topography 
    z[:] = b + s
    
    #now remove the soil from a set of cells, to lower the FS and make a barrier
    # to landslide propagation
    #remove the soil from one column and see if the landslide still propagates up
    alterednodes = [41, 28, 15] # top, middle, bottom, with one row of closed nodes above and below
    s[alterednodes] = [0.1, 0.1, 0]   # new soil depths (m), one per node, in the same order
    # 0.1 should make a stable cell. and 0 should test no soil edge cases

    z[:] = b + s                 # keep topography = bedrock + soil

    
    #run the flow router
    fd = PriorityFloodFlowRouter(mg, separate_hill_flow=True, suppress_out=True)
    fd.run_one_step()
    
    #now run the landlsider
    # note the FS should be unstable on the selected cell
    #the stable stopper value is set to 2, which should
    
    sls = SoilLandsliderGeo(
        mg,
        angle_int_frict=0.8,
        landslides_on_boundary_nodes=False,
        critical_sliding_nodes=[30],
        verbose_landslides=True,
        stable_stopper_factor=2,
    )
    sls.run_one_step(dt=1)
    
    # The test, did node 27 fail, it should be above the barrier
    test_cell_erosion=mg.at_node["landslide__erosion"][27]

    err_msg = "FS stable stopper didn't stop the slide"
    testing.assert_equal(test_cell_erosion, 0, err_msg=err_msg)
    






