import pytest
import numpy as np
from landlab import RasterModelGrid
from landlab.components import LandslideProbability

@pytest.fixture
def spacing():
    spacing = (10e0, 10e0)
    return spacing

@pytest.fixture
def shape():
    shape = (5, 4)
    return shape
@pytest.fixture
def dimensions():
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    return dimensions
    
@pytest.fixture
def ls_prob_zeros():
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    shape = dimensions[0]
    spacing = dimensions[1]
    coordinates = dimensions[2]
    
    grid = RasterModelGrid(shape=shape,xy_spacing=spacing,xy_of_lower_left=coordinates)
    grid.add_zeros("topographic__slope", at="node", dtype=float)
    grid.add_zeros("topographic__specific_contributing_area", at="node", dtype=int)
    grid.add_zeros("soil__transmissivity", at="node", dtype=float)
    grid.add_zeros("soil__saturated_hydraulic_conductivity", at="node", dtype=float)
    grid.add_zeros("soil__mode_total_cohesion", at="node", dtype=int)
    grid.add_zeros("soil__minimum_total_cohesion", at="node", dtype=int)
    grid.add_zeros("soil__maximum_total_cohesion", at="node", dtype=int)
    grid.add_zeros("soil__internal_friction_angle", at="node", dtype=int)
    grid.add_zeros("soil__density", at="node", dtype=int)
    grid.add_zeros("soil__thickness", at="node", dtype=float)

    return LandslideProbability(grid)

@pytest.fixture
def unit_default_value():

    unit_default_value={}
    unit_default_value["topographic__slope"]=(0.1,0.8)
    unit_default_value["topographic__specific_contributing_area"]=(10,100)
    unit_default_value["soil__transmissivity"]=(0.1, 10)
    unit_default_value["soil__mode_total_cohesion"]=(30, 900)
    mode = unit_default_value["soil__mode_total_cohesion"]
    # code used to run the model uses a random scatter <cohesion_scatter = np.random.randint(1,10)>
    scatter = np.random.randint(1,10)  #smallest min is the largest scatter, minus the mode cohesion
    unit_default_value["soil__minimum_total_cohesion"]=((mode[0]-scatter),mode[1]-scatter)
    #large max is the largest scatter, plus the mode cohesion
    unit_default_value["soil__maximum_total_cohesion"]=((mode[0]+scatter),mode[1]+scatter)
    unit_default_value["soil__internal_friction_angle"]=(26, 37)
    unit_default_value["soil__thickness"]=(0.1,1) #functionally setting to uniform uniform 1 meter default; coded for the option of a range of values
    unit_default_value["soil__density"]=(1999,2001) #functionally setting to uniform 2000
    #Note on Ksat default range calculation
    #reverse order: (see above) lowest Transmissivity value are in upper right of synthethic grid,
    #not reversed: (see above) highest soil thickness in upper right of synthethic grid,
    #caution creating Ksat since parameterization of T and hs will have inverse relationship on Ksat
    #Kmin  = Tmin/hs_max   -slow- think of a deep thick clay layer 
    kmin = unit_default_value['soil__transmissivity'][0]/unit_default_value['soil__thickness'][1] 
    #Kmax = Tmax/hs_min    -fast- think of a sieve or a parking lot
    kmax = unit_default_value['soil__transmissivity'][1]/unit_default_value['soil__thickness'][0] 
    unit_default_value["soil__saturated_hydraulic_conductivity"]= (kmin,kmax)

    #Output value range used in tests will change as a function of input default values
    unit_default_value["landslide__probability_of_failure"]= (0,1)
    unit_default_value["soil__mean_relative_wetness"]=(0,1)
    unit_default_value["soil__mean_watertable_depth"]=(0,1000)
    unit_default_value["soil__mean_recharge"]=(0,1000)
    unit_default_value["soil__probability_of_saturation"]=(0,1)
    return unit_default_value

@pytest.fixture
def ls_prob():
    dimensions = ((5, 4), (10e0, 10e0), (0.0, 0.0))
    shape = dimensions[0]
    spacing = dimensions[1]
    coordinates = dimensions[2]
    grid = RasterModelGrid(shape=shape,xy_spacing=spacing,xy_of_lower_left=coordinates)
    gridnum = grid.number_of_nodes
    np.random.seed(seed=5)
    scatter = np.random.randint(1,10)  #smallest min is the largest scatter, minus the mode cohesion
    
    unit_default_value={}
    unit_default_value["topographic__slope"]=(0.1,0.8)
    unit_default_value["topographic__specific_contributing_area"]=(10,100)
    unit_default_value["soil__transmissivity"]=(0.1, 10)
    unit_default_value["soil__mode_total_cohesion"]=(30, 900)
    mode = unit_default_value["soil__mode_total_cohesion"]
    # code used to run the model uses a random scatter <cohesion_scatter = np.random.randint(1,10)>
    scatter = np.random.randint(1,10)  #smallest min is the largest scatter, minus the mode cohesion
    unit_default_value["soil__minimum_total_cohesion"]=((mode[0]-scatter),mode[1]-scatter)
    #large max is the largest scatter, plus the mode cohesion
    unit_default_value["soil__maximum_total_cohesion"]=((mode[0]+scatter),mode[1]+scatter)
    unit_default_value["soil__internal_friction_angle"]=(26, 37)
    unit_default_value["soil__thickness"]=(0.1,1) #functionally setting to uniform uniform 1 meter default; coded for the option of a range of values
    unit_default_value["soil__density"]=(1999,2001) #functionally setting to uniform 2000
    #Note on Ksat default range calculation
    #reverse order: (see above) lowest Transmissivity value are in upper right of synthethic grid,
    #not reversed: (see above) highest soil thickness in upper right of synthethic grid,
    #caution creating Ksat since parameterization of T and hs will have inverse relationship on Ksat
    #Kmin  = Tmin/hs_max   -slow- think of a deep thick clay layer 
    kmin = unit_default_value['soil__transmissivity'][0]/unit_default_value['soil__thickness'][1] 
    #Kmax = Tmax/hs_min    -fast- think of a sieve or a parking lot
    kmax = unit_default_value['soil__transmissivity'][1]/unit_default_value['soil__thickness'][0] 
    unit_default_value["soil__saturated_hydraulic_conductivity"]= (kmin,kmax)

    #Output value range used in tests will change as a function of input default values
    unit_default_value["landslide__probability_of_failure"]= (0,1)
    unit_default_value["soil__mean_relative_wetness"]=(0,1)
    unit_default_value["soil__mean_watertable_depth"]=(0,1000)
    unit_default_value["soil__mean_recharge"]=(0,1000)
    unit_default_value["soil__probability_of_saturation"]=(0,1)
  
    parameter_range = unit_default_value

    #Slope
    grid.at_node["topographic__slope"]= np.sort(np.random.uniform(
        parameter_range["topographic__slope"][0], 
        parameter_range["topographic__slope"][1], 
        gridnum).astype(float)
    )    
    #reverse order: lowest slopes are in upper right of synthethic grid
    grid['node']['topographic__slope'] = grid['node']['topographic__slope'][::-1]  

    #Area
    grid.at_node["topographic__specific_contributing_area"] = np.sort(np.random.randint(
        parameter_range["topographic__specific_contributing_area"][0],
        parameter_range["topographic__specific_contributing_area"][1], 
        gridnum).astype(int)
    #no reverse order: highest contributing areas are in upper right of synthethic grid
    )
    
    #Transmissivity
    grid.at_node["soil__transmissivity"] = np.sort(np.random.uniform(
        parameter_range["soil__transmissivity"][0], 
        parameter_range["soil__transmissivity"][1], 
        gridnum).astype(float)
    )    
    #reverse order: lowest Transmissivity value are in upper right of synthethic grid
    grid['node']['soil__transmissivity'] = grid['node']['soil__transmissivity'][::-1]  
       
    #Cohesion
    grid.at_node["soil__mode_total_cohesion"] = np.random.randint(
        parameter_range["soil__mode_total_cohesion"][0], 
        parameter_range["soil__mode_total_cohesion"][1], 
        gridnum).astype(int)
      
    grid.at_node["soil__minimum_total_cohesion"] = (
        grid.at_node["soil__mode_total_cohesion"] - scatter 
        ).astype(int)
              
    grid.at_node["soil__maximum_total_cohesion"] = (
        grid.at_node["soil__mode_total_cohesion"] + scatter 
        ).astype(int)
                                                           
    #Internal angle of friction
    grid.at_node["soil__internal_friction_angle"] = np.sort(np.random.randint(
        parameter_range["soil__internal_friction_angle"][0], 
        parameter_range["soil__internal_friction_angle"][1], 
        gridnum).astype(int)
    )            

    #soil thickness
    grid.at_node["soil__thickness"] = np.sort(np.random.uniform(
        parameter_range["soil__thickness"][0],
        parameter_range["soil__thickness"][1], 
        gridnum).astype(float)
    #no reverse order: greatest soil thickness is in upper right of synthethic grid
    )
    #soil density
    grid.at_node["soil__density"] = np.sort(np.random.randint(
        parameter_range["soil__density"][0],
        parameter_range["soil__density"][1], 
        gridnum).astype(int)
    #no reverse order: greatest density is in upper right of synthethic grid
    )
    
    #K sat - #calculate T in driver with K * D; also usable as calibration multiplier value for T        
    grid['node']['soil__saturated_hydraulic_conductivity']= \
    grid['node']['soil__transmissivity']/grid['node']['soil__thickness']
    #reverse order: (see above) lowest Transmissivity value are in upper right of synthethic grid,
    #not reversed: (see above) highest soil thickness in upper right of synthethic grid,
    #caution creating Ksat since parameterization of T and hs will have inverse relationship on Ksat
    
    return LandslideProbability(grid)