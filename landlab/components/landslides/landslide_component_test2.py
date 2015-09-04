#!/usr/bin/env python
"""
Purpose is to calculate a factor-of-safety for nodes or grid cells
based on topograhic and soils inputs and recharge from a hydrologic model.
Uses Monte Carlo simulations to estimate probability of failure at each node.

@author: rstrauch - Univerity of Washington
Created on Thu Aug 20 16:47:11 2015
Last edit Aug 31, 2015
"""
#%% Import
import numpy as np
from math import sin,cos,radians
from scipy import interpolate

#%% Instantiate

class Factor_of_Safety():
    """
    This component is designed to calculate a stability index (Factor of Safety) and
    probability of failure at the grid node based on the infinite slope stability model.
    
    The driving force for failure is provided by the user in the form of groundwater
    recharge generated from an hydrology model.  The model also uses topographic 
    and soils characteristics input in the landslide_driver.
    
    A factor-of-safey calcuation function provides the user with the mean soil
    wetness index, mean factor-of-safety, and probabilty of failure at each node.
    """

# component name
    _name = 'Factor_of_Safety'

# component requires these values to do its calculation, get from driver
    _input_var_names = set([
        'number_of_interations',
        'contributing_area',
        'slope',
        'soil_transmissivity__mean',
        'combined_cohesion__mean',
        'soil_internal_angle_friction__mean',
        'soil_density',
        'soil_thickness__mean',        
        'daily_recharge',
        'daily_recharge_probability',
        'number_days',
    ])
    
#  component creates these values   
    _output_var_names = set([
        'Wetness_index',        
        'Factor_of_Safety',
        'Probability_of_failure',
    ])
    
# units for each field
    _var_units = {
        'contributing_area': 'm',
        'slope': 'degrees',
        'soil_transmissivity__mean': 'm2/day',
        'combined_cohesion__mean': 'Pa or kg/m-s2',
        'soil_internal_angle_friction__mean': 'degrees',
        'soil_density': 'kg/m3',
        'soil_thickness__mean': 'm',        
        'daily_recharge': 'mm/day',
    }
# grid centering of each field
    _var_mapping = {
        'contributing_area': 'node',
        'slope': 'node',
        'soil_transmissivity__mean': 'node',
        'combined_cohesion__mean': 'node',
        'soil_internal_angle_friction__mean': 'node',
        'soil_density': 'node',
        'soil_thickness__mean': 'node',        
        'daily_recharge': 'node',
        'daily_recharge_probability': 'node',
        'Wetness_index': 'node',        
        'Factor_of_Safety': 'node',
        'Probability_of_failure': 'node',
    }

# short description of each field
    _var_defs = {                                                               # may need to be more specific on these?
        'contributing_area': 'specific contributing area (upslope area/cell face length) that drains to node',
        'slope': 'slope of surface at cell',
        'soil_transmissivity__mean': 'mean rate of water transmitted through a unit width of saturated soil',
        'combined_cohesion__mean': 'mean combined root and soil cohesion at node',
        'soil_internal_angle_friction__mean': 'The critical angle just before failure due to friction between particles',
        'soil_density': 'wet bulk density of soil',
        'soil_thickness__mean': 'soil depth to restrictive layer',        
        'daily_recharge': 'distribution of baseflow and runoff from hydrologic model at node; currently annual daily maximums',
        'daily_recharge_probability': 'probabilty associated with daily recharge distribution generated from ECDF', 
        'Wetness_index': 'Indicator of soil wetness; relative depth perched water table within the soil layer',        
        'Factor_of_Safety': '(FS) dimensionless index of stability based on infinite slope stabiliity model',
        'Probability_of_failure': 'number of times FS is <1 out of number of interations user selected',        
    }
  
# Run Component

    def __init__(self, number_of_interations, contributing_area, slope,
        soil_transmissivity__mean, combined_cohesion__mean, 
        soil_internal_angle_friction__mean, soil_density, soil_thickness__mean,
        daily_recharge, daily_recharge__probability, number_of_data):        
        self.n = number_of_interations                                         # [] No. of Monte Carlo simulations
        self.a = contributing_area                                             # [m]
        self.theta = slope                                                     # [degrees]
        self.Tmean = soil_transmissivity__mean                                 # [m2/day]
        self.Cmean = combined_cohesion__mean                                   # [Pa or kg/m-s2]
        self.phi_mean = soil_internal_angle_friction__mean                     # [degrees]
        self.rho = soil_density		                                            # density of soil [kg/m3]
        self.hs_mean = soil_thickness__mean                                     # [m]
        self.daily_recharge = daily_recharge                                   # [mm/d] (presorted)
        self.daily_recharge__probability = daily_recharge__probability           # [] probability
        self.size = number_of_data                                                # [] No. of recharge data
        self.g = 9.81 			                                                # gravity [m/s2]
        self.initialize()  #run these methods
        self.calculate_factor_of_safety()

    def initialize(self):
        # generate distributions to sample from
    ## in future, we'll try triangle distribution with np.random.triangular()
        # Transmissivity (T)
        Tstdev = 0.5*self.Tmean      ## standard deviation for distribution         !!PLACE HOLDER!!
        self.T = np.random.normal(self.Tmean, Tstdev, self.n)                  # normal probability distribution; nx1 array
        # Cohesion
        Cstdev = 0.5*self.Cmean
        self.C = np.random.normal(self.Cmean, Cstdev, self.n)
        # phi - internal angle of friction
        phi_stdev = 0.5*self.phi_mean
        self.phi = np.random.normal(self.phi_mean, phi_stdev, self.n)
        # soil thickness
        hs_stdev = 0.5*self.hs_mean
        self.hs = np.random.normal(self.hs_mean, hs_stdev, self.n)
        # recharge distribution
        Yrand = np.sort(np.random.rand(self.n))                                # n random numbers (0 to 1) in a column
        f = interpolate.interp1d(self.daily_recharge__probability,
            self.daily_recharge, bounds_error=False, fill_value=min(self.daily_recharge)) # interpolate function based on recharge data & probabilty
        daily_recharge_interpolated = f(Yrand)                                 # Recharge interpolated from Yrand probabilities (n count)
        self.Re = np.array(daily_recharge_interpolated)

    def calculate_factor_of_safety(self):  #!!! need to check theta, phi if should be np.arctan!!!
        # calculate components of FS equation
        self.X = self.C/(self.hs*self.rho*self.g)                              # demensionless cohesion
        self.WI = ((self.Re/1000.0)/self.T)*(self.a/np.sin(
            np.radians(self.theta)))                                           # wetness index
        np.place(self.WI, self.WI > 1, 1.0)                                    # make maximum WI = 1.0
        self.WI__mean = np.mean(self.WI)                                        # mean wetnes index
        self.Y = np.tan(np.radians(self.phi))*(1 - self.WI*0.5)

        # calculate Factor-of-safety
        self.FS = self.X/sin(radians(self.theta)) + cos(radians(self.theta))*(
            self.Y/sin(radians(self.theta)))                                   # factor of safety
        self.FS__mean = np.mean(self.FS)		                              # mean factor of safety
        count = 0
        for val in self.FS:                                                    # find how many FS values < 1
            if val < 1.0:
                count = count + 1
        self.FS_L1 = float(count)
#        self.FS_L1 = float(sum(self.FS < 1))			                      # number with unstable FS values (<1)
        self.Probability_of_failure = self.FS_L1/self.n 	   # probability: No. unstable values/total No. of values (n)

#%% Finalize and output

#%% check the following
# units in radians everywhere for slopes/angles?