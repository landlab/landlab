########## rainfall_driver.py ###############
##
## This component generates rainfall  
## events based on statistical distributions.
##
## No particular units must be used, but it was
## written with the storm units in hours (hr)
## and depth units in millimeters (mm)
##
##
##
## Written by Jordan Adams, 2013.
###########################################

import os

import numpy as np 
import random
from landlab.model_parameter_dictionary import ModelParameterDictionary


_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'preciptest.in')

class PrecipitationDistribution:

    def __init__(self):
        '''All initial values are set to zero until initalized
        using the initialize() method which reads in data using
        the ModelParameter Dictionary and sets the random variables
        according to their respective distribution'''
        self.mean_storm = 0.0
        self.mean_intensity = 0.0
        self.mean_interstorm = 0.0
        self.mean_storm_depth = 0.0

        '''Storm_duration, interstorm_duration, storm_depth and intensity
        are not read in, rather, they are generated through initialize() and
        updated using the update() method'''
        self.storm_duration = 0.0
        self.interstorm_duration = 0.0
        self.storm_depth = 0.0
        self.intensity = 0.0

    def initialize(self, input_file=None):
        MPD = ModelParameterDictionary()
        '''
        We imported methods from ModelParameterDictionary
        to read the parameters from the input file.
 
        Necessary parameters to run this code include:
    
            mean_storm (type: float)
            mean_intensity (type: float)
            mean_interstorm (type: float)
    
        mean_storm_depth is calculated using mean storm and intensity data

        The random variables of storm_duration, interstorm_duration, storm_depth
        and intensity are found using methods declared later in the class
        '''
        if input_file is None:
            input_file = _DEFAULT_INPUT_FILE
            ###input_file = raw_input('Enter the input file location: ')
        MPD.read_from_file(input_file)
        
        
        
        self.mean_storm = MPD.read_float( 'mean_storm')
        self.mean_intensity = MPD.read_float( 'mean_intensity')
        self.mean_interstorm = MPD.read_float( 'mean_interstorm' )
        self.mean_storm_depth = self.mean_intensity * self.mean_storm

        self.storm_duration = self.get_precipitation_event_duration()
        self.interstorm_duration = self.get_interstorm_event_duration()
        self.storm_depth = self.get_storm_depth()
        self.intensity = self.get_storm_intensity()


    def update(self):
        
        '''This method updates the variables of storm_duration, interstorm_duration,
        storm_depth and intensity that were found in initialize() if needed.'''
        
        self.storm_duration = self.get_precipitation_event_duration()
        self.interstorm_duration = self.get_interstorm_event_duration()
        self.storm_depth = self.get_storm_depth()
        self.intensity = self.get_storm_intensity()
        
    def get_precipitation_event_duration(self):
        '''
    This method is the storm generator.
    
    This method has one argument: the mean_storm parameter.
    (In Eagleson (1978), this parameter was called Tr.)
    
    It finds a random storm_duration value
    based on the poisson distribution about the mean.
    This is accomplished using the expovariate function
    from the 'random' standard library.
    Additionally, it is rounded to contain 4 significant figures, 
    for neatness. 
    
    The if-else statement is very important here. Values of 0
    can exist in the Poission distribution, but it does not make
    sense to have 0 duration storms, so to avoid that,
    we return a storm duration IF it is greater than 0,
    otherwise, recursion is employed to re-call the storm
    generator function and get a new value.
    '''
        storm = round(random.expovariate(1/self.mean_storm),2)
        while storm == 0: 
            storm = round(random.expovariate(1/self.mean_storm),2)
        self.storm_duration = storm
        return self.storm_duration

    
    
    def get_interstorm_event_duration(self):
        ''' 
    This method is the interstorm duration generator
    
    This method takes one argument, the mean_interstorm parameter.
    (In Eagleson (1978), this parameter was called Tb.)
    
    This method is modeled identically to get_precipitation_event_duration()
    
    This method finds a random value for interstorm_duration
    based on the poisson distribution about the mean.
    This is accomplished using the expovariate function
    from the 'random' standard library.
    Additionally, it is rounded to contain 4 significant figures, for neatness. 
    
    The if-else statement is very important here. Values of 0
    can exist in the Poission distribution, but it does not make
    sense to have 0 hour interstorm durations (DOES IT NOT?)
    To avoid 0 hour interstorm durations, we return a 
    interstorm duration IF it is greater than 0,
    otherwise, recursion is employed to re-call the interstorm
    duration generator function and get a new value. '''
    
        interstorm = round(random.expovariate(1/self.mean_interstorm),2)
        while interstorm == 0:
            interstorm = round(random.expovariate(1/self.mean_interstorm),2)
        self.interstorm_duration = interstorm
        return self.interstorm_duration
     
                            
    def get_storm_depth(self):
        '''  
        This method is the storm depth generator.
    Storm depth is used to generate a realistic 
    intensity for different storm events.
    
    (In Eagleson (1978) this parameter was called 'h')
    
    This method requires storm_duration, mean_storm duration
    and the mean_storm_depth. Storm_duration is generated through
    the initialize() or update() method. mean_storm and mean_storm_depth
    are read in using the ModelParameterDictionary.
    
    Numpy has a random number generator to get values
    from a given Gamma distribution. It takes two arguments,
    alpha (or the shape parameter), which is the generated over the mean event
    and beta (or the scale parameter), which is the mean value
    These are all arguments in the function, which returns storm depth.
    '''
        
        shape_parameter = (self.storm_duration/self.mean_storm)
        scale_parameter = (self.mean_storm_depth)
        self.storm_depth = np.random.gamma(shape_parameter, scale_parameter)
        return self.storm_depth

    
    def get_storm_intensity(self):
        '''  
        This method draws storm intensity out of the storm depth
    generated by get_storm_depth. 
    
    This method requires the storm_depth and storm_duration
    and is the same as the parameter ('i') in Eagleson (1978), but instead of
    being drawn from Poission, this is drawn from the Gamma distribution
    of ('h'), as h = i*Tr. 
    
                            '''
        self.intensity = self.storm_depth / self.storm_duration
        return self.intensity
