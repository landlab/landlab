########## radiation.py ###############
##
##  This code calculates incident Shortwave radiation on a 
##  surface taking into account location, slope and aspect
##  of the surface.
##
##  ref: Javier H. Flores Cervantes et al. 2012 "A geomorphic .." Ecohydrol
##
## Written by Sai Nudurupati & Erkan Istanbulluoglu, 2013.
###########################################

import os
from numpy import *
from math import *
import landlab
from landlab import RasterModelGrid  # add landlab. as prefix if have landlab installed
from landlab import ModelParameterDictionary # add landlab. as prefix
from matplotlib.pyplot import *
from random import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from landlab.plot import imshow_grid, imshow_active_cells

# Input file 
_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                 'radiation_input.txt')


class Radiation():

    def __init__( self ):

        self._latitude = 0.0       # Latitude in deg
        self._t = 0.0              # Local Time
        self._dt = 0.0             # Local Time Step 
        self._beta = 0.0           # Aspect
        self._slope = 0.0          # Slope
        self._julian = 0.0         # Julian Day
        self._current_dt = 0.0     # Time Step
        self._current_time = 0.0   # Time
        self._phi = 0.0       # Latitude in Radians
        self._delta = 0.0     # Declination Angle of Sun
        self._tau = 0.0       # Hour Angle
        self._alpha = 0.0     # Solar Altitude
        self._phisun = 0.0    # Sun's Azhimuth
        self._flat = 0.0      # Flat Surface Reference Cos(Solar Incidence)
        #self._sig = 5.6e-8   # Stefan Boltzmann's Constant
        self._n = 2.0         # Clear Sky Turbidity Factor 
        self._m = 0.0         # Optical Air Mass
        self._N = 0.0         # Cloudiness
        self._A = 0.0         # Albedo
        self._Io = 1353.0     # Solar Constant
        
        
    def initialize( self ):
        
        MPD = ModelParameterDictionary()
        MPD.read_from_file( _DEFAULT_INPUT_FILE )

        # Reading Input Parameters
        self._N = MPD.read_float( 'CLOUDINESS' )        
        self._latitude = MPD.read_float( 'LATITUDE' )
        self._A = MPD.read_float( 'ALBEDO' )

        
    def update( self, RMG ):
        
        self._Si = RMG.empty(centering = 'cell')  # Cos(Solar Angle of Incidence)
        self._BETA = RMG.empty(centering = 'cell') # Aspect (Clckwise from N)
        self._Rs = RMG.empty(centering = 'cell')   # Incoming Shortwave Radn            
        self._slope_ = RMG.empty(centering = 'cell') # Slope - Rise/Run
        self._slope = RMG.empty(centering = 'cell')  # Slope - in Radians
        self._angles = RMG.empty(centering = 'cell') # Angles - Clckwise from N      
                       
        self._t = 12                                     # Assuming noon
        self._julian = floor( (self._current_time - floor( self._current_time))\
                                                        * 365 ) # Julian day

        print 'Julian Day = ', self._julian
                
        self._phi = pi/180.0 * self._latitude             # Latitude in Radians
        
        self._delta = 23.45 * pi/180.0 * cos(2*pi/365 * (172 - self._julian))  \
                                                            # Declination angle

        self._tau = (self._t + 12) * pi/12.0                   # Hour Angle

        self._alpha = asin(sin(self._delta) * sin(self._phi) + cos(self._delta)\
                        *cos(self._phi)*cos(self._tau))       # Solar Altitude

        if self._alpha <= 0.25*pi/180:          # If altitude is -ve, 
            self._alpha = 0.25*pi/180           # sun is beyond the horizon

        self._Rgl = (1 - self._A) * (1 - 0.65 * pow((self._N),2)) *            \
                 (self._Io*exp((-1) * self._n * (0.128 - 0.054 *               \
                    log10(1/sin(self._alpha)))*(1/sin(self._alpha))))
                    # Counting for Albedo, Cloudiness and Atmospheric turbidity

        self._phisun = atan(-sin(self._tau)/(tan(self._delta)*cos(self._phi)   \
                        - sin(self._phi)*cos(self._tau)))     # Sun's Azhimuth

        if ( self._phisun >= 0 and -sin(self._tau) <= 0 ):
            self._phisun = self._phisun + pi
            
        elif ( self._phisun <= 0 and -sin(self._tau) >= 0 ):
            self._phisun = self._phisun + pi
            

        self._flat = cos(atan(0)) * sin(self._alpha) + sin(atan(0)) *          \
                        cos(self._alpha) * cos(self._phisun - 0)
                                                       # flat surface reference 
             
        for i in range(0, RMG.number_of_cells):

            self._slope_[i], self._angles[i] =                                 \
                RMG.calculate_max_gradient_across_node( self._Z,               \
                  RMG.cell_node[i])

            self._slope[i] = atan(self._slope_[i])   
            self._BETA[i] = self._angles[i]*pi/180  
            
            self._Si[i] = (cos(self._slope[i]) * sin(self._alpha) +            \
                                sin(self._slope[i]) * cos(self._alpha)         \
                                  * cos(self._phisun - self._BETA[i]))
            
            if self._Si[i] <= 0:
                self._Si[i] = 0.0

            self._Rs[i] = self._Rgl * self._Si[i]    # Incoming Shortwave Radn
            

    def plott( self, RMG ):
        
        figure(1)        
        imshow_grid(RMG,self._Z, values_at = 'node')
        savefig('figure1')

        figure(2)
        self._slope_[self._slope_ >= 1] = 1
        imshow_active_cells(RMG, self._slope_, 'Slope', None, ('m','m'))
        savefig('figure2')

        figure(3)
        imshow_active_cells(RMG,self._Si,                                    \
                   'Cosine of Solar Angle of Incidence', None, ('m','m'))
        savefig('figure3') 
        
        figure(4)
        imshow_active_cells(RMG,self._BETA,                                    \
            'Aspect (clockwise from North)', 'radians', ('m','m'))
        savefig('figure4')
                
        figure(5)
        imshow_active_cells(RMG,self._Rs,                                    \
                   'Incoming Shortwave Radiation', 'W/m^2', ('m','m'))     
        savefig('figure5')
    
    def set_current_time( self, current_time ):
        
        self._current_time = current_time
