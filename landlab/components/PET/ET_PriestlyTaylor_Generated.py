########## ET_PriestlyTaylor.py ###############
##
##  This component calculates and updates Priestly Taylor Evapotranspiration.
##
##  Equations from ACSE EWRI report (Jan 2005) are used
##  Inputs are read from the file 'PET_input.txt'
##  
##  
##
## Written by Sai Nudurupati & E.I. Nov 2013.
###########################################

import os
from numpy import * 
from landlab import ModelParameterDictionary
from math import *
from matplotlib.pyplot import *

""" Default input file """
_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'PET_input.txt')


class ET():

    def __init__( self ):
        pass

    def initialize( self ):

        """ Reading data from input file """
        MPD = ModelParameterDictionary()
        MPD.read_from_file( _DEFAULT_INPUT_FILE )

        """ Priestly Taylor Constant """
        self._alpha = MPD.read_float( 'alpha' )
        """ Surface Albedo """
        self._a = MPD.read_float( 'a' )
        """ Latent Heat of Vaporization """
        self._pwhv = MPD.read_float( 'pwhv' )
        """ Psychometric Constant """
        self._y = MPD.read_float( 'y' )
        """ Stefan Boltzmann Constant """
        self._sigma = MPD.read_float( 'sigma' )
        """ Solar Constant """
        self._Gsc = MPD.read_float( 'Gsc' )
        """ Latitude of Location under Study """
        self._lat = MPD.read_float( 'lat' )
        """ Elevation of Location under Study """
        self._z = MPD.read_float( 'z' )
        """ Adjustment Coefficient """
        self._Krs = MPD.read_float( 'Krs' )        

        """ Initialize Arrays """

        """ Net Short Wave Radiation """
        self._Rns = 0
        """ Net Radiation """
        self._Rn = 0
        """ Clear Sky Solar Radiation """
        self._Rso = 0
        """ Short Wave Radiation """
        self._Rs = 0
        """ Relative cloudiness factor """
        self._u = 0
        """ Priestly Taylor evapotranspiration """
        self._ETp = 0
        """ Extraterrestrial Radiation """
        self._Ra = 0
        """ Julian Day """
        self._J = 0
        """ Sunset hour angle """
        self._ws = 0
        self._x = 0
        """ Solar declination angle """
        self._sdecl = 0
        """ Inverse relative distance factor """
        self._dr = 0
        """ Cloudiness function """
        self._fcd = 0
        """ Latitude in radians """
        self._phi = (3.14/180)*self._lat

           
    def update( self, TIME ):

        """
            Julian Day - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (25)
        """        
        self._J = floor( (TIME - floor( TIME)) * 365 ) 

        """
            Solar Declination Angle - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (51)
        """
        self._sdecl = 0.409*sin((((2.0*3.14)/365.0)*self._J)-1.39)

        """
            Inverse Relative Distance Factor - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (50)
        """
        self._dr = 1 + (0.033*cos((2.0*3.14/365.0)*self._J))

        """
            To calculate ws - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (61)
        """
        self._x = 1.0-((pow((tan(self._phi)),2.0))*(pow(tan(self._sdecl),2.0)))  
        if self._x <= 0:
            self._x = 0.00001;
            
        """
            Sunset Hour Angle - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (60)
        """
        self._ws = (3.14/2.0)-atan((-1*tan(self._phi)*tan(self._sdecl))/(pow(self._x,2.0)))

        """
            Extraterrestrial radmodel.docx - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (48)
        """
        self._Ra = 11.57*(24.0/3.14)*4.92*self._dr*((self._ws*sin(self._phi)*sin(self._sdecl))+ \
                        (cos(self._phi)*cos(self._sdecl)*(sin(self._ws))))

        """
            Clear-sky Solar Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (47)
        """
        self._Rs = (0.75+((2.0*pow(10,-5.0))*self._z))*self._Ra

        """
            Net Short Wave Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (43)
        """
        self._Rns = self._Rs*(1-self._a)

        """
            Relative Cloudiness - ASCE-EWRI Task Committee Report, Jan-2005 - Page 35
        """
        self._u = 0.3
        
        """
            Cloudiness Function - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (45)
        """
        self._fcd = (1.35*self._u)-0.35

        """
            Net Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (42)
        """                                         
        self._Rn = self._Rns 

        """
            Priestly Taylor Evapotranspiration in mm/day - Priestly et. al. 1972
        """
        self._ETp = max(self._alpha*(self._delta/(self._delta+self._y) \
                        )*(self._Rn/self._pwhv), 0)                                                                                                       


    def plott( self ):

        figure(1)
        self._px = xrange( 1, self._length+1)
        plot( self._px, self._Rs, label = 'Modeled' )
        hold(True)
        plot( self._px, self._RsObs, '*', label = 'Observed' )
        hold(False)
        xlabel( 'Days' )
        ylabel( 'Short Wave Radiation' )
        title( 'Shortwave Radiation' )
        legend()
        

        figure(2)
        plot( self._px, self._ETp )
        xlabel( 'Days' )
        ylabel( 'Evapotranspiration PET (in mm/day)' )
        title( 'Evapotranspiration' )
        show()

