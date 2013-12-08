########## ET_PriestlyTaylor.py ###############
##
##  This component calculates and updates Priestly Taylor Evapotranspiration.
##
##  Equations from ACSE EWRI report (Jan 2005) are used
##  Inputs are read from the file 'PET_input.txt'
##  Input observed data from .mat file 'Bondville_IL.mat' is also
##  read for comparisions with modeled data.
##
##
## Written by Erkan Istanbulluoglu & Sai Nudurupati, 2013.
###########################################

import os
from numpy import * 
from landlab import ModelParameterDictionary
from math import *
from matplotlib.pyplot import *

_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'PET_input.txt')

_DEFAULT_INPUT_FILE_1 = os.path.join(os.path.dirname(__file__),
                                  'Bondville_IL.mat')

class ET:

    def __init__( self ):
        pass

    def initialize( self ):

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

        """ Reading .mat file for comparision of simulations with observed data.
            Observed data is included in .mat file """
        
        MPD.read_from_mat( _DEFAULT_INPUT_FILE_1 )
        """ Day of the year """
        self._Dm = MPD.data['Dm']
        """ Month of the year """
        self._M = MPD.data['M']
        """ Year """
        self._Y = MPD.data['Y']
        """ Maximum temperature of the day """
        self._Tmax = MPD.data['Tmax']
        """ Minimum temperature of the day """
        self._Tmin = MPD.data['Tmin']        
        """ Observed Short Wave Radiation """
        self._RsObs = MPD.data['RsObs']
        
        self._length = len(self._RsObs)

        """ Initialize Arrays """
        """ Saturation Vapor Pressure Gradient """
        self._delta = zeros( self._length, dtype = float )
        """ Saturation Vapor Pressure """
        self._es = zeros( self._length, dtype = float )
        """ Actual Vapor Pressure """
        self._ea = zeros( self._length, dtype = float )
        """ Net Long Wave Radiation """
        self._Rnl = zeros( self._length, dtype = float )
        """ Net Short Wave Radiation """
        self._Rns = zeros( self._length, dtype = float )
        """ Net Radiation """
        self._Rn = zeros( self._length, dtype = float )
        """ Clear Sky Solar Radiation """
        self._Rso = zeros( self._length, dtype = float )
        """ Short Wave Radiation """
        self._Rs = zeros( self._length, dtype = float )        
        self._u = zeros( self._length, dtype = float )
        self._ETp = zeros( self._length, dtype = float )
        """ Extraterrestrial Radiation """
        self._Ra = zeros( self._length, dtype = float )
        """ Julian Day """
        self._J = zeros( self._length, dtype = float )
        self._ws = zeros( self._length, dtype = float )
        self._x = zeros( self._length, dtype = float )
        """ Solar Declination Angle """
        self._sdecl = zeros( self._length, dtype = float )
        self._dr = zeros( self._length, dtype = float )
        self._Tavg = zeros( self._length, dtype = float )
        self._fcd = zeros( self._length, dtype = float )
        self._phi = (3.14/180)*self._lat

           
    def update( self ):

        for i in xrange( 0, self._length ):
            """
                Julian Day - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (25)
            """
            self._J[i] = self._Dm[i] - 32.0 + round(275.0*(self._M[i]/9.0)) + \
                        (2*round(3/(self._M[i]+1))) + round((self._M[i]/100.0)- \
                        (self._Y[i]%4.0)/4.0+0.975)
            """
                Saturation Vapor Pressure - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (37)
            """
            self._es[i] = 0.6108*exp((17.27*self._Tavg[i])/(237.7+self._Tavg[i]))
            """
                Actual Vapor Pressure - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (38)
            """
            self._ea[i] = 0.6108*exp((17.27*self._Tmin[i])/(237.7+self._Tmin[i]))
            """
                Slope of Saturation Vapor Pressure - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (36)
            """
            self._delta[i] = (4098.0*self._es[i])/(pow((237.3+self._Tavg[i]),2.0))
            """
                Solar Declination Angle - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (51)
            """
            self._sdecl[i] = 0.409*sin((((2.0*3.14)/365.0)*self._J[i])-1.39)
            """
                Inverse Relative Distance Factor - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (50)
            """
            self._dr[i] = 1 + (0.033*cos((2.0*3.14/365.0)*self._J[i]))
            """
                To calculate ws - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (61)
            """
            self._x[i] = 1.0-((pow((tan(self._phi)),2.0))*(pow(tan(self._sdecl[i]),2.0)))  
            if self._x[i] <= 0:
                self._x[i] = 0.00001;
                
            """
                Sunset Hour Angle - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (60)
            """
            self._ws[i] = (3.14/2.0)-atan((-1*tan(self._phi)*tan(self._sdecl[i]))/(pow(self._x[i],2.0)))
            """
                Extraterrestrial radmodel.docx - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (48)
            """
            self._Ra[i] = 11.57*(24.0/3.14)*4.92*self._dr[i]*((self._ws[i]*sin(self._phi)*sin(self._sdecl[i]))+ \
                          (cos(self._phi)*cos(self._sdecl[i])*(sin(self._ws[i]))))
            """
                Clear-sky Solar Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (47)
            """
            self._Rso[i] = (0.75+((2.0*pow(10,-5.0))*self._z))*self._Ra[i]
            self._Rs[i] = min(self._Krs*self._Ra[i]*sqrt(self._Tmax[i]-self._Tmin[i]), self._Rso[i])
            """
                Net Short Wave Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (43)
            """
            self._Rns[i] = self._Rs[i]*(1-self._a)

            """
                Relative Cloudiness - ASCE-EWRI Task Committee Report, Jan-2005 - Page 35
            """
            if self._Rso[i] > 0:
                self._u[i] = self._Rs[i]/self._Rso[i]
            else:
                self._u[i] = 0

            if self._u[i] < 0.3:
                self._u[i] = 0.3
            elif self._u[i] > 1:
                self._u[i] = 1.0

            """
                Cloudiness Function - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (45)
            """
            self._fcd[i] = (1.35*self._u[i])-0.35
            """
                Net Long Wave Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (44)
            """
            self._Rnl[i] = self._sigma*self._fcd[i]*(0.34-(0.14*sqrt(self._ea[i]))* \
                           ((pow((self._Tmax[i]+273.16),4.0)+ \
                            pow((self._Tmin[i]+273.16),4.0))/2.0))

            """
                Net Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (42)
            """                                         
            self._Rn[i] = self._Rns[i] - self._Rnl[i]

            self._ETp[i] = max(self._alpha*(self._delta[i]/(self._delta[i]+self._y) \
                           )*(self._Rn[i]/self._pwhv), 0)                                                                                                       


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
        ylabel( 'Evapotranspiration PET' )
        title( 'Evapotranspiration' )
        show()

