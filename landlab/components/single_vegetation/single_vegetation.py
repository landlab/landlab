######### single_vegetation.py ###########
##
## 
## This code is an 'In Progress' application
## of vegetation dynamics at a point.
##
##  04 Aug 2013 - EI & SN
#############################################

import os
from numpy import *
from math import *
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from soil_moisture import SoilMoisture
from landlab import ModelParameterDictionary
from matplotlib.pyplot import *

_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'single_vegetation.in')

class SingleVegetation ():

    def __init__( self ):

        """ List of model input parameters initialized """

        """ vegetation water use efficiency """
        self._WUE = 0.0
        """ maximum vegetation LAI sum live and dead biomass """
        self._LAI_max = 0.0
        """ specific leaf area for live biomass """
        self._cb = 0.0
        """ specific leaf area for dead biomass """
        self._cd = 0.0
        """ live biomass senescence coefficient """
        self._ksg = 0.0
        """ dead biomass decay coefficient """
        self._kdd = 0.0
        
        """ List of model varaibles from other components initialized """
    
        """ Tb - Time between storms """
        Tb = 0.0
        """ ETA - Evapotranspiration """
        ETA = 0.0

        """ List of model states initialized """
        """ Net Primary Productivity """
        NPP = 0.0
        """ Initial Live biomass """
        Blive_ini = 0.0
        """ Initial Dead biomass """
        Bdead_ini = 0.0
        """ Live biomass """
        Blive = 0.0
        """ Dead biomass """
        Bdead = 0.0
        """ Maximum biomass """
        Bmax = 0.0 
        """ Live LAI """
        LAIlive = 0.0
        """ Dead LAI """
        LAIdead = 0.0
        """ Constant for analytical solution """
        Yconst = 0.0
        """ Total vegeattion cover """
        Vt = 0.0

        self._LAIlive = 0.0
        self._LAIdead = 0.0
        self._Vt = 0.0
        self._Blive_ini = 0.0
        self._Bdead_ini = 0.0
    
        
        """ Model iterations based on the number of storms """
        self._iterate_storm = 0

    def initialize( self ):

        MPD = ModelParameterDictionary()
        MPD.read_from_file(_DEFAULT_INPUT_FILE)
        
        self._vegcover = MPD.read_float( 'VEG_COV' )
        self._WUE = MPD.read_float( 'WUE' )
        self._LAI_max = MPD.read_float( 'LAI_MAX' )
        self._cb = MPD.read_float( 'CB' )
        self._cd = MPD.read_float( 'CD' )
        self._ksg = MPD.read_float( 'KSG' )
        self._kdd = MPD.read_float( 'KDD' )
        self._Blive_ini = MPD.read_float( 'BLIVE_INI' )
        self._Bdead_ini = MPD.read_float( 'BDEAD_INI' )
        

    def update( self, P, Tb, Tr, ETT ): # Please note that an extra variable ETT
        # to accomodate data from SM.get_ET() is added as extra input

        LAIdead = self._LAIdead # retrieving data from previous call
        LAIlive = self._LAIlive # Retrieving data from previous call
        
        NPP = ETT*self._WUE*24*0.55*1000
        Bmax = (self._LAI_max - LAIdead)/self._cb # LAId = LAIdead ???
        Yconst = (1/((1/Bmax)+(self._ksg/NPP)))
        Blive = (self._Blive_ini - Yconst) *exp(-(NPP/Yconst)*(Tb/24)) + Yconst
        Bdead = (self._Bdead_ini + (Blive - max(Blive * exp(-self._ksg * Tb/24),0.0001)))*exp(-self._kdd * Tb/24)

        LAIlive =  min (self._cb * Blive, self._LAI_max)
        LAIdead =  min (self._cd * Bdead, self._LAI_max - LAIlive)
        Vt = 1 - exp(-0.75*(-0.75 * (LAIlive + LAIdead)))

        self._LAIlive = LAIlive
        self._LAIdead = LAIdead
        self._Vt = Vt

        self._Blive_ini = Blive
        self._Bdead_ini = Bdead

        
 
    """ Return vegetation variables """
    def get_LAIlive( self ):
        return self._LAIlive

    def get_LAIdead( self ):
        return self._LAIdead

    def get_VegCov( self ):
        return self._Vt
