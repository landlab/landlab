########## soil_moisture.py ###############
##
##  This component calculates and updates soil moisture after each storm.
##  Soil moisture is represented as a single bucket. Rainfall depth fills the bucket
##  and the soil moisture decays as a result of leakage and ET following the analytical solution of
##  Laio et al., (2001) 
##  This component can operate on any raster grid
##  Input file is named soilmoisture_input.txt and is temporarily placed under landlab.components.
##  All input values can be altered using set_*** functions
##  All output arrays can be retrieved using get_*** functions
##
##  Storms are considered to be instantaneous events.
##  Storm duration, depth and interstorm duration are obtained from RainfallDriver.py
##  Storm depth is obtained in mm and storm duration and interstorm duration are obtained in hours
##
##  Sai Nudurupati and E.I. 11/17/2013.
###########################################

import os
from numpy import *
from math import *
import landlab
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from landlab.plot import imshow_grid, imshow_active_cells
from matplotlib.pyplot import *


_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'soilmoisture_input.in')

class SoilMoisture():

    def __init__( self ):

        """ vegcover - vegetation cover fraction """
        self._vegcover = 0.0
        """ interception_cap - Maximum interception capacity """
        self._interception_cap = 0.0
        """ zr - root depth """
        self._zr = 0.0
        """ runon - runon from upstream """
        self._runon = 0.0
        """ PET - potential evapotranspiration """
        self._PET = 0.0
        """ fbare - scaling coefficient for bare soil Potential Evapotranspiration """
        self._fbare = 0.0
        """ Ib - Bare soil infiltration capacity [mm/h] """
        self._soil_Ib = 0.0
        """ Iv - Vegetated surface infiltration capacity [mm/h] """
        self._soil_Iv = 0.0
        """ Ew - Residual Evaporation after wilting """
        self._soil_Ew = 0.0
        """ pc - porosity """
        self._soil_pc = 0.0
        """ fc - field capacity """
        self._soil_fc = 0.0
        """ sc - stomata closure moisture fraction """
        self._soil_sc = 0.0
        """ wp - wilting point """ 
        self._soil_wp = 0.0
        """ hgw - Hygroscopic water """
        self._soil_hgw = 0.0
        """ beta - Deep percolation rate constant """
        self._soil_beta = 0.0
     
        
        """ initialize """
        P = 0
        Tb = 0
        Tr = 0
        ETA = 0.0
        Water_Stress = 0.0
        
        """ sini - soil moisture with storm at the end of time intervals """
        sini = 0.0
        """ s - final soil moisture """
        s = 0.0
        """ Dd - Drainage from rootzone """
        Dd = 0.0
        """ Runoff - runoff after storm """
        Runoff = 0.0
        
        """ Iterations - Number of iterations/storms """
        self._iterate_storm = 0
        """ so - Soil moisture of the ground before storm """
        so = 0.0

        """ Initialize output variables  """
        
        self._Peff = 0.0        
        self._time = 0.0
        
       
    def initialize( self ):

        """ Read Input File """        
        MPD = ModelParameterDictionary()
        MPD.read_from_file(_DEFAULT_INPUT_FILE)
        
        """ Read Input Parameters """        
        self._iterate_storm = MPD.read_int( 'N_STORMS' )
        self._vegcover = MPD.read_float( 'VEG_COV' )
        self._interception_cap = MPD.read_float( 'INTERCEPT_CAP' )
        self._zr = MPD.read_float( 'ZR' )
        self._runon = MPD.read_float( 'RUNON' )
        self._fbare = MPD.read_float( 'F_BARE' )
        self._soil_Ib = MPD.read_float( 'I_B' )
        self._soil_Iv = MPD.read_float( 'I_V' )
        self._soil_Ew = MPD.read_float( 'EW' )
        self._soil_pc = MPD.read_float( 'PC' )
        self._soil_fc = MPD.read_float( 'FC' )
        self._soil_sc = MPD.read_float( 'SC' )
        self._soil_wp = MPD.read_float( 'WP' )
        self._soil_hgw = MPD.read_float( 'HGW' )
        self._soil_beta = MPD.read_float( 'BETA' )        
        

    def update( self, P, Tb, Tr, fveg, SO ):

        """
            Reading in 'P' storm depth, 'Tb' Interstorm duration, 'Tr' Storm
            duration, 'fveg' Vegetation factor, 'SO' Initial Soil Moisture
        """
        
        """
            Initializing parameters that are constant in space for a given Storm
        """
        
        V = fveg     # total veg cover used for interception
        fbare = self._fbare
        ZR = self._zr
        pc = self._soil_pc
        fc = self._soil_fc
        sc = self._soil_sc
        wp = self._soil_wp
        hgw = self._soil_hgw
        beta = self._soil_beta

        """
            Calculating parameters that are constant over space
        """
        
        Inf_cap = self._soil_Ib*(1-V) + self._soil_Iv*V # Infiltration capacity
        Int_cap = min(V*self._interception_cap, P)  # Interception capacity
        Peff = max(P-Int_cap, 0.0)             # Effective precipitation depth
        mu = (Int_cap/1000.0)/(pc*ZR*(exp(beta*(1-fc))-1))        
        
        """ 
            Initializing variables for spatial run for a given storm
        """
        
        (mm,nn) = SO.shape
        
        self._ETA = zeros((SO.shape))
        self._S = zeros((SO.shape))
        self._D = zeros((SO.shape))
        self._Ro = zeros((SO.shape))
        self._WS = zeros((SO.shape))

        """
            Looping over space for the time period (Tr+Tb)
        """
        
        """ Rows """
                
        for i in range(0,mm):
            
            """ Columns """    
            
            for j in range(0,nn):
                
                """
                    Calculating parameters that are variable over space
                """
                
                pet = self._PET[i][j]           # Potential Evapotranspiration
                Ep = max((pet*V+fbare*pet*(1-V)) - Int_cap, 0.0)  #
                nu = ((Ep/24.0)/1000.0)/(pc*ZR) # Loss function parameter
                nuw = ((Ep*0.1/24)/1000.0)/(pc*ZR) # Loss function parameter
                so = SO[i][j]               # Initial Soil Moisture
                sini = so + (Peff/(pc*ZR*1000.0))+self._runon       
        
                """
                    Evaluating Runoff and accordingly adjusting initial soil
                    moisture
                """          
                      
                if sini>1:
                    runoff = (sini-1)*pc*ZR*1000
                    sini = 1
                else:
                    runoff = 0

                """
                    Evaluating the current 'regime' of soil moisture (Laio et al.
                    Fig. 5) and calculating the evolution time and subsequently
                    calculating final soil moisture, drainage and actual 
                    evapotranspiration (Analytical solution, Laio et al.)
                """
        
                if sini>=fc:
                    tfc = (1.0/(beta*(mu-nu)))*(beta*(fc-sini)+                \
                             log((nu-mu+mu*exp(beta*(sini-fc)))/nu))

                    tsc = ((fc-sc)/nu)+tfc
                    twp = ((sc-wp)/(nu-nuw))*log(nu/nuw)+tsc
        
                    if Tb<tfc:
                        s = abs(sini-(1/beta)*log(((nu-mu+mu*                  \
                             exp(beta*(sini-fc)))*exp(beta*(nu-mu)*Tb)         \
                             -mu*exp(beta*(sini-fc)))/(nu-mu)))

                        Dd = ((pc*ZR*1000)*(sini-s))-(Tb*(Ep/24))
                        ETA = (Tb*(Ep/24))
        
                    elif Tb>=tfc and Tb<tsc:
                        s = fc-(nu*(Tb-tfc))
                        Dd = ((pc*ZR*1000)*(sini-fc))-((tfc)*(Ep/24))
                        ETA = (Tb*(Ep/24))
        
                    elif Tb>=tsc and Tb<twp:
                        s = wp+(sc-wp)*((nu/(nu-nuw))*exp((-1)*((nu-nuw)/(sc-wp))*(Tb-tsc))-(nuw/(nu-nuw)))
                        Dd = ((pc*ZR*1000)*(sini-fc))-(tfc*Ep/24)
                        ETA = (1000*ZR*pc*(sini-s))-Dd
        
                    else:
                        s = hgw+(wp-hgw)*exp((-1)*(nuw/(wp-hgw))*max(Tb-twp,0))
                        Dd = ((pc*ZR*1000)*(sini-fc))-(tfc*Ep/24)
                        ETA = (1000*ZR*pc*(sini-s))-Dd
        
                elif sini<fc and sini>=sc:
                    tfc = 0
                    tsc = (sini-sc)/nu
                    twp = ((sc-wp)/(nu-nuw))*log(nu/nuw)+tsc
        
                    if Tb<tsc:
                        s = sini - nu*Tb
                        Dd = 0
                        ETA = 1000*ZR*pc*(sini-s)
        
                    elif Tb>=tsc and Tb<twp:
                        s = wp+(sc-wp)*((nu/(nu-nuw))*exp((-1)*((nu-nuw)/(sc-wp))*(Tb-tsc))-(nuw/(nu-nuw)))
                        Dd = 0
                        ETA = (1000*ZR*pc*(sini-s))
        
                    else:
                        s = hgw+(wp-hgw)*exp((-1)*(nuw/(wp-hgw))*(Tb-twp))
                        Dd = 0
                        ETA = (1000*ZR*pc*(sini-s))
                        
                elif sini<sc and sini>=wp:
                    tfc = 0
                    tsc = 0
                    twp = ((sc-wp)/(nu-nuw))*log(1+(nu-nuw)*(sini-wp)/(nuw*(sc-wp)))
        
                    if Tb<twp:
                        s = wp+((sc-wp)/(nu-nuw))*((exp((-1)*((nu-nuw)/(sc-wp))*Tb))*(nuw+((nu-nuw)/(sc-wp))*(sini-wp))-nuw)
                        Dd = 0
                        ETA = (1000*ZR*pc*(sini-s))
        
                    else:
                        s = hgw+(wp-hgw)*exp((-1)*(nuw/(wp-hgw))*(Tb-twp))
                        Dd = 0
                        ETA = (1000*ZR*pc*(sini-s))
        
                else:
                    tfc = 0
                    tsc = 0
                    twp = 0
        
                    s = hgw+(sini-hgw)*exp((-1)*(nuw/(wp-hgw))*Tb)
                    Dd = 0
                    ETA = (1000*ZR*pc*(sini-s))
        
                Water_Stress = min(max(pow(((sc - ((s+sini)/2)) / (sc - wp)),4),0.0),1.0) 
            
                """ varaibles for get/set functions """
                self._Peff = Peff
                self._ETA[i][j] = ETA
                self._S[i][j] = s
                self._D[i][j] = Dd
                self._Ro[i][j] = runoff                
                self._WS[i][j] = Water_Stress
                #self._so[i][j] = s 
        self._time = self._time + (Tb + Tr)/(24*365.4)  

    """ Return Effective precipittion mm """
    def get_Peff( self ):
        return self._Peff
    
    """ Return Current Time """  
    def get_time( self ):
        return self._time
 
    """ Return Interstorm ET mm """
    def get_ET( self ):
        return self._ETA

    """ Return Soil Moisture  """
    def get_S( self ):
        return self._S

    """ Return Drainage from Root Zone mm """
    def get_D( self ):
        return self._D

    """ Return Runoff mm """
    def get_Ro( self ):
        return self._Ro
    
    """ Return Water Stress """
    def get_WaterStress( self ):
        return self._WS

    """ Set Number of Storms """
    def set_strms( self, strms ):
        self._iterate_storm = strms        
        
    """ Set Potential Evapotranspiration Rate """
    def set_PET( self, PET ):
        self._PET = PET


    def plott( self, P_, S_, ETA_, D_, Ro_, Time_, RMG ):

        figure(1)
        plot( Time_, P_ )
        ylabel( 'Storm Depth [mm]' )
        title( 'Storm generation' )
        xlabel( 'Time:Year' )        

        figure(2)
        imshow_grid(RMG, S_, values_at = 'cell') 

        show()
