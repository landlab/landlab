########## soil_moisture.py ###############
##
##  This component calculates and updates soil moisture after each storm.
##  Multiple storms can be considered which is fed as an input.
##  This code implements Laio's(2001) solution for soil moisture calculation.
##  Input file is named soilmoisture_input.txt and is temporarily placed under landlab.components.
##
##  Storms are considered to be instantaneous events.
##  Storm duration, depth and interstorm duration are obtained from RainfallDriver.py
##  Storm depth is obtained in mm and storm duration and interstorm duration are obtained in hours
##
##  After every storm, Laio's solution for soil moisture drying is implemented. 
##
##
## Written by Erkan Istanbulluoglu & Sai Nudurupati, 2013.
###########################################

import os
from numpy import *
from math import *
import landlab
from landlab.components.RainfallDriver import PrecipitationDistribution
from landlab.model_parameter_dictionary import ModelParameterDictionary
from matplotlib.pyplot import *


_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'soilmoisture_input.txt')

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

        """ soil_type - Soil type to be selected """
        self._soil_type = 1
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
        
        """ P - Precipitation """
        P = 0.0
        """ ETA - Evapotranspiration """
        ETA = 0.0
        """ sini - soil moisture with storm at the end of time intervals """
        sini = 0.0
        """ s - final soil moisture """
        s = 0.0
        """ Dd - Drainage from rootzone """
        Dd = 0.0
        """ Tb - Time between storms """
        Tb = 0.0
        """ Tr - Storm duration """
        Tr = 0.0
        """ Runoff - runoff after storm """
        Runoff = 0.0

        """ Iterations - Number of iterations/storms """
        self._iterate_storm = 0
        """ so - Soil moisture of the ground before storm """
        so = 0.0

    def initialize( self ):

        MPD = ModelParameterDictionary()
        MPD.read_from_file(_DEFAULT_INPUT_FILE)

        self._iterate_storm = MPD.read_int( 'strms' )
        print '\nNumber of storms: ', self._iterate_storm
        self._vegcover = MPD.read_float( 'VegCover' )
        print '\n Veg cover: ', self._vegcover
        self._interception_cap = MPD.read_float( 'InterceptionCap' )
        print '\n Interception cap', self._interception_cap
        self._zr = MPD.read_float( 'ZR' )
        print '\n Root depth', self._zr
        self._runon = MPD.read_float( 'Runon' )
        print '\n Initial runoff:', self._runon
        self._PET = MPD.read_float( 'PET' )
        print '\n Fixed PET: ', self._PET
        self._fbare = MPD.read_float( 'F_Bare' )
        print '\n F_bare: ', self._fbare
        self._so = MPD.read_float( 'so' )        
        self._si = self._so
        print '\n Initial soil moisture: ', self._so
        
        self._soil_type = MPD.read_float( 'SoilType' )
        self._soil_Ib = MPD.read_float( 'Ib' )
        print '\n Ib: ', self._soil_Ib
        self._soil_Iv = MPD.read_float( 'Iv' )
        print '\n Iv: ', self._soil_Iv
        self._soil_Ew = MPD.read_float( 'Ew' )
        print '\n Ew: ', self._soil_Ew
        self._soil_pc = MPD.read_float( 'pc' )
        print '\n pc: ', self._soil_pc
        self._soil_fc = MPD.read_float( 'fc' )
        print '\n fc: ', self._soil_fc
        self._soil_sc = MPD.read_float( 'sc' )
        print '\n sc: ', self._soil_sc
        self._soil_wp = MPD.read_float( 'wp' )
        print '\n wp: ', self._soil_wp
        self._soil_hgw = MPD.read_float( 'hgw' )
        print '\n hgw: ', self._soil_hgw
        self._soil_beta = MPD.read_float( 'Beta' )
        print '\n Beta: ', self._soil_beta
        
        self._length = self._iterate_storm
        """ Initializing arrays to store values of variables """
        self._Pin = zeros( self._length, dtype = float )
        self._Pa = zeros( self._length, dtype = float )
        self._P = zeros( self._length, dtype = float )
        self._ETA = zeros( self._length, dtype = float )
        self._Sini = zeros( self._length, dtype = float )
        self._S = zeros( self._length, dtype = float )
        self._D = zeros( self._length, dtype = float )
        self._Tb = zeros( self._length, dtype = float )
        self._Tr = zeros( self._length, dtype = float )
        self._Ro = zeros( self._length, dtype = float )
        self._tfc = zeros( self._length, dtype = float )
        self._tsc = zeros( self._length, dtype = float )
        self._twp = zeros( self._length, dtype = float )

        self._Total = zeros( self._length, dtype = float )
        self._percentage = zeros( self._length, dtype = float )

    def update( self ):

        V = self._vegcover
        fb = self._fbare
        ZR = self._zr
        fbare = self._fbare
        PET = self._PET
        pc = self._soil_pc
        fc = self._soil_fc
        sc = self._soil_sc
        wp = self._soil_wp
        hgw = self._soil_hgw
        beta = self._soil_beta
           
        """ Create an instance of Storm Class """
        PD = PrecipitationDistribution()
        """ Initialize Storm - Create first storm """
        PD.initialize()
        
        for i in range( 0, self._iterate_storm ):
            
           if i != 0:
               """ Create a new storm """ 
               PD.update()

           """ Soil moisture dynamics """
           P = PD.storm_depth
           Pin = PD.intensity
           Tb = PD.interstorm_duration
           Tr = PD.storm_duration                      
           
           Inf_cap = self._soil_Ib*(1-V) + self._soil_Iv*V
           Int_cap = min(V*self._interception_cap, P)
           Peff = max(P-Int_cap, 0.0)
           mu = (Int_cap/1000.0)/(pc*ZR*(exp(beta*(1-fc))-1))
           Ep = max((PET*V+fbare*PET*(1-V)) - Int_cap, 0.0)
           nu = ((Ep/24.0)/1000.0)/(pc*ZR)
           nuw = ((Ep*0.1/24)/1000.0)/(pc*ZR)

           sini = self._so + (Peff/(pc*ZR*1000.0))+self._runon

           if sini>1:
               runoff = (sini-1)*pc*ZR*1000
               sini = 1
           else:
               runoff = 0

           if sini>=fc:
                tfc = (1.0/(beta*(mu-nu)))*(beta*(fc-sini)+log((nu-mu+mu*exp(beta*(sini-fc)))/nu))
                tsc = ((fc-sc)/nu)+tfc
                twp = ((sc-wp)/(nu-nuw))*log(nu/nuw)+tsc

                if Tb<tfc:
                    s = abs(sini-(1/beta)*log(((nu-mu+mu*exp(beta*(sini-fc)))*exp(beta*(nu-mu)*Tb)-mu*exp(beta*(sini-fc)))/(nu-mu)))
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

           self._Pa[i] = P
           self._Pin[i] = Pin
           self._P[i] = Peff
           self._ETA[i] = ETA
           self._Sini[i] = sini
           self._S[i] = s
           self._D[i] = Dd
           self._Tb[i] = Tb
           self._Tr[i] = Tr
           self._Ro[i] = runoff
           self._tfc[i] = tfc
           self._tsc[i] = tsc
           self._twp[i] = twp

           self._so = s

    def plott( self ):

        figure(1)
        self._px = xrange( 0, self._iterate_storm )
        plot( self._px, self._P )
        ylabel( 'Effective Precipitation Depth in mm' )
        title( 'Precipitation Record' )
        xlabel( 'Storms' )

        figure(2)
        plot( self._px, self._S )
        ylabel( 'Soil Moisture - fraction' )
        title( 'Soil Moisture Dynamics' )
        xlabel( 'Storms' )

        figure(3)
        plot( self._px, self._ETA, '+', label = 'Evapotranspiration' )
        hold(True)
        plot( self._px, self._D, '-', label = 'Drainage from Root Zone' )
        hold(True)
        plot( self._px, self._Ro, '*', label = 'Runoff' )
        hold(False)
        xlabel( 'Storms' )
        ylabel( 'Depth in mm' )
        legend()
        title( 'Leakage/Losses' )
                
        show()



        
