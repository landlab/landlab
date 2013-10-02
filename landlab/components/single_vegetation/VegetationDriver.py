######### VegetationDriver.py ###########
##
## 04 Aug 2013 - EI & SN
## 
## This code is an example for implementation of 
## single_vegetation
##
#############################################

# vegetation driver 

from single_vegetation import *
from soil_moisture import *


""" Create a Soil Moisture Object """
SM = SoilMoisture()
""" Initialize the Soil Moisture """
SM.initialize()

""" Create a Vegetation Object """
VEG = SingleVegetation()
""" Initialize Vegetation """
VEG.initialize()


SM._length = SM._iterate_storm

""" Create an instance of Storm Class """
PD = PrecipitationDistribution()
""" Initialize Storm - Create first storm """
PD.initialize()

#######

""" Loop through the storms """
for i in range( 0, SM._iterate_storm ):
    if i != 0:
       """ Create a new storm """ 
       PD.update()
    """ Soil moisture dynamics """
    P = PD.storm_depth
    Tb = PD.interstorm_duration
    Tr = PD.storm_duration

    """ Update Soil Moisture """
    SM.update( P, Tb, Tr )

    ETT = SM.get_ET()  # Since we are not sending the object SM to VEG.update,
    # we need to send the information of  SM.get_ET() as another variable - here ETT
    
    VEG.update(P, Tb, Tr, ETT ) # Note that ETT (variable - denoting
    # data from SM.get_ET()) 
