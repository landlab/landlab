######### soil_moisture_driver.py ###########
##
## This component can be considered as an example for implementation
## of soil_moisture.py
##
##
##
#############################################

# Soil Moisture Driver

from soil_moisture import *

""" Create an Soil Moisture Object """
SM = SoilMoisture()

""" Initialize Soil Moisture """
SM.initialize()

SM._length = SM._iterate_storm

""" Create an instance of Storm Class """
PD = PrecipitationDistribution()
""" Initialize Storm - Create first storm """
PD.initialize()

""" Initializing OPTIONAL arrays to store values of variables and plotting """
P_ = zeros( SM._length, dtype = float )
Peff_ = zeros( SM._length, dtype = float )
ETA_ = zeros( SM._length, dtype = float )
Time_ = zeros( SM._length, dtype = float )
S_ = zeros( SM._length, dtype = float )
D_ = zeros( SM._length, dtype = float )
Ro_ = zeros( SM._length, dtype = float )

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
    Fveg = 0.5 #self._soil_Fveg

    """ Update Soil Moisture """
    SM.update( P, Tb, Tr, Fveg)
    
    P_[i] = P
    Peff_[i] = SM.get_Peff()  # Example of using a get function
    S_[i] = SM.get_S()
    ETA_[i] = SM.get_ET()     # SM._ETA  Another example of using data stored in the object
    Ro_[i] = SM.get_Ro()
    D_[i] = SM._D             #SM.get_D() also works
    Time_[i] = SM.get_time()    

#######

""" Calculate RHS of Water Balance """
Total = (S_[SM._iterate_storm-1]-SM._si)*1000*SM._soil_pc*SM._zr+sum(ETA_)+sum(D_)+sum(Ro_)
""" Calculate percentage error in LHS and RHS of water balance equation """
Percent = (sum(Peff_)-Total)/(sum(Peff_))*100

print '\nPercent Error: ', Percent, '\n'

""" Plot output """
SM.plott( P_, S_, ETA_, D_, Ro_, Time_ )
