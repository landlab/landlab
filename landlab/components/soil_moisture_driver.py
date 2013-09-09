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


def main():
    """ Create an Soil Moisture Object """
    SM = SoilMoisture()
    """ Initialize Soil Moisture """
    SM.initialize()
    SM.set_strms(50000)   # Example of using a set function

    SM._length = SM._iterate_storm
    """ Initializing arrays to store values of variables """
    Pa_ = zeros( SM._length, dtype = float )
    P_ = zeros( SM._length, dtype = float )
    ETA_ = zeros( SM._length, dtype = float )
    Sini_ = zeros( SM._length, dtype = float )
    S_ = zeros( SM._length, dtype = float )
    D_ = zeros( SM._length, dtype = float )
    Ro_ = zeros( SM._length, dtype = float )

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
        Pin = PD.intensity
        Tb = PD.interstorm_duration
        Tr = PD.storm_duration

        """ Update Soil Moisture """
        SM.update( P, Pin, Tb, Tr )

        Pa_[i] = Pin
        P_[i] = SM.get_Peff() # Example of using a get function
        Sini_[i] = SM._Sini   # Example of using data stored in the object
        S_[i] = SM.get_S()
        ETA_[i] = SM._ETA     # Another example of using data stored in the object
        Ro_[i] = SM.get_Ro()
        D_[i] = SM.get_D()
        
    #######

    """ Calculate RHS of Water Balance """ 
    Total = (S_[SM._iterate_storm-1]-SM._si)*1000*SM._soil_pc*SM._zr+sum(ETA_)+sum(D_)+sum(Ro_)
    """ Calculate percentage error in LHS and RHS of water balance equation """
    Percent = (sum(P_)-Total)/(sum(P_))*100
    print '\nPercent Error: ', Percent, '\n'


    """ Plot output """
    SM.plott( P_, S_, ETA_, D_, Ro_ )

if __name__ == '__main__':
    main()
