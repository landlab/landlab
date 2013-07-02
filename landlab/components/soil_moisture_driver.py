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
""" Update Soil Moisture """
SM.update()


""" Calculate RHS of Water Balance """ 
Total = (SM._S[SM._iterate_storm-1]-SM._si)*1000*SM._soil_pc*SM._zr+sum(SM._ETA)+sum(SM._D)+sum(SM._Ro)
""" Calculate percentage error in LHS and RHS of water balance equation """
Percent = (sum(SM._P)-Total)/(sum(SM._P))*100
print '\nPercent Error: ', Percent, '\n'


""" Plot output """
SM.plott()
