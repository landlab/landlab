import numpy as np
from photosynthesis import *
import pandas as pd

day = 150
latitude = 39*np.pi/180 # in radians
PAR = 1124.98 # day 100 in 2010 on POWER_Point_Daily_20091231_20151230_039d81N_074d09W_LST.csv
increment_hour = 12 # or can use the below for loop
last_biomass = pd.DataFrame({'root_biomass':[3.0], 'stem_biomass':[2.0], 'leaf_biomass':[1.5]})
lai=3

p = Photosynthesis(latitude, 100) 
# Solar Variables
p.update_solar_variables(day)
solar_elevation = p.calculate_solar_elevation(increment_hour)
print("Declination="+str(p._solar_declination))
print("Sunrise="+str(p._sunrise))
print("Sunset="+str(p._sunset))
print("Daylength="+str(p._sunlit_increment))
print("Solar elevation at "+str(increment_hour)+"="+str(solar_elevation))
print("Incremental PAR="+str(p.calculate_incremental_PAR(increment_hour, PAR))) 
print("Hourly direct light extinction="+str(p.calculate_hourly_direct_light_extinction(solar_elevation)))
print("Hourly diffuse light extinction="+str(p.calculate_hourly_diffuse_light_extinction(lai)))
print("Sunlit/Shaded LAI="+str(p.calculate_sunlit_shaded_LAI_proportion(solar_elevation,lai)))
print("Absorbed Incremental PAR="+str(p.calculate_absorbed_incremental_PAR(increment_hour, solar_elevation, PAR, lai)))
print("Photosynthesis="+str(p.photosynthesize(PAR, last_biomass, lai, day)))

#for day_increment in p.gauss_integration_params:
#    (abscissa,weight)=day_increment
#    increment_hour = p._sunrise+abscissa*p._sunlit_increment
#    print(increment_hour)