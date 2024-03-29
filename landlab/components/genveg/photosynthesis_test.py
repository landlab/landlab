import numpy as np
from photosynthesis import Photosynthesis
import pandas as pd

day = 150
latitude = 39 * np.pi / 180  # in radians
PAR = (
    np.ones(1) * 1124.98
)  # day 100 in 2010 on POWER_Point_Daily_20091231_20151230_039d81N_074d09W_LST.csv
increment_hour = 12  # or can use the below for loop
last_biomass = pd.DataFrame(
    {
        "root_biomass": [3.0],
        "stem_biomass": [2.0],
        "leaf_biomass": [1.5],
        "shoot_sys_width": [0.5],
    }
)
lai = 3
min_temp = (
    1.85  # day 100 in POWER_Point_Daily_20100101_20191231_045d67N_0123d94W_LST.csv
)
max_temp = (
    10.51  # day 100 in POWER_Point_Daily_20100101_20191231_045d67N_0123d94W_LST.csv
)
hour_temp = 10.06  # temperature at 12

p = Photosynthesis(latitude, 100)
# Solar Variables
p.update_solar_variables(day)
solar_elevation = p.calculate_solar_elevation(increment_hour)
print("Declination=" + str(p._solar_declination))
print("Sunrise=" + str(p._sunrise))
print("Sunset=" + str(p._sunset))
print("Daylength=" + str(p._sunlit_increment))
print("Solar elevation at " + str(increment_hour) + "=" + str(solar_elevation))
print(
    "Incremental PAR="
    + str(p.calculate_incremental_PAR(increment_hour, solar_elevation, PAR, day))
)
print(
    "Hourly direct light extinction="
    + str(p.calculate_hourly_direct_light_extinction(solar_elevation))
)
print(
    "Hourly diffuse light extinction="
    + str(p.calculate_hourly_diffuse_light_extinction(lai))
)
print(
    "Sunlit/Shaded LAI="
    + str(p.calculate_sunlit_shaded_LAI_proportion(solar_elevation, lai))
)
print(
    "Absorbed Incremental PAR="
    + str(
        p.calculate_absorbed_incremental_PAR(
            increment_hour, solar_elevation, PAR, lai, day
        )
    )
)
print("ET Rad=" + str(p.calculate_hourly_ET_rad(solar_elevation, day)))
print("Light limits=" + str(p.calculate_light_limits(PAR, hour_temp)))
print("Assimilation limits=" + str(p.calculate_assimilation_limits()))
print("Rubisco limits=" + str(p.get_rubsico_limits(hour_temp)))
print("Sink limits=" + str(p.get_sink_limits(hour_temp)))
print("CO2 comp=" + str(p.get_CO2_comp(hour_temp)))
print(
    "Leaf assimilation="
    + str(p.calculate_leaf_assimilation(increment_hour, PAR, min_temp, max_temp))
)
print(
    "Photosynthesis="
    + str(p.photosynthesize(PAR, min_temp, max_temp, last_biomass, lai, day))
)

#  for day_increment in p.gauss_integration_params:
#    (abscissa,weight)=day_increment
#    increment_hour = p._sunrise+abscissa*p._sunlit_increment
#    print(increment_hour)
