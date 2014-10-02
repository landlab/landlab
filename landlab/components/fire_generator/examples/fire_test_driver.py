""" fire_test_driver.py

This is a sample driver simply showing the histogram of fire recurrence
values drawn from the Weibull distribution.

This uses the default file in the generate_fire.py component. Because we are
drawing from a stochastic, statistical distribution, the output may change
after the component is reloaded. 

The default file is currently set using parameters from central Colorado,
given by Cannon et al., (2008) and Moody and Martin (2001). This sample driver
should not be applied at other sites without careful consideration.

"""

from landlab.components.fire_generator.generate_fire import FireGenerator
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


# Initializing the FireGenerator() class from fire_generator.py
FG = FireGenerator()
FG.initialize()

# Finding unknown scale parameter given the mean fire recurrence value
FG.get_scale_parameter()

# Finding a time series based on the Weibull distribution
FG.generate_fire_time_series()

# Sorting the fire series to test the distribution...
FG.firelist.sort()

# Plotting the histogram to show the distribution.
plt.hist(FG.firelist)
plt.title('Weibull Distribution of Fire Recurrence Values')
plt.xlabel('Histogram of Fire Recurrence (years)', fontsize=14)
plt.ylabel('Frequency of Fire Recurrece Values', fontsize=14)
plt.show()
