""" fire_test_series_example.py

This is a sample driver simply showing the histogram of fire recurrence
values drawn from the Weibull distribution.

This uses the default file in the generate_fire.py component. Because we are
drawing from a stochastic, statistical distribution, the output may change
after the component is reloaded.

The default file is currently set using parameters from central Colorado,
given by Cannon et al., (2008) and Moody and Martin (2001). This sample driver
should not be applied at other sites without careful consideration.

"""

from landlab.components.fire_generator import FireGenerator
from matplotlib import pyplot as plt
import numpy as np

mean_fire = 20.0 # Twenty years between fires
shape_param = 5.5 # Right skewed curve


# Initializing the FireGenerator() class from fire_generator.py
FG = FireGenerator(mean_fire_recurrence = mean_fire,
                   shape_parameter = shape_param)

# Finding a time series based on the Weibull distribution

# First, a list of each fire event as it occurs in time. Time = 0 is the start
# of the model run.
fire_events = [0]
time  = 0

# This will store the randomly generated data points, to be plotted as a
# histogram.
sampled_time_to_fire = []

# How long the time series will be.
time_to_run = 50000.0

# Helper iterator.
i = 0

# Loop to generate time series
while time <= time_to_run:
    # Generate event
    new_event = FG.generate_fire_recurrence()

    # Adjust the time so it fits in with our time series.
    time_elapsed = new_event + fire_events[i]

    # Append our new events to the fire time series.
    fire_events.append(time_elapsed)

    # Append our sampled time to fire to create a historgram.
    sampled_time_to_fire.append(new_event)
    i+= 1
    time = time_elapsed



# Plotting the histogram to show the distribution.
plt.figure(1)
plt.hist(sampled_time_to_fire)
plt.title('Weibull Distribution of Fire Recurrence Values')
plt.xlabel('Histogram of Fire Recurrence (years)', fontsize=14)
plt.ylabel('Frequency of Fire Recurrece Values', fontsize=14)

# Plot the first 100 years, dots represent fire events in time.
plt.figure(2)
ones = np.ones(len(fire_events))
plt.plot(fire_events, ones, 'yo')
plt.xlim(0, 100)
plt.title('First 100 years, dots are fire events')
plt.xlabel('Time (years)', fontsize=14)
plt.tick_params(axis='y', which='both', left='off', right='off',
                labelleft='off')
