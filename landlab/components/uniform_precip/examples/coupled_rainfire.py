""" coupled_rainfire.py

This driver evaluates rainfall and fire time series looking for an overlap
which could signify the potential for a large erosion event.

This code runs the fire and rainfall generators for 100 years baed on
input values in the default input file 'fireraininput.txt'.

This simulation parameters are based on central Colorado mean values, based on
values from Cannon et al., 2008 and Moody and Martin, 2001. It is not advised to
use this sample driver and parameters to model fire and precipitation-driven
erosion events in other study areas.

Written by Jordan Marie Adams, 2014.

"""

import os
from matplotlib import pyplot as plt
from landlab.components.uniform_precip import PrecipitationDistribution
from landlab.components.fire_generator import FireGenerator
import numpy as np
from math import ceil

# Input text file name and location
filename = os.path.join(os.path.dirname(__file__), 'fireraininput.txt')

# Initializing the PrecipitationDistribution class using the default file
# and getting the time series needed for comparison against the fire time series.

Rain = PrecipitationDistribution(filename)
Rain.get_storm_time_series()
storm= Rain.storm_time_series

# Initializing the FireGenerator class using the default file and getting the
# time series needed for comparison against the precipitation time series.

# As an additional step, we should find the scale parameter and set it.
# The default value is set to 0.

Fire = FireGenerator(filename)
Fire.get_scale_parameter()
Fire.generate_fire_time_series()
fires = Fire.fire_events

## Methods used to find these potentially erosion-inducing events.

def set_threshold(arr):

    # Gets threshold using the intensity-duration threshold
    # relationship described in Cannon et al., 2008, for the
    # Coal Seam Fire, East of Denver, Colorado.

    # Threshold value changes for each storm event.

    for event in arr:
        indx = arr.index(event)
        D = arr[indx][1] - arr[indx][0]
        I = 6.5*D**(-0.7)
        arr[indx].append(I)
    return arr

def get_max_intensity(storm):

    # This finds the maximum storm intensity
    # and is used simply for plotting purposes.

    intensity_list = []
    for intensity in storm:
        indx = storm.index(intensity)
        intensity_list.append(storm[indx][2])
    max_int = max(intensity_list)
    max_intensity = int(ceil(max_int/100.0))*100
    return max_intensity

def threshold_exceed(storm_arr):

    # This finds storm events where the
    # storm duration and intensity exceed
    # the threshold described by Cannon et al., 2008.

    threshold_exceeded = []
    for each in storm_arr:
        index_each = storm_arr.index(each)
        if storm_arr[index_each][2] >= storm_arr[index_each][3]:
            threshold_exceeded.append(storm_arr[index_each])
    return threshold_exceeded

def find_overlap(fire_arr, rain_arr):

    # This function takes the fire array and threshold-exceeding
    # events to identify where in the timeline there is an overlap
    # between fire and potentially erosion-inducing storm.

    overlap = []
    for item in fire_arr:
        i=0
        indx = fire_arr.index(item)
        start_fire = fire_arr[indx][0]
        end_fire = fire_arr[indx][1]
        while i < len(rain_arr):
            if rain_arr[i][0] < start_fire and rain_arr[i][1] > start_fire and rain_arr[i][1] < end_fire:
                overlap.append(rain_arr[i])
            if rain_arr[i][0] > start_fire and rain_arr[i][1] < end_fire:
                overlap.append(rain_arr[i])
            if rain_arr[i][0] < end_fire and rain_arr[i][1] > end_fire:
                overlap.append(rain_arr[i])
            i+=1
    return overlap

def create_thresh_plot(fire_arr,storm_arr, exceed_arr, max_intense=1000):

    # This function takes the different arrays of fire, storm and threshold-
    # exceeding events and plots them to show potentially erosion-inducing storms.

    plt.figure(1)
    plt.xlabel('Time (years)', fontsize=16) ## x axis label
    plt.ylabel('Rainfall Intensity (mm/day)', fontsize=16) ## y axis label
    plt.title('Randomly Generated Rainfall Time Series', fontsize=18) ## chart title
    ax = plt.gca()
    tick_locations=[0, 3652.42, 7304.84, 10957.26, 14609.68, 18262.1, 21914.52, 25566.96, 29219.36, 32871.78, 36524.2]#labels = range(ticks.size)
    tick_labels=[0,10,20,30,40,50,60,70,80,90,100]
    plt.xticks(tick_locations, tick_labels)
    ax.tick_params(labelsize=16)
    plt.ylim(ymin=0, ymax=max_intense)
    plt.xlim(0, 36524.2)
    for f in fire_arr:
        y = fire_arr.index(f)
        start = fire_arr[y][0]
        end = fire_arr[y][0] + 1.0
        plt.broken_barh([(start, 1)], ((max_intense-200),max_intense), label='Fire', color='orange')

    for s in storm_arr:
        x = storm_arr.index(s) ##creates rainfall graph. x=length of storm. y=intensity
        start = storm_arr[x][0]
        end = storm_arr[x][1] - storm_arr[x][0]
        plt.broken_barh([(start, end)], (0,storm_arr[x][2]), label='Rain', color = 'blue') ## for each Tr period

    for t in exceed_arr:
        z = exceed_arr.index(t)
        start = exceed_arr[z][0]
        end = exceed_arr[z][1] - exceed_arr[z][0]
        plt.broken_barh([(start, end)], (0,exceed_arr[z][2]), label='Threshold Exceeded', color='red')
    plt.show()

# Now to call the functions which do all of the work.

set_threshold(storm)
max_intense = get_max_intensity(storm)
threshold_exceeded = threshold_exceed(storm)
thresh = find_overlap(fires, threshold_exceeded)
create_thresh_plot(fires, storm, thresh, max_intense)
