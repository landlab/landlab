## Evaluating rainfall and fire time series for overlap 
##
## Runs the fire and rainfall generators for 100 years based
## on input values selected by the user. Right now, this model
## is designed to run a simulation for central Colorado, based on
## threshold information calculated by Cannon et al., 2008 for the
## Coal Seam Fire and fire recurrence from Moody and Martin, 2001 
## for the Buffalo Creek Fire.
##
## Current design not necessarily applicable to other study areas.
##
## Written by Jordan Marie Adams
## Last updated: October 29, 2013

import os
from matplotlib import pyplot as plt
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from landlab.components.fire_generator.generate_fire import FireGenerator
import numpy as np
from math import ceil

## INPUT TXT FILE WITH NECESSARY PARAMETERS ##
filename = os.path.join(os.path.dirname(__file__), 'fireraininput.txt')

## INITIALIZING THE CLASSES IN LANDLAB ##

Rain = PrecipitationDistribution()
Rain.initialize(filename)
Rain.get_storm_time_series() ## UNITS IN DAYS
storm= Rain.storm_time_series

Fire = FireGenerator()
Fire.initialize(filename)
Fire.get_scale_parameter() 
Fire.generate_fire_time_series()
fires = Fire.fire_events

## FUNCTIONS TO GET POTENTIAL EROSION EVENTS

## set_threshold() ##
## 
## GETS THRESHOLD BASED ON
## CANNON ET AL., 2008 RELATIONSHIP
## FOR THE COAL SEAM FIRE, EAST OF
## DENVER, COLORADO

def set_threshold(arr): 
    for event in arr:
        indx = arr.index(event)
        D = arr[indx][1] - arr[indx][0]
        I = 6.5*D**(-0.7)
        arr[indx].append(I)
    return arr
    
## get_max_intensity() ##
## 
## GETS THE MAXIMUM INTENSITY
## NEEDED FOR PLOTTING. 

def get_max_intensity(storm):
    intensity_list = []
    for intensity in storm:
        indx = storm.index(intensity)
        intensity_list.append(storm[indx][2])
    max_int = max(intensity_list)
    max_intensity = int(ceil(max_int/100.0))*100
    return max_intensity

## threshold_exceed()
##
## FINDS STORM EVENTS WHERE THE
## ACTUAL STORM DURATION AND INTENSITY
## EXCEEDS THE THRESHOLD 

def threshold_exceed(storm_arr):
    threshold_exceeded = []
    for each in storm_arr:
        index_each = storm_arr.index(each)
        if storm_arr[index_each][2] >= storm_arr[index_each][3]:
            threshold_exceeded.append(storm_arr[index_each])
    return threshold_exceeded

##  find_overlap()
## 
## THIS FUNCTION TAKES THE FIRE_ARRAY
## AND THRESHOLD ARRAY AND COMPARES THE TIMELINES
## TO FIND WHERE THERE IS A OVERLAP AND A 
## POTENTIALLY EROSION-INDUCING STORM
    
def find_overlap(fire_arr, rain_arr):
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
    
## write_to_file() 
##
## Takes any arrays generated prior and 
## writes values to .txt files so that
## additional plots can be made elsewhere if needed.
    
def write_to_file(fire_arr, storm_arr, exceed_arr, firetxt, stormtxt, exceedtxt):
    np.savetxt(firetxt, fire_arr)
    np.savetxt(stormtxt, storm_arr)
    np.savetxt(exceedtxt, exceed_arr)
    
## CALLING THE FUNCTIONS TO DO ALL THE WORK....
    
set_threshold(storm)
max_intense = get_max_intensity(storm)         
threshold_exceeded = threshold_exceed(storm)
thresh = find_overlap(fires, threshold_exceeded)

## create_thresh_plot()
## 
## TAKES THE DIFFERENT ARRAYS (FIRES, STORMS
## AND THRESHOLD EXCEEDENCE) AND PLOTS THEM 
## TO SHOW POTENTIALLY EROSION-INDUCING STORMS

def create_thresh_plot(fire_arr,storm_arr, exceed_arr, max_intense=1000):
    plt.figure(1)
    plt.xlabel('Time (years)', fontsize=20) ## x axis label
    plt.ylabel('Rainfall Intensity (mm/day)', fontsize=20) ## y axis label
    plt.title('Randomly Generated Rainfall Time Series', fontsize=24) ## chart title
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
        
create_thresh_plot(fires, storm, thresh, max_intense)
##write_to_file(fires, storm, thresh, fire_output.txt, storm_output.txt, thresh_output.txt)
