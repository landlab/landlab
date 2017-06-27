""" uniform_precip_test.py

Sample driver file for the PrecipitationDistribution component.

This file will generate storm characteristics, update those characteristics,
and finally generate a storm time series across 100 years. After the storm
time series is generated, a plot will be drawn with storm intensity on the y-axis
and the storm duration along the x - axis.

Keep in mind - because these are drawn from a statistical distribution in
a stochastic manner, each run of this driver will differ from the last run.

Written by Jordan Marie Adams, 2013.

"""
from __future__ import print_function

from landlab.components.uniform_precip import PrecipitationDistribution
from matplotlib import pyplot as plt


def create_precip_plot(storm_arr):
    # Plotting precipitation distribution time series

    # Creating a new figure instance
    plt.figure(1)

    # Labeling the x and y axes
    plt.xlabel('Time (years)', fontsize=14)
    plt.ylabel('Rainfall Intensity (mm/day)’, fontsize=14)

    # Setting the plot title
    plt.title('Randomly Generated Rainfall Time Series', fontsize=16)

    # Now to plot the axes the way we'd like them...
    ax = plt.gca()

    # Manually set the ten times (in hours) to plot across 100 years
    tick_locations = ([0, 3652.42, 7304.84, 10957.26, 14609.68, 18262.1,
                       21914.52, 25566.96, 29219.36, 32871.78, 36524.2])

    # This next list will actually replace the hours with time in years.
    tick_labels=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # Swapping out the two tick labels to plot intensity against years.
    plt.xticks(tick_locations, tick_labels)

    # Setting tick label size.
    ax.tick_params(labelsize=14)

    # Setting the limits for the x and y
    plt.xlim(0, 36524.2)
    plt.ylim(ymin=0, ymax=20)

    # Looping through the storm array to plot the intensity as the height of each bar plot
    # and the width will correspond to the storm duration.

    for s in storm_arr:
        x = storm_arr.index(s)
        start = storm_arr[x][0]
        end = storm_arr[x][1] - storm_arr[x][0]
        plt.broken_barh([(start, end)], (0,storm_arr[x][2]), label='Rain',
                                                                color = 'blue')

    plt.show()


def main():
    # First we will create an instance of PrecipitationDistribution
    PD = PrecipitationDistribution(mean_storm_duration = 2.0,
                                   mean_interstorm_duration = 50.0,
                                   mean_storm_depth = 0.05, total_t = 37000.)

    # Because the values for storm duration, interstorm duration, storm
    # depth and intensity are set stochastically in the initialization
    # phase, we should see that they seem reasonable.

    print("Mean storm duration is: ", PD.mean_storm_duration, " hours, while",
             "the value from the Poisson distribution is: ", PD.storm_duration)
    print("Mean interstorm Duration is: ", PD.mean_interstorm_duration,
          'hours, while the value from the Poisson distribution is: ',
          PD.interstorm_duration)
    print("Mean storm depth is: ", PD.mean_storm_depth, "mm, while the value",
          "from the Poisson distribution is: ", PD.storm_depth)
    print("Mean intensity is: ", PD.mean_intensity, "mm/hr, while the value",
          "from the Poisson distribution is: ", PD.intensity)
    print('\n')


    # If we generate a time series we can plot a precipitation distribution
    PD.get_storm_time_series()

    # And get the storm array from the component..
    storm_arr = PD.storm_time_series

    # And now to call the plotting method.
    create_precip_plot(storm_arr)

if __name__ == '__main__':
    main()

