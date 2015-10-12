""" generate_fire.py

This component generates random numbers using the Weibull distribution
(Weibull, 1951). No particular units must be used, but it was written with
the fire recurrence units in time (yrs).

Using the Weibull Distribution assumes two things: All elements within the study
area have the same fire regime. Each element must have (on average) a constant
fire regime during the time span of the study.

As of Sept. 2013, fires are considered instantaneous events independent of
other fire events in the time series.

Written by Jordan Marie Adams, 2013.

"""

import os
from random import weibullvariate
from scipy import special
from landlab import ModelParameterDictionary

_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__), 'fire.txt')

class FireGenerator:

    def __init__(self, input_file=None):
        """All initial values are set to zero until initalized
        using the initialize() method which reads in data using
        the ModelParameterDictionary and sets the random variables
        according to their respective distribution.

        REQUIRED PARAMETERS
        ------------------

        Shape Parameter: Describes the skew of the Weibull distribution.
        If shape < 3.5, data skews left.
        If shape == 3.5, data is normal.
        If shape > 3.5, data skews right.

        Scale Parameter: Describes the peak of the Weibull distribution,
        located at 63.5% value of the cumulative distribution function. If unknown,
        it can be found using mean fire recurrence value and the get_scale_parameter()
        method described later.

        Mean Fire Recurrence : Average recurrence for a given area, elevation, veg type, etc.

        Total Run Time : Total model run time

        Delta T : Model time step.


        -------
        Time to Next Fire: Value generated from the random.weibullvariate() function based
        on the scale and shape parameters
        """
        # We import methods from ModelParameterDictionary
        # to read the parameters from the input file.

        #If the scale parameter is unknown, it can be found using the mean
        # fire recurrence value, which MUST be known or estimated to run the
        # get_scale_parameter() method."""

        MPD = ModelParameterDictionary()

        if input_file is None:
            input_file = _DEFAULT_INPUT_FILE
        MPD.read_from_file(input_file)

        self.shape_parameter = MPD.read_float("SHAPE_PARAMETER")
        self.scale_parameter = MPD.read_float("SCALE_PARAMETER")
        self.mean_fire_recurrence = MPD.read_float("MEAN_FIRE_RECURRENCE")
        self.total_run_time = MPD.read_float("RUN_TIME")
        self.delta_t = MPD.read_int("DELTA_T")
        self.time_to_next_fire = 0.0

    def get_scale_parameter(self):

        """ If the scale factor is unknown, we can draw it from the mean
        fire recurrence interval, given the following equation:

        This ONLY works if we have the shape parameter. Generally, the shape
        parameter will be greater than 3.5, creating a distribution that is skewed
        to the higher values (skewed to the right).

        Mean_fire_recurrence = scale_parameter*(gamma_function(1+(1/shape)))

        :returns: scale_parameter as a float"""

        if self.scale_parameter == 0.0:
            shape_in_gamma_func = float(1+(1/self.shape_parameter))
            gamma_func = special.gamma(shape_in_gamma_func)
            self.scale_parameter = (self.mean_fire_recurrence/gamma_func)
            return self.scale_parameter
        else:
            return self.scale_parameter

    def generate_fire_recurrence(self):

        """ Finds the time to next fire (fire recurrence) based on the scale parameter (63.5% of
        fire Weibull distribution) and the shape parameter (describes the skew of the histogram, shape = 3.5
        represents a normal distribution).

        Rounds the time to next fire to 4 significant figures, for neatness.

        :returns: time_to_next_fire as a float"""

        self.time_to_next_fire = round(weibullvariate(self.scale_parameter, self.shape_parameter),2)
        return self.time_to_next_fire

    def generate_fire_time_series(self):

        """Allows for a series of fire events to be generated given
        a total time.

        Created for situations where total run time is definite, and number
        of fires can change across different runs.

        :creates: array with several fire events, all values are float"""

        self.fire_events =[]
        event = self.generate_fire_recurrence()
        end_event = event + 365.0
        self.fire_events.append([event, end_event])
        t = 0
        i = 0
        while t <= self.total_run_time:
            fire = self.generate_fire_recurrence()
            start_fire = self.fire_events[i][0] + (fire)
            end_fire = start_fire + (365.0)
            self.fire_events.append([start_fire, end_fire])
            t += end_fire
            i+=1

    def update(self):
        """Update function allows us to update the value of "time_to_next_fire" by
        re-calling the generate_fire_reccurence() function.

        Created for instances when a definite number of fires need to
        be generated.

        :returns: updated value for time_to_next_fire as a float"""

        self.time_to_next_fire = self.generate_fire_recurrence()
        return self.time_to_next_fire
