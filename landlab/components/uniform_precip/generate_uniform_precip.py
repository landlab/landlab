""" generate_uniform_precip.py
 This component generates rainfall
 events based on statistical distributions.

 No particular units must be used, but it was
 written with the storm units in hours (hr)
and depth units in millimeters (mm)


 Written by Jordan Adams, 2013.
"""

import os
import numpy as np
import random
from landlab import Component, ModelParameterDictionary
from landlab.core.model_parameter_dictionary import MissingKeyError


_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'preciptest.in')

class PrecipitationDistribution(Component):
    """Landlab component that generates precipitation events
    using the rectangular Poisson pulse model described in
    Eagleson (1978).

    This component can generate a random storm duration, interstorm
    duration, precipitation intensity or storm depth from a Poisson
    distribution when given a mean value.

    Default input file is named 'preciptest.in' and can be found in
    the landlab.components.uniform_precip folder.

        Inputs
        ------
        input_file : Contains necessary inputs. If not given, default input file is used.
            - MEAN_STORM: (type : float) the mean storm duration if not provided in initialization
            - MEAN_DEPTH: (type : float) the mean storm depth if not provided in initizalization
            - MEAN_INTERSTORM: (type : float) the mean interstorm duration if not provided in initialization
            - RUN_TIME: (type : float) total model run time if not provided in initialization
            - DELTA_T: (type : int) external time step increment if not provided in initialization
                (it is not obligtory to provide a DELTA_T)

    So, without an input file (selecting the default), we can call this component like...

    >>> from landlab.components.uniform_precip import PrecipitationDistribution
    >>> precip = PrecipitationDistribution()

    To use hard-coded values for mean storm, mean interstorm, mean depth, model run time and delta t...
    Say we use 1.5 for mean storm, 15 for mean interstorm, 0.5 for mean depth, 100 for model run time and 1 for delta t...

    >>> precip = PrecipitationDistribution(input_file=None,
    ...     mean_storm=1.5, mean_interstorm=15.0, mean_storm_depth=0.5,
    ...     total_t=100.0, delta_t=1)
    """

    def __init__(self, input_file=None, mean_storm=None, mean_interstorm=None, mean_storm_depth=None, total_t=None, delta_t=None):
        """ This reads in information from the input_file (either default or user
            assigned, and creates an instantaneous storm event drawn from the Poisson distribution
        """

        # First we create an instance of the Model Parameter Dictionary
        MPD = ModelParameterDictionary()

        # If no input_file is given,the default file is used
        if input_file is None:
            input_file = _DEFAULT_INPUT_FILE

        # This reads in the file information
        MPD.read_from_file(input_file)

        # And now we set our different parameters
        # using the model parameter dictionary...
        if mean_storm == None:
            self.mean_storm = MPD.read_float( 'MEAN_STORM')
        else:
            self.mean_storm = mean_storm

        if mean_interstorm == None:
            self.mean_interstorm = MPD.read_float( 'MEAN_INTERSTORM')
        else:
            self.mean_interstorm =mean_interstorm

        if mean_storm_depth== None:
            self.mean_storm_depth = MPD.read_float( 'MEAN_DEPTH')
        else:
            self.mean_storm_depth =mean_storm_depth

        if total_t== None:
            self.run_time = MPD.read_float( 'RUN_TIME')
        else:
            self.run_time =total_t

        if delta_t== None:
            if input_file != _DEFAULT_INPUT_FILE:
                try:
                    self.delta_t = MPD.read_float( 'DELTA_T')
                except MissingKeyError:
                    self.delta_t = None
            else:
                self.delta_t = None
        else:
            self.delta_t =delta_t

        # Mean_intensity is not set by the MPD, but can be drawn from
        # the mean storm depth and mean storm duration.
        self.mean_intensity = self.mean_storm_depth / self.mean_storm

        # If a time series is created later, this blank list will be used.
        self.storm_time_series =[]

        # Given the mean values assigned above using either the model
        # parameter dictionary or the init function, we can call the
        # different methods to assign values from the Poisson distribution.

        self.storm_duration = self.get_precipitation_event_duration()
        self.interstorm_duration = self.get_interstorm_event_duration()
        self.storm_depth = self.get_storm_depth()
        self.intensity = self.get_storm_intensity()
        self._elapsed_time = 0.


    def update(self):
        """If new values for storm duration, interstorm duration, storm depth
        and intensity are needed, this method can be used to update those values
        one time.

        >>> from landlab.components.uniform_precip import (
        ...     PrecipitationDistribution)
        >>> PD = PrecipitationDistribution()
        >>> PD.update()

        Additionally, if we wanted to update several times, a loop could be
        utilized to accomplish this. Say we want 5 storm_durations; this
        pseudo-code represents a way to accomplish this...

        >>> PD = PrecipitationDistribution()
        >>> storm_duration_list=[]
        >>> i = 0
        >>> while i < 4:
        ...     storm_duration_list.append(PD.storm_duration)
        ...     PD.update()
        ...     i+=1
        """

        self.storm_duration = self.get_precipitation_event_duration()
        self.interstorm_duration = self.get_interstorm_event_duration()
        self.storm_depth = self.get_storm_depth()
        self.intensity = self.get_storm_intensity()

    def get_precipitation_event_duration(self):
        """This method is the storm generator.

    This method has one argument: the mean_storm parameter.
    (In Eagleson (1978), this parameter was called Tr.)

    It finds a random storm_duration value
    based on the poisson distribution about the mean.
    This is accomplished using the expovariate function
    from the "random" standard library.
    Additionally, it is rounded to contain 4 significant figures,
    for neatness.

    The if-else statement is very important here. Values of 0
    can exist in the Poission distribution, but it does not make
    sense to have 0 duration storms, so to avoid that,
    we return a storm duration IF it is greater than 0,
    otherwise, recursion is employed to re-call the storm
    generator function and get a new value.

    :returns: storm_duration as a float
    """
        storm = round(random.expovariate(1/self.mean_storm),2)
        while storm == 0:
            storm = round(random.expovariate(1/self.mean_storm),2)
        self.storm_duration = storm
        return self.storm_duration



    def get_interstorm_event_duration(self):
        """ This method is the interstorm duration generator

    This method takes one argument, the mean_interstorm parameter.
    (In Eagleson (1978), this parameter was called Tb.)

    This method is modeled identically to get_precipitation_event_duration()

    This method finds a random value for interstorm_duration
    based on the poisson distribution about the mean.
    This is accomplished using the expovariate function
    from the "random" standard library.
    Additionally, it is rounded to contain 4 significant figures, for neatness.

    The if-else statement is very important here. Values of 0
    can exist in the Poission distribution, but it does not make
    sense to have 0 hour interstorm durations.
    To avoid 0 hour interstorm durations, we return a
    interstorm duration IF it is greater than 0,
    otherwise, recursion is employed to re-call the interstorm
    duration generator function and get a new value.

    :returns: interstorm_duration as a float"""

        interstorm = round(random.expovariate(1/self.mean_interstorm),2)
        while interstorm == 0:
            interstorm = round(random.expovariate(1/self.mean_interstorm),2)
        self.interstorm_duration = interstorm
        return self.interstorm_duration


    def get_storm_depth(self):
        """  This method is the storm depth generator.
    Storm depth is used to generate a realistic
    intensity for different storm events.

    (In Eagleson (1978) this parameter was called "h")

    This method requires storm_duration, mean_storm duration
    and the mean_storm_depth. Storm_duration is generated through
    the initialize() or update() method. mean_storm and mean_storm_depth
    are read in using the ModelParameterDictionary.

    Numpy has a random number generator to get values
    from a given Gamma distribution. It takes two arguments,
    alpha (or the shape parameter), which is the generated over the mean event
    and beta (or the scale parameter), which is the mean value
    These are all arguments in the function, which returns storm depth.

    :returns: storm_depth as a float
    """

        shape_parameter = (self.storm_duration/self.mean_storm)
        scale_parameter = (self.mean_storm_depth)
        self.storm_depth = np.random.gamma(shape_parameter, scale_parameter)
        return self.storm_depth


    def get_storm_intensity(self):
        """   This method draws storm intensity out of the storm depth
    generated by get_storm_depth.

    This method requires the storm_depth and storm_duration
    and is the same as the parameter ("i") in Eagleson (1978), but instead of
    being drawn from Poission, this is drawn from the Gamma distribution
    of ("h"), as h = i*Tr.

    :returns: storm_intensity as a float

                            """
        self.intensity = self.storm_depth / self.storm_duration
        return self.intensity


    def get_storm_time_series(self):
        """
        This method creates a time series of storms based on storm_duration, and
        interstorm_duration. From these values it will calculate a complete
        time series.

        The storm_time_series returned by this method is made up of sublists, each comprising of three
        sub-parts (e.g. [[x,y,z], [a,b,c]]) where x and a are the beginning times of a precipitation
        event, y and b are the ending times of the precipitation event and z and c represent the
        average intensity (mm/hr) of the storm lasting from x to y and a to be, respectively.
        :returns: array containing several sub-arrays of events [start, finish, intensity]


        There is a helper and an iterator. Open to suggestions on how to make this
        sleeker. helper keeps track of storm events so that we can properly adjust the time for the
        storm_time_series list.

        Storm iterator makes sure that the run goes through to the specified time, either read in as "run_time" by
        the ModelParameterDictionary or specified the fire time this method is called."""

        storm = self.get_precipitation_event_duration()
        self.get_storm_depth()
        intensity = self.get_storm_intensity()
        self.storm_time_series.append([0, storm, intensity])

        storm_helper = storm
        storm_iterator = storm
        while storm_iterator <= self.run_time:
            next_storm_start  = storm_helper + round(self.get_interstorm_event_duration(),2)
            next_storm_end = next_storm_start + round(self.get_precipitation_event_duration(),2)
            intensity = round(self.get_storm_intensity(),2)
            self.get_storm_depth()
            self.storm_time_series.append([next_storm_start, next_storm_end, intensity])
            storm_iterator = storm_helper
            storm_helper = next_storm_end
            storm_iterator = storm_helper
        return self.storm_time_series

    def yield_storm_interstorm_duration_intensity(self, subdivide_interstorms=False):
        """
        This method is intended to be equivalent to get_storm_time_series,
        but instead offers a generator functionality. This will be useful in
        cases where the whole sequence of storms and interstorms doesn't need
        to be stored, where we can save memory this way.

        The method keeps track of the DELTA_T such that if a storm needs to be
        generated longer than this supplied model timestep, the generator will
        return the storm in "chunks", until there is no more storm duration.
        e.g.,
        storm of intensity 1. is 4.5 long, the DELTA_T is 2., the generator
        yields (2.,1.) -> (2.,1.) -> (0.5,1.) -> ...

        If DELTA_T is None or not supplied, no subdivision occurs.

        Once a storm has been generated, this method will follow it with the
        next interstorm, yielded as (interstorm_duration, 0.). Note that the
        interstorm will NOT be subdivided according to DELTA_T unless you set
        the flag *subdivide_interstorms* to True.

        The method will keep yielding until it reaches the RUN_TIME, where it
        will terminate.

        YIELDS:
            - a tuple, (interval_duration, rainfall_rate_in_interval)

        One recommended procedure is to instantiate the generator, then call
        instance.next() repeatedly to get the sequence.

        Added DEJH, Dec 2014
        """
        delta_t = self.delta_t
        if delta_t == None:
            assert subdivide_interstorms == False, 'You specified you wanted storm subdivision, but did not provide a DELTA_T to allow this!'
        self._elapsed_time = 0.
        while self._elapsed_time<self.run_time:
            storm_duration = self.get_precipitation_event_duration()
            step_time=0.
            self.get_storm_depth()
            intensity = self.get_storm_intensity() #this is a VELOCITY, i.e., a rainfall rate
            if self._elapsed_time+storm_duration>self.run_time:
                storm_duration = self.run_time-self._elapsed_time
            while delta_t!=None and storm_duration-step_time>delta_t:
                yield (delta_t, intensity)
                step_time+=delta_t
            yield (storm_duration-step_time, intensity)
            self._elapsed_time += storm_duration

            interstorm_duration = self.get_interstorm_event_duration()
            if self._elapsed_time+interstorm_duration>self.run_time:
                interstorm_duration = self.run_time-self._elapsed_time
            if subdivide_interstorms:
                step_time=0.
                while interstorm_duration-step_time>delta_t:
                    yield (delta_t, 0.)
                    step_time += delta_t
                yield (interstorm_duration-step_time, 0.)
            else:
                yield (interstorm_duration, 0.)
            self._elapsed_time += interstorm_duration

    @property
    def elapsed_time(self):
        """
        Return the elapsed time recorded by the module.
        This will be particularly useful in the midst of a yield loop.
        """
        return self._elapsed_time
