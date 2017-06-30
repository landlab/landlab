"""Generate precipitation using the Poisson pulse model.

Landlab component that generates precipitation events using the rectangular
Poisson pulse model described in Eagleson (1978, Water Resources Research).


No particular units must be used, but it was written with the storm units in
hours (hr) and depth units in millimeters (mm)

Written by Jordan Adams, 2013, updated May 2016
"""


import random
import numpy as np
from landlab import Component


class PrecipitationDistribution(Component):

    """Generate precipitation events.

    This component can generate a random storm duration, interstorm
    duration, precipitation intensity or storm depth from a Poisson
    distribution when given a mean value.

    Construction::

        PrecipitationDistribution(mean_storm_duration=0.0,
                                  mean_interstorm_duration=0.0,
                                  mean_storm_depth=0.0, total_t=0.0,
                                  delta_t=0.0, random_seed=0)

    Parameters
    ----------
    mean_storm_duration : float
        Average duration of a precipitation event.
    mean_interstorm_duration : float
        Average duration between precipitation events.
    mean_storm_depth : float
        Average depth of precipitation events.
    total_t : float, optional
        If generating a time series, the total amount of time .
    delta_t : float, optional
        If you want to break up storms into determined subsections using
        yield_storm_interstorm_duration_intensity, a delta_t is needed.
    random_seed : int or float, optional
        Seed value for random-number generator.

    Examples
    --------
    >>> from landlab.components.uniform_precip import PrecipitationDistribution
    >>> import numpy as np
    >>> np.random.seed(np.arange(10))

    To use hard-coded values for mean storm, mean interstorm, mean depth,
    model run time and delta t...  Say we use 1.5 for mean storm, 15 for mean
    interstorm, 0.5 for mean depth, 100 for model run time and 1 for delta t...

    >>> precip = PrecipitationDistribution(mean_storm_duration = 1.5,
    ...     mean_interstorm_duration = 15.0, mean_storm_depth = 0.5,
    ...     total_t = 100.0, delta_t = 1.)
    >>> for (dt, rate) in precip.yield_storm_interstorm_duration_intensity():
    ...     pass  # and so on
    """

    _name = 'PrecipitationDistribution'

    _input_var_names = tuple()

    _output_var_names = tuple()

    _var_units = dict()

    _var_mapping = dict()

    _var_doc = dict()

    def __init__(self, mean_storm_duration=0.0, mean_interstorm_duration=0.0,
                 mean_storm_depth=0.0, total_t=0.0, delta_t=None,
                 random_seed=0, **kwds):
        """Create the storm generator.

        Parameters
        ----------
        mean_storm_duration : float
            Average duration of a precipitation event.
        mean_interstorm_duration : float
            Average duration between precipitation events.
        mean_storm_depth : float
            Average depth of precipitation events.
        total_t : float, optional
            If generating a time series, the total amount of time.
        delta_t : float or None, optional
            If you want to break up storms into determined subsections using
            yield_storm_interstorm_duration_intensity, a delta_t is needed.
        random_seed : int or float, optional
            Seed value for random-number generator.
        """

        self.mean_storm_duration = mean_storm_duration

        self.mean_interstorm_duration = mean_interstorm_duration

        self.mean_storm_depth = mean_storm_depth

        self.run_time = total_t

        self.delta_t = delta_t

        if self.delta_t == 0.:
            self.delta_t = None

        # Mean_intensity is not set by the MPD, but can be drawn from
        # the mean storm depth and mean storm duration.
        self.mean_intensity = self.mean_storm_depth / self.mean_storm_duration

        # If a time series is created later, this blank list will be used.
        self.storm_time_series = []

        # Seed the random-number generator
        self.seed_generator(random_seed)

        # Given the mean values assigned above using either the model
        # parameter dictionary or the init function, we can call the
        # different methods to assign values from the Poisson distribution.

        self.storm_duration = self.get_precipitation_event_duration()
        self.interstorm_duration = self.get_interstorm_event_duration()
        self.storm_depth = self.get_storm_depth()
        self.intensity = self.get_storm_intensity()
        self._elapsed_time = 0.

    def update(self):
        """Update the storm values.

        If new values for storm duration, interstorm duration, storm depth
        and intensity are needed, this method can be used to update those
        values one time.

        Examples
        --------
        >>> from landlab.components import PrecipitationDistribution
        >>> precip = PrecipitationDistribution(mean_storm_duration=1.5,
        ...     mean_interstorm_duration=15.0, mean_storm_depth=0.5,
        ...     total_t=100.0, delta_t=1)

        Additionally, if we wanted to update several times, a loop could be
        utilized to accomplish this. Say we want 5 storm_durations; this
        pseudo-code represents a way to accomplish this...

        >>> storm_duration_list = []
        >>> i = 0
        >>> while i < 4:
        ...     storm_duration_list.append(precip.storm_duration)
        ...     precip.update()
        ...     i += 1
        """
        self.storm_duration = self.get_precipitation_event_duration()
        self.interstorm_duration = self.get_interstorm_event_duration()
        self.storm_depth = self.get_storm_depth()
        self.intensity = self.get_storm_intensity()

    def get_precipitation_event_duration(self):
        """This method is the storm generator.

        This method has one argument: the mean_storm_duration parameter.
        (In Eagleson (1978), this parameter was called Tr.)

        It finds a random storm_duration value
        based on the poisson distribution about the mean.
        This is accomplished using the expovariate function
        from the "random" standard library.

        Returns
        -------
        float
            The storm duration.
        """
        return random.expovariate(1.0 / self.mean_storm_duration)

    def get_interstorm_event_duration(self):
        """Generate interstorm events.

        This method takes one argument, the mean_interstorm_duration parameter.
        (In Eagleson (1978), this parameter was called Tb.)

        This method is modeled identically to
        get_precipitation_event_duration()

        This method finds a random value for interstorm_duration
        based on the poisson distribution about the mean.
        This is accomplished using the expovariate function
        from the "random" standard library.

        Returns
        -------
        float
            The interstorm duration.
        """
        return random.expovariate(1.0 / self.mean_interstorm_duration)

    def get_storm_depth(self):
        """Generate storm depth.

        Storm depth is used to generate a realistic
        intensity for different storm events.

        (In Eagleson (1978) this parameter was called "h")

        This method requires storm_duration, mean_storm_duration
        and the mean_storm_depth. Storm_duration is generated through
        the initialize() or update() method.

        Numpy has a random number generator to get values
        from a given Gamma distribution. It takes two arguments,
        alpha (or the shape parameter), which is the generated over the mean
        event and beta (or the scale parameter), which is the mean value
        These are all arguments in the function, which returns storm depth.

        Returns
        -------
        float
            The storm depth.
        """

        shape_parameter = (self.storm_duration / self.mean_storm_duration)
        scale_parameter = (self.mean_storm_depth)
        self.storm_depth = np.random.gamma(shape_parameter, scale_parameter)
        return self.storm_depth

    def get_storm_intensity(self):
        """Get the storm intensity.

        This method draws storm intensity out of the storm depth generated by
        get_storm_depth.

        This method requires the storm_depth and storm_duration
        and is the same as the parameter ("i") in Eagleson (1978), but instead
        of being drawn from Poission, this is drawn from the Gamma distribution
        of (*h*), as :math:`h = i * Tr`.

        Returns
        -------
        float
            The storm intensity.
        """
        self.intensity = self.storm_depth / self.storm_duration
        return self.intensity

    def get_storm_time_series(self):
        """Get a time series of storms.

        This method creates a time series of storms based on storm_duration,
        and interstorm_duration. From these values it will calculate a complete
        time series.

        The storm_time_series returned by this method is made up of sublists,
        each comprising of three sub-parts (e.g. [[x,y,z], [a,b,c]]) where x
        and a are the beginning times of a precipitation event, y and b are the
        ending times of the precipitation event and z and c represent the
        average intensity (mm/hr) of the storm lasting from x to y and a to be,
        respectively.

        Returns
        -------
        array
            containing several sub-arrays of events [start, finish, intensity]
        """

        storm = self.get_precipitation_event_duration()
        self.get_storm_depth()
        intensity = self.get_storm_intensity()
        self.storm_time_series.append([0, storm, intensity])

        storm_helper = storm
        storm_iterator = storm
        while storm_iterator <= self.run_time:
            next_storm_start = storm_helper + (round(
                                    self.get_interstorm_event_duration(), 2))
            next_storm_end = next_storm_start + (round(
                                self.get_precipitation_event_duration(), 2))
            intensity = round(self.get_storm_intensity(), 2)
            self.get_storm_depth()
            self.storm_time_series.append(
                            [next_storm_start, next_storm_end, intensity])
            storm_iterator = storm_helper
            storm_helper = next_storm_end
            storm_iterator = storm_helper
        return self.storm_time_series

    def yield_storm_interstorm_duration_intensity(self,
                                                  subdivide_interstorms=False):
        """Iterator for a time series of storms.

        This method is intended to be equivalent to get_storm_time_series,
        but instead offers a generator functionality. This will be useful in
        cases where the whole sequence of storms and interstorms doesn't need
        to be stored, where we can save memory this way.

        The method keeps track of the delta_t such that if a storm needs to be
        generated longer than this supplied model timestep, the generator will
        return the storm in "chunks", until there is no more storm duration.
        e.g.,
        storm of intensity 1. is 4.5 long, the delta_t is 2., the generator
        yields (2.,1.) -> (2.,1.) -> (0.5,1.) -> ...

        If delta_t is None or not supplied, no subdivision occurs.

        Once a storm has been generated, this method will follow it with the
        next interstorm, yielded as (interstorm_duration, 0.). Note that the
        interstorm will NOT be subdivided according to delta_t unless you set
        the flag *subdivide_interstorms* to True.

        The method will keep yielding until it reaches the RUN_TIME, where it
        will terminate.

        Yields
        ------
        tuple of float
            (interval_duration, rainfall_rate_in_interval)

        Notes
        -----
        One recommended procedure is to instantiate the generator, then call
        instance.next() repeatedly to get the sequence.
        """
        # Added DEJH, Dec 2014
        delta_t = self.delta_t
        if delta_t is None:
            assert subdivide_interstorms is False, (
                'You specified you wanted storm subdivision, but did not ' +
                'provide a delta_t to allow this!')
        self._elapsed_time = 0.
        while self._elapsed_time < self.run_time:
            storm_duration = self.get_precipitation_event_duration()
            step_time = 0.
            self.get_storm_depth()
            intensity = self.get_storm_intensity()  # this is a rainfall rate
            if self._elapsed_time + storm_duration > self.run_time:
                storm_duration = self.run_time - self._elapsed_time
            while delta_t is not None and storm_duration - step_time > delta_t:
                yield (delta_t, intensity)
                step_time += delta_t
            yield (storm_duration - step_time, intensity)
            self._elapsed_time += storm_duration

            # If the last storm did not use up all our elapsed time, generate
            # an inter-storm period.
            if self._elapsed_time < self.run_time:
                interstorm_duration = self.get_interstorm_event_duration()
                if self._elapsed_time + interstorm_duration > self.run_time:
                    interstorm_duration = self.run_time - self._elapsed_time
                if subdivide_interstorms:
                    step_time = 0.
                    while interstorm_duration-step_time > delta_t:
                        yield (delta_t, 0.)
                        step_time += delta_t
                    yield (interstorm_duration - step_time, 0.)
                else:
                    yield (interstorm_duration, 0.)
                self._elapsed_time += interstorm_duration

    def generate_from_stretched_exponential(self, scale, shape):
        """Generate and return a random variable from a stretched exponential
        distribution with given scale and shape.
        
        Examples
        --------
        >>> np.random.seed(0)
        >>> np.round(np.random.rand(3), 6)  # these are our 3 rand #s to test
        array([ 0.548814,  0.715189,  0.602763])
        >>> from landlab.components import PrecipitationDistribution
        >>> pd = PrecipitationDistribution(mean_storm_duration=1.0,
        ...                                mean_interstorm_duration=1.0,
        ...                                mean_storm_depth=1.0)
        >>> np.random.seed(0)  # re-set seed so we get the same 3 #s
        >>> np.round(1000 * pd.generate_from_stretched_exponential(2.0, 0.5))
        720.0
        >>> np.round(1000 * pd.generate_from_stretched_exponential(2.0, 0.5))
        225.0
        >>> np.round(1000 * pd.generate_from_stretched_exponential(2.0, 0.5))
        513.0
        """
        return scale * ((-np.log(np.random.rand())) ** (1.0 / shape))

    def seed_generator(self, seedval=0):
        """Seed the random-number generator.

        The examples illustrate:
        (1) that we can get the same sequence again by re-seeding with the
            same value (the default is zero)
        (2) when we use a value other than the default, we get a different
            sequence

        Examples
        --------
        >>> precip = PrecipitationDistribution(mean_storm_duration = 1.5,
        ...     mean_interstorm_duration = 15.0, mean_storm_depth = 0.5,
        ...     total_t = 100.0, delta_t = 1.)
        >>> round(precip.storm_duration, 2)
        2.79
        >>> round(precip.interstorm_duration, 2)
        21.28
        >>> round(precip.storm_depth, 2)
        2.45
        >>> round(precip.intensity, 2)
        0.88
        >>> precip.seed_generator() # re-seed and get same sequence again
        >>> round(precip.get_precipitation_event_duration(), 2)
        2.79
        >>> round(precip.get_interstorm_event_duration(), 2)
        21.28
        >>> round(precip.get_storm_depth(), 2)
        2.45
        >>> round(precip.get_storm_intensity(), 2)
        0.88
        >>> precip = PrecipitationDistribution(mean_storm_duration = 1.5,
        ...     mean_interstorm_duration = 15.0, mean_storm_depth = 0.5,
        ...     total_t = 100.0, delta_t = 1., random_seed=1)
        >>> round(precip.storm_duration, 2) # diff't vals with diff't seed
        0.22
        >>> round(precip.interstorm_duration, 2)
        28.2
        >>> round(precip.storm_depth, 4)
        0.0012
        >>> round(precip.intensity, 4)
        0.0054
        """
        random.seed(seedval)
        np.random.seed(seedval)

    @property
    def elapsed_time(self):
        """Get the elapsed time recorded by the module.

        This will be particularly useful in the midst of a yield loop.
        """
        return self._elapsed_time
