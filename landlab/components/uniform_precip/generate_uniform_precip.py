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
from landlab import ModelGrid


class PrecipitationDistribution(Component):
    """Generate precipitation events.

    This component can generate a random storm duration, interstorm
    duration, precipitation intensity or storm depth from a Poisson
    distribution when given a mean value.

    Examples
    --------
    >>> from landlab.components.uniform_precip import PrecipitationDistribution
    >>> import numpy as np
    >>> np.random.seed(np.arange(10))

    To use hard-coded values for mean storm, mean interstorm, mean depth,
    model run time and delta t...  Say we use 1.5 for mean storm, 15 for mean
    interstorm, 0.5 for mean depth, 100 for model run time and 1 for delta t...

    >>> precip = PrecipitationDistribution(
    ...     mean_storm_duration=1.5,
    ...     mean_interstorm_duration=15.0,
    ...     mean_storm_depth=0.5,
    ...     total_t=100.0,
    ...     delta_t=1.0,
    ... )
    >>> for dt, rate in precip.yield_storm_interstorm_duration_intensity():
    ...     pass  # and so on
    ...

    Alternatively, we can pass a grid to the component, and call yield_storms()
    to generate storm-interstorm float pairs while the intensity data is stored
    in the grid scalar field 'rainfall__flux':

    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((4, 5))
    >>> precip = PrecipitationDistribution(
    ...     mg,
    ...     mean_storm_duration=1.5,
    ...     mean_interstorm_duration=15.0,
    ...     mean_storm_depth=0.5,
    ...     total_t=46.0,
    ... )
    >>> storm_dts = []
    >>> interstorm_dts = []
    >>> intensities = []
    >>> precip.seed_generator(seedval=1)
    >>> for storm_dt, interstorm_dt in precip.yield_storms():
    ...     storm_dts.append(storm_dt)
    ...     interstorm_dts.append(interstorm_dt)
    ...     intensities.append(mg.at_grid["rainfall__flux"])
    ...
    >>> len(storm_dts) == 4  # 4 storms in the simulation
    True
    >>> len(interstorm_dts) == len(storm_dts)
    True
    >>> rf_intensities_to_test = np.array(
    ...     [
    ...         0.8138257984406472,
    ...         0.15929112025199238,
    ...         0.17254519305000884,
    ...         0.09817611240558813,
    ...     ]
    ... )
    >>> np.allclose(intensities, rf_intensities_to_test)
    True
    >>> np.isclose(sum(storm_dts) + sum(interstorm_dts), 46.0)  # test total_t
    True
    >>> np.isclose(interstorm_dts[-1], 0.0)  # sequence truncated as necessary
    True

    We can also turn the generator straight into a list, like this:

    >>> precip.seed_generator(seedval=1)
    >>> steps = [(dt + istorm_dt) for (dt, istorm_dt) in precip.yield_storms()]
    >>> np.isclose(sum(steps), 46.0)
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Eagleson, P. (1978). Climate, soil, and vegetation: 2. The distribution of
    annual precipitation derived from observed storm sequences. Water Resources
    Research  14(5), 713-721. https://dx.doi.org/10.1029/wr014i005p00713

    """

    _name = "PrecipitationDistribution"

    _unit_agnostic = True

    _info = {
        "rainfall__flux": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "[depth unit]/[time unit]",
            "mapping": "grid",
            "doc": "Depth of water delivered per unit time in each storm",
        }
    }

    def __init__(
        self,
        grid=None,
        mean_storm_duration=1.0,
        mean_interstorm_duration=1.0,
        mean_storm_depth=1.0,
        total_t=0.0,
        delta_t=None,
        random_seed=0,
    ):
        """Create the storm generator.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab grid (optional). If provided, storm intensities will be
            stored as a grid scalar field as the component simulates storms.
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
        super().__init__(grid)

        self._mean_storm_duration = mean_storm_duration

        self._mean_interstorm_duration = mean_interstorm_duration

        self._mean_storm_depth = mean_storm_depth

        self._run_time = total_t

        self._delta_t = delta_t

        if self._delta_t == 0.0:
            self._delta_t = None

        # Mean_intensity is not set by the MPD, but can be drawn from
        # the mean storm depth and mean storm duration.
        self._mean_intensity = self._mean_storm_depth / self._mean_storm_duration

        # If a time series is created later, this blank list will be used.
        self._storm_time_series = []

        # Seed the random-number generator
        self.seed_generator(random_seed)

        # Given the mean values assigned above using either the model
        # parameter dictionary or the init function, we can call the
        # different methods to assign values from the Poisson distribution.

        self._storm_duration = self.get_precipitation_event_duration()
        self._interstorm_duration = self.get_interstorm_event_duration()
        self._storm_depth = self.get_storm_depth()
        self._elapsed_time = 0.0

        # Test if we got a grid. If we did, then assign it to _grid, and we
        # are able to use the at_grid field. If not, that's cool too.
        if grid is not None:
            assert isinstance(grid, ModelGrid)  # must be a grid

        # build LL fields, if a grid is supplied:
        if grid is not None:
            self.grid.add_field("rainfall__flux", 0.0, at="grid")
            self._gridupdate = True
        else:
            self._gridupdate = False

        self._intensity = self.get_storm_intensity()

    @property
    def storm_duration(self):
        """Duration of storm.

        [T]
        """
        return self._storm_duration

    @property
    def interstorm_duration(self):
        """Interstorm duration.

        [T]
        """
        return self._interstorm_duration

    @property
    def storm_depth(self):
        """Depth of water in the storm.

        [L]
        """
        return self._storm_depth

    def update(self):
        """Update the storm values.

        If new values for storm duration, interstorm duration, storm depth
        and intensity are needed, this method can be used to update those
        values one time.

        Examples
        --------
        >>> from landlab.components import PrecipitationDistribution
        >>> precip = PrecipitationDistribution(
        ...     mean_storm_duration=1.5,
        ...     mean_interstorm_duration=15.0,
        ...     mean_storm_depth=0.5,
        ...     total_t=100.0,
        ...     delta_t=1,
        ... )

        Additionally, if we wanted to update several times, a loop could be
        utilized to accomplish this. Say we want 5 storm_durations; this
        pseudo-code represents a way to accomplish this...

        >>> storm_duration_list = []
        >>> i = 0
        >>> while i < 4:
        ...     storm_duration_list.append(precip.storm_duration)
        ...     precip.update()
        ...     i += 1
        ...

        Note though that alternatively we could also do this, avoiding the
        method entirely...


        >>> # ^^this lets you "manually" get the next item from the iterator
        >>> precip = PrecipitationDistribution(
        ...     mean_storm_duration=1.5,
        ...     mean_interstorm_duration=15.0,
        ...     mean_storm_depth=0.5,
        ...     total_t=46.0,
        ... )
        >>> storm_duration_list = []
        >>> for i in range(5):
        ...     storm_duration_list.append(
        ...         next(precip.yield_storm_interstorm_duration_intensity())[0]
        ...     )
        ...

        Notice that doing this will *not* automatically stop the iteration,
        however - it will continue ad infinitum.
        """
        self._storm_duration = self.get_precipitation_event_duration()
        self._interstorm_duration = self.get_interstorm_event_duration()
        self._storm_depth = self.get_storm_depth()
        self._intensity = self.get_storm_intensity()

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
        return random.expovariate(1.0 / self._mean_storm_duration)

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
        return random.expovariate(1.0 / self._mean_interstorm_duration)

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

        shape_parameter = self._storm_duration / self._mean_storm_duration
        scale_parameter = self._mean_storm_depth
        self._storm_depth = np.random.gamma(shape_parameter, scale_parameter)
        return self._storm_depth

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
        self._intensity = self._storm_depth / self._storm_duration
        if self._gridupdate:
            self._grid.at_grid["rainfall__flux"] = self._intensity
        return self._intensity

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

        Even if a grid was passed to the component at instantiation, calling
        this method does not update the grid fields.

        Returns
        -------
        array
            containing several sub-arrays of events [start, finish, intensity]
        """

        storm = self.get_precipitation_event_duration()
        self.get_storm_depth()
        intensity = self.get_storm_intensity()
        self._storm_time_series.append([0, storm, intensity])

        storm_helper = storm
        storm_iterator = storm
        while storm_iterator <= self._run_time:
            next_storm_start = storm_helper + (
                round(self.get_interstorm_event_duration(), 2)
            )
            next_storm_end = next_storm_start + (
                round(self.get_precipitation_event_duration(), 2)
            )
            intensity = round(self.get_storm_intensity(), 2)
            self.get_storm_depth()
            self._storm_time_series.append(
                [next_storm_start, next_storm_end, intensity]
            )
            storm_iterator = storm_helper
            storm_helper = next_storm_end
            storm_iterator = storm_helper
        return self._storm_time_series

    def yield_storm_interstorm_duration_intensity(self, subdivide_interstorms=False):
        """Iterator for a time series of storms interspersed with interstorms.

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
        instance.next() (in Python 2) or next(instance) (in Python 3)
        repeatedly to get the sequence.
        """
        # Added DEJH, Dec 2014
        # Modified to use an optional output field, DEJH 1/8/17

        delta_t = self._delta_t
        if delta_t is None:
            assert subdivide_interstorms is False, (
                "You specified you wanted storm subdivision, but did not "
                + "provide a delta_t to allow this!"
            )
        self._elapsed_time = 0.0
        while self._elapsed_time < self._run_time:
            storm_duration = self.get_precipitation_event_duration()
            step_time = 0.0
            self.get_storm_depth()
            self._intensity = self.get_storm_intensity()  # this is a rate
            # ^ this updates the grid field, if needed
            if self._elapsed_time + storm_duration > self._run_time:
                storm_duration = self._run_time - self._elapsed_time
            while delta_t is not None and storm_duration - step_time > delta_t:
                yield (delta_t, self._intensity)
                step_time += delta_t
            yield (storm_duration - step_time, self._intensity)
            self._elapsed_time += storm_duration

            # If the last storm did not use up all our elapsed time, generate
            # an inter-storm period.
            if self._elapsed_time < self._run_time:
                interstorm_duration = self.get_interstorm_event_duration()
                if self._elapsed_time + interstorm_duration > self._run_time:
                    interstorm_duration = self._run_time - self._elapsed_time
                self._intensity = 0.0
                if self._gridupdate:
                    self._grid.at_grid["rainfall__flux"] = 0.0
                if subdivide_interstorms:
                    step_time = 0.0
                    while interstorm_duration - step_time > delta_t:
                        yield (delta_t, 0.0)
                        step_time += delta_t
                    yield (interstorm_duration - step_time, 0.0)
                else:
                    yield (interstorm_duration, 0.0)
                self._elapsed_time += interstorm_duration

    def yield_storms(self):
        """Iterator for a time series of storm-interstorm pairs.

        This method is very similar to this component's other generator,
        yield_storm_interstorm_duration_intensity(), but the way it yields is
        slightly different. Instead of yielding (interval_duration, rf_rate),
        with interstorms represented as intervals with rf_rate = 0, it yields:

            (storm_duration, interstorm_duration)

        When each tuple pair is yielded, the grid scalar field 'rainfall__flux'
        is updated with the rainfall rate occuring during storm_duration.

        This generator method is designed for direct equivalence with the
        spatially resolved generators found elsewhere in Landlab.

        This generator will be useful in cases where the whole sequence of
        storms and interstorms doesn't need to be stored, where we can save
        memory this way.

        This method does not attempt to subdivide timesteps. If you want that,
        provide a delta_t at instantiation, and use
        yield_storm_interstorm_duration_intensity(subdivide_interstorms=True).

        The method will keep yielding until it reaches the RUN_TIME, where it
        will terminate. If it terminates during a storm, the final tuple will
        be (truncated_storm_duration, 0.). Otherwise it will be
        (storm_duration, truncated_interstorm_duration.)

        Yields
        ------
        tuple of float
            (storm_duration, interstorm_duration)

        Notes
        -----
        One recommended procedure is to instantiate the generator, then call
        next(instance) (in Python 3) repeatedly to get the sequence (See
        Examples, below).

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((4, 5))
        >>> precip = PrecipitationDistribution(
        ...     mg,
        ...     mean_storm_duration=1.5,
        ...     mean_interstorm_duration=15.0,
        ...     mean_storm_depth=0.5,
        ...     total_t=46.0,
        ... )
        >>> storm_dts = []
        >>> interstorm_dts = []
        >>> intensities = []
        >>> precip.seed_generator(seedval=1)
        >>> for storm_dt, interstorm_dt in precip.yield_storms():
        ...     storm_dts.append(storm_dt)
        ...     interstorm_dts.append(interstorm_dt)
        ...     intensities.append(mg.at_grid["rainfall__flux"])
        ...
        >>> len(storm_dts) == 4  # 4 storms in the simulation
        True
        >>> len(interstorm_dts) == len(storm_dts)
        True
        >>> rf_intensities_to_test = np.array(
        ...     [
        ...         0.8138257984406472,
        ...         0.15929112025199238,
        ...         0.17254519305000884,
        ...         0.09817611240558813,
        ...     ]
        ... )
        >>> np.allclose(intensities, rf_intensities_to_test)
        True
        >>> np.isclose(sum(storm_dts) + sum(interstorm_dts), 46.0)  # total_t
        True
        >>> np.isclose(interstorm_dts[-1], 0.0)  # sequence truncated
        True

        An alternative way to use the generator might be:

        >>> # ^^this lets you "manually" get the next item from the iterator
        >>> flux = mg.at_grid.pop("rainfall__flux")  # remove the existing field
        >>> precip = PrecipitationDistribution(
        ...     mg,
        ...     mean_storm_duration=1.5,
        ...     mean_interstorm_duration=15.0,
        ...     mean_storm_depth=0.5,
        ...     total_t=46.0,
        ... )
        >>> precip.seed_generator(seedval=1)
        >>> mystorm_generator = precip.yield_storms()
        >>> my_list_of_storms = []
        >>> for i in range(4):
        ...     my_list_of_storms.append(next(mystorm_generator))
        ...

        Note that an exception is thrown if you go too far:

        >>> my_list_of_storms.append(next(mystorm_generator))  # the 5th iter
        Traceback (most recent call last):
          ...
        StopIteration

        Also note that the generator won't terminate if you try this without
        first instantiating the generator:

        >>> allmytimes = []
        >>> for i in range(20):  # this will run just fine
        ...     allmytimes.append(next(precip.yield_storms()))
        ...
        >>> total_t = sum([sum(storm) for storm in allmytimes])
        >>> total_t > 46.0
        True
        """
        # we must have instantiated with a grid, so check:
        assert hasattr(self, "_grid")

        # now exploit the existing generator to make this easier & less
        # redundant:
        delta_t = self._delta_t
        self._delta_t = None  # this is necessary to suppress chunking behaviour
        # in the other generator
        othergen = self.yield_storm_interstorm_duration_intensity()
        # enter a loop, to break as needed:
        tobreak = False
        while not tobreak:
            # we always start with a storm, so:
            try:
                (storm_dur, storm_int) = next(othergen)
            except StopIteration:
                break  # stop dead. We terminated at a good place
            try:
                (interstorm_dur, _) = next(othergen)
            except StopIteration:
                tobreak = True
                interstorm_dur = 0.0
            # reset the rainfall__flux field, that got overstamped in the
            # interstorm iter:
            self._grid.at_grid["rainfall__flux"] = storm_int
            yield (storm_dur, interstorm_dur)
        # now, just in case, restore self._delta_t:
        self._delta_t = delta_t

    def generate_from_stretched_exponential(self, scale, shape):
        """Generate and return a random variable from a stretched exponential
        distribution with given scale and shape.

        Examples
        --------
        >>> np.random.seed(0)
        >>> np.round(np.random.rand(3), 6)  # these are our 3 rand #s to test
        array([0.548814, 0.715189, 0.602763])
        >>> from landlab.components import PrecipitationDistribution
        >>> pd = PrecipitationDistribution(
        ...     mean_storm_duration=1.0,
        ...     mean_interstorm_duration=1.0,
        ...     mean_storm_depth=1.0,
        ... )
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

        1. That we can get the same sequence again by re-seeding with the
           same value (the default is zero)
        2. When we use a value other than the default, we get a different
           sequence

        Examples
        --------
        >>> precip = PrecipitationDistribution(
        ...     mean_storm_duration=1.5,
        ...     mean_interstorm_duration=15.0,
        ...     mean_storm_depth=0.5,
        ...     total_t=100.0,
        ...     delta_t=1.0,
        ... )
        >>> round(precip.storm_duration, 2)
        2.79
        >>> round(precip.interstorm_duration, 2)
        21.28
        >>> round(precip.storm_depth, 2)
        2.45
        >>> round(precip.intensity, 2)
        0.88
        >>> precip.seed_generator()  # re-seed and get same sequence again
        >>> round(precip.get_precipitation_event_duration(), 2)
        2.79
        >>> round(precip.get_interstorm_event_duration(), 2)
        21.28
        >>> round(precip.get_storm_depth(), 2)
        2.45
        >>> round(precip.get_storm_intensity(), 2)
        0.88
        >>> precip = PrecipitationDistribution(
        ...     mean_storm_duration=1.5,
        ...     mean_interstorm_duration=15.0,
        ...     mean_storm_depth=0.5,
        ...     total_t=100.0,
        ...     delta_t=1.0,
        ...     random_seed=1,
        ... )
        >>> round(precip.storm_duration, 2)  # diff't vals with diff't seed
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

    @property
    def intensity(self):
        """Get the intensity of the most recent storm simulated."""
        return self.get_storm_intensity()
