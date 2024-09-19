"""Landlab component that generates a random fire event in time.

This component generates a random fire event or fire time series from the
Weibull statistical distribution.

.. codeauthor:: Jordan Adams

This component generates random numbers using the Weibull distribution
(Weibull, 1951). No particular units must be used, but it was written with
the fire recurrence units in time (yrs).

Using the Weibull Distribution assumes two things: All elements within the
study area have the same fire regime. Each element must have (on average) a
constant fire regime during the time span of the study.

As of Sept. 2013, fires are considered instantaneous events independent of
other fire events in the time series.

Written by Jordan M. Adams, 2013. Updated April 2016.


Examples
--------
>>> from landlab.components.fire_generator import FireGenerator
>>> from landlab import RasterModelGrid

Create an instance of the FireGenerator component

>>> mg = RasterModelGrid((10, 10))
>>> fg = FireGenerator(mg, mean_fire_recurrence=15.0, shape_parameter=4.5)

This creates an instance of the component that has a mean_fire_recurrence, or
average interval between fires of 15 years. We gave it a shape parameter of
4.5, suggesting the cumulative distribution function that is skewed right.

Since we didn't pass it a scale parameter, the component calculates it for you.
Testing this...

>>> fg.scale_parameter
16.437036931437866

To get a time to next fire:

>>> fg.generate_fire_recurrence()  # doctest: +SKIP
10.68

References
----------
**Required Software Citation(s) Specific to this Component**

None Listed

**Additional References**

Polakow, D., Dunne, T. (1999). Modelling fire-return interval T: stochasticity
and censoring in the two-parameter Weibull model Ecological Modelling  121(1),
79-102. https://dx.doi.org/10.1016/s0304-3800(99)00074-5

"""

from random import weibullvariate

from scipy import special

from landlab import Component


class FireGenerator(Component):
    """Generate a random fire event or time series.

    Parameters
    ----------
    mean_fire_recurrence : float
        Average time between fires for a given location
    shape_parameter : float
        Describes the skew of the Weibull distribution.
        If shape < 3.5, data skews left.
        If shape == 3.5, data is normal.
        If shape > 3.5, data skews right.
        To approximate a normal bell curve, use a value of 3.5
    scale_parameter : float, optional
        Describes the peak of the Weibull distribution, located at the
        63.5% value of the cumulative distribution function. If unknown,
        it can be found using mean fire recurrence value and the
        get_scale_parameter() method described later.
    """

    _name = "FireGenerator"

    _unit_agnostic = True

    _info = {}

    def __init__(
        self, grid, mean_fire_recurrence=1.0, shape_parameter=3.5, scale_parameter=None
    ):
        """Generate a random fire event in time.

        Parameters
        ----------
        grid: landlab model grid
        mean_fire_recurrence : float
            Average time between fires for a given location
        shape_parameter : float
            Describes the skew of the Weibull distribution.
            If shape < 3.5, data skews left.
            If shape == 3.5, data is normal.
            If shape > 3.5, data skews right.
        scale_parameter : float, optional
            Describes the peak of the Weibull distribution, located at the
            63.5% value of the cumulative distribution function. If unknown,
            it can be found using mean fire recurrence value and the
            get_scale_parameter().
        """
        super().__init__(grid)
        self._mean_fire_recurrence = mean_fire_recurrence

        self._shape_parameter = shape_parameter

        if scale_parameter is None:
            self.get_scale_parameter()

        else:
            self._scale_parameter = scale_parameter

    @property
    def scale_parameter(self):
        """Scale parameter for the random distribution."""
        return self._scale_parameter

    def get_scale_parameter(self):
        """Get the scale parameter.

        ::
            mean_fire_recurrence = (scale_parameter * (
                                    special.gamma(1 + (1 / shape))))

        sets the scale parameter.
        """

        shape_in_gamma_func = float(1 + (1 / self._shape_parameter))
        gamma_func = special.gamma(shape_in_gamma_func)
        self._scale_parameter = self._mean_fire_recurrence / gamma_func

    def generate_fire_recurrence(self):
        """Get time to next fire.

        Finds the time to next fire (fire recurrence) based on the scale
        parameter (63.5% of fire Weibull distribution) and the shape parameter
        (describes the skew of the histogram, shape = 3.5
        represents a normal distribution).

        Rounds the time to next fire to 4 significant figures, for neatness.

        Returns
        -------
        float
            Updated value for the time to next fire.
        """
        self._time_to_next_fire = round(
            weibullvariate(self._scale_parameter, self._shape_parameter), 2
        )
        return self._time_to_next_fire
