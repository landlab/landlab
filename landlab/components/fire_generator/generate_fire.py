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

Create an instance of the FireGenerator component


>>> fg = FireGenerator(mean_fire_recurrence = 15.0, shape_parameter = 4.5)

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

"""

from random import weibullvariate
from scipy import special


class FireGenerator(object):

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

    def __init__(self, mean_fire_recurrence=0.0, shape_parameter=0.0,
                 scale_parameter=None):
        """Generate a random fire event in time.

        Parameters
        ----------
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

        self.mean_fire_recurrence = mean_fire_recurrence

        self.shape_parameter = shape_parameter

        if scale_parameter is None:
            self.get_scale_parameter()

        else:
            self.scale_parameter = scale_parameter


    def get_scale_parameter(self):
        """Get the scale parameter.

        ::
            mean_fire_recurrence = (scale_parameter * (
                                    special.gamma(1 + (1 / shape))))

        sets the scale parameter.
        """

        shape_in_gamma_func = float( 1+ (1 / self.shape_parameter))
        gamma_func = special.gamma(shape_in_gamma_func)
        self.scale_parameter = (self.mean_fire_recurrence / gamma_func)


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
        self.time_to_next_fire = round(weibullvariate(self.scale_parameter,
                                                      self.shape_parameter), 2)
        return self.time_to_next_fire
