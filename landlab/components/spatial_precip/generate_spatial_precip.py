
import numpy as np
import os
import inspect
from six.moves import range
from matplotlib.pyplot import figure
from scipy.stats import genextreme, fisk
from landlab import RasterModelGrid, CLOSED_BOUNDARY, Component


class PrecipitationDistribution(Component):

    _name = 'PrecipitationDistribution'

    _input_var_names = (
        'topographic__elevation',
    )

    _output_var_names = (
        'rainfall__flux',
        'rainfall__total_depth_per_year',
    )

    _var_units = {
        'topographic__elevation': 'm',
        'rainfall__flux': 'mm/hr',
        'rainfall__total_depth_per_year': 'mm/yr',
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'rainfall__flux': 'node',
        'rainfall__total_depth_per_year': 'node',
    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'rainfall__flux':
            'Depth of water delivered per unit time in each storm',
        'rainfall__total_depth_per_year':
            'Depth of water delivered in total in each model year',
    }

    def __init__(self, grid, number_of_years=1, orographic_scenario=None):
        """
        A component to generate a sequence of spatially resolved storms over a
        grid, following a lightly modified version (see below) of the
        stochastic methods of Singer & Michaelides, Env Res Lett 12, 104011,
        2017, & Singer et al., Geosci. Model Dev., submitted.

        The method is heavily stochastic, and at the present time is intimately
        calibrated against the conditions at Walnut Gulch, described in those
        papers. In particular, assumptions around intensity-duration
        calibration and orographic rainfall are "burned in" for now, and are
        not accessible to the user. The various probability distributions
        supplied to the various run methods default to WG values, but are
        easily modified.

        The component has two ways of simulating a "year". This choice is
        controlled by the 'limit' parameter of the yield methods. If limit==
        'total_rainfall', the component will continue to run until the total
        rainfall for the season and/or year exceeds a stochastically generated
        value. This method is directly comparable to the Singer & Michaelides
        method, but will almost always result in years which are not one
        calendar year long, unless the input distributions are very carefully
        recalibrated for each use case. If limit=='total_time', the component
        will terminate a season and/or year once the elapsed time exceeds one
        year. In this case, the total rainfall will not correspond to the
        stochastically generated total. You can access the actual total for the
        last season using the property `(median_)total_rainfall_last_season`.

        Key methods are:
        yield_storms
            Generate a timeseries of storm:interstorm duration pairs, alongside
            a field that describes the spatial distribution of rain during that
            storm.
        yield_years
            Generate a timeseries of ints giving number of storms per year,
            alongside a field that describes the spatial distribution of total
            rainfall across that year.
        yield_seasons
            Generate a timeseries of ints giving number of storms per season,
            alongside a field that describes the spatial distribution of total
            rainfall across that season.
        calc_annual_rainfall
            Produce a timeseries of tuples giving total rainfall each season,
            without resolving the storms spatially (i.e., fast!).

        A large number of properties are available to access storm properties
        during generation:
        - current_year
        - current_season
        - storm_depth_last_storm
        - storm_recession_value_last_storm
        - storm_duration_last_storm
        - storm_area_last_storm
        - storm_intensity_last_storm
        - total_rainfall_this_season
        - total_rainfall_this_year
        - total_rainfall_last_season
        - total_rainfall_last_year
        - median_total_rainfall_this_season
        - median_total_rainfall_this_year
        - median_total_rainfall_last_season
        - median_total_rainfall_last_year
        - number_of_nodes_under_storm
        - nodes_under_storm

        Note that becuase these are medians not means,
        median_total_rainfall_last_season + median_total_rainfall_this_season
        != median_total_rainfall_this_year.

        Significant differences between this component and the Singer code are:
        - The component does not model evapotranspiration. Use a separate
            Landlab component for this.
        - The component runs only over a LL grid; there is no such thing as a
            validation or simulation run.
        - It produces "fuzz" around intensity values using a continuous
            distribution; Singer does this with integer steps.
        - Step changes mid-run cannot be explicitly modelled. Instead, run the
            component for a fixed duration, make the change to the
            distribution input parameter, then run it again.
        - Storms can be centred at any point, not just over nodes.
        - Edge buffering is now dynamic; i.e., big storms have a bigger edge
            buffer than smaller storms. Storms can be centred off the grid
            edges.
        - Storms are never discarded - once a storm is drawn, it must hit the
            catchment, and positions are repeatedly selected until this can
            happen. Singer's method would discard such a storm and draw a new
            one.
        - Durations are not rescaled to ensure both total duration and total
            precip are both satisfied at the same time, as in Singer's method.
            Instead, the component either matches a year's duration, *or*
            exactly a year's worth of rain. This choice is dictated by the
            `limit` parameter in the yield methods.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab model grid of any type.
        number_of_years : int
            The number of years over which to generate storms.
        orographic_scenario : {None, 'Singer'}
            Whether to use no orographic rule, or to adopt S&M's 2017
            calibration for Walnut Gulch.

        """
        self._grid = grid
        gaugecount = (grid.status_at_node != CLOSED_BOUNDARY).sum()
        self._gauge_dist_km = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots1 = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots2 = np.zeros(gaugecount, dtype='float')
        self._numyrs = number_of_years

        self._max_numstorms = 5000
        # This is for initializing matrices. Trailing zeros are deleted from
        # matrixes at the end of the code.

        assert orographic_scenario in (None, 'Singer')
        self._orographic_scenario = orographic_scenario

        # build LL fields:
        self.initialize_output_fields()
        # bind the field to the internal variable:
        self._rain_int_gauge = self.grid.at_node['rainfall__flux']
        self._total_rf_year = self.grid.at_node[
            'rainfall__total_depth_per_year']

        # store some info on the open node grid extent:
        open_nodes = self.grid.status_at_node != CLOSED_BOUNDARY
        self._minx = self.grid.node_x[open_nodes].min()
        self._maxx = self.grid.node_x[open_nodes].max()
        self._miny = self.grid.node_y[open_nodes].min()
        self._maxy = self.grid.node_y[open_nodes].max()
        self._widthx = self._maxx - self._minx
        self._widthy = self._maxy - self._miny
        self._running_total_rainfall_this_year = self.grid.zeros('node')
        self._running_total_rainfall_this_season = self.grid.zeros('node')

    def yield_storms(self, limit='total_time', style='whole_year',
                     total_rf_trend=0., storminess_trend=0.,
                     monsoon_fraction_of_year=0.42,
                     monsoon_total_rf_gaussian={
                         'sigma': 64., 'mu': 207.},
                     monsoon_storm_duration_GEV={
                         'shape': -0.570252, 'sigma': 35.7389, 'mu': 34.1409,
                         'trunc_interval': (1., 1040.)},
                     monsoon_storm_area_GEV={
                         'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                         'trunc_interval': (5.e+06, 3.e+08)},
                     monsoon_storm_interarrival_GEV={
                         'shape': -0.807971, 'sigma': 9.49574, 'mu': 10.6108,
                         'trunc_interval': (0., 120.)},
                     monsoon_storm_radial_weakening_gaussian={
                         'sigma': 0.08, 'mu': 0.25,
                         'trunc_interval': (0.15, 0.67)},
                     winter_total_rf_gaussian={
                         'sigma': 52., 'mu': 1.65},
                     winter_storm_duration_fisk={
                         'c': 1.0821, 'scale': 68.4703,
                         'trunc_interval': (1., 5000.)},
                     winter_storm_area_GEV={
                         'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                         'trunc_interval': (5.e+06, 3.e+08)},
                     winter_storm_interarrival_GEV={
                         'shape': 1.1131, 'sigma': 53.2671, 'mu': 47.4944,
                         'trunc_interval': (1., 720.)},
                     winter_storm_radial_weakening_gaussian={
                         'sigma': 0.08, 'mu': 0.25,
                         'trunc_interval': (0.15, 0.67)}):
        """
        Yield a timeseries giving the number if storms occurring each year in
        a rainfall simulation.

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, 2018 & Singer et al, submitted.

        Parameters
        ----------
        limit : ('total_time', 'total_rainfall')
            Controls whether a season is defined based on its total rainfall
            (and can be any length), or by its duration (and can have any
            amount of rainfall). If 'total_time', monsoon_fraction_of_year
            sets the fraction of a year occupied by the monsoon.
        style : ('whole_year', 'monsoonal', 'winter')
            Controls whether the component seeks to simulate a western US-
            style "monsoonal" climate, a western US-style winter climate,
            or a full year combining both. These distributions are by default
            based on Singer et al.'s calibrations. Note if 'monsoonal',
            the total duration of a "year" will appear to be only
            `monsoon_fraction_of_year`, and the opposite for `winter`.
        total_rf_trend : float
            Controls if a drift is applied to the total rainfall distribution
            through time. If 0., no trend. If positive, rainfall totals
            increase gradually through time. If negative, they fall through
            time. S&M recommend +/- 0.07 for a realistic climate change driven
            drift at Walnut Gulch.
        storminess_trend : float
            Controls if a drift is applied to the expected intensity of
            individual storms through time. If 0., no trend. If positive,
            storms get more intense through time, if negative, less so. S&M
            recommend +/- 0.01 for a realistic climate change driven drift at
            Walnut Gulch.
        monsoon_fraction_of_year : float
            If limit == 'total_time', sets the fraction of one year occupied
            by the monsoon season. If not, ignored. Singer's monsoon runs from
            May to September, inclusive, and the default reflects this.

       monsoon_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm. In MIN.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling. In MIN.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        winter_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        Yields
        ------
        (storm_t, interval_t) : (float, float)
            Tuple pair of duration of a single storm, then the interstorm
            interval that follows it. In hrs. The rainfall__flux field
            describes the rainfall rate during the interval storm_t as the
            tuple is yielded. In HRS.
            Note that the rainfall__total_depth_per_year field gives the total
            accumulated rainfall depth during the *last completed* model year,
            not the year to the point of yield. For the latter, use the
            property `total_rainfall_this_year`.
        """
        return self._run_the_process(
            yield_storms=True, yield_years=False, yield_seasons=False,
            limit=limit, style=style,
            monsoon_fraction_of_year=monsoon_fraction_of_year,
            total_rf_trend=total_rf_trend, storminess_trend=storminess_trend,
            monsoon_total_rf_gaussian=monsoon_total_rf_gaussian,
            monsoon_storm_duration_GEV=monsoon_storm_duration_GEV,
            monsoon_storm_area_GEV=monsoon_storm_area_GEV,
            monsoon_storm_interarrival_GEV=monsoon_storm_interarrival_GEV,
            monsoon_storm_radial_weakening_gaussian=
                                       monsoon_storm_radial_weakening_gaussian,
            winter_total_rf_gaussian=winter_total_rf_gaussian,
            winter_storm_duration_fisk=winter_storm_duration_fisk,
            winter_storm_area_GEV=winter_storm_area_GEV,
            winter_storm_interarrival_GEV=winter_storm_interarrival_GEV,
            winter_storm_radial_weakening_gaussian=
                                        winter_storm_radial_weakening_gaussian)

    def yield_years(self, limit='total_time', style='whole_year',
                    total_rf_trend=0., storminess_trend=0.,
                    monsoon_fraction_of_year=0.42,
                    monsoon_total_rf_gaussian={
                        'sigma': 64., 'mu': 207.},
                    monsoon_storm_duration_GEV={
                        'shape': -0.570252, 'sigma': 35.7389, 'mu': 34.1409,
                        'trunc_interval': (1., 1040.)},
                    monsoon_storm_area_GEV={
                        'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                        'trunc_interval': (5.e+06, 3.e+08)},
                    monsoon_storm_interarrival_GEV={
                        'shape': -0.807971, 'sigma': 9.49574, 'mu': 10.6108,
                        'trunc_interval': (0., 120.)},
                    monsoon_storm_radial_weakening_gaussian={
                        'sigma': 0.08, 'mu': 0.25,
                        'trunc_interval': (0.15, 0.67)},
                    winter_total_rf_gaussian={
                        'sigma': 52., 'mu': 1.65},
                    winter_storm_duration_fisk={
                        'c': 1.0821, 'scale': 68.4703,
                        'trunc_interval': (1., 5000.)},
                    winter_storm_area_GEV={
                        'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                        'trunc_interval': (5.e+06, 3.e+08)},
                    winter_storm_interarrival_GEV={
                        'shape': 1.1131, 'sigma': 53.2671, 'mu': 47.4944,
                        'trunc_interval': (1., 720.)},
                    winter_storm_radial_weakening_gaussian={
                        'sigma': 0.08, 'mu': 0.25,
                        'trunc_interval': (0.15, 0.67)}):
        """
        Yield a timeseries giving the number if storms occurring each year in
        a rainfall simulation.

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, 2018 & Singer et al, submitted.

        Parameters
        ----------
        limit : ('total_time', 'total_rainfall')
            Controls whether a season is defined based on its total rainfall
            (and can be any length), or by its duration (and can have any
            amount of rainfall). If 'total_time', monsoon_fraction_of_year
            sets the fraction of a year occupied by the monsoon.
        style : ('whole_year', 'monsoonal', 'winter')
            Controls whether the component seeks to simulate a western US-
            style "monsoonal" climate, a western US-style winter climate,
            or a full year combining both. These distributions are by default
            based on Singer et al.'s calibrations. Note if 'monsoonal',
            the total duration of a "year" will appear to be only
            `monsoon_fraction_of_year`, and the opposite for `winter`.
        total_rf_trend : float
            Controls if a drift is applied to the total rainfall distribution
            through time. If 0., no trend. If positive, rainfall totals
            increase gradually through time. If negative, they fall through
            time. S&M recommend +/- 0.07 for a realistic climate chage driven
            drift at Walnut Gulch.
        storminess_trend : float
            Controls if a drift is applied to the expected intensity of
            individual storms through time. If 0., no trend. If positive,
            storms get more intense through time, if negative, less so. S&M
            recommend +/- 0.01 for a realistic climate change driven drift at
            Walnut Gulch.
        monsoon_fraction_of_year : float
            If limit == 'total_time', sets the fraction of one year occupied
            by the monsoon season. If not, ignored. Singer's monsoon runs from
            May to September, inclusive, and the default reflects this.

        monsoon_total_rf_gaussian is a normal distribution controlling the
            total rainfall expected in each year. S&M use 'mu' in {143., 271.}
            for step changes up/down in rainfall totals. In mm.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm. In MIN.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling. In MIN.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        winter_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        Yields
        ------
        number_of_storms_per_year : float
            Float that gives the number of storms simulated in the year that
            elapsed since the last yield. The rainfall__total_depth_per_year
            field gives the total accumulated rainfall depth during the year
            preceding the yield. rainfall__flux gives the rainfall intensity of
            the last storm in that year.
        """
        return self._run_the_process(
            yield_storms=False, yield_years=True, yield_seasons=False,
            limit=limit, style=style,
            total_rf_trend=total_rf_trend, storminess_trend=storminess_trend,
            monsoon_fraction_of_year=monsoon_fraction_of_year,
            monsoon_total_rf_gaussian=monsoon_total_rf_gaussian,
            monsoon_storm_duration_GEV=monsoon_storm_duration_GEV,
            monsoon_storm_area_GEV=monsoon_storm_area_GEV,
            monsoon_storm_interarrival_GEV=monsoon_storm_interarrival_GEV,
            monsoon_storm_radial_weakening_gaussian=
                                       monsoon_storm_radial_weakening_gaussian,
            winter_total_rf_gaussian=winter_total_rf_gaussian,
            winter_storm_duration_fisk=winter_storm_duration_fisk,
            winter_storm_area_GEV=winter_storm_area_GEV,
            winter_storm_interarrival_GEV=winter_storm_interarrival_GEV,
            winter_storm_radial_weakening_gaussian=
                                        winter_storm_radial_weakening_gaussian)

    def yield_seasons(self, limit='total_time', style='whole_year',
                      total_rf_trend=0., storminess_trend=0.,
                      monsoon_fraction_of_year=0.42,
                      monsoon_total_rf_gaussian={
                          'sigma': 64., 'mu': 207.},
                      monsoon_storm_duration_GEV={
                          'shape': -0.570252, 'sigma': 35.7389, 'mu': 34.1409,
                          'trunc_interval': (1., 1040.)},
                      monsoon_storm_area_GEV={
                          'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                          'trunc_interval': (5.e+06, 3.e+08)},
                      monsoon_storm_interarrival_GEV={
                          'shape': -0.807971, 'sigma': 9.49574, 'mu': 10.6108,
                          'trunc_interval': (0., 120.)},
                      monsoon_storm_radial_weakening_gaussian={
                          'sigma': 0.08, 'mu': 0.25,
                          'trunc_interval': (0.15, 0.67)},
                      winter_total_rf_gaussian={
                          'sigma': 52., 'mu': 1.65},
                      winter_storm_duration_fisk={
                          'c': 1.0821, 'scale': 68.4703,
                          'trunc_interval': (1., 5000.)},
                      winter_storm_area_GEV={
                          'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                          'trunc_interval': (5.e+06, 3.e+08)},
                      winter_storm_interarrival_GEV={
                          'shape': 1.1131, 'sigma': 53.2671, 'mu': 47.4944,
                          'trunc_interval': (1., 720.)},
                      winter_storm_radial_weakening_gaussian={
                          'sigma': 0.08, 'mu': 0.25,
                          'trunc_interval': (0.15, 0.67)}):
        """
        Yield a timeseries giving the number if storms occurring each season in
        a rainfall simulation. Only meaningfully different from yield_years if
        style=='whole_year'.

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, 2018 & Singer et al, submitted.

        Parameters
        ----------
        limit : ('total_time', 'total_rainfall')
            Controls whether a season is defined based on its total rainfall
            (and can be any length), or by its duration (and can have any
            amount of rainfall). If 'total_time', monsoon_fraction_of_year
            sets the fraction of a year occupied by the monsoon.
        style : ('whole_year', 'monsoonal', 'winter')
            Controls whether the component seeks to simulate a western US-
            style "monsoonal" climate, a western US-style winter climate,
            or a full year combining both. These distributions are by default
            based on Singer et al.'s calibrations. Note if 'monsoonal',
            the total duration of a "year" will appear to be only
            `monsoon_fraction_of_year`, and the opposite for `winter`.
        total_rf_trend : float
            Controls if a drift is applied to the total rainfall distribution
            through time. If 0., no trend. If positive, rainfall totals
            increase gradually through time. If negative, they fall through
            time. S&M recommend +/- 0.07 for a realistic climate chage driven
            drift at Walnut Gulch.
        storminess_trend : float
            Controls if a drift is applied to the expected intensity of
            individual storms through time. If 0., no trend. If positive,
            storms get more intense through time, if negative, less so. S&M
            recommend +/- 0.01 for a realistic climate change driven drift at
            Walnut Gulch.
        monsoon_fraction_of_year : float
            If limit == 'total_time', sets the fraction of one year occupied
            by the monsoon season. If not, ignored. Singer's monsoon runs from
            May to September, inclusive, and the default reflects this.

        monsoon_total_rf_gaussian is a normal distribution controlling the
            total rainfall expected in each year. S&M use 'mu' in {143., 271.}
            for step changes up/down in rainfall totals. In mm.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm. In MIN.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling. In MIN.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        winter_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        Yields
        ------
        number_of_storms_per_season : float
            Float that gives the number of storms simulated in the season that
            elapsed since the last yield. The rainfall__total_depth_per_year
            field gives the total accumulated rainfall depth during the *year*
            preceding the yield, *so far*. rainfall__flux gives the rainfall
            intensity of the last storm in that year.
        NB: Use the component property total_rainfall_last_season to access
            the *actual* amount of rainfall in the season that has the number
            of storms that the method generates.
        """
        return self._run_the_process(
            yield_storms=False, yield_years=False, yield_seasons=True,
            limit=limit, style=style,
            total_rf_trend=total_rf_trend, storminess_trend=storminess_trend,
            monsoon_fraction_of_year=monsoon_fraction_of_year,
            monsoon_total_rf_gaussian=monsoon_total_rf_gaussian,
            monsoon_storm_duration_GEV=monsoon_storm_duration_GEV,
            monsoon_storm_area_GEV=monsoon_storm_area_GEV,
            monsoon_storm_interarrival_GEV=monsoon_storm_interarrival_GEV,
            monsoon_storm_radial_weakening_gaussian=
                                       monsoon_storm_radial_weakening_gaussian,
            winter_total_rf_gaussian=winter_total_rf_gaussian,
            winter_storm_duration_fisk=winter_storm_duration_fisk,
            winter_storm_area_GEV=winter_storm_area_GEV,
            winter_storm_interarrival_GEV=winter_storm_interarrival_GEV,
            winter_storm_radial_weakening_gaussian=
                                        winter_storm_radial_weakening_gaussian)

    def _run_the_process(self, yield_storms=True, yield_years=False,
                         yield_seasons=False, limit='total_time',
                         style='whole_year', monsoon_fraction_of_year=0.42,
                         total_rf_trend=0., storminess_trend=0.,
                         monsoon_total_rf_gaussian={
                             'sigma': 64., 'mu': 207.},
                         monsoon_storm_duration_GEV={
                             'shape': -0.570252, 'sigma': 35.7389,
                             'mu': 34.1409, 'trunc_interval': (1., 1040.)},
                         monsoon_storm_area_GEV={
                             'shape': 0., 'sigma': 2.83876e+07,
                             'mu': 1.22419e+08,
                             'trunc_interval': (5.e+06, 3.e+08)},
                         monsoon_storm_interarrival_GEV={
                             'shape': -0.807971, 'sigma': 9.49574,
                             'mu': 10.6108, 'trunc_interval': (0., 120.)},
                         monsoon_storm_radial_weakening_gaussian={
                             'sigma': 0.08, 'mu': 0.25,
                             'trunc_interval': (0.15, 0.67)},
                         winter_total_rf_gaussian={
                             'sigma': 52., 'mu': 1.65},
                         winter_storm_duration_fisk={
                             'c': 1.0821, 'scale': 68.4703,
                             'trunc_interval': (1., 5000.)},
                         winter_storm_area_GEV={
                             'shape': 0., 'sigma': 2.83876e+07,
                             'mu': 1.22419e+08,
                             'trunc_interval': (5.e+06, 3.e+08)},
                         winter_storm_interarrival_GEV={
                             'shape': 1.1131, 'sigma': 53.2671, 'mu': 47.4944,
                             'trunc_interval': (1., 720.)},
                         winter_storm_radial_weakening_gaussian={
                             'sigma': 0.08, 'mu': 0.25,
                             'trunc_interval': (0.15, 0.67)}):
        """
        This is the underlying process that runs the component, but it should
        be run by a user through the yield_storms and yield_years methods.

        Fuzz to the chosen values is now selected from a continuous
        distribution, not from integer values.

        total_rf_trend controls if a drift is applied to the total rainfall
        distribution through time. If 0., no trend. If positive, rainfall
        totals increase gradually through time. If negative, they fall through
        time. S&M recommend +/- 0.07 for a realistic climate chage driven drift
        at Walnut Gulch.

        storminess_trend controls if a drift is applied to the expected
        intensity of individual storms through time. If 0., no trend. If
        positive, storms get more intense through time, if negative, less so.
        S&M recommend +/- 0.01 for a realistic climate change driven drift at
        Walnut Gulch.

        All default distributions reflect values for Walnut Gulch, see Singer &
        Michaelides, submitted:

        monsoon_total_rf_gaussian is a normal distribution controlling the
            total rainfall expected in each year. S&M use 'mu' in {143., 271.}
            for step changes up/down in rainfall totals. In mm.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm. In MIN.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.
        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling. In MIN.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS.
        winter_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        """
        FUZZMETHOD = 'DEJH'
        FUZZWIDTH = 5.  # if DEJH
        self._phantom_storm_count = 0
        # ^this property tracks the number of storms in the run that received
        # zero intensity (and thus didn't really exist)
        self._opennodes = self.grid.status_at_node != CLOSED_BOUNDARY
        self._total_rainfall_last_season = self.grid.zeros('node')
        
        # safety check for init conds:
        if yield_storms:
            assert yield_years is False
            assert yield_seasons is False
        if yield_years:
            assert yield_storms is False
            assert yield_seasons is False
        if yield_seasons:
            assert yield_storms is False
            assert yield_years is False

        # add variable for number of simulations of simyears
        simyears = self._numyrs  # number of years to simulate
        numcurves = 11  # number of intensity-duration curves (see below for
        # curve equations)
        hrsinyr = 24.*365.
        hrsinmonsoon = monsoon_fraction_of_year * hrsinyr
        hrsinwinter = (1. - monsoon_fraction_of_year) * hrsinyr

        assert limit in ('total_rainfall', 'total_time')

        assert style in ('whole_year', 'monsoonal', 'winter')
        if style == 'whole_year':
            reps = 2
        else:
            reps = 1

        opennodes = self._opennodes
        num_opennodes = np.sum(opennodes)
        IDs_open = np.where(opennodes)[0]  # need this later
        X1 = self.grid.node_x
        Y1 = self.grid.node_y
        Xin = X1[opennodes]
        Yin = Y1[opennodes]
        Zz = self.grid.at_node['topographic__elevation'][opennodes]
        numgauges = Xin.size  # number of rain gauges in the basin.
        # NOTE: In this version this produces output on a grid, rather than at
        # real gauge locations.

        if FUZZMETHOD == 'MS':
            raise NameError(
                'The Singer method for fuzz is no longer supported')

        # NOTE this is overly specific
        if self._orographic_scenario == 'Singer':
            # These are elevation ranges for the 3 orographic groups
            OroGrp1 = np.arange(int(np.round(Zz.min())), 1350)
            OroGrp2 = np.arange(1350, 1500)
            OroGrp3 = np.arange(1500, int(np.round(Zz.max())))

        # lambda_, kappa, and C are parameters of the intensity-duration curves
        # of the form: intensity =
        # lambda*exp(-0.508*duration)+kappa*exp(-0.008*duration)+C
        lambda_ = [642.2, 578.0, 513.8, 449.5, 385.3, 321.1, 256.9, 192.7,
                   128.4, 64.1, 21.0]
        kappa = [93.1, 83.8, 74.5, 65.2, 55.9, 46.6, 37.2, 27.9, 18.6, 9.3,
                 0.9]
        C = [4.5, 4., 3.5, 3., 2.5, 2., 1.5, 1., 0.5, 0.25, 0.05]

        # Unlike MS's original implementation, we no longer pull ET values, as
        # this should be a different component.

        self._Ptot_ann_global = np.zeros(simyears)
        self._Ptot_monsoon_global = np.zeros(simyears)

        Intensity_local_all = 0  # initialize all variables (concatenated
        # matrices of generated output)
        Storm_totals_all = 0
        Duration_local_all = 0

        master_storm_count = 0
        last_year_count = 0
        storm_trend = 0

        for syear in range(simyears):
            self._year = syear
            year_time = 0.  # tracks simulation time per year in hours
            storm_trend += storminess_trend
            Ptotal = 0.
            year_storm_count = 0
            breaker = False
            Storm_total_local_year = np.zeros(
                (self._max_numstorms, num_opennodes))
            self._storm_running_sum_of_seasons = np.zeros(num_opennodes)
            self._storm_running_sum_1st_seas = np.zeros(num_opennodes)

            storms_yr_so_far = 0
            for seas in range(reps):
                seas_time = 0.  # tracks elapsed season time in hours
                Storm_running_sum_seas = np.zeros((2, num_opennodes))
                # ^ 1st col is running total, 2nd is data to add to it
                if seas == 0 and not style == 'winter':
                    self._current_season = 'M'
                    # This is the pdf fitted to all available station precip
                    # data (normal dist). It will be sampled below.
                    Ptot_pdf_norm = monsoon_total_rf_gaussian
                    # This step change case can now be handled direct from
                    # input flags
                    # if self._ptot_scenario == 'ptot+':
                    #     Ptot_pdf_norm = {'sigma': 64., 'mu': 271.}
                    # elif self._ptot_scenario == 'ptot-':
                    #     Ptot_pdf_norm = {'sigma': 64., 'mu': 143.}
                    # else:
                    #     Ptot_pdf_norm = {'sigma': 64., 'mu': 207.}
                    # the trending cases need to be handled in the loop

                    # This is the pdf fitted to all available station duration
                    # data (GEV dist). It will be sampled below.
                    # #### matlab's GEV is (shape_param, scale(sigma), pos(mu))
                    # note that in Scipy, we must add a minus to the shape
                    # param for a GEV to match Matlab's implementation
                    Duration_pdf = monsoon_storm_duration_GEV
                    # This is the pdf fitted to all available station area
                    # data (EV dist). It will be sampled below.
                    # #### matlab's EV is (mu, sigma)
                    Area_pdf_EV = monsoon_storm_area_GEV
                    # This is the pdf fitted to all available station area
                    # data (GEV dist). It will be sampled below.
                    Int_arr_pdf_GEV = monsoon_storm_interarrival_GEV
                    # This is the pdf of storm gradient recession coefficients
                    # from Morin et al, 2005 (normal dist). It will be sampled
                    # below.
                    Recess_pdf_norm = monsoon_storm_radial_weakening_gaussian
                    seas_total = hrsinmonsoon
                else:
                    self._current_season = 'W'
                    Ptot_pdf_norm = winter_total_rf_gaussian
                    Duration_pdf = winter_storm_duration_fisk
                    Area_pdf_EV = winter_storm_area_GEV
                    Int_arr_pdf_GEV = winter_storm_interarrival_GEV
                    Recess_pdf_norm = winter_storm_radial_weakening_gaussian
                    seas_total = hrsinwinter

                if not np.isclose(total_rf_trend, 0.):
                    mu = Ptot_pdf_norm.pop('mu')
                    mu += mu * total_rf_trend
                    Ptot_pdf_norm['mu'] = mu
                # sample from normal distribution and saves global value of
                # Ptot (that must be equalled or exceeded) for each year
                season_rf_limit = self.calc_annual_rainfall(
                    style=style, monsoon_total_rf_gaussian=Ptot_pdf_norm)[seas]
                self._Ptot_ann_global[syear] += season_rf_limit
                if seas == 0 and not style == 'winter':
                    self._Ptot_monsoon_global[syear] = season_rf_limit
                Storm_total_local_seas = np.zeros(
                    (self._max_numstorms, num_opennodes))
                seas_cum_Ptot_gauge = np.zeros(numgauges)
                self._entries = 0
                seas_storm_count = 0

                for storm in range(self._max_numstorms):
                    self._rain_int_gauge.fill(0.)
                    int_arr_val = genextreme.rvs(
                        c=Int_arr_pdf_GEV['shape'], loc=Int_arr_pdf_GEV['mu'],
                        scale=Int_arr_pdf_GEV['sigma'])
                    try:
                        int_arr_val = np.clip(
                            int_arr_val, Int_arr_pdf_GEV['trunc_interval'][0],
                            Int_arr_pdf_GEV['trunc_interval'][1])
                    except KeyError:
                        # ...just in case
                        if int_arr_val < 0.:
                            int_arr_val = 0.
                    self._int_arr_val = int_arr_val
                    # ^Samples from distribution of interarrival times (hr).
                    # This can be used to develop STORM output for use in
                    # rainfall-runoff models or any water balance application.
                    # sample uniformly from storm center matrix from grid w
                    # 10 m spacings covering basin:

                    area_val = genextreme.rvs(c=Area_pdf_EV['shape'],
                                              loc=Area_pdf_EV['mu'],
                                              scale=Area_pdf_EV['sigma'])
                    try:
                        area_val = np.clip(
                            area_val, Area_pdf_EV['trunc_interval'][0],
                            Area_pdf_EV['trunc_interval'][1])
                    except KeyError:
                        # ...just in case
                        if area_val < 0.:
                            area_val = 0.
                    self._area_val = area_val
                    # ^Samples from distribution of storm areas

                    r = np.sqrt(area_val/np.pi)  # value here shd be selected
                    rsq = r**2
                    # based on area above in meters to match the UTM values

                    # This way of handling storm locations is really quite
                    # different to MS's. He uses a fixed buffer width, and
                    # throws away any storm that doesn't intersect. We
                    # instead retain all storms, and *make sure* the storm
                    # intersects using a dynamic buffer. MS's method will
                    # preferentially sample larger storms, though unclear
                    # what that would mean in practice.
                    # MS also snaps his storms onto the grid. This seems
                    # unnecessary, and we don't do it here.
                    while 1:
                        cx, cy = self._locate_storm(r)
                        # Determine which gauges are hit by Euclidean geometry:
                        gdist = (Xin-cx)**2 + (Yin-cy)**2
                        mask_name = (gdist <= rsq)  # this is defacto MS's aa
                        # this short circuits the storm loop in the case that
                        # the storm does not affect any 'gauging' location
                        if np.any(np.equal(mask_name, True)):
                            break

                    self._x = cx
                    self._y = cy
                    year_storm_count += 1
                    seas_storm_count += 1
                    master_storm_count += 1
                    gauges_hit = np.where(mask_name)[0]
                    num_gauges_hit = gauges_hit.size
                    # this routine below allows for orography in precip by
                    # first determining the closest gauge and then determining
                    # its orographic grouping
                    cc = np.argmin(gdist)
                    closest_gauge = np.round(Zz[cc])  # this will be compared
                    # against orographic gauge groupings to determine the
                    # appropriate set of intensity-duration curves
                    ######

                    # This routine below determines to which orographic group
                    # the closest gauge to the storm center belongs to, and
                    # censors the number of curves accordingly
                    # missing top curve in GR1, top and bottom curves for GR2,
                    # and bottom curve for GR3
                    # new version of orography compares local 'gauge' elevation
                    # to elevation bands called OroGrp, defined above
                    # NOTE again, DEJH thinks this could be generalised a lot
                    if self._orographic_scenario == 'Singer':
                        if closest_gauge in OroGrp1:
                            baa = 'a'
                        elif closest_gauge in OroGrp2:
                            baa = 'b'
                        elif closest_gauge in OroGrp3:
                            baa = 'c'
                        else:
                            raise ValueError(
                                'closest_gauge not found in curve lists!')
                    elif self._orographic_scenario is None:
                        baa = None
                    if seas == 0 and not style == 'winter':
                        duration_val = genextreme.rvs(
                            c=Duration_pdf['shape'], loc=Duration_pdf['mu'],
                            scale=Duration_pdf['sigma'])
                    else:
                        duration_val = fisk.rvs(
                            c=Duration_pdf['c'], scale=Duration_pdf['scale'])
                    # hacky fix to prevent occasional < 0 values:
                    # (I think because Matlab is able to set limits manually)
                    try:
                        duration_val = np.clip(
                            duration_val, Duration_pdf['trunc_interval'][0],
                            Duration_pdf['trunc_interval'][1])
                    except KeyError:
                        # ...just in case
                        if duration_val < 0.:
                            duration_val = 0.
                    durationhrs = duration_val / 60.
                    self._durationhrs = durationhrs
                    year_time += durationhrs
                    seas_time += durationhrs
                    # we will always do the next storm, even if it exceeds the
                    # specified "total" time

                    # original curve# probs for 30%-20%-10%: [0.0636, 0.0727,
                    # 0.0819, 0.0909, 0.0909, 0.0909, 0.0909, 0.0909, 0.1001,
                    # 0.1090, 0.1182]
                    # original curve# probs are modified as below
                    # add weights to reflect reasonable probabilities that
                    # favor lower curves:
                    if baa == 'a':
                        wgts = [0.0318, 0.0759, 0.0851, 0.0941, 0.0941, 0.0941,
                                0.0941, 0.0941, 0.1033, 0.1121, 0.1213]
                    elif baa == 'b':
                        wgts = [0.0478, 0.0778, 0.0869, 0.0959, 0.0959, 0.0959,
                                0.0959, 0.0959, 0.1051, 0.1141, 0.0888]
                    elif baa == 'c':
                        wgts = [0.0696, 0.0786, 0.0878, 0.0968, 0.0968, 0.0968,
                                0.0968, 0.0968, 0.1060, 0.1149, 0.0591]
                    elif baa is None:
                        wgts = [0.0636, 0.0727, 0.0819, 0.0909, 0.0909, 0.0909,
                                0.0909, 0.0909, 0.1001, 0.1090, 0.1182]
                    # which curve did we pick?:
                    int_dur_curve_val = np.random.choice(numcurves, p=wgts)

                    intensity_val = (lambda_[int_dur_curve_val] *
                                     np.exp(-0.508 * duration_val) +
                                     kappa[int_dur_curve_val] *
                                     np.exp(-0.008*duration_val) +
                                     C[int_dur_curve_val])
                    # ...these curves are based on empirical data from WG

                    if FUZZMETHOD == 'MS':
                        fuzz_int_val = np.random.choice(fuzz)
                    elif FUZZMETHOD == 'DEJH':
                        # this dist should look identical, w/o discretisation
                        fuzz_int_val = FUZZWIDTH*2.*(np.random.rand()-0.5)
                    else:
                        raise NameError
                    intensity_val += fuzz_int_val
                    # ^this allows for specified fuzzy tolerance around
                    # selected intensity (but it can go -ve)
                    # formerly, here MS used a rounding and threshold to
                    # prevent storms with a value < 1. We're going to remove
                    # the rounding and threshold at zero instead. (below)

                    # This scales the storm center intensity upward, so the
                    # values at each gauge are realistic once the gradient is
                    # applied.
                    intensity_val += intensity_val * storm_trend
                    # storminess trend is applied and its effect rises each
                    # year of simulation
                    # DEJH has removed the rounding
                    # Note that is is now possible for intensity_val to go
                    # negative, so:
                    if intensity_val < 0.:
                        intensity_val = 0.
                        self._phantom_storm_count += 1
                    # note storms of zero intensity are now permitted (though
                    # should hopefully remain pretty rare.)
                    self._intensity_val = intensity_val

                    # area to determine which gauges are hit:
                    recess_val = np.random.normal(
                        loc=Recess_pdf_norm['mu'],
                        scale=Recess_pdf_norm['sigma'])
                    try:
                        recess_val = np.clip(
                            recess_val, Recess_pdf_norm['trunc_interval'][0],
                            Recess_pdf_norm['trunc_interval'][1])
                    except KeyError:
                        pass  # this one is OK <0., I think
                    self._recess_val = recess_val
                    # this pdf of recession coefficients determines how
                    # intensity declines with distance from storm center (see
                    # below)
                    # determine cartesian distances to all hit gauges and
                    # associated intensity values at each gauge hit by the
                    # storm
                    # This is a data storage solution to avoid issues that can
                    # arise with slicing grid areas with heavy tailed sizes
                    self._entries = np.sum(mask_name)  # only open nodes
                    entries = self._entries
                    # NOTE _gauge_dist_km only contains nodes under the storm!
                    # The remaining entries are garbage
                    # Xin -> only the open nodes, note
                    self._gauge_dist_km[:entries] = np.sqrt(
                        gdist[mask_name]) / 1000.
                    self._temp_dataslots2[:entries] = gdist[mask_name] / 1.e6
                    self._temp_dataslots2[:entries] *= -2. * recess_val**2
                    np.exp(self._temp_dataslots2[:entries],
                           out=self._temp_dataslots2[:entries])
                    self._temp_dataslots2[:entries] *= intensity_val
                    mask_incl_closed = IDs_open[mask_name]
                    self._nodes_hit = mask_incl_closed
                    # ^note this is by ID, not bool
                    self._rain_int_gauge[
                        mask_incl_closed] = self._temp_dataslots2[:entries]
                    # calc of _rain_int_gauge follows Rodriguez-Iturbe et al.,
                    # 1986; Morin et al., 2005 but sampled from a distribution
                    # only need to add the bit that got rained on, so:
                    self._temp_dataslots2[:entries] *= duration_val / 60.
                    seas_cum_Ptot_gauge[mask_name] += self._temp_dataslots2[
                        :entries]
                    # collect storm totals for all gauges into rows by storm
                    Storm_total_local_seas[storm, :] = (
                        self._rain_int_gauge[opennodes] * duration_val / 60.)
                    Storm_total_local_year[(storm+storms_yr_so_far), :] = \
                        Storm_total_local_seas[storm, :]
                    self._max_storm_depth = Storm_total_local_seas[
                        storm, :].max()

                    self._Storm_total_local_seas = Storm_total_local_seas
                    self._Storm_total_local_year = Storm_total_local_year
                    Storm_running_sum_seas[
                        1, :] = Storm_total_local_seas[storm, :]
                    np.nansum(Storm_running_sum_seas, axis=0,
                              out=Storm_running_sum_seas[0, :])
                    if np.any(Storm_total_local_seas < 0.):
                        raise ValueError(syear, storm)
                    self._median_seas_rf_total = np.nanmedian(
                        Storm_running_sum_seas[0, :])
                    self._Storm_running_sum_seas = Storm_running_sum_seas[0, :]

                    if limit == 'total_time':
                        if seas_time + int_arr_val > seas_total:
                            int_arr_val = (seas_total - seas_time).clip(0.)
                            breaker = True
                    else:
                        if self._median_seas_rf_total > season_rf_limit:
                            breaker = True
                    if yield_storms is True:
                        yield (durationhrs, int_arr_val)
                    seas_time += int_arr_val
                    year_time += int_arr_val
                    if breaker:
                        # Don't create Ptotal_local per MS... just
                        breaker = False
                        break
                    if storm + 1 == self._max_numstorms:
                        raise ValueError(
                            '_max_numstorms set too low for this run')
                storms_yr_so_far = seas_storm_count
                self._storm_running_sum_of_seasons += Storm_running_sum_seas[
                    0, :]
                self._total_rainfall_last_season[
                    self._opennodes] = Storm_running_sum_seas[0, :]
                self._storm_running_sum_1st_seas += Storm_running_sum_seas[
                    0, :]
                if yield_seasons is True:
                    yield seas_storm_count

            self._total_rf_year[opennodes] = self._storm_running_sum_of_seasons
            if yield_years is True and yield_seasons is False:
                yield year_storm_count

    def calc_annual_rainfall(self, style='whole_year',
                             monsoon_total_rf_gaussian={
                                'sigma': 64., 'mu': 207.},
                             winter_total_rf_gaussian={
                                 'sigma': 52., 'mu': 1.65}):
        """
        Return a tuple of rainfall totals (mm) for the year, with entries
        subdividing the yearly total into seasons as appropriate.

        Parameters
        ----------
        style : ('whole_year', 'monsoonal', 'winter')
            Whether to simulate 2 seasons, or a single season.
        monsoon_total_rf_gaussian : dict of sigma and mu for the summer
            distribution (if used). Defaults to Walnut Gulch.
        winter_total_rf_gaussian : dict of sigma and mu for the summer
            distribution (if used). Defaults to Walnut Gulch.

        Returns
        -------
        tuple : (first_season_total, [second_season_total])
            If style=='monsoonal' or 'winter', a len(1) tuple of the total rf.
            If style=='whole_year', a len(2) tuple of (monsoon, winter) totals.

        Examples
        --------
        >>> mg = RasterModelGrid((10, 10), 500.)
        >>> z = mg.add_zeros('node', 'topographic__elevation')
        >>> rain = PrecipitationDistribution(mg)
        >>> mytotals = []
        >>> for yr in range(5):
        ...     mytotals.append(rain.calc_annual_rainfall(style='whole_year'))
        >>> [len(x) == 2 for x in mytotals]
        [True, True, True, True, True]
        >>> mytotals = []
        >>> for yr in range(3):
        ...     mytotals.append(rain.calc_annual_rainfall(style='monsoonal'))
        >>> [len(x) == 1 for x in mytotals]
        [True, True, True]
        """
        assert style in ('whole_year', 'monsoonal', 'winter')
        if style in ('whole_year', 'monsoonal'):
            # sample from normal distribution and saves global value of Ptot
            # (that must be equalled or exceeded) for each year
            summer_rf_limit = np.random.normal(
                loc=monsoon_total_rf_gaussian['mu'],
                scale=monsoon_total_rf_gaussian['sigma'])
            try:
                summer_rf_limit = np.clip(
                    summer_rf_limit,
                    monsoon_total_rf_gaussian['trunc_interval'][0],
                    monsoon_total_rf_gaussian['trunc_interval'][1])
            except KeyError:
                # ...just in case
                if summer_rf_limit < 0.:
                    summer_rf_limit = 0.
        if style in ('whole_year', 'winter'):
            # sample from normal distribution and saves global value of Ptot
            # (that must be equalled or exceeded) for each year
            winter_rf_limit = np.random.normal(
                loc=winter_total_rf_gaussian['mu'],
                scale=winter_total_rf_gaussian['sigma'])
            try:
                winter_rf_limit = np.clip(
                    winter_rf_limit,
                    winter_total_rf_gaussian['trunc_interval'][0],
                    winter_total_rf_gaussian['trunc_interval'][1])
            except KeyError:
                # ...just in case
                if winter_rf_limit < 0.:
                    winter_rf_limit = 0.
        if style == 'monsoonal':
            return (summer_rf_limit, )
        elif style == 'winter':
            return (winter_rf_limit, )
        else:
            return (summer_rf_limit, winter_rf_limit)

    def _locate_storm(self, storm_radius):
        """
        Because of the way the stats fall out, any simulated storm from the
        distribution must intersect the catchment somewhere. Trying to write
        this in a grid-agnostic fashion...
        """
        stormposx = np.random.rand()*(self._widthx + 2.*storm_radius)
        stormposy = np.random.rand()*(self._widthy + 2.*storm_radius)
        stormx = self._minx - storm_radius + stormposx
        stormy = self._miny - storm_radius + stormposy
        return stormx, stormy

    @property
    def current_year(self):
        """Get the current year as an int."""
        return self._year

    @property
    def current_season(self):
        """Get the current season. 'M' is monsoon, 'W' is winter."""
        return self._current_season

    @property
    def storm_depth_last_storm(self):
        """Get the maximum storm depth during the last storm (mm)."""
        return self._max_storm_depth

    @property
    def storm_recession_value_last_storm(self):
        """Get the recession parameter (radial die-off) for the last storm."""
        return self._recess_val

    @property
    def storm_duration_last_storm(self):
        """Get the duration (in hrs) of the last storm."""
        return self._durationhrs

    @property
    def storm_area_last_storm(self):
        """Get the area (in m**2) of the last storm."""
        return self._area_val

    @property
    def storm_intensity_last_storm(self):
        """
        Get the intensity (mm/hr) of the last storm, averaged under the storm
        footprint. Note that duration * intensity != storm max depth.
        """
        return self._intensity_val

    @property
    def total_rainfall_last_season(self):
        """
        Get the total recorded rainfall over the last (completed) simulated
        season, spatially resolved (mm).
        """
        return self._total_rainfall_last_season

    @property
    def total_rainfall_last_year(self):
        """
        Get the total recorded rainfall over the last (completed) simulated
        year, spatially resolved (mm). Equivalent to the field
        'rainfall__total_depth_per_year'.
        """
        return self._total_rf_year

    @property
    def total_rainfall_this_season(self):
        """
        Get the accumulated, spatially resolved total rainfall over the
        grid for the season so far (mm).
        """
        self._running_total_rainfall_this_season[
            self._opennodes] = self._Storm_running_sum_seas
        return self._running_total_rainfall_this_season

    @property
    def total_rainfall_this_year(self):
        """
        Get the accumulated, spatially resolved total rainfall over the
        grid for the year so far (mm).
        """
        self._running_total_rainfall_this_year[
            self._opennodes] = (self._storm_running_sum_1st_seas +
                                self._Storm_running_sum_seas)
        return self._running_total_rainfall_this_year

    @property
    def median_total_rainfall_last_season(self):
        """
        Get the median total rainfall recorded over the open nodes of the grid
        during the last (completed) simulated season (mm).
        """
        return np.nanmedian(self._total_rainfall_last_season[self._opennodes])

    @property
    def median_total_rainfall_last_year(self):
        """
        Get the median total rainfall recorded over the open nodes of the grid
        during the last (completed) simulated year (mm).
        """
        return np.nanmedian(self.total_rainfall_last_year[self._opennodes])

    @property
    def median_total_rainfall_this_season(self):
        """
        Get the accumulated median total rainfall over the open nodes of the
        grid so far this season (mm).
        """
        return self._median_seas_rf_total

    @property
    def median_total_rainfall_this_year(self):
        """
        Get the accumulated median total rainfall over the open nodes of the
        grid so far this year (mm).
        """
        return np.nanmedian(self.total_rainfall_this_year[self._opennodes])

    @property
    def number_of_nodes_under_storm(self):
        """Get the number of nodes under the last storm."""
        return self._entries

    @property
    def nodes_under_storm(self):
        """Get the IDs of the nodes under the last storm."""
        return self._nodes_hit

    @property
    def coordinates_of_last_storm_center(self):
        """
        Get the coordinates of the center of the last storm as (x, y).
        """
        return (self._x, self._y)


if __name__ == "__main__":
    from landlab.plot import imshow_grid_at_node
    from matplotlib.pyplot import show

    nx = 40
    ny = 40
    dx = 250.
    mg = RasterModelGrid((nx, ny), dx)
    # closed_nodes = np.zeros((nx, ny), dtype=bool)
    # closed_nodes[:, :10] = True
    # closed_nodes[:, 30:] = True
    # closed_nodes[30:, :] = True
    # mg.status_at_node[closed_nodes.flatten()] = CLOSED_BOUNDARY
    # imshow_grid_at_node(mg, mg.status_at_node)
    # show()

    z = mg.add_zeros('node', 'topographic__elevation')
    z += 1400.
    rain = PrecipitationDistribution(mg, number_of_years=2)
    count = 0
    total_t = 0.
    for dt, interval_t in rain.yield_storms(style='whole_year',
                                            limit='total_time'):
        count += 1
        total_t += dt + interval_t
        # print(rain.median_total_rainfall_this_year)
        if count % 50 == 0:
            # imshow_grid_at_node(mg, rain.total_rainfall_this_year,
            #                     cmap='Blues')
            print('Season:', rain.current_season, 'of yr', rain.current_year)
            print('Current storm:', count)
            print('MEANS')
            print(rain.total_rainfall_this_season[rain._opennodes].mean())
            print(rain.total_rainfall_last_season[rain._opennodes].mean())
            print(rain.total_rainfall_this_year[rain._opennodes].mean())
            print(rain.total_rainfall_last_year[rain._opennodes].mean())
            print('-----')
            print('MEDIANS')
            print(rain.median_total_rainfall_this_season)
            print(rain.median_total_rainfall_last_season)
            print(rain.median_total_rainfall_this_year)
            print(rain.median_total_rainfall_last_year)
            print('*****')
            show()
    print("Effective total years:")
    print(total_t/24./365.)
    print("Storms simulated:")
    print(count)

    # mg = RasterModelGrid((100, 100), 500.)
    # # mg.status_at_node[closed_nodes.flatten()] = CLOSED_BOUNDARY
    # # imshow_grid_at_node(mg, mg.status_at_node)
    # # show()
    # z = mg.add_zeros('node', 'topographic__elevation')
    # z += 1000.
    # rain = PrecipitationDistribution(mg, number_of_years=2)
    # count = 0
    # total_storms = 0.
    # for storms_in_year in rain.yield_years():
    #     count += 1
    #     total_storms += storms_in_year
    #     print(storms_in_year)
    #     imshow_grid_at_node(mg, 'rainfall__total_depth_per_year', cmap='jet')
    #     show()
    # 
    # _ = mg.at_node.pop('rainfall__total_depth_per_year')
    # _ = mg.at_node.pop('rainfall__flux')
    # rain = PrecipitationDistribution(mg, number_of_years=1)
    # for storms_in_season in rain.yield_seasons():
    #     print(storms_in_season)
    #     imshow_grid_at_node(mg, rain.total_rainfall_last_season, cmap='jet')
    #     show()
    # 
    # for yr in range(30):
    #     print(rain.calc_annual_rainfall(style='whole_year'))
    # 
    # from landlab import VoronoiDelaunayGrid
    # 
    # x = np.random.rand(2000)*50000.
    # y = np.random.rand(2000)*50000.
    # vdg = VoronoiDelaunayGrid(x, y)
    # vdg.add_zeros('node', 'topographic__elevation')
    # rain = PrecipitationDistribution(vdg, number_of_years=2)
    # count = 0
    # total_storms = 0.
    # for storms_in_year in rain.yield_years(limit='total_rainfall'):
    #     count += 1
    #     total_storms += storms_in_year
    #     print(storms_in_year)
    #     imshow_grid_at_node(vdg, 'rainfall__total_depth_per_year',
    #     cmap='jet')
    #     print('Stochastically-set expected rainfall last season')
    #     print()
    #     print('Actual rainfall last season')
    #     print()
    #     show()
