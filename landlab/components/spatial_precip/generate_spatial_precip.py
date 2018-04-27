# %----------------------------------------------------------------------------
# %STOchastic Rainfall Model (STORM). A rainstorm generator in this case based
# on empirical data (Walnut Gulch, AZ) # %Code name: WG_storms_v3_01  %this
# version allows for simulations that are much longer than the input series in
# order to compare the distributions of storm characteristics to the original.
# %It also includes the 'mean gauge approach' to determine when a simulation
# year stops. This involves summing storm totals at each gauge for each year
# until the mean for all gauges exceeds the selected annual total precip value.
# It also allows for fuzzy selection of intensity at the storm center based on
# a fixed value of duration and %incorporates orographic effects, wherein there
# are separate intensity-duration curves derived for three intervals of basin
# elevation (1200-1350m, 1351-1500m, 1501-1650m) %Current version also includes
# interarrival times between storms, allowing for output to drive other model
# frameworks (rainfall-ruonff, water balance,LEMs) %This version will also
# include output at each grid location, rather than only at gauge locations.
# %Author: Michael Singer 2017 %Date created: 2015-6
# ----------------------------------------------------------------------------

# Major diffs to Singer code:
# No ET modelled
# No such thing as validation vs simulation
# It produces fuzz by a continuous distribution
# Step changes are not explicitly included
# Storms can be centred at any point, not just nodes
# Edge buffering is now dynamic; i.e., big storms have a bigger edge buffer
# Storms are never discarded - once a storm is drawn, it must hit the catchment

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

    def __init__(self, grid,
                 number_of_years=1, buffer_width=5000,
                 orographic_scenario='Singer', save_outputs=False,
                 path_to_input_files=None):
        """
        It's on the user to ensure the grid is big enough to permit the buffer.
        save_outputs : if not None, str path to save
        orographic_scenario : {None, 'Singer'}

        The Storm_matrix is:
        0  : master storm count
        1  : storm area (m**2)
        2  : storm duration (rounded to nearest min, in min)
        3  : which curve is selected to draw area/duration from? (0-10)
        4  : peak storm intensity (mm/hr)
        5  : number of gauges/nodes hit
        6  : recession value (how rapidly does intensity wane from centre?)
        7  : accumulated rf during storm (mm)
        8  : storm centre x coordinate (grid unit)
        9  : storm centre y coordinate (grid unit)
        10 : year in the simulation
        11 : the season. 0 is monsoon, 1 is winter
        """
        self._grid = grid
        gaugecount = (grid.status_at_node != CLOSED_BOUNDARY).sum()
        self._gauge_dist_km = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots1 = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots2 = np.zeros(gaugecount, dtype='float')
        self._numyrs = number_of_years
        self._buffer_width = buffer_width
        if save_outputs is not None:
            assert type(save_outputs) in (bool, str)
        self._savedir = save_outputs

        self._max_numstorms = 5000
        # This is for initializing matrices. Trailing zeros are deleted from
        # the matrix at the end of the code.

        if self._savedir is not False:
            if type(self._savedir) is str:
                dimension = self._max_numstorms
            else:
                dimension = self._max_numstorms*number_of_years
            self._Storm_matrix = np.zeros((dimension, 12))

        assert orographic_scenario in {None, 'Singer'}
        self._orographic_scenario = orographic_scenario

        if path_to_input_files is None:
            self._path = os.path.join(
                os.path.dirname(os.path.realpath(__file__)), 'model_input')
        else:
            assert type(path_to_input_files) is str
            self._path = path_to_input_files

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

    def yield_storms(self, style='whole_year',
                     total_rf_trend=0., storminess_trend=0.,
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
        a rainfall simulation. As

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, 2018 & Singer et al, submitted.

        Parameters
        ----------
        style : ('whole_year', 'monsoonal', 'winter')
            Controls whether the component seeks to simulate a western US-
            style "monsoonal" climate, a western US-style winter climate,
            or a full year combining both. These distributions are by default
            based on Singer et al.'s calibrations.
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

       monsoon_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
        winter_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        Yields
        ------
        (storm_t, interval_t) : (float, float)
            Tuple pair of duration of a single storm, then the interstorm
            interval that follows it. The rainfall__flux field describes the
            rainfall rate during the interval storm_t as the tuple is yielded.
            Note that the rainfall__total_depth_per_year field gives the total
            accumulated rainfall depth during the *last completed* model year,
            not the year to the point of yield.
        """
        return self._run_the_process(
            yield_storms=True, yield_years=False, yield_seasons=False,
            style=style,
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

    def yield_years(self, style='whole_year',
                    total_rf_trend=0., storminess_trend=0.,
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
        style : ('whole_year', 'monsoonal', 'winter')
            Controls whether the component seeks to simulate a western US-
            style "monsoonal" climate, a western US-style winter climate,
            or a full year combining both. These distributions are by default
            based on Singer et al.'s calibrations.
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

        monsoon_total_rf_gaussian is a normal distribution controlling the
            total rainfall expected in each year. S&M use 'mu' in {143., 271.}
            for step changes up/down in rainfall totals. In mm.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
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
            style=style,
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

    def yield_seasons(self, style='whole_year',
                      total_rf_trend=0., storminess_trend=0.,
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
        style : ('whole_year', 'monsoonal', 'winter')
            Controls whether the component seeks to simulate a western US-
            style "monsoonal" climate, a western US-style winter climate,
            or a full year combining both. These distributions are by default
            based on Singer et al.'s calibrations.
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

        monsoon_total_rf_gaussian is a normal distribution controlling the
            total rainfall expected in each year. S&M use 'mu' in {143., 271.}
            for step changes up/down in rainfall totals. In mm.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
        winter_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.

        Yields
        ------
        number_of_storms_per_season : float
            Float that gives the number of storms simulated in the season that
            elapsed since the last yield. The rainfall__total_depth_per_year
            field gives the total accumulated rainfall depth during the year
            preceding the yield, *so far*. rainfall__flux gives the rainfall
            intensity of the last storm in that year.
        NB: Use the component property total_rainfall_last_season to access
            the *actual* amount of rainfall in the season that has the number
            of storms that the method generates.
        """
        return self._run_the_process(
            yield_storms=False, yield_years=False, yield_seasons=True,
            style=style,
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

    def _run_the_process(self, yield_storms=True, yield_years=False,
                         yield_seasons=False, style='whole_year',
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

        Use the hardwired FUZZMETHOD switch {'MS', 'DEJH'} to control whether
        the fuzz is applied in a discretised fashion from input file (MS), or
        as a continuous random variable (DEJH).

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

        NOTE still an issue with how we do intensity stepchange.

        All default distributions reflect values for Walnut Gulch, see Singer &
        Michaelides, submitted:

        monsoon_total_rf_gaussian is a normal distribution controlling the
            total rainfall expected in each year. S&M use 'mu' in {143., 271.}
            for step changes up/down in rainfall totals. In mm.
        monsoon_storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm.
        monsoon_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
        monsoon_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.
        winter_total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk is a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV is a generalised extreme value
            distribution controlling the interarrival time between each storm.
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
        # what's the dir of this component?
        # this is a nasty hacky way for now
        thisdir = self._path
        # Initialize variables for annual rainfall total (mm/h)
        # storm center location (RG1-RG85), etc.

        # add variable for number of simulations of simyears
        simyears = self._numyrs  # number of years to simulate
        numcurves = 11  # number of intensity-duration curves (see below for
        # curve equations)

        assert style in ('whole_year', 'monsoonal', 'winter')
        if style == 'whole_year':
            reps = 2
        else:
            reps = 1

        opennodes = self.grid.status_at_node != CLOSED_BOUNDARY
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
            # a vector of fuzzy tolerance values for intensity selection
            # width is +/-5, discretised every 1
            fuzz = np.loadtxt(os.path.join(thisdir, 'fuzz.csv'))
            fuzz = fuzz.astype(float)

        # now build a buffered target area of nodes:
        target_area_nodes = self.grid.zeros('node', dtype=bool)
        # which are within buffer_width of the perimeter? Try to do this
        # in a memory efficient fashion.
        # True catchment edges must have link statuses that are CLOSED:
        closed_links = self.grid.status_at_link == CLOSED_BOUNDARY
        # one of their end nodes must be not CLOSED:
        edge_link_head_open = self.grid.status_at_node[
            self.grid.node_at_link_head][closed_links] != CLOSED_BOUNDARY
        head_open_node_IDs = self.grid.node_at_link_head[closed_links][
            edge_link_head_open]
        tail_open_node_IDs = self.grid.node_at_link_tail[closed_links][
            np.logical_not(edge_link_head_open)]
        # Together, this is a list of the IDs of all the nodes on the catchment
        # perimeter. So:
        for node_list in (head_open_node_IDs, tail_open_node_IDs):
            for edgenode in node_list:
                edgenode_x = self.grid.x_of_node[edgenode]
                edgenode_y = self.grid.y_of_node[edgenode]
                dists_to_edgenode = self.grid.calc_distances_of_nodes_to_point(
                    (edgenode_x, edgenode_y))
                target_area_nodes[
                    dists_to_edgenode <= self._buffer_width] = True
        # finish off by stamping the core nodes over the top:
        target_area_nodes[opennodes] = True

        Xxin = self.grid.x_of_node[target_area_nodes]
        Yyin = self.grid.y_of_node[target_area_nodes]

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
            if type(self._savedir) is str:
                self._Storm_matrix.fill(0.)
            calendar_time = 0  # tracks simulation time per year in hours
            storm_trend += storminess_trend
            Ptotal = 0.
            year_storm_count = 0
            Storm_total_local_year = np.zeros(
                (self._max_numstorms, num_opennodes))
            Storm_running_sum_year = np.zeros(num_opennodes)

            storms_yr_so_far = 0
            for seas in range(reps):
                Storm_running_sum_seas = np.zeros((2, num_opennodes))
                # ^ 1st col is running total, 2nd is data to add to it
                if seas == 0 and not style == 'winter':
                    print('MONSOON')
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
                else:
                    print('WINTER')
                    Ptot_pdf_norm = winter_total_rf_gaussian
                    Duration_pdf = winter_storm_duration_fisk
                    Area_pdf_EV = winter_storm_area_GEV
                    Int_arr_pdf_GEV = winter_storm_interarrival_GEV
                    Recess_pdf_norm = winter_storm_radial_weakening_gaussian

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
                    # ^Samples from distribution of storm areas
                    r = np.sqrt(area_val/np.pi)  # value here shd be selected
                    # based on area above in meters to match the UTM values
                    while 1:
                        cx, cy = self._locate_storm(r)
                        # Determine which gauges are hit by Euclidean geometry:
                        gdist = (Xin-cx)**2 + (Yin-cy)**2
                        mask_name = (gdist <= r**2)  # this is defacto MS's aa
                        # this short circuits the storm loop in the case that
                        # the storm does not affect any 'gauging' location
                        if np.any(np.equal(mask_name, True)):
                            break
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
                    # round to nearest minute for consistency w measured data:
                    duration_val = round(duration_val)
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
                    # we're not going to store the calendar time (DEJH change)

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
                    # this pdf of recession coefficients determines how
                    # intensity declines with distance from storm center (see
                    # below)
                    # determine cartesian distances to all hit gauges and
                    # associated intensity values at each gauge hit by the
                    # storm
                    # This is a data storage solution to avoid issues that can
                    # arise with slicing grid areas with heavy tailed sizes
                    xcoords = Xin
                    ycoords = Yin
                    self._entries = np.sum(mask_name)
                    entries = self._entries
                    # NOTE _gauge_dist_km only contains nodes under the storm!
                    # The remaining entries are garbage
                    self._gauge_dist_km[:entries] = xcoords[mask_name]
                    self._gauge_dist_km[:entries] -= cx
                    np.square(self._gauge_dist_km[:entries],
                              out=self._gauge_dist_km[:entries])
                    self._temp_dataslots1[:entries] = ycoords[mask_name]
                    self._temp_dataslots1[:entries] -= cy
                    np.square(self._temp_dataslots1[:entries],
                              out=self._temp_dataslots1[:entries])
                    self._gauge_dist_km[:entries] += self._temp_dataslots1[
                        :entries]
                    if np.any(self._gauge_dist_km[:entries] < 0.):
                        raise ValueError()
                    np.sqrt(self._gauge_dist_km[:entries],
                            out=self._gauge_dist_km[:entries])
                    self._gauge_dist_km[:entries] /= 1000.
                    # _rain_int_gauge has been zeroed earlier in loop, so
                    self._temp_dataslots2[:entries] = self._gauge_dist_km[
                        :entries]
                    # this copy would not be necessary if we didn't want to
                    # preserve self._gauge_dist_km
                    np.square(self._temp_dataslots2[:entries],
                              out=self._temp_dataslots2[:entries])
                    self._temp_dataslots2[:entries] *= -2. * recess_val**2
                    np.exp(self._temp_dataslots2[:entries],
                           out=self._temp_dataslots2[:entries])
                    self._temp_dataslots2[:entries] *= intensity_val
                    mask_incl_closed = IDs_open[mask_name]
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

                    self._Storm_total_local_seas = Storm_total_local_seas
                    self._Storm_total_local_year = Storm_total_local_year
                    Storm_running_sum_seas[
                        1, :] = Storm_total_local_seas[storm, :]
                    self._Storm_running_sum_seas = Storm_running_sum_seas
                    np.nansum(Storm_running_sum_seas, axis=0,
                              out=Storm_running_sum_seas[0, :])
                    if np.any(Storm_total_local_seas < 0.):
                        raise ValueError(syear, storm)
                    self._median_seas_rf_total = np.nanmedian(
                        Storm_running_sum_seas[0, :])

                    if self._savedir is not False:
                        if type(self._savedir) is str:
                            rowID = year_storm_count - 1
                        else:
                            rowID = master_storm_count - 1
                        # save some properties:
                        Storm_matrix = self._Storm_matrix
                        Storm_matrix[rowID, 0] = master_storm_count
                        Storm_matrix[rowID, 1] = area_val
                        Storm_matrix[rowID, 2] = duration_val
                        Storm_matrix[rowID, 3] = int_dur_curve_val
                        Storm_matrix[rowID, 4] = intensity_val
                        Storm_matrix[rowID, 5] = num_gauges_hit
                        Storm_matrix[rowID, 6] = recess_val
                        Storm_matrix[rowID, 7] = (intensity_val *
                                                  duration_val / 60.)
                        Storm_matrix[rowID, 8] = cx
                        Storm_matrix[rowID, 9] = cy
                        Storm_matrix[rowID, 10] = syear
                        if seas == 0 and not style == 'winter':
                            Storm_matrix[rowID, 11] = 0
                        else:
                            Storm_matrix[rowID, 11] = 1

                    if yield_storms is True:
                        yield (duration_val, int_arr_val)
                    # now blank the field for the interstorm period
                    self._rain_int_gauge.fill(0.)
                    if self._median_seas_rf_total > season_rf_limit:
                        # Don't create Ptotal_local per MS... just
                        break
                    if storm + 1 == self._max_numstorms:
                        raise ValueError(
                            '_max_numstorms set too low for this run')
                storms_yr_so_far = seas_storm_count
                Storm_running_sum_year += Storm_running_sum_seas[0, :]
                self._total_rainfall_last_season = Storm_running_sum_seas[0, :]
                if yield_seasons is True:
                    yield seas_storm_count

            if type(self._savedir) is str:
                # crop it down to save:
                [maxrow, maxcol] = np.argwhere(self._Storm_matrix).max(axis=0)
                np.savetxt(os.path.join(
                    self._savedir, 'Storm_matrix_y' + str(syear) + '.csv'),
                    self._Storm_matrix[:(maxrow+1), :(maxcol+1)])
            self._total_rf_year[opennodes] = Storm_running_sum_year
            if yield_years is True and yield_seasons is False:
                yield year_storm_count

        if self._savedir is True:
            # crop it down for convenience:
            [maxrow, maxcol] = np.argwhere(self._Storm_matrix).max(axis=0)
            self._Storm_matrix = self._Storm_matrix[:(maxrow+1), :(maxcol+1)]

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
    def total_rainfall_last_season(self):
        """Get the total recorded rainfall over the last simulated season."""
        return self._total_rainfall_last_season


from landlab.plot import imshow_grid_at_node
from matplotlib.pyplot import show
mg = RasterModelGrid((100, 100), 500.)
closed_nodes = np.zeros((100, 100), dtype=bool)
closed_nodes[:, :30] = True
closed_nodes[:, 70:] = True
closed_nodes[70:, :] = True

mg.status_at_node[closed_nodes.flatten()] = CLOSED_BOUNDARY
# imshow_grid_at_node(mg, mg.status_at_node)
# show()
z = mg.add_zeros('node', 'topographic__elevation')
z += 1000.
rain = PrecipitationDistribution(mg, number_of_years=2, save_outputs=True)
count = 0
total_t = 0.
for dt, interval_t in rain.yield_storms():
    count += 1
    total_t += dt + interval_t
    print rain._median_seas_rf_total
    if count % 100 == 0:
        imshow_grid_at_node(mg, 'rainfall__flux', cmap='Blues')
        show()
print("Effective total years:")
print(total_t/24./365.)
#
# print('*****')
#
# mg = RasterModelGrid((100, 100), 500.)
# mg.status_at_node[closed_nodes.flatten()] = CLOSED_BOUNDARY
# # imshow_grid_at_node(mg, mg.status_at_node)
# # show()
# z = mg.add_zeros('node', 'topographic__elevation')
# z += 1000.
# rain = PrecipitationDistribution(mg, number_of_years=3)
# count = 0
# total_t = 0.
# for dt, interval_t in rain.yield_storms():
#     count += 1
#     total_t += dt + interval_t
#     print rain._median_seas_rf_total
#     if count % 100 == 0:
#         imshow_grid_at_node(mg, 'rainfall__flux')
#         show()
# print("Effective total years:")
# print(total_t/24./365.)
# print('*****')
mg = RasterModelGrid((100, 100), 500.)
# mg.status_at_node[closed_nodes.flatten()] = CLOSED_BOUNDARY
# imshow_grid_at_node(mg, mg.status_at_node)
# show()
z = mg.add_zeros('node', 'topographic__elevation')
z += 1000.
rain = PrecipitationDistribution(mg, number_of_years=3)
count = 0
total_storms = 0.
for storms_in_year in rain.yield_years():
    count += 1
    total_storms += storms_in_year
    print(storms_in_year)
    imshow_grid_at_node(mg, 'rainfall__total_depth_per_year', cmap='jet')
    show()

_ = mg.at_node.pop('rainfall__total_depth_per_year')
_ = mg.at_node.pop('rainfall__flux')
rain = PrecipitationDistribution(mg, number_of_years=2)
for storms_in_season in rain.yield_seasons():
    print(storms_in_season)
    imshow_grid_at_node(mg, rain.total_rainfall_last_season, cmap='jet')
    show()

for yr in range(30):
    print(rain.calc_annual_rainfall(style='whole_year'))
