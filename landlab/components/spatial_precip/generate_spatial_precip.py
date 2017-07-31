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

import numpy as np
import os
import inspect
from six.moves import range
from matplotlib.pyplot import figure
from scipy.stats import genextreme
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

    def __init__(self, grid, mode='simulation', number_of_simulations=1,
                 number_of_years=1, buffer_width=5000,
                 orographic_scenario='Singer', save_outputs=False,
                 path_to_input_files='/Users/daniel/development/landlab/landlab/components/spatial_precip'):
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
        """
        self._grid = grid
        assert mode in ('simulation', 'validation')
        self._mode = mode
        if mode == 'simulation':
            gaugecount = (grid.status_at_node != CLOSED_BOUNDARY).sum()
        else:
            Eastings = np.loadtxt(os.path.join(thisdir, 'Easting.csv'))
            gaugecount = Eastings.size
        self._gauge_dist_km = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots1 = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots2 = np.zeros(gaugecount, dtype='float')
        self._numsims = number_of_simulations
        self._numyrs = number_of_years
        self._buffer_width = buffer_width
        if save_outputs is not None:
            assert type(save_outputs) in (bool, str)
        self._savedir = save_outputs

        self._max_numstorms = 2000
        # This is for initializing matrices. Trailing zeros are deleted from
        # the matrix at the end of the code.

        if self._savedir is not False:
            if type(self._savedir) is str:
                dimension = self._max_numstorms
            else:
                dimension = self._max_numstorms*number_of_years
            self._Storm_matrix = np.zeros((dimension, 11))

        assert orographic_scenario in {None, 'Singer'}
        self._orographic_scenario = orographic_scenario

        self._path = path_to_input_files

        # build LL fields:
        self.initialize_output_fields()
        # bind the field to the internal variable:
        self._rain_int_gauge = self.grid.at_node['rainfall__flux']
        self._total_rf_year = self.grid.at_node[
            'rainfall__total_depth_per_year']

    def yield_storms(self, total_rf_trend=0., storminess_trend=0.,
                     total_rf_gaussian={
                         'sigma': 64., 'mu': 207.},
                     storm_duration_GEV={
                         'shape': -0.570252, 'sigma': 35.7389, 'mu': 34.1409,
                         'trunc_interval': (1., 1040.)},
                     storm_area_GEV={
                         'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                         'trunc_interval': (5.e+06, 3.e+08)},
                     storm_interarrival_GEV={
                         'shape': -0.807971, 'sigma': 9.49574, 'mu': 10.6108,
                         'trunc_interval': (0., 120.)},
                     storm_radial_weakening_gaussian={
                         'sigma': 0.08, 'mu': 0.25,
                         'trunc_interval': (0.15, 0.67)}):
        """
        Yield a timeseries giving the number if storms occurring each year in
        a rainfall simulation. As

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, submitted.

        Parameters
        ----------
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
        NOTE still an issue with how we do intensity stepchange.

        total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm.
        storm_area_GEV is a generalised extreme value distribution controlling
            the plan view area of each storm. S&M use 'shape': 0., which
            collapses the distribution to a plain extreme value distribution.
        storm_interarrival_GEV is a generalised extreme value distribution
            controlling the interarrival time between each storm.
        storm_radial_weakening_gaussian is a normal distribution controlling
            the rate of intensity decline with distance from storm center. For
            more detail see Rodriguez-Iturbe et al., 1986; Morin et al., 2005.

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
            yield_storms=True, yield_years=False,
            total_rf_trend=total_rf_trend, storminess_trend=storminess_trend,
            total_rf_gaussian=total_rf_gaussian,
            storm_duration_GEV=storm_duration_GEV,
            storm_area_GEV=storm_area_GEV,
            storm_interarrival_GEV=storm_interarrival_GEV,
            storm_radial_weakening_gaussian=storm_radial_weakening_gaussian)

    def yield_years(self, total_rf_trend=0., storminess_trend=0.,
                    total_rf_gaussian={
                        'sigma': 64., 'mu': 207.},
                    storm_duration_GEV={
                        'shape': -0.570252, 'sigma': 35.7389, 'mu': 34.1409,
                        'trunc_interval': (1., 1040.)},
                    storm_area_GEV={
                        'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                        'trunc_interval': (5.e+06, 3.e+08)},
                    storm_interarrival_GEV={
                        'shape': -0.807971, 'sigma': 9.49574, 'mu': 10.6108,
                        'trunc_interval': (0., 120.)},
                    storm_radial_weakening_gaussian={
                        'sigma': 0.08, 'mu': 0.25,
                        'trunc_interval': (0.15, 0.67)}):
        """
        Yield a timeseries giving the number if storms occurring each year in
        a rainfall simulation. As

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, submitted.

        Parameters
        ----------
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
        NOTE still an issue with how we do intensity stepchange.

        total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm.
        storm_area_GEV is a generalised extreme value distribution controlling
            the plan view area of each storm. S&M use 'shape': 0., which
            collapses the distribution to a plain extreme value distribution.
        storm_interarrival_GEV is a generalised extreme value distribution
            controlling the interarrival time between each storm.
        storm_radial_weakening_gaussian is a normal distribution controlling
            the rate of intensity decline with distance from storm center. For
            more detail see Rodriguez-Iturbe et al., 1986; Morin et al., 2005.

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
            yield_storms=False, yield_years=True,
            total_rf_trend=total_rf_trend, storminess_trend=storminess_trend,
            total_rf_gaussian=total_rf_gaussian,
            storm_duration_GEV=storm_duration_GEV,
            storm_area_GEV=storm_area_GEV,
            storm_interarrival_GEV=storm_interarrival_GEV,
            storm_radial_weakening_gaussian=storm_radial_weakening_gaussian)

    def _run_the_process(self, yield_storms=True, yield_years=False,
                         total_rf_trend=0., storminess_trend=0.,
                         total_rf_gaussian={
                             'sigma': 64., 'mu': 207.},
                         storm_duration_GEV={
                             'shape': -0.570252, 'sigma': 35.7389,
                             'mu': 34.1409, 'trunc_interval': (1., 1040.)},
                         storm_area_GEV={
                             'shape': 0., 'sigma': 2.83876e+07,
                             'mu': 1.22419e+08,
                             'trunc_interval': (5.e+06, 3.e+08)},
                         storm_interarrival_GEV={
                             'shape': -0.807971, 'sigma': 9.49574,
                             'mu': 10.6108, 'trunc_interval': (0., 120.)},
                         storm_radial_weakening_gaussian={
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

        total_rf_gaussian is a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        storm_duration_GEV is a generalised extreme value distribution
            controlling the duration of each storm.
        storm_area_GEV is a generalised extreme value distribution controlling
            the plan view area of each storm. S&M use 'shape': 0., which
            collapses the distribution to a plain extreme value distribution.
        storm_interarrival_GEV is a generalised extreme value distribution
            controlling the interarrival time between each storm.
        storm_radial_weakening_gaussian is a normal distribution controlling
            the rate of intensity decline with distance from storm center. For
            more detail see Rodriguez-Iturbe et al., 1986; Morin et al., 2005.

        """
        FUZZMETHOD = 'DEJH'
        FUZZWIDTH = 5.  # if DEJH
        self._phantom_storm_count = 0
        # ^this property tracks the number of storms in the run that received
        # zero intensity (and thus didn't really exist)
        # safety check for init conds:
        if yield_storms:
            assert yield_years is False
        if yield_years:
            assert yield_storms is False
        # what's the dir of this component?
        # this is a nasty hacky way for now
        thisdir = self._path
        # Initialize variables for annual rainfall total (mm/h)
        # storm center location (RG1-RG85), etc.

# NOTE this one is currently redundant:
        # This scalar specifies the value of fractional step change in
        # intensity when storms+ or storms- are applied in STORMINESS_SCENARIO
        storm_stepchange = 0.25

        # add variable for number of simulations of simyears
        simyears = self._numyrs  # number of years to simulate
        numcurves = 11  # number of intensity-duration curves (see below for
        # curve equations)

        storm_scaling = 1.  # No storm scaling, as problem appears to be fixed
        # with smaller grid spacing.
        # Formerly, this scaled the storm center intensity upward, so the
        # values at each gauge are realistic once the gradient is applied.

        # This is the pdf fitted to all available station precip data (normal
        # dist). It will be sampled below.
# NOTE right now we ignore all poss scenarios, i.e., use the C cases (? Check)
        Ptot_pdf_norm = total_rf_gaussian
        # This step change case can now be handled direct from input flags
        # if self._ptot_scenario == 'ptot+':
        #     Ptot_pdf_norm = {'sigma': 64., 'mu': 271.}
        # elif self._ptot_scenario == 'ptot-':
        #     Ptot_pdf_norm = {'sigma': 64., 'mu': 143.}
        # else:
        #     Ptot_pdf_norm = {'sigma': 64., 'mu': 207.}
        # the trending cases need to be handled in the loop

        # This is the pdf fitted to all available station duration data
        # (GEV dist). It will be sampled below.
        # #### matlab's GEV is (shape_param, scale(sigma), pos(mu))
        # note that in Scipy, we must add a minus to the shape param for a GEV
        # to match Matlab's implementation
        Duration_pdf_GEV = storm_duration_GEV
        # This is the pdf fitted to all available station area data (EV dist).
        # It will be sampled below.
        # #### matlab's EV is (mu, sigma)
        Area_pdf_EV = storm_area_GEV
        # This is the pdf fitted to all available station area data (GEV dist).
        # It will be sampled below.
        Int_arr_pdf_GEV = storm_interarrival_GEV
        # This is the pdf of storm gradient recession coefficiencts from Morin
        # et al, 2005 (normal dist). It will be sampled below.
        Recess_pdf_norm = storm_radial_weakening_gaussian

        opennodes = self.grid.status_at_node != CLOSED_BOUNDARY
        num_opennodes = np.sum(opennodes)
        IDs_open = np.where(opennodes)[0]  # need this later
        if self._mode == 'validation':
            raise ValueError('Not currently supported!')
            Easting = np.loadtxt(os.path.join(thisdir, 'Easting.csv'))
            # This is the Longitudinal data for each gauge.
            Northing = np.loadtxt(os.path.join(thisdir, 'Northing.csv'))
            # This is the Latitudinal data for each gauge.
            gauges = np.loadtxt(os.path.join(thisdir, 'gauges.csv'))
            # This is the list of gauge numbers. It will be sampled below.
            gauge_elev = np.loadtxt(os.path.join(thisdir, 'gauge_elev.csv'))
            numgauges = gauges.size
        else:
            X1 = self.grid.node_x
            Y1 = self.grid.node_y
            Xin = X1[opennodes]
            Yin = Y1[opennodes]
            Zz = self.grid.at_node['topographic__elevation'][opennodes]
            numgauges = Xin.size  # number of rain gauges in the basin.
        # NOTE: In this version this produces output on a grid, rather than at
        # real gauge locations.

# NOTE this block is for validation & shouldn't live here
        # # This is the storm depth data for use in model evaluation.
        # Storm_depth_data = np.loadtxt(os.path.join(thisdir,
        #                                            'Storm_depth_data.csv'))
        # # This is the intensity data for use in model evaluation.
        # Intensity_data = np.loadtxt(os.path.join(thisdir,
        #                                          'Intensity_data.csv'))
        # # This is the duration data for use in model evaluation.
        # Duration_data = np.loadtxt(os.path.join(thisdir,
        #                                         'Duration_data.csv'))
        if FUZZMETHOD == 'MS':
            # a vector of fuzzy tolerance values for intensity selection
            # width is +/-5, discretised every 1
            fuzz = np.loadtxt(os.path.join(thisdir, 'fuzz.csv'))
            fuzz = fuzz.astype(float)
        ET_monthly_day = np.loadtxt(os.path.join(thisdir,
                                                 'ET_monthly_day.txt'))
        ET_monthly_night = np.loadtxt(os.path.join(thisdir,
                                                   'ET_monthly_night.txt'))
        # This are matrices of averaged day/nighttime values of ET grouped as
        # one column per month.

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
            Ptotal = 0
            if not np.isclose(total_rf_trend, 0.):
                mu = Ptot_pdf_norm.pop('mu')
                mu += mu * total_rf_trend
                Ptot_pdf_norm['mu'] = mu
            # sample from normal distribution and saves global value of Ptot
            # (that must be equalled or exceeded) for each year
            annual_limit = np.random.normal(
                loc=Ptot_pdf_norm['mu'], scale=Ptot_pdf_norm['sigma'])
            try:
                annual_limit = np.clip(
                    annual_limit, Ptot_pdf_norm['trunc_interval'][0],
                    Ptot_pdf_norm['trunc_interval'][1])
            except KeyError:
                # ...just in case
                if annual_limit < 0.:
                    annual_limit = 0.
            self._Ptot_ann_global[syear] = annual_limit
            Storm_total_local_year = np.zeros(
                (self._max_numstorms, num_opennodes))
            Storm_running_sum = np.zeros((2, num_opennodes))
            # ^ 1st col is running total, 2nd is data to add to it
            ann_cum_Ptot_gauge = np.zeros(numgauges)
            self._entries = 0
            storm_count = 0
            for storm in range(self._max_numstorms):
                int_arr_val = genextreme.rvs(c=Int_arr_pdf_GEV['shape'],
                                             loc=Int_arr_pdf_GEV['mu'],
                                             scale=Int_arr_pdf_GEV['sigma'])
                try:
                    int_arr_val = np.clip(
                        int_arr_val, Int_arr_pdf_GEV['trunc_interval'][0],
                        Int_arr_pdf_GEV['trunc_interval'][1])
                except KeyError:
                    # ...just in case
                    if int_arr_val < 0.:
                        int_arr_val = 0.
                # ^Samples from distribution of interarrival times (hr). This
                # can be used to develop STORM output for use in rainfall-
                # runoff models or any water balance application.
                # sample uniformly from storm center matrix from grid with 10 m
                # spacings covering basin:
# NOTE DEJH believes this should be a true random spatial
                # sample in a next iteration
                center_val_X = np.random.choice(Xxin)
                center_val_Y = np.random.choice(Yyin)
                # ^sample uniformly from storm center matrix from grid with
                # even spacings within a specified buffer around the basin.
                North = center_val_Y
                East = center_val_X

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
                # value of coord should be set to storm center selected
                # (same below)
                cx = East
                cy = North
                r = np.sqrt(area_val/np.pi)  # value here should be selected
                # based on area above in meters to match the UTM values in
                # North and East vectors.
                # Determine which gauges are hit by Euclidean geometry:
                if self._mode == 'simulation':
                    gdist = (Xin-cx)**2 + (Yin-cy)**2
                elif self._mode == 'validation':
                    gdist = (Easting-cx)**2 + (Northing-cy)**2
                mask_name = (gdist <= r**2)  # this is defacto Mike's aa
                # this short circuits the storm loop in the case that the storm
                # does not affect any 'gauging' location
                if np.all(np.equal(mask_name, False)):
                    continue
                storm_count += 1
                master_storm_count += 1
                gauges_hit = np.where(mask_name)[0]
                num_gauges_hit = gauges_hit.size
                # this routine below allows for orography in precip by first
                # determining the closest gauge and then determining its
                # orographic grouping
                cc = np.argmin(gdist)
                closest_gauge = np.round(Zz[cc])  # this will be compared
                # against orographic gauge groupings to determine the
                # appropriate set of intensity-duration curves
                ######

                # This routine below determines to which orographic group the
                # closest gauge to the storm center belongs to, and censors the
                # number of curves accordingly
                # missing top curve in GR1, top and bottom curves for GR2, and
                # bottom curve for GR3
                # new version of orography compares local 'gauge' elevation to
                # elevation bands called OroGrp, defined above
#### NOTE again, DEJH thinks this could be simplified a lot
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
                duration_val = genextreme.rvs(c=Duration_pdf_GEV['shape'],
                                              loc=Duration_pdf_GEV['mu'],
                                              scale=Duration_pdf_GEV['sigma'])
                # round to nearest minute for consistency with measured data:
                duration_val = round(duration_val)
                # hacky fix to prevent occasional < 0 values:
                # (I think because Matlab is able to set limits manually)
                try:
                    duration_val = np.clip(
                        duration_val, Duration_pdf_GEV['trunc_interval'][0],
                        Duration_pdf_GEV['trunc_interval'][1])
                except KeyError:
                    # ...just in case
                    if duration_val < 0.:
                        duration_val = 0.
                # we're not going to store the calendar time (DEJH change)

                # original curve# probs for 30%-20%-10%: [0.0636, 0.0727,
                # 0.0819, 0.0909, 0.0909, 0.0909, 0.0909, 0.0909, 0.1001,
                # 0.1090, 0.1182]
                # original curve# probs are modified as below
                # add weights to reflect reasonable probabilities that favor
                # lower curves:
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
                # ...these curves are based on empirical data from Walnut Gulch
# NOTE DEJH wants to know exactly how these are defined

                if FUZZMETHOD == 'MS':
                    fuzz_int_val = np.random.choice(fuzz)
                elif FUZZMETHOD == 'DEJH':
                    # this dist should look identical, w/o discretisation
                    fuzz_int_val = FUZZWIDTH * 2. * (np.random.rand() - 0.5)
                else:
                    raise NameError
                intensity_val += fuzz_int_val
                # ^this allows for specified fuzzy tolerance around selected
                # intensity (but it can go -ve)
                # formerly, here MS used a rounding and threshold to prevent
                # storms with a value < 1. We're going to remove the rounding
                # and threshold at zero instead. (below)

                # This scales the storm center intensity upward, so the values
                # at each gauge are realistic once the gradient is applied.
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
                    loc=Recess_pdf_norm['mu'], scale=Recess_pdf_norm['sigma'])
                try:
                    recess_val = np.clip(
                        recess_val, Recess_pdf_norm['trunc_interval'][0],
                        Recess_pdf_norm['trunc_interval'][1])
                except KeyError:
                    pass  # this one is OK <0., I think
                # this pdf of recession coefficients determines how intensity
                # declines with distance from storm center (see below)
                # determine cartesian distances to all hit gauges and
                # associated intensity values at each gauge hit by the storm
                # This is a data storage solution to avoid issues that can
                # arise with slicing grid areas with heavy tailed sizes
                if self._mode == 'validation':
                    xcoords = Easting
                    ycoords = Northing
                else:
                    xcoords = Xin
                    ycoords = Yin
                self._entries = np.sum(mask_name)
                entries = self._entries
                # NOTE _gauge_dist_km only contains the nodes under the storm!
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
                # _rain_int_gauge has already been zeroed earlier in loop, so
                self._temp_dataslots2[:entries] = self._gauge_dist_km[:entries]
                # this copy would not be necessary if we didn't want to
                # preserve self._gauge_dist_km
                np.square(self._temp_dataslots2[:entries],
                          out=self._temp_dataslots2[:entries])
                self._temp_dataslots2[:entries] *= -2. * recess_val**2
                np.exp(self._temp_dataslots2[:entries],
                       out=self._temp_dataslots2[:entries])
                self._temp_dataslots2[:entries] *= intensity_val
                mask_incl_closed = IDs_open[mask_name]
                self._rain_int_gauge[mask_incl_closed] = self._temp_dataslots2[
                    :entries]
                # calc of _rain_int_gauge follows Rodriguez-Iturbe et al.,
                # 1986; Morin et al., 2005 but sampled from a distribution
                # only need to add the bit that got rained on, so:
                self._temp_dataslots2[:entries] *= duration_val / 60.
                ann_cum_Ptot_gauge[mask_name] += self._temp_dataslots2[
                    :entries]
                # collect storm total data for all gauges into rows by storm
                Storm_total_local_year[storm, :] = (
                    self._rain_int_gauge[opennodes] * duration_val / 60.)

                self._Storm_total_local_year = Storm_total_local_year
                Storm_running_sum[1, :] = Storm_total_local_year[storm, :]
                self._Storm_running_sum = Storm_running_sum
                np.nansum(Storm_running_sum, axis=0,
                          out=Storm_running_sum[0, :])
                if np.any(Storm_total_local_year < 0.):
                    raise ValueError(syear, storm)
                self._median_rf_total = np.nanmedian(Storm_running_sum[0, :])

                if self._savedir is not False:
                    if type(self._savedir) is str:
                        rowID = storm_count - 1
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

                if yield_storms is True:
                    yield (duration_val, int_arr_val)
                # now blank the field for the interstorm period
                self._rain_int_gauge.fill(0.)
                if self._median_rf_total > self._Ptot_ann_global[syear]:
                    # we're not going to create Ptotal_local for now... just
                    break
                if storm + 1 == self._max_numstorms:
                    raise ValueError('_max_numstorms set too low for this run')

            if type(self._savedir) is str:
                # crop it down to save:
                [maxrow, maxcol] = np.argwhere(self._Storm_matrix).max(axis=0)
                np.savetxt(os.path.join(
                    self._savedir, 'Storm_matrix_y' + str(syear) + '.csv'),
                    self._Storm_matrix[:(maxrow+1), :(maxcol+1)])
            self._total_rf_year[opennodes] = Storm_running_sum[0, :]
            if yield_years is True:
                yield storm_count

        if self._savedir is True:
            # crop it down for convenience:
            [maxrow, maxcol] = np.argwhere(self._Storm_matrix).max(axis=0)
            self._Storm_matrix = self._Storm_matrix[:(maxrow+1), :(maxcol+1)]



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
for dt, interval_t in rain.yield_storms(storm_radial_weakening_gaussian={
                                            'sigma': 0.08, 'mu': 0.25,
                                            'trunc_interval': (0.15, 0.67)}):
    count += 1
    total_t += dt + interval_t
    print rain._median_rf_total
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
# total_storms = 0.
# for storms_in_year in rain.yield_years():
#     count += 1
#     total_storms += storms_in_year
#     print(storms_in_year)
#     imshow_grid_at_node(mg, 'rainfall__total_depth_per_year', cmap='jet')
#     show()
