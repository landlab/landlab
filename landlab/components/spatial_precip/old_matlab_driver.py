# %----------------------------------------------------------------------------
# %STOchastic Rainfall Model (STORM). A rainstorm generator in this case based on empirical data (Walnut Gulch, AZ)
#
# %Code name: WG_storms_v3_01  %this version allows for simulations that are much longer than the input series in order to compare the distributions of storm characteristics to the original.
# %It also includes the 'mean gauge approach' to determine when a simulation year stops. This involves summing storm totals at each gauge for each year until the mean for all gauges exceeds the
# %selected annual total precip value. It also allows for fuzzy selection of intensity at the storm center based on a fixed value of duration and
# %incorporates orographic effects, wherein there are separate intensity-duration curves derived for three intervals of basin elevation (1200-1350m, 1351-1500m, 1501-1650m)
# %Current version also includes interarrival times between storms, allowing for output to drive other model frameworks (rainfall-ruonff, water balance,LEMs)
# %This version will also include output at each grid location, rather than only at gauge locations.
# %Author: Michael Singer 2017
# %Date created: 2015-6
# %----------------------------------------------------------------------------

import numpy as np
import os
import inspect
from six.moves import range
from matplotlib.pyplot import figure
from scipy.stats import genextreme
from landlab import RasterModelGrid, CLOSED_BOUNDARY, Component


class PrecipitationDistribution(Component):

    def __init__(self, grid, mode='simulation', number_of_simulations=1,
                 number_of_years=1, buffer_width=5000, ptot_scenario='ptotC',
                 storminess_scenario='stormsC', save_outputs=None,
                 path_to_input_files='/Users/daniel/development/landlab/landlab/components/spatial_precip'):
        """
        It's on the user to ensure the grid is big enough to permit the buffer.
        save_outputs : if not None, str path to save

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
        self._rain_int_gauge = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots1 = np.zeros(gaugecount, dtype='float')
        self._temp_dataslots2 = np.zeros(gaugecount, dtype='float')
        self._numsims = number_of_simulations
        self._numyrs = number_of_years
        assert ptot_scenario in ('ptotC', 'ptot+', 'ptot-', 'ptotT+', 'ptotT-')
        self._ptot_scenario = ptot_scenario
        assert storminess_scenario in ('stormsC', 'storms+', 'storms-',
                                       'stormsT+', 'stormsT-')
        self._storms_scenario = storminess_scenario
        self._buffer_width = buffer_width
        self._savedir = save_outputs

        self._max_numstorms = 1000
        # This is for initializing matrices. Trailing zeros are deleted from
        # the matrix at the end of the code.

        self._Storm_matrix = np.zeros(
            (self._max_numstorms*number_of_years, 11))

        self._path = path_to_input_files

    def simple_run(self):
        # what's the dir of this component?
        # this is a nasty hacky way for now
        thisdir = self._path
        # unnecessary as related to output file gen & documentation?
        # t0 = time()
        # t1 = [datestr(floor(now)) '_' datestr(rem(now,1))];
        # t2=regexprep(t1,'[^\w'']',''); %for naming output directories and files by current date/time
        # mkdir('C:\bliss\sacbay\papers\WG_Rainfall_Model\model_output\',t2)

        # Initialize variables for annual rainfall total (mm/h) storm center location (RG1-RG85), etc.

        # This scalar specifies the fractional change in intensity per year
        # when storm_trend is applied in STORMINESS_SCENARIO
# NOTE this needs to be set dynamically
        ptot_scaling_factor = 0.07
        # This scalar specifies the fractional change in intensity per year
        # when storm_trend is applied in STORMINESS_SCENARIO
        storminess_scaling_factor = 0.01
        # This scalar specifies the value of fractional step change in
        # intensity when storms+ or storms- are applied in STORMINESS_SCENARIO
        storm_stepchange = 0.25

        # add variable for number of simulations of simyears
        simyears = self._numyrs  # number of years to simulate
        numcurves = 11  # number of intensity-duration curves (see below for curve equations)

        storm_scaling = 1.  # No storm scaling, as problem appears to be fixed with smaller grid spacing.
        #%storm_scaling = 1.15  # This scales the storm center intensity upward, so the values at each gauge are realistic once the gradient is applied.

# NOTE risk here that matlab & Python present mirrored GEV & EV dists
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Ptot_pdf % This is the pdf fitted to all available station precip data (normal dist). It will be sampled below.
        # #### This to be replaced by a Ptot_mu and Ptot_sigma
# NOTE right now we ignore all poss scenarios, i.e., use the C cases (? Check)
        if self._ptot_scenario == 'ptot+':
            Ptot_pdf_norm = {'sigma': 63.9894, 'mu': 271.,
                             'trunc_interval': (np.nan, np.nan)}
        elif self._ptot_scenario == 'ptot-':
            Ptot_pdf_norm = {'sigma': 63.9894, 'mu': 143.,
                             'trunc_interval': (np.nan, np.nan)}
        else:
            Ptot_pdf_norm = {'sigma': 63.9894, 'mu': 207.489,
                             'trunc_interval': (1., 460.)}
        # the trending cases need to be handled in the loop

        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Duration_pdf % This is the pdf fitted to all available station duration data (GEV dist). It will be sampled below.
        # #### matlab's GEV is (shape_param, scale(sigma), pos(mu))
        # note that in Scipy, we must add a minus to the shape param for a GEV
        # to match Matlab's implementation
        Duration_pdf_GEV = {'shape': -0.570252, 'sigma': 35.7389, 'mu': 34.1409,
                            'trunc_interval': (1., 1040.)}
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Area_pdf % This is the pdf fitted to all available station area data (EV dist). It will be sampled below.
        # #### matlab's EV is (mu, sigma)
        Area_pdf_EV = {'shape': 0., 'sigma': 2.83876e+07, 'mu': 1.22419e+08,
                       'trunc_interval': (5.e+06, 3.e+08)}
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Int_arr_pdf % This is the pdf fitted to all available station area data (GEV dist). It will be sampled below.
        Int_arr_pdf_GEV = {'shape': -0.807971, 'sigma': 9.49574, 'mu': 10.6108,
                           'trunc_interval': (0., 120.)}
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Recess_pdf % This is the pdf of storm gradient recession coefficiencts from Morin et al, 2005 (normal dist). It will be sampled below.
        Recess_pdf_norm = {'sigma': 0.08, 'mu': 0.25,
                           'trunc_interval': (0.15, 0.67)}

        opennodes = self.grid.status_at_node != CLOSED_BOUNDARY
        if self._mode == 'validation':
            Easting = np.loadtxt(os.path.join(thisdir, 'Easting.csv'))  # This is the Longitudinal data for each gauge.
            Northing = np.loadtxt(os.path.join(thisdir, 'Northing.csv'))  # This is the Latitudinal data for each gauge. It will be sampled below.
            gauges = np.loadtxt(os.path.join(thisdir, 'gauges.csv'))  # This is the list of gauge numbers. It will be sampled below.
            gauge_elev = np.loadtxt(os.path.join(thisdir, 'gauge_elev.csv'))
            numgauges = gauges.size
        else:
            X1 = self.grid.node_x
            Y1 = self.grid.node_y
            Xin = X1[opennodes]
            Yin = Y1[opennodes]
            Zz = self.grid.at_node['topographic__elevation'][opennodes]
            numgauges = Xin.size  # number of rain gauges in the basin. NOTE: In this version this produces output on a grid, rather than at real gauge locations.

        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\X % This is the Longitudinal data for each grid point. It will be sampled below to determine storm center location.
        X = np.loadtxt(os.path.join(thisdir, 'X.csv'))
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Y % This is the Latitudinal data for each grid point. It will be sampled below to determine storm center location.
        Y = np.loadtxt(os.path.join(thisdir, 'Y.csv'))
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Storm_depth_data % This is the storm depth data for use in model evaluation.
        Storm_depth_data = np.loadtxt(os.path.join(thisdir,
                                                   'Storm_depth_data.csv'))
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Intensity_data % This is the intensity data for use in model evaluation.
        Intensity_data = np.loadtxt(os.path.join(thisdir,
                                                 'Intensity_data.csv'))
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Duration_data % This is the duration data for use in model evaluation.
        Duration_data = np.loadtxt(os.path.join(thisdir,
                                                'Duration_data.csv'))
        # These three are in mm and hr
        # # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Gauge_GR1 % This is the lowest elevation gauge grouping used for orography (1200-1350m).
        # # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Gauge_GR2 % This is the middle elevation gauge grouping used for orography (1351-1500m).
        # # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\Gauge_GR3 % This is the highest elevation gauge grouping used for orography (1501-1650m).
        # Gauge_GR1 = np.loadtxt(os.path.join(thisdir, 'Gauge_GR1.csv'))
        # Gauge_GR2 = np.loadtxt(os.path.join(thisdir, 'Gauge_GR2.csv'))
        # Gauge_GR3 = np.loadtxt(os.path.join(thisdir, 'Gauge_GR3.csv'))
        # load C:\bliss\sacbay\papers\WG_Rainfall_Model\model_input\fuzz % This is a vector of fuzzy tolerace values for intensity selection.
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

        # Unlike MS's original implementation, we get our ET rates from the
        # generator fn, below

        Storm_matrix = self._Storm_matrix
        Ptot_ann_global = np.zeros(simyears)
# NOTE we're trying to avoid needing to use Gauge_matrix_... structures

        Intensity_local_all = 0  # initialize all variables (concatenated matrices of generated output)
        #%Intensity_local_all = zeros(85*simyears,1); %initialize local_all variables (concatenated vector of generated output at each gauge location)
        Storm_totals_all = 0
        Duration_local_all = 0

        storm_count = 0
        master_storm_count = 0
        storm_trend = 0

        for syear in range(simyears):
            calendar_time = 0  # tracks simulation time per year in hours
            storm_trend += storminess_scaling_factor
            Ptotal = 0
            if self._ptot_scenario == 'ptotT+':
                mu = Ptot_pdf_norm.pop('mu')
                mu += mu * ptot_scaling_factor
                Ptot_pdf_norm['mu'] = mu
            elif self._ptot_scenario == 'ptotT-':
                mu = Ptot * ptot_pdf_norm.pop('mu')
                mu -= mu * ptot_scaling_factor
                Ptot_pdf_norm['mu'] = mu
            # sample from normal distribution and saves global value of Ptot
            # (that must be equalled or exceeded) for each year
            Ptot_ann_global[syear] = np.random.normal(
                loc=Ptot_pdf_norm['mu'], scale=Ptot_pdf_norm['sigma'])
            Storm_total_local_year = np.zeros((self._max_numstorms, numgauges))
            Storm_running_sum = np.zeros((2, numgauges))
            # ^ 1st col is running total, 2nd is data to add to it
            ann_cum_Ptot_gauge = np.zeros(numgauges)
            self._entries = 0
            for storm in range(self._max_numstorms):
                int_arr_val = genextreme.rvs(c=Int_arr_pdf_GEV['shape'],
                                             loc=Int_arr_pdf_GEV['mu'],
                                             scale=Int_arr_pdf_GEV['sigma'])
                # ^Samples from distribution of interarrival times (hr). This
                # can be used to develop STORM output for use in rainfall-
                # runoff models or any water balance application.
                self._rain_int_gauge.fill(0.)
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
                # ^Samples from distribution of storm areas
                # value of coord should be set to storm center selected
                # (same below)
                cx = East
                cy = North
                r = np.sqrt(area_val/np.pi)  # value here should be selected based
                # on area above in meters to match the UTM values in North and
                # East vectors.
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
                # save some properties:
                Storm_matrix[master_storm_count, 0] = master_storm_count
                Storm_matrix[master_storm_count, 1] = area_val
                Storm_matrix[master_storm_count, 8] = cx
                Storm_matrix[master_storm_count, 9] = cy
                Storm_matrix[master_storm_count, 10] = syear
                # this routine below allows for orography in precip by first
                # determining the closest gauge and then determining its
                # orographic grouping
                cc = np.argmin(gdist)
                closest_gauge = np.round(Zz[cc])  # this will be compared
                # against orographic gauge groupings to determine the
                # appropriate set of intensity-duration curves
                ######
                Storm_matrix[master_storm_count, 5] = num_gauges_hit
# NOTE used to calc Gauges_hit all. seems redundant. Replace w anal of mask_name as needed (Flag??)


                # This routine below determines to which orographic group the
                # closest gauge to the storm center belongs to, and censors the
                # number of curves accordingly
                # missing top curve in GR1, top and bottom curves for GR2, and
                # bottom curve for GR3
                # new version of orography compares local 'gauge' elevation to
                # elevation bands called OroGrp, defined above
#### NOTE again, DEJH thinks this could be simplified a lot
                if closest_gauge in OroGrp1:
                    # %int_dur_curve_num = (2:numcurves)'; % these were empirically determined based on data from WG (monsoon rainfall only)-lowest orographic group.
                    baa = 'a'
                elif closest_gauge in OroGrp2:
                    # %int_dur_curve_num = (2:numcurves-1)'; % these were empirically determined based on data from WG (monsoon rainfall only)-middle orographic group.
                    baa = 'b'
                elif closest_gauge in OroGrp3:
                    # %int_dur_curve_num = (1:numcurves-1)'; % these were empirically determined based on data from WG (monsoon rainfall only)-highest orographic group.
                    baa = 'c'
                else:
                    raise ValueError('closest_gauge not found in curve lists!')
                duration_val = genextreme.rvs(c=Duration_pdf_GEV['shape'],
                                              loc=Duration_pdf_GEV['mu'],
                                              scale=Duration_pdf_GEV['sigma'])
                # round to nearest minute for consistency with measured data:
                duration_val = round(duration_val)
                # %Duration_global(storm,year) = duration_val;
                Storm_matrix[master_storm_count, 2] = duration_val
                # we're not going to store the calendar time (DEJH change)

                # original curve# probs for 30%-20%-10%: [0.0636 0.0727 0.0819
                # 0.0909 0.0909 0.0909 0.0909 0.0909 0.1001 0.1090 0.1182]
                # original curve# probs are modified as below
                # add weights to reflect reasonable probabilities that favor
                # lower curves:
                if baa == 'a':
                    wgts = [0.0318, 0.0759, 0.0851, 0.0941, 0.0941, 0.0941,
                            0.0941, 0.0941, 0.1033, 0.1121, 0.1213]
                elif baa == 'b':
                    wgts = [0.0478, 0.0778, 0.0869, 0.0959, 0.0959, 0.0959,
                            0.0959, 0.0959, 0.1051, 0.1141, 0.0888]
                else:
                    wgts = [0.0696, 0.0786, 0.0878, 0.0968, 0.0968, 0.0968,
                            0.0968, 0.0968, 0.1060, 0.1149, 0.0591]
                # which curve did we pick?:
                int_dur_curve_val = np.random.choice(numcurves, p=wgts)
                Storm_matrix[master_storm_count, 3] = int_dur_curve_val

                intensity_val = (lambda_[int_dur_curve_val] *
                                 np.exp(-0.508 * duration_val) +
                                 kappa[int_dur_curve_val] *
                                 np.exp(-0.008*duration_val) +
                                 C[int_dur_curve_val])
                # ...these curves are based on empirical data from Walnut Gulch

                fuzz_int_val = np.random.choice(fuzz)
                intensity_val2 = intensity_val + fuzz_int_val
                # ^this allows for specified fuzzy tolerance around selected
                # intensity
# NOTE DEJH believes this is pretty sketchy:
                # if intensity_val2 < 1.:  # cannot have zero or -ve intensity
                #     intensity_val = 1.
                # else:
                #     intensity_val = intensity_val2
                intensity_val = intensity_val2.clip(1.)
                intensity_val *= storm_scaling
                # This scales the storm center intensity upward, so the values
                # at each gauge are realistic once the gradient is applied.
                # Note the clip to 1 (not 0); < 1mm rain is forbidden
                if self._storms_scenario == 'stormsT+':
                    intensity_val += intensity_val * storm_trend
                    # storminess trend is applied and its effect rises each
                    # year of simulation
                    # DEJH has removed the rounding as he believes this will
                    # prevent trending in the intensity_val for low storm_trend
                if self._storms_scenario == 'stormsT-':
                    intensity_val -= intensity_val * storm_trend
# NOTE this needs clarification, as doesn't seem right in MS's code
                if self._storms_scenario == 'storms+':
                    raise("this doesn't work right yet")
                    # intensity_val += intensity_val * storm_stepchange
                if self._storms_scenario == 'storms-':
                    raise("this doesn't work right yet")
                Storm_matrix[master_storm_count, 4] = intensity_val

                # area to determine which gauges are hit:
                recess_val = np.random.normal(
                    loc=Recess_pdf_norm['mu'], scale=Recess_pdf_norm['sigma'])
                # this pdf of recession coefficients determines how intensity
                # declines with distance from storm center (see below)
                Storm_matrix[master_storm_count, 6] = recess_val
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
                # self._gauge_dist_km[mask_name] = xcoords[mask_name]
                # self._gauge_dist_km[mask_name] -= cx
                # np.square(self._gauge_dist_km[mask_name],
                #           out=self._gauge_dist_km[mask_name])
                # self._temp_dataslots[mask_name] = ycoords[mask_name]
                # self._temp_dataslots[mask_name] -= cy
                # np.square(self._temp_dataslots[mask_name],
                #           out=self._temp_dataslots[mask_name])
                # self._gauge_dist_km[mask_name] += self._temp_dataslots[
                #     mask_name]
                # if np.any(self._gauge_dist_km[mask_name] < 0.):
                #     raise ValueError(self._test_save, )
                # np.sqrt(self._gauge_dist_km[mask_name],
                #         out=self._gauge_dist_km[mask_name])
                # self._gauge_dist_km[mask_name] /= 1000.
                # # _rain_int_gauge has already been zeroed earlier in loop, so
                # self._rain_int_gauge[mask_name] = self._gauge_dist_km[
                #     mask_name]
                # np.square(self._rain_int_gauge[mask_name],
                #           out=self._rain_int_gauge[mask_name])
                # self._rain_int_gauge[mask_name] *= -2. * recess_val**2
                # np.exp(self._rain_int_gauge[mask_name],
                #        out=self._rain_int_gauge[mask_name])
                # self._rain_int_gauge[mask_name] *= intensity_val
                # # calc of _rain_int_gauge follows Rodriguez-Iturbe et al.,
                # # 1986; Morin et al., 2005 but sampled from a distribution
                # # only need to add the bit that got rained on, so:
                # self._temp_dataslots[mask_name] = self._rain_int_gauge[
                #     mask_name]
                # self._temp_dataslots[mask_name] *= duration_val / 60.
                # ann_cum_Ptot_gauge[mask_name] += self._temp_dataslots[
                #     mask_name]
                self._entries = np.sum(mask_name)
                entries = self._entries
                # NOTE _gauge_dist_km only contains the nodes under the storm!
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
                self._rain_int_gauge[mask_name] = self._temp_dataslots2[
                    :entries]
                # calc of _rain_int_gauge follows Rodriguez-Iturbe et al.,
                # 1986; Morin et al., 2005 but sampled from a distribution
                # only need to add the bit that got rained on, so:
                self._temp_dataslots2[:entries] *= duration_val / 60.
                ann_cum_Ptot_gauge[mask_name] += self._temp_dataslots2[
                    :entries]
# TESTER
                self._ann_cum_Ptot_gauge = ann_cum_Ptot_gauge
# NOTE DEJH thinks this is going to get unworkable real fast for big grids.
# Need a much better data storage solution. Come back to this when purpose
# gets more obvious.
#             for jj = 1:numgauges
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,1) = syear;']) %year
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,2) = master_storm_count;']) %storm #
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,3) = rain_int_gauge(jj);'])
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,4) = duration_val;'])
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,5) = rain_int_gauge(jj)*duration_val/60;']) %storm total
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,6) = ann_cum_Ptot_gauge(jj);']) %ann cum total (Ptot)
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,7) = int_arr_val;']) %interarrival time in hours
#                 eval(['Gauge_matrix_',num2str(jj),'(master_storm_count,8) = calendar_time;']) %simulation time per year in hours
#             end
#             Intensity_local_all = [Intensity_local_all,rain_int_gauge]; %#ok<AGROW> %collect into vector of all simulated intensities (at all gauges)
#             dur_step(1:numgauges) = duration_val;
#             Duration_local_all = [Duration_local_all,dur_step]; %#ok<AGROW>
                # collect storm total data for all gauges into rows by storm
                Storm_total_local_year[storm, :] = (self._rain_int_gauge *
                                                    duration_val / 60.)
    #             Storm_totals_all = [Storm_totals_all,rain_int_gauge.*duration_val/60]; %#ok<AGROW>
                Storm_matrix[master_storm_count, 7] = (intensity_val *
                                                       duration_val / 60.)
                self._Storm_total_local_year = Storm_total_local_year
    #             Ptotal = zeros(1,numgauges);
    #             for k = 1:numgauges
    #                 Ptotal(k) = nansum(Storm_total_local_year(:,k)); % %sum the annual storm total at each gauge
    #             end
                # This bit replaced with:
                Storm_running_sum[1, :] = Storm_total_local_year[storm, :]
                self._Storm_running_sum = Storm_running_sum
                np.nansum(Storm_running_sum, axis=0, out=Storm_running_sum[0, :])

    #             Ptotal_test = find(nanmedian(Ptotal) > Ptot_ann_global(syear), 1);    %once the median of all gauges exceeds the selected annual storm total, a new simulation year beginsPtotal_test = find(Ptotal > Ptot_ann_global(syear));    %once the first gauge exceeds the selected annual storm total, a new simulation year begins
                # if ~isempty(Ptotal_test)
                #     %eval(['Ptotal_local_',num2str(syear),'(1:numgauges) = Ptotal;'])
                #     %Ptotal_local(syear,1:numgauges) = Ptotal;
                #     for l = 1:numgauges
                #         Ptotal_local(syear,l) = Ptotal(l); %#ok
                #     end
                #     break %end storm lopp and start a new simulation year
                # end

#                yield self._rain_int_gauge

                if np.nanmedian(Storm_running_sum[0, :]) > Ptot_ann_global[syear]:
                    # we're not going to create Ptotal_local for now... just
                    break

#             eval(['Storm_total_local_year_',num2str(syear),'(storm,1:numgauges) = Storm_total_local_year(storm,1:numgauges);']) %collect all local annual storm totals for each gauge.



#
#     Storm_matrix = Storm_matrix(Storm_matrix(:,9)>0,:); %#ok %gets rid of trailing zeros from initialized matrix
#     AB = find(Gauge_matrix_1(:,2)>0); %#ok
#     for CC = 1:numgauges
#         eval(['Gauge_matrix_',num2str(CC),' = Gauge_matrix_',num2str(CC),'(AB,:);'])
#     end
#     Gauges_hit_all(Gauges_hit_all == 0) = NaN; %#ok
#     GZ = find(Intensity_local_all>0);
#     Intensity_all = Intensity_local_all(GZ); %#ok
#     Duration_all = Duration_local_all(GZ); %#ok
#     Storm_totals_all = Storm_totals_all(GZ); %#ok
#     Ptot_ann_global = Ptot_ann_global(2:length(Ptot_ann_global)); %#ok
#     Gauges_hit_all = Gauges_hit_all(2:length(Gauges_hit_all)); %#ok
#
# %     eval(['save model_output\',tx0,'\',t2,'\Ptot_ann_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_global_',t2,' Ptot_ann_global'])
# %     eval(['save model_output\',tx0,'\',t2,'\Storm_matrix_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' Storm_matrix'])
# %     eval(['save model_output\',tx0,'\',t2,'\Gauges_hit_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_all_',t2,' Gauges_hit_all'])
# %     eval(['save model_output\',tx0,'\',t2,'\Storm_totals_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_all_',t2,' Storm_totals_all'])
# %     eval(['save model_output\',tx0,'\',t2,'\Intensity_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_selected_',t2,' Intensity_all'])
# %     eval(['save model_output\',tx0,'\',t2,'\Duration_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_selected_',t2,' Duration_all'])
# %
# %     eval(['save model_output\',tx0,'\',t2,'\ET_matrix_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' ET_matrix'])
#
#     eval(['save ',tx2,'\Ptot_ann_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_global_',t2,' Ptot_ann_global'])
#     eval(['save ',tx2,'\Storm_matrix_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' Storm_matrix'])
#     eval(['save ',tx2,'\Gauges_hit_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_all_',t2,' Gauges_hit_all'])
#     eval(['save ',tx2,'\Storm_totals_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_all_',t2,' Storm_totals_all'])
#     eval(['save ',tx2,'\Intensity_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_selected_',t2,' Intensity_all'])
#     eval(['save ',tx2,'\Duration_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_selected_',t2,' Duration_all'])
#
#     eval(['save ',tx2,'\ET_matrix_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' ET_matrix'])
#
#     for kk = 1:numgauges
#         %eval(['save model_output\',t2,'\Gauge_matrix_',num2str(kk),'_y_',t2,'.txt Gauge_matrix_',num2str(kk), ' -ASCII';])
#         eval(['save ',tx2,'\Gauge_matrix',num2str(kk),'_',num2str(NUMSIMS),'sims_',num2str(NUMSIMYRS),'y_',t2,' Gauge_matrix_',num2str(kk);])
#     end
#
#     boo = etime(clock,t0);
#     runtime_seconds = boo; %#ok
#     runtime_minutes = boo/60 %#ok
#
# end
