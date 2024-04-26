import contextlib

import numpy as np
from scipy.stats import fisk
from scipy.stats import genextreme

from landlab import Component
from landlab import RasterModelGrid


class SpatialPrecipitationDistribution(Component):
    """Generate spatially resolved precipitation events.

    A component to generate a sequence of spatially resolved storms over a
    grid, following a lightly modified version (see below) of the
    stochastic methods of Singer & Michaelides, Env Res Lett 12, 104011,
    2017, & Singer et al., Geosci. Model Dev., accepted, 10.5194/gmd-2018-86.

    The method is heavily stochastic, and at the present time is intimately
    calibrated against the conditions at Walnut Gulch, described in those
    papers. In particular, assumptions around intensity-duration
    calibration and orographic rainfall are "burned in" for now, and are
    not accessible to the user. The various probability distributions
    supplied to the various run methods default to WG values, but are
    easily modified.  This calibration reflects a US desert southwest
    "monsoonal" climate, and the component distinguishes (optionally)
    between two seasons, "monsoonal" and "winter". The intensity-duration
    relationship is shared between the seasons, and so may prove useful in
    a variety of storm-dominated contexts.

    The default is to disable the orographic rainfall functionality of the
    component. However, if orographic_scenario == 'Singer', the component
    requires a 'topographic__elevation' field to already exist on the grid
    at the time of instantiation.

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

    Note that this component cannot simulate the occurrence of more than one
    storm at the same time. Storms that should be synchronous will instead
    occur sequentially, with no interstorm time. This limitation means that
    if enough storms occur in a year that numstorms*mean_storm_duration
    exceeds one year, the number of simulated storms will saturate. This
    limitation may be relaxed in the future.

    The component offers the option to modify the maximum number of storms
    simulated per year. If you find simulations encountering this limit too
    often, you may need to raise this limit. Conversely, it could be lowered
    to reduce memory usage over small grids. However, in increasing the value,
    beware - the component maintains two limit*nnodes arrays, which will chew
    through memory if the limit gets too high. The default will happily
    simulate grids up to around 50 km * 50 km using the default probability
    distributions.

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
        - target_median_total_rainfall_this_season

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
        - Storms can be centred at any spatial coordinate, not just over nodes.
        - Edge buffering is now dynamic; i.e., big storms have a bigger edge
            buffer than smaller storms. Storms can be centered off the grid
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

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid, VoronoiDelaunayGrid
    >>> mg = RasterModelGrid((10, 10), xy_spacing=1000.0)
    >>> rain = SpatialPrecipitationDistribution(mg)

    Calling yield_storms will produce storm-interstorm duration (hr) pairs
    until the model runtime has elapsed.

    >>> np.random.seed(1)
    >>> total_t_each_step = [
    ...     (storm + interstorm) for (storm, interstorm) in rain.yield_storms()
    ... ]
    >>> len(total_t_each_step)
    41
    >>> np.isclose(sum(total_t_each_step) / 24.0, 365.0)
    True

    The actual rainfall intensities during that interval are accessible in the
    'rainfall__flux' field (mm/hr). The storm centre does not have to be over
    the grid, but in this case, it was for the last simulated storm:

    >>> mg.at_node["rainfall__flux"].argmax()
    80

    We can also run the component for only one season (i.e., only using one
    of the pdf sets describing the storm properties):

    >>> for field in ("rainfall__flux", "rainfall__total_depth_per_year"):
    ...     _ = mg.at_node.pop(field)  # clear out the existing fields
    ...
    >>> rain = SpatialPrecipitationDistribution(mg, number_of_years=2)
    >>> np.random.seed(5)
    >>> total_t_each_step = [
    ...     (storm + interstorm)
    ...     for (storm, interstorm) in rain.yield_storms(
    ...         style="monsoonal", monsoon_fraction_of_year=0.35
    ...     )
    ... ]
    >>> np.isclose(sum(total_t_each_step) / 24.0 / 365.0 / 2.0, 0.35)
    True

    Note this behaviour can be stopped by upping monsoon_fraction_of_year:

    >>> np.random.seed(5)
    >>> total_t_each_step = [
    ...     (storm + interstorm)
    ...     for (storm, interstorm) in rain.yield_storms(
    ...         style="monsoonal", monsoon_fraction_of_year=1.0
    ...     )
    ... ]
    >>> np.isclose(round(sum(total_t_each_step) / 24.0 / 365.0 / 2.0, 2), 1.0)
    True

    yield_years yields the number of storms in the last whole year.
    Use 'rainfall__total_depth_per_year' to access the rainfall map for the
    last fully elapsed year, or equivalently, the total_rainfall_last_year
    property. Note the component seamlessly handles non-raster grid types:

    >>> vdg = VoronoiDelaunayGrid(
    ...     np.random.rand(100) * 1000.0, np.random.rand(100) * 1000.0
    ... )
    >>> np.random.seed(3)
    >>> rain = SpatialPrecipitationDistribution(vdg, number_of_years=3)
    >>> storms_each_year = []
    >>> for total_storms in rain.yield_years(
    ...     style="monsoonal", total_rf_trend=0.05, storminess_trend=-0.02
    ... ):
    ...     storms_each_year.append(total_storms)
    ...     assert np.all(
    ...         np.equal(
    ...             vdg.at_node["rainfall__total_depth_per_year"],
    ...             rain.total_rainfall_last_year,
    ...         )
    ...     )
    >>> sum(storms_each_year)
    11

    yield_seasons yields rainfall statistics for individual seasons. Access
    these using the various provided component properties. Note that we can
    get the component to yield a total rainfall that is calibrated to the
    supplied total_rf_gaussians if we set limit to 'total__rainfall' rather
    than 'total_time' (at the cost of exactly matching the season length):

    >>> for field in ("rainfall__flux", "rainfall__total_depth_per_year"):
    ...     _ = mg.at_node.pop(field)  # clear out the existing fields
    ...
    >>> rain = SpatialPrecipitationDistribution(mg, number_of_years=2)
    >>> np.random.seed(5)
    >>> season_list = []
    >>> theoretical_median_rf_season = []
    >>> median_rf_season = []
    >>> median_rf_last_year = []
    >>> mean_rf_season = []
    >>> mean_rf_last_year = []
    >>> for storm_number in rain.yield_seasons(limit="total_rainfall"):
    ...     season_list.append(rain.current_season)
    ...     theoretical_median_rf_season.append(
    ...         rain.target_median_total_rainfall_this_season
    ...     )
    ...     median_rf_season.append(rain.median_total_rainfall_this_season)
    ...     median_rf_last_year.append(rain.median_total_rainfall_last_year)
    ...     mean_rf_season.append(rain.total_rainfall_this_season.mean())
    ...     mean_rf_last_year.append(rain.total_rainfall_last_year.mean())
    ...
    >>> season_list == ["M", "W", "M", "W"]
    True
    >>> [
    ...     meas > sim
    ...     for (meas, sim) in zip(median_rf_season, theoretical_median_rf_season)
    ... ]  # must exceed
    [True, True, True, True]
    >>> np.isclose(median_rf_last_year[0], 0.0)
    True
    >>> for season in (0, 2):  # this property must be the same in both seasons
    ...     np.isclose(median_rf_last_year[season], median_rf_last_year[season + 1])
    ...
    True
    True

    Note that because we work here with medians, the seasonal medians don't sum
    to the year median, but the means do:

    >>> np.isclose(median_rf_last_year[2], median_rf_season[0] + median_rf_season[1])
    False
    >>> np.isclose(mean_rf_last_year[2], mean_rf_season[0] + mean_rf_season[1])
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Singer, M., Michaelides, K., Hobley, D. (2018). STORM 1.0: a simple,
    flexible, and parsimonious stochastic rainfall generator for simulating
    climate and climate change. Geoscientific Model Development  11(9),
    3713-3726. https://dx.doi.org/10.5194/gmd-11-3713-2018

    **Additional References**

    None Listed

    """

    _name = "SpatialPrecipitationDistribution"

    _unit_agnostic = False

    _cite_as = """@Article{gmd-2018-86,
        title={STORM: A simple, flexible, and parsimonious stochastic rainfall
               generator for simulating climate and climate change},
        author={Singer, M. B. and Michaelides, K. and Hobley, D. E. J.},
        journal={Geoscientific Model Development Discussions},
        volume={2018},
        pages={1--25},
        year={2018},
        url={https://www.geosci-model-dev-discuss.net/gmd-2018-86/},
        doi={10.5194/gmd-2018-86}
        }"""

    _info = {
        "rainfall__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm/hr",
            "mapping": "node",
            "doc": "Depth of water delivered per unit time in each storm",
        },
        "rainfall__total_depth_per_year": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm/yr",
            "mapping": "node",
            "doc": "Depth of water delivered in total in each model year",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self, grid, number_of_years=1, orographic_scenario=None, max_numstorms=5000
    ):
        """Create the SpatialPrecipitationDistribution generator component.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab model grid of any type.
        number_of_years : int
            The number of years over which to generate storms.
        orographic_scenario : {None, 'Singer', func}
            Whether to use no orographic rule, or to adopt S&M's 2017
            calibration for Walnut Gulch. Alternatively, provide a function
            here that turns the provided elevation of the storm center into
            a length-11 curve weighting to select which orographic scenario
            to apply.
        """
        super().__init__(grid)

        gaugecount = (grid.status_at_node != grid.BC_NODE_IS_CLOSED).sum()
        self._gauge_dist_km = np.zeros(gaugecount, dtype="float")
        self._temp_dataslots1 = np.zeros(gaugecount, dtype="float")
        self._temp_dataslots2 = np.zeros(gaugecount, dtype="float")
        self._numyrs = number_of_years

        self._max_numstorms = max_numstorms
        # This is for initializing matrices. Trailing zeros are deleted from
        # matrixes at the end of the code.

        assert orographic_scenario in (None, "Singer")
        self._orographic_scenario = orographic_scenario

        # build LL fields:
        self.initialize_output_fields()
        # bind the field to the internal variable:
        self._rain_int_gauge = self._grid.at_node["rainfall__flux"]
        self._total_rf_year = self._grid.at_node["rainfall__total_depth_per_year"]

        # store some info on the open node grid extent:
        open_nodes = self._grid.status_at_node != self._grid.BC_NODE_IS_CLOSED
        self._minx = self._grid.node_x[open_nodes].min()
        self._maxx = self._grid.node_x[open_nodes].max()
        self._miny = self._grid.node_y[open_nodes].min()
        self._maxy = self._grid.node_y[open_nodes].max()
        self._widthx = self._maxx - self._minx
        self._widthy = self._maxy - self._miny
        self._running_total_rainfall_this_year = self._grid.zeros("node")
        self._running_total_rainfall_this_season = self._grid.zeros("node")

        self._open_area = self._grid.cell_area_at_node[open_nodes].sum()
        self._scaling_to_WG = self._open_area / 275710702.0
        # ^ this is the relative size of the catchment compared to WG

    def yield_storms(
        self,
        limit="total_time",
        style="whole_year",
        total_rf_trend=0.0,
        storminess_trend=0.0,
        monsoon_fraction_of_year=0.42,
        monsoon_total_rf_gaussian=(("sigma", 64.0), ("mu", 207.0)),
        monsoon_storm_duration_GEV=(
            ("shape", -0.570252),
            ("sigma", 35.7389),
            ("mu", 34.1409),
            ("trunc_interval", (0.0, 1040.0)),
        ),
        monsoon_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        monsoon_storm_interarrival_GEV=(
            ("shape", -0.807971),
            ("sigma", 9.4957),
            ("mu", 10.6108),
            ("trunc_interval", (0.0, 720.0)),
        ),
        monsoon_storm_radial_weakening_gaussian=(("sigma", 0.08), ("mu", 0.25)),
        winter_total_rf_gaussian=(("sigma", 52.0), ("mu", 1.65)),
        winter_storm_duration_fisk=(
            ("c", 1.0821),
            ("scale", 68.4703),
            ("trunc_interval", (0.0, 5000.0)),
        ),
        winter_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        winter_storm_interarrival_GEV=(
            ("shape", 1.1131),
            ("sigma", 53.2671),
            ("mu", 47.4944),
            ("trunc_interval", (0.0, 720.0)),
        ),
        winter_storm_radial_weakening_gaussian=(
            ("sigma", 0.08),
            ("mu", 0.25),
            ("trunc_interval", (0.15, 0.67)),
        ),
    ):
        """Yield a timeseries giving the number of storms occurring each year
        in a rainfall simulation.

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, 2017 & Singer et al, submitted.

        Parameters
        ----------
        limit : str
            Controls whether a season is defined based on its total rainfall
            (and can be any length), or by its duration (and can have any
            amount of rainfall). One of 'total_time' or 'total_rainfall'.
            If 'total_time', monsoon_fraction_of_year
            sets the fraction of a year occupied by the monsoon.
        style : str
            Controls whether the component seeks to simulate a western US-
            style "monsoonal" climate, a western US-style winter climate,
            or a full year combining both. One of 'whole_year', 'monsoonal',
            or 'winter' These distributions are by default
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
        monsoon_total_rf_gaussian : dict
            Parameters defining the normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        monsoon_storm_duration_GEV : dict
            Parameters defining a generalised extreme value distribution
            controlling the duration of each storm. In minutes.
        monsoon_storm_area_GEV : dict
            Parameters defining a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        monsoon_storm_interarrival_GEV : dict
            Parameters defining a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS. Note that this calibration is specifically to Walnut Gulch,
            which has an area of 275 km**2. The generator directly scales this
            resulting distribution to the area ratio of Walnut Gulch to the
            open cells of the grid. This crudely accounts for the fact that
            bigger catchments will have more storms, but note that the heavy
            tail on this distribution means the default distribution shape
            will not be trustworthy for catchments with big differences in
            size from Walnut Gulch.
        monsoon_storm_radial_weakening_gaussian : dict
            Parameters defining a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.
        winter_total_rf_gaussian : dict
            Parameters defining a normal distribution controlling the total
            rainfall expected in each year. S&M use 'mu' in {143., 271.} for
            step changes up/down in rainfall totals.
        winter_storm_duration_fisk : dict
            Parameters defining a Fisk (i.e., log-logistic) distribution
            controlling the duration of each storm. Note this differs from the
            summer scaling. In Minutes.
        winter_storm_area_GEV is a generalised extreme value distribution
            controlling the plan view area of each storm. S&M use 'shape': 0.,
            which collapses the distribution to a plain extreme value
            distribution.
        winter_storm_interarrival_GEV : dict
            Parameters defining a generalised extreme value
            distribution controlling the interarrival time between each storm.
            In HRS. The same considerations apply here as for the monsoonal
            interstorm equivalent.
        winter_storm_radial_weakening_gaussian : dict
            Parameters defining a normal distribution
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
            yield_storms=True,
            yield_years=False,
            yield_seasons=False,
            limit=limit,
            style=style,
            monsoon_fraction_of_year=monsoon_fraction_of_year,
            total_rf_trend=total_rf_trend,
            storminess_trend=storminess_trend,
            monsoon_total_rf_gaussian=monsoon_total_rf_gaussian,
            monsoon_storm_duration_GEV=monsoon_storm_duration_GEV,
            monsoon_storm_area_GEV=monsoon_storm_area_GEV,
            monsoon_storm_interarrival_GEV=monsoon_storm_interarrival_GEV,
            monsoon_storm_radial_weakening_gaussian=monsoon_storm_radial_weakening_gaussian,
            winter_total_rf_gaussian=winter_total_rf_gaussian,
            winter_storm_duration_fisk=winter_storm_duration_fisk,
            winter_storm_area_GEV=winter_storm_area_GEV,
            winter_storm_interarrival_GEV=winter_storm_interarrival_GEV,
            winter_storm_radial_weakening_gaussian=winter_storm_radial_weakening_gaussian,
        )

    def yield_years(
        self,
        limit="total_time",
        style="whole_year",
        total_rf_trend=0.0,
        storminess_trend=0.0,
        monsoon_fraction_of_year=0.42,
        monsoon_total_rf_gaussian=(("sigma", 64.0), ("mu", 207.0)),
        monsoon_storm_duration_GEV=(
            ("shape", -0.570252),
            ("sigma", 35.7389),
            ("mu", 34.1409),
            ("trunc_interval", (1.0, 1040.0)),
        ),
        monsoon_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        monsoon_storm_interarrival_GEV=(
            ("shape", -0.807971),
            ("sigma", 9.4957),
            ("mu", 10.6108),
            ("trunc_interval", (0.0, 720.0)),
        ),
        monsoon_storm_radial_weakening_gaussian=(
            ("sigma", 0.08),
            ("mu", 0.25),
            ("trunc_interval", (0.15, 0.67)),
        ),
        winter_total_rf_gaussian=(("sigma", 52.0), ("mu", 1.65)),
        winter_storm_duration_fisk=(
            ("c", 1.0821),
            ("scale", 68.4703),
            ("trunc_interval", (1.0, 5000.0)),
        ),
        winter_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        winter_storm_interarrival_GEV=(
            ("shape", 1.1131),
            ("sigma", 53.2671),
            ("mu", 47.4944),
            ("trunc_interval", (0.0, 720.0)),
        ),
        winter_storm_radial_weakening_gaussian=(
            ("sigma", 0.08),
            ("mu", 0.25),
            ("trunc_interval", (0.15, 0.67)),
        ),
    ):
        """Yield a timeseries giving the number if storms occurring each year
        in a rainfall simulation.

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, 2017 & Singer et al, submitted.

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
            In HRS. Note that this calibration is specifically to Walnut Gulch,
            which has an area of 275 km**2. The generator directly scales this
            resulting distribution to the area ratio of Walnut Gulch to the
            open cells of the grid. This crudely accounts for the fact that
            bigger catchments will have more storms, but note that the heavy
            tail on this distribution means the default distribution shape
            will not be trustworthy for catchments with big differences in
            size from Walnut Gulch.
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
            In HRS. The same considerations apply here as for the monsoonal
            interstorm equivalent.
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
            yield_storms=False,
            yield_years=True,
            yield_seasons=False,
            limit=limit,
            style=style,
            total_rf_trend=total_rf_trend,
            storminess_trend=storminess_trend,
            monsoon_fraction_of_year=monsoon_fraction_of_year,
            monsoon_total_rf_gaussian=monsoon_total_rf_gaussian,
            monsoon_storm_duration_GEV=monsoon_storm_duration_GEV,
            monsoon_storm_area_GEV=monsoon_storm_area_GEV,
            monsoon_storm_interarrival_GEV=monsoon_storm_interarrival_GEV,
            monsoon_storm_radial_weakening_gaussian=monsoon_storm_radial_weakening_gaussian,
            winter_total_rf_gaussian=winter_total_rf_gaussian,
            winter_storm_duration_fisk=winter_storm_duration_fisk,
            winter_storm_area_GEV=winter_storm_area_GEV,
            winter_storm_interarrival_GEV=winter_storm_interarrival_GEV,
            winter_storm_radial_weakening_gaussian=winter_storm_radial_weakening_gaussian,
        )

    def yield_seasons(
        self,
        limit="total_time",
        style="whole_year",
        total_rf_trend=0.0,
        storminess_trend=0.0,
        monsoon_fraction_of_year=0.42,
        monsoon_total_rf_gaussian=(("sigma", 64.0), ("mu", 207.0)),
        monsoon_storm_duration_GEV=(
            ("shape", -0.570252),
            ("sigma", 35.7389),
            ("mu", 34.1409),
            ("trunc_interval", (1.0, 1040.0)),
        ),
        monsoon_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        monsoon_storm_interarrival_GEV=(
            ("shape", -0.807971),
            ("sigma", 9.4957),
            ("mu", 10.6108),
            ("trunc_interval", (0.0, 720.0)),
        ),
        monsoon_storm_radial_weakening_gaussian=(
            ("sigma", 0.08),
            ("mu", 0.25),
            ("trunc_interval", (0.15, 0.67)),
        ),
        winter_total_rf_gaussian=(("sigma", 52.0), ("mu", 1.65)),
        winter_storm_duration_fisk=(
            ("c", 1.0821),
            ("scale", 68.4703),
            ("trunc_interval", (1.0, 5000.0)),
        ),
        winter_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        winter_storm_interarrival_GEV=(
            ("shape", 1.1131),
            ("sigma", 53.2671),
            ("mu", 47.4944),
            ("trunc_interval", (0.0, 720.0)),
        ),
        winter_storm_radial_weakening_gaussian=(
            ("sigma", 0.08),
            ("mu", 0.25),
            ("trunc_interval", (0.15, 0.67)),
        ),
    ):
        """Yield a timeseries giving the number if storms occurring each season
        in a rainfall simulation. Only meaningfully different from yield_years
        if style=='whole_year'.

        All default distributions specified as parameters reflect values for
        Walnut Gulch, see Singer & Michaelides, 2017 & Singer et al, submitted.

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
            In HRS. Note that this calibration is specifically to Walnut Gulch,
            which has an area of 275 km**2. The generator directly scales this
            resulting distribution to the area ratio of Walnut Gulch to the
            open cells of the grid. This crudely accounts for the fact that
            bigger catchments will have more storms, but note that the heavy
            tail on this distribution means the default distribution shape
            will not be trustworthy for catchments with big differences in
            size from Walnut Gulch.
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
            In HRS. The same considerations apply here as for the monsoonal
            interstorm equivalent.
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
            yield_storms=False,
            yield_years=False,
            yield_seasons=True,
            limit=limit,
            style=style,
            total_rf_trend=total_rf_trend,
            storminess_trend=storminess_trend,
            monsoon_fraction_of_year=monsoon_fraction_of_year,
            monsoon_total_rf_gaussian=monsoon_total_rf_gaussian,
            monsoon_storm_duration_GEV=monsoon_storm_duration_GEV,
            monsoon_storm_area_GEV=monsoon_storm_area_GEV,
            monsoon_storm_interarrival_GEV=monsoon_storm_interarrival_GEV,
            monsoon_storm_radial_weakening_gaussian=monsoon_storm_radial_weakening_gaussian,
            winter_total_rf_gaussian=winter_total_rf_gaussian,
            winter_storm_duration_fisk=winter_storm_duration_fisk,
            winter_storm_area_GEV=winter_storm_area_GEV,
            winter_storm_interarrival_GEV=winter_storm_interarrival_GEV,
            winter_storm_radial_weakening_gaussian=winter_storm_radial_weakening_gaussian,
        )

    def _run_the_process(
        self,
        yield_storms=True,
        yield_years=False,
        yield_seasons=False,
        limit="total_time",
        style="whole_year",
        monsoon_fraction_of_year=0.42,
        total_rf_trend=0.0,
        storminess_trend=0.0,
        monsoon_total_rf_gaussian=(("sigma", 64.0), ("mu", 207.0)),
        monsoon_storm_duration_GEV=(
            ("shape", -0.570252),
            ("sigma", 35.7389),
            ("mu", 34.1409),
            ("trunc_interval", (1.0, 1040.0)),
        ),
        monsoon_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        monsoon_storm_interarrival_GEV=(
            ("shape", -0.807971),
            ("sigma", 9.4957),
            ("mu", 10.6108),
            ("trunc_interval", (0.0, 720.0)),
        ),
        monsoon_storm_radial_weakening_gaussian=(
            ("sigma", 0.08),
            ("mu", 0.25),
            ("trunc_interval", (0.15, 0.67)),
        ),
        winter_total_rf_gaussian=(("sigma", 52.0), ("mu", 1.65)),
        winter_storm_duration_fisk=(
            ("c", 1.0821),
            ("scale", 68.4703),
            ("trunc_interval", (1.0, 5000.0)),
        ),
        winter_storm_area_GEV=(
            ("shape", 0.0),
            ("sigma", 2.83876e07),
            ("mu", 1.22419e08),
            ("trunc_interval", (5.0e06, 3.0e08)),
        ),
        winter_storm_interarrival_GEV=(
            ("shape", 1.1131),
            ("sigma", 53.2671),
            ("mu", 47.4944),
            ("trunc_interval", (0.0, 720.0)),
        ),
        winter_storm_radial_weakening_gaussian=(
            ("sigma", 0.08),
            ("mu", 0.25),
            ("trunc_interval", (0.15, 0.67)),
        ),
    ):
        """This is the underlying process that runs the component, but it
        should be run by a user through the yield_storms and yield_years
        methods.

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
            In HRS. Note that this calibration is specifically to Walnut Gulch,
            which has an area of 275 km**2. The generator directly scales this
            resulting distribution to the area ratio of Walnut Gulch to the
            open cells of the grid. This crudely accounts for the fact that
            bigger catchments will have more storms, but note that the heavy
            tail on this distribution means the default distribution shape
            will not be trustworthy for catchments with big differences in
            size from Walnut Gulch.
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
            In HRS. The same considerations apply here as for the monsoonal
            interstorm equivalent.
        winter_storm_radial_weakening_gaussian is a normal distribution
            controlling the rate of intensity decline with distance from storm
            center. For more detail see Rodriguez-Iturbe et al., 1986; Morin
            et al., 2005.
        """
        monsoon_total_rf_gaussian = dict(monsoon_total_rf_gaussian)
        monsoon_storm_duration_GEV = dict(monsoon_storm_duration_GEV)
        monsoon_storm_area_GEV = dict(monsoon_storm_area_GEV)
        monsoon_storm_interarrival_GEV = dict(monsoon_storm_interarrival_GEV)
        monsoon_storm_radial_weakening_gaussian = dict(
            monsoon_storm_radial_weakening_gaussian
        )
        winter_total_rf_gaussian = dict(winter_total_rf_gaussian)
        winter_storm_duration_fisk = dict(winter_storm_duration_fisk)
        winter_storm_area_GEV = dict(winter_storm_area_GEV)
        winter_storm_interarrival_GEV = dict(winter_storm_interarrival_GEV)
        winter_storm_radial_weakening_gaussian = dict(
            winter_storm_radial_weakening_gaussian
        )

        FUZZMETHOD = "DEJH"
        FUZZWIDTH = 5.0  # if DEJH
        self._phantom_storm_count = 0
        # ^this property tracks the number of storms in the run that received
        # zero intensity (and thus didn't really exist)
        self._opennodes = self._grid.status_at_node != self._grid.BC_NODE_IS_CLOSED
        self._total_rainfall_last_season = self._grid.zeros("node")

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
        hrsinyr = 24.0 * 365.0
        hrsinmonsoon = monsoon_fraction_of_year * hrsinyr
        hrsinwinter = (1.0 - monsoon_fraction_of_year) * hrsinyr

        assert limit in ("total_rainfall", "total_time")

        assert style in ("whole_year", "monsoonal", "winter")
        if style == "whole_year":
            reps = 2
        else:
            reps = 1

        opennodes = self._opennodes
        num_opennodes = np.sum(opennodes)
        IDs_open = np.where(opennodes)[0]  # need this later
        X1 = self._grid.node_x
        Y1 = self._grid.node_y
        Xin = X1[opennodes]
        Yin = Y1[opennodes]
        try:
            Zz = self._grid.at_node["topographic__elevation"][opennodes]
        except KeyError:
            assert self._orographic_scenario is None
        numgauges = Xin.size  # number of rain gauges in the basin.
        # NOTE: In this version this produces output on a grid, rather than at
        # real gauge locations.

        assert FUZZMETHOD == "DEJH", "The Singer method for fuzz is no longer supported"

        # lambda_, kappa, and C are parameters of the intensity-duration curves
        # of the form: intensity =
        # lambda*exp(-0.508*duration)+kappa*exp(-0.008*duration)+C
        lambda_ = [
            642.2,
            578.0,
            513.8,
            449.5,
            385.3,
            321.1,
            256.9,
            192.7,
            128.4,
            64.1,
            21.0,
        ]
        kappa = [93.1, 83.8, 74.5, 65.2, 55.9, 46.6, 37.2, 27.9, 18.6, 9.3, 0.9]
        C = [4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.25, 0.05]

        # Unlike MS's original implementation, we no longer pull ET values, as
        # this should be a different component.

        self._Ptot_ann_global = np.zeros(simyears)
        self._Ptot_monsoon_global = np.zeros(simyears)

        master_storm_count = 0
        storm_trend = 0

        for syear in range(simyears):
            self._year = syear
            year_time = 0.0  # tracks simulation time per year in hours
            storm_trend += storminess_trend
            year_storm_count = 0
            breaker = False
            Storm_total_local_year = np.zeros((self._max_numstorms, num_opennodes))
            self._storm_running_sum_of_seasons = np.zeros(num_opennodes)
            self._storm_running_sum_1st_seas = np.zeros(num_opennodes)

            storms_yr_so_far = 0
            for seas in range(reps):
                seas_time = 0.0  # tracks elapsed season time in hours
                Storm_running_sum_seas = np.zeros((2, num_opennodes))
                # ^ 1st col is running total, 2nd is data to add to it
                if seas == 0 and style != "winter":
                    self._current_season = "M"
                    # This is the pdf fitted to all available station precip
                    # data (normal dist). It will be sampled below.
                    Ptot_pdf_norm = monsoon_total_rf_gaussian

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
                    self._current_season = "W"
                    Ptot_pdf_norm = winter_total_rf_gaussian
                    Duration_pdf = winter_storm_duration_fisk
                    Area_pdf_EV = winter_storm_area_GEV
                    Int_arr_pdf_GEV = winter_storm_interarrival_GEV
                    Recess_pdf_norm = winter_storm_radial_weakening_gaussian
                    seas_total = hrsinwinter

                if not np.isclose(total_rf_trend, 0.0):
                    mu = Ptot_pdf_norm.pop("mu")
                    mu += mu * total_rf_trend
                    Ptot_pdf_norm["mu"] = mu
                # sample from normal distribution and saves global value of
                # Ptot (that must be equalled or exceeded) for each year
                season_rf_limit = self.calc_annual_rainfall(
                    style=style, monsoon_total_rf_gaussian=Ptot_pdf_norm
                )[seas]
                self._season_rf_limit = season_rf_limit
                self._Ptot_ann_global[syear] += season_rf_limit
                if seas == 0 and style != "winter":
                    self._Ptot_monsoon_global[syear] = season_rf_limit
                Storm_total_local_seas = np.zeros((self._max_numstorms, num_opennodes))
                seas_cum_Ptot_gauge = np.zeros(numgauges)
                self._entries = 0

                for seas_storm_count, storm in enumerate(range(self._max_numstorms)):
                    self._rain_int_gauge.fill(0.0)
                    int_arr_val = genextreme.rvs(
                        c=Int_arr_pdf_GEV["shape"],
                        loc=Int_arr_pdf_GEV["mu"],
                        scale=Int_arr_pdf_GEV["sigma"],
                    )
                    try:
                        int_arr_val = np.clip(
                            int_arr_val,
                            Int_arr_pdf_GEV["trunc_interval"][0],
                            Int_arr_pdf_GEV["trunc_interval"][1],
                        )
                    except KeyError:
                        # ...just in case
                        if int_arr_val < 0.0:
                            int_arr_val = 0.0
                    # now, correct the scaling relative to WG
                    int_arr_val /= self._scaling_to_WG
                    self._int_arr_val = int_arr_val
                    # ^Samples from distribution of interarrival times (hr).
                    # This can be used to develop STORM output for use in
                    # rainfall-runoff models or any water balance application.
                    # sample uniformly from storm center matrix from grid w
                    # 10 m spacings covering basin:

                    area_val = genextreme.rvs(
                        c=Area_pdf_EV["shape"],
                        loc=Area_pdf_EV["mu"],
                        scale=Area_pdf_EV["sigma"],
                    )
                    try:
                        area_val = np.clip(
                            area_val,
                            Area_pdf_EV["trunc_interval"][0],
                            Area_pdf_EV["trunc_interval"][1],
                        )
                    except KeyError:
                        # ...just in case
                        if area_val < 0.0:
                            area_val = 0.0
                    self._area_val = area_val
                    # ^Samples from distribution of storm areas

                    r = np.sqrt(area_val / np.pi)  # value here shd be selected
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
                        gdist = (Xin - cx) ** 2 + (Yin - cy) ** 2
                        mask_name = gdist <= rsq  # this is defacto MS's aa
                        # this short circuits the storm loop in the case that
                        # the storm does not affect any 'gauging' location
                        if np.any(np.equal(mask_name, True)):
                            break

                    self._x = cx
                    self._y = cy
                    year_storm_count += 1
                    seas_storm_count += 1
                    master_storm_count += 1

                    # This routine below determines to which orographic group
                    # the closest gauge to the storm center belongs to, and
                    # censors the number of curves accordingly
                    # missing top curve in GR1, top and bottom curves for GR2,
                    # and bottom curve for GR3
                    # NOTE again, DEJH thinks this could be generalised a lot

                    # original curve# probs for 30%-20%-10%: [0.0636, 0.0727,
                    # 0.0819, 0.0909, 0.0909, 0.0909, 0.0909, 0.0909, 0.1001,
                    # 0.1090, 0.1182]
                    # original curve# probs are modified as below
                    # add weights to reflect reasonable probabilities that
                    # favor lower curves:
                    if self._orographic_scenario is not None:
                        # this routine below allows for orography in precip by
                        # first determining the closest gauge and then
                        # determining its orographic grouping
                        cc = np.argmin(gdist)
                        closest_gauge_z = Zz[cc]  # this will be
                        # compared against orographic gauge groupings to
                        # determine the appropriate set of intensity-duration
                        # curves
                        if self._orographic_scenario == "Singer":
                            wgts = Singer_orographic_rainfall(closest_gauge_z)
                        else:
                            wgts = self._orographic_scenario(closest_gauge_z)
                    elif self._orographic_scenario is None:
                        wgts = [
                            0.0636,
                            0.0727,
                            0.0819,
                            0.0909,
                            0.0909,
                            0.0909,
                            0.0909,
                            0.0909,
                            0.1001,
                            0.1090,
                            0.1182,
                        ]
                    if seas == 0 and style != "winter":
                        duration_val = genextreme.rvs(
                            c=Duration_pdf["shape"],
                            loc=Duration_pdf["mu"],
                            scale=Duration_pdf["sigma"],
                        )
                    else:
                        duration_val = fisk.rvs(
                            c=Duration_pdf["c"], scale=Duration_pdf["scale"]
                        )
                    # hacky fix to prevent occasional < 0 values:
                    # (I think because Matlab is able to set limits manually)
                    try:
                        duration_val = np.clip(
                            duration_val,
                            Duration_pdf["trunc_interval"][0],
                            Duration_pdf["trunc_interval"][1],
                        )
                    except KeyError:
                        # ...just in case
                        if duration_val < 0.0:
                            duration_val = 0.0
                    durationhrs = duration_val / 60.0
                    self._durationhrs = durationhrs
                    year_time += durationhrs
                    seas_time += durationhrs
                    # we will always do the next storm, even if it exceeds the
                    # specified "total" time

                    # which curve did we pick?:
                    int_dur_curve_val = np.random.choice(numcurves, p=wgts)

                    intensity_val = (
                        lambda_[int_dur_curve_val] * np.exp(-0.508 * duration_val)
                        + kappa[int_dur_curve_val] * np.exp(-0.008 * duration_val)
                        + C[int_dur_curve_val]
                    )
                    # ...these curves are based on empirical data from WG

                    # this dist should look identical, w/o discretisation
                    fuzz_int_val = FUZZWIDTH * 2.0 * (np.random.rand() - 0.5)

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
                    if intensity_val < 0.0:
                        intensity_val = 0.0
                        self._phantom_storm_count += 1
                    # note storms of zero intensity are now permitted (though
                    # should hopefully remain pretty rare.)
                    self._intensity_val = intensity_val

                    # area to determine which gauges are hit:
                    recess_val = np.random.normal(
                        loc=Recess_pdf_norm["mu"], scale=Recess_pdf_norm["sigma"]
                    )
                    with contextlib.suppress(KeyError):
                        recess_val = np.clip(
                            recess_val,
                            Recess_pdf_norm["trunc_interval"][0],
                            Recess_pdf_norm["trunc_interval"][1],
                        )
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
                    self._gauge_dist_km[:entries] = np.sqrt(gdist[mask_name]) / 1000.0
                    self._temp_dataslots2[:entries] = gdist[mask_name] / 1.0e6
                    self._temp_dataslots2[:entries] *= -2.0 * recess_val**2
                    np.exp(
                        self._temp_dataslots2[:entries],
                        out=self._temp_dataslots2[:entries],
                    )
                    self._temp_dataslots2[:entries] *= intensity_val
                    mask_incl_closed = IDs_open[mask_name]
                    self._nodes_hit = mask_incl_closed
                    # ^note this is by ID, not bool
                    self._rain_int_gauge[mask_incl_closed] = self._temp_dataslots2[
                        :entries
                    ]
                    # calc of _rain_int_gauge follows Rodriguez-Iturbe et al.,
                    # 1986; Morin et al., 2005 but sampled from a distribution
                    # only need to add the bit that got rained on, so:
                    self._temp_dataslots2[:entries] *= duration_val / 60.0
                    seas_cum_Ptot_gauge[mask_name] += self._temp_dataslots2[:entries]
                    # collect storm totals for all gauges into rows by storm
                    Storm_total_local_seas[storm, :] = (
                        self._rain_int_gauge[opennodes] * duration_val / 60.0
                    )
                    Storm_total_local_year[(storm + storms_yr_so_far), :] = (
                        Storm_total_local_seas[storm, :]
                    )
                    self._max_storm_depth = Storm_total_local_seas[storm, :].max()

                    self._Storm_total_local_seas = Storm_total_local_seas
                    self._Storm_total_local_year = Storm_total_local_year
                    Storm_running_sum_seas[1, :] = Storm_total_local_seas[storm, :]
                    np.nansum(
                        Storm_running_sum_seas, axis=0, out=Storm_running_sum_seas[0, :]
                    )
                    if np.any(Storm_total_local_seas < 0.0):
                        raise ValueError(syear, storm)
                    self._median_seas_rf_total = np.nanmedian(
                        Storm_running_sum_seas[0, :]
                    )
                    self._Storm_running_sum_seas = Storm_running_sum_seas[0, :]

                    if limit == "total_time":
                        if seas_time + int_arr_val > seas_total:
                            int_arr_val = (seas_total - seas_time).clip(0.0)
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
                        raise ValueError("_max_numstorms set too low for this run")
                storms_yr_so_far = seas_storm_count
                self._storm_running_sum_of_seasons += Storm_running_sum_seas[0, :]
                self._total_rainfall_last_season[self._opennodes] = (
                    Storm_running_sum_seas[0, :]
                )
                self._storm_running_sum_1st_seas += Storm_running_sum_seas[0, :]
                if yield_seasons:
                    yield seas_storm_count

            self._total_rf_year[opennodes] = self._storm_running_sum_of_seasons
            if yield_years is True and yield_seasons is False:
                yield year_storm_count

    def calc_annual_rainfall(
        self,
        style="whole_year",
        monsoon_total_rf_gaussian=(("sigma", 64.0), ("mu", 207.0)),
        winter_total_rf_gaussian=(("sigma", 52.0), ("mu", 1.65)),
    ):
        """Return a tuple of rainfall totals (mm) for the year, with entries
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
        >>> mg = RasterModelGrid((10, 10), xy_spacing=500.0)
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> rain = SpatialPrecipitationDistribution(mg)
        >>> mytotals = []
        >>> for yr in range(5):
        ...     mytotals.append(rain.calc_annual_rainfall(style="whole_year"))
        ...
        >>> [len(x) == 2 for x in mytotals]
        [True, True, True, True, True]
        >>> mytotals = []
        >>> for yr in range(3):
        ...     mytotals.append(rain.calc_annual_rainfall(style="monsoonal"))
        ...
        >>> [len(x) == 1 for x in mytotals]
        [True, True, True]
        """
        monsoon_total_rf_gaussian = dict(monsoon_total_rf_gaussian)
        winter_total_rf_gaussian = dict(winter_total_rf_gaussian)

        assert style in ("whole_year", "monsoonal", "winter")
        if style in ("whole_year", "monsoonal"):
            # sample from normal distribution and saves global value of Ptot
            # (that must be equalled or exceeded) for each year
            summer_rf_limit = np.random.normal(
                loc=monsoon_total_rf_gaussian["mu"],
                scale=monsoon_total_rf_gaussian["sigma"],
            )
            try:
                summer_rf_limit = np.clip(
                    summer_rf_limit,
                    monsoon_total_rf_gaussian["trunc_interval"][0],
                    monsoon_total_rf_gaussian["trunc_interval"][1],
                )
            except KeyError:
                # ...just in case
                if summer_rf_limit < 0.0:
                    summer_rf_limit = 0.0
        if style in ("whole_year", "winter"):
            # sample from normal distribution and saves global value of Ptot
            # (that must be equalled or exceeded) for each year
            winter_rf_limit = np.random.normal(
                loc=winter_total_rf_gaussian["mu"],
                scale=winter_total_rf_gaussian["sigma"],
            )
            try:
                winter_rf_limit = np.clip(
                    winter_rf_limit,
                    winter_total_rf_gaussian["trunc_interval"][0],
                    winter_total_rf_gaussian["trunc_interval"][1],
                )
            except KeyError:
                # ...just in case
                if winter_rf_limit < 0.0:
                    winter_rf_limit = 0.0
        if style == "monsoonal":
            return (summer_rf_limit,)
        elif style == "winter":
            return (winter_rf_limit,)
        else:
            return (summer_rf_limit, winter_rf_limit)

    def _locate_storm(self, storm_radius):
        """Because of the way the stats fall out, any simulated storm from the
        distribution must intersect the catchment somewhere.

        Note written in a grid-agnostic fashion.
        """
        stormposx = np.random.rand() * (self._widthx + 2.0 * storm_radius)
        stormposy = np.random.rand() * (self._widthy + 2.0 * storm_radius)
        stormx = self._minx - storm_radius + stormposx
        stormy = self._miny - storm_radius + stormposy
        return stormx, stormy

    @property
    def current_year(self):
        """Get the current year as an int."""
        return self._year

    @property
    def current_season(self):
        """Get the current season.

        'M' is monsoon, 'W' is winter.
        """
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
        """Get the intensity (mm/hr) of the last storm, averaged under the
        storm.

        footprint. Note that duration * intensity != storm max depth.
        """
        return self._intensity_val

    @property
    def total_rainfall_last_season(self):
        """Get the total recorded rainfall over the last (completed) simulated
        season, spatially resolved (mm)."""
        return self._total_rainfall_last_season

    @property
    def total_rainfall_last_year(self):
        """Get the total recorded rainfall over the last (completed) simulated
        year, spatially resolved (mm).

        Equivalent to the field 'rainfall__total_depth_per_year'.
        """
        return self._total_rf_year

    @property
    def total_rainfall_this_season(self):
        """Get the accumulated, spatially resolved total rainfall over the grid
        for the season so far (mm)."""
        self._running_total_rainfall_this_season[self._opennodes] = (
            self._Storm_running_sum_seas
        )
        return self._running_total_rainfall_this_season

    @property
    def total_rainfall_this_year(self):
        """Get the accumulated, spatially resolved total rainfall over the grid
        for the year so far (mm)."""
        self._running_total_rainfall_this_year[self._opennodes] = (
            self._storm_running_sum_1st_seas + self._Storm_running_sum_seas
        )
        return self._running_total_rainfall_this_year

    @property
    def median_total_rainfall_last_season(self):
        """Get the median total rainfall recorded over the open nodes of the
        grid during the last (completed) simulated season (mm)."""
        return np.nanmedian(self._total_rainfall_last_season[self._opennodes])

    @property
    def median_total_rainfall_last_year(self):
        """Get the median total rainfall recorded over the open nodes of the
        grid during the last (completed) simulated year (mm)."""
        return np.nanmedian(self.total_rainfall_last_year[self._opennodes])

    @property
    def median_total_rainfall_this_season(self):
        """Get the accumulated median total rainfall over the open nodes of the
        grid so far this season (mm)."""
        return self._median_seas_rf_total

    @property
    def median_total_rainfall_this_year(self):
        """Get the accumulated median total rainfall over the open nodes of the
        grid so far this year (mm)."""
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
        """Get the coordinates of the center of the last storm as (x, y)."""
        return (self._x, self._y)

    @property
    def target_median_total_rainfall_this_season(self):
        """Get the stochastically generated "target" average total rainfall
        amount over the catchment for the current season.

        If limit == 'total_rainfall', this will be very close to
        median_total_rainfall_last_season. If 'total_time', it will
        diverge from this value.
        """
        return self._season_rf_limit


def Singer_orographic_rainfall(z_closest_node_to_center):
    """Return a set of curve weights for a provided z, assuming an orographic
    rule following that presented in Singer & Michaelides 2017 & Singer et al.
    2018 and applicable specifically to Walnut Gulch. i.e., there are three
    orographic divisions, divided at 1350 m and 1500 m.

    Parameters
    ----------
    z_closest_node_to_center : float
        The elevation of the node closest to the storm center.

    Returns
    -------
    wgts : length 11 list
        The weighting parameters to use in selecting a storm distribution
        curve.
    """
    if z_closest_node_to_center <= 1350:
        wgts = [
            0.0318,
            0.0759,
            0.0851,
            0.0941,
            0.0941,
            0.0941,
            0.0941,
            0.0941,
            0.1033,
            0.1121,
            0.1213,
        ]
    elif 1350 < z_closest_node_to_center <= 1500:
        wgts = [
            0.0478,
            0.0778,
            0.0869,
            0.0959,
            0.0959,
            0.0959,
            0.0959,
            0.0959,
            0.1051,
            0.1141,
            0.0888,
        ]
    elif z_closest_node_to_center > 1500:
        wgts = [
            0.0696,
            0.0786,
            0.0878,
            0.0968,
            0.0968,
            0.0968,
            0.0968,
            0.0968,
            0.1060,
            0.1149,
            0.0591,
        ]
    return wgts


if __name__ == "__main__":
    from matplotlib.pyplot import show

    nx = 17
    ny = 8
    dx = 1000.0
    mg = RasterModelGrid((nx, ny), xy_spacing=dx)

    z = mg.add_zeros("topographic__elevation", at="node")
    z += 1400.0
    rain = SpatialPrecipitationDistribution(mg, number_of_years=1)
    total_t = 0.0
    for count, dt, interval_t in enumerate(
        rain.yield_storms(style="whole_year", limit="total_time")
    ):
        total_t += dt + interval_t
        print(dt, interval_t)
        if count % 100 == 0:
            print("Season:", rain.current_season, "of yr", rain.current_year)
            print("Current storm:", count)
            show()
    print("Effective total years:")
    print(total_t / 24.0 / 365.0)
    print("Storms simulated:")
    print(count)
