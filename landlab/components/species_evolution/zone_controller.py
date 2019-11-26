#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ZoneController of SpeciesEvolver."""
import numpy as np
from scipy.ndimage.measurements import label

from . import ZoneSpecies, Zone
from .record import Record
from .zone import _update_zones


class ZoneController(object):
    """Controls zones and populates them with species.

    This object manages 'zones' that are used to evaluate the spatial aspect
    of species. A zone represents a portion of a model grid. It is made up of
    spatially continuous grid nodes.

    This controller creates zones using the parameter, ``zone_function`` set at
    initialization. This function identifies all of the grid nodes where zones
    are to be created. A zone is created for each cluster of spatially
    continuous nodes. Zones are updated also using this function when the
    ``run_one_step`` method of this controller is called.

    The structure of an example model grid is diagrammed below to demonstrate
    how zones are created. The grid contains six columns and six rows. In this
    example, the function returns True where node values are greater than 0.
    Nodes marked with an ``*`` are nodes where a zone exists and the mask is
    True. All other nodes are marked with a ``·``. A zone is created for each
    cluster of continuous nodes where the mask is True.

    values         mask returned
    evaluated      by zone function
    0 0 0 0 0 0    · · · · · ·
    0 0 4 5 4 0    · · * * * ·
    0 0 6 3 0 0    · · * * · ·
    0 0 0 0 0 0    · · · · · ·
    0 6 4 0 4 0    · * * · * ·
    2 0 0 0 0 0    * · · · · ·

    The above example is continued in the following four diagrams that
    demonstrate how individual zones are identified. Each zone is marked with a
    ``x``, ``o``, ``+``, or ``@``. Clusters can be identified using ``D8``
    where diagonal neighbors are included or ``D4`` where diagonal neighbors
    are excluded. A minimum zone area can be enforced with the ``minimum_area``
    parameter.

    D8             D4             D8             D4
    min area = 0   min area = 0   min area = 2   min area = 2
    · · · · · ·    · · · · · ·    · · · · · ·    · · · · · ·
    · · + + + ·    · · + + + ·    · · + + + ·    · · + + + ·
    · · + + · ·    · · + + · ·    · · + + · ·    · · + + · ·
    · · · · · ·    · · · · · ·    · · · · · ·    · · · · · ·
    · x x · o ·    · @ @ · o ·    · x x · · ·    · x x · · ·
    x · · · · ·    x · · · · ·    x · · · · ·    · · · · · ·

    By default, ``ZoneSpecies`` are used with this controller, and the
    following paragraphs make that assumption. Speciation occurs when the
    species exists in more than one zone once the allopatric wait time has been
    exceeded in that zone. See ZoneSpecies documentation for more about
    allopatric wait time.

    Here, a different example grid demonstrates the temporal connectivity of
    zones. The grid represents the time, ``T0`` with the nodes of a zone
    marked with ``x``. The following examples will use D8 neighborhoods and a
    minimum zone area of 0.

    T0
    · · · · · ·
    · · · · · ·
    · x x x x ·
    · x x x x ·
    · · · · · ·
    · · · · · ·

    Below are variations of the grid at a later time, ``T1`` in four
    variations where each contains one zone. In ``T1a``, ``T1b``, and ``T1c``
    the zone stayed the same, moved, and changed size, respectively. Species
    migrate with the zone when at least one zone node overlaps between the two
    time steps. However, in ``T1d``, no nodes overlaps, therefore species
    do not disperse from the zone in T0 to the zone in T1d.

    T1a            T1b            T1c            T1d
    · · · · · ·    · · · · · ·    · · · · · ·    · · · · · ·
    · · · · · ·    · · · · · ·    · · · * · *    · · · · · ·
    · * * * * ·    · · · · · ·    · · * * * *    · · · · · ·
    · * * * * ·    * * * * · ·    · · · · · ·    · · · · · ·
    · · · · · ·    * * * * · ·    · · · · · ·    * * * * · ·
    · · · · · ·    · · · · · ·    · · · · · ·    * * * * · ·

    Another ``T1`` variation, now demonstrating two zones, ``*`` and ``x``.
    Multiple zones overlapping a zone in the prior time step can be interpreted
    as a zone that fragmented, which may affect resident species. The number of
    zone fragmentations can be viewed in the ``record_data_frame`` attribute.

    T1e
    · · · · · ·
    · · · · · ·
    · · · * * *
    x x · * * ·
    x x · · · ·
    x x · · · ·

    The controller had to decide which of the two clusters of continuous nodes
    in T1 should be designated as the same zone in T0. This decision must be
    made in general when multiple clusters overlap the same zone in the prior
    time step. The zone in the current time step that overlaps the prior time
    step zone the most becomes the same zone in the earlier time step. In this
    example, the cluster to the right overlapped four nodes and the left
    cluster overlapped only one node, therefore the right cluster became the
    star zone. However, this is merely for creating new zones objects.

    The grid diagrammed below continues from T1e. The continuous nodes
    overlapped two zones in T1e. When multiple zones overlap, one zone is
    assumed to be the prior zone and the others are considered captured zones.
    The number of zone captures can be viewed in the ``record_data_frame``
    attribute.

    T2
    · · · · · ·
    · · · · · ·
    · · · · · ·
    · x x x · ·
    · x x x · ·
    · · · · · ·

    The controller had to again decide which of the two clusters of continuous
    nodes in T1e should be designated as the same zone in T2. This decision
    must be made in general when multiple clusters in the prior time step
    overlap a zone in the current time step. The zone in the prior time step
    that overlaps the current time step zone the most becomes the zone in the
    earlier time step. In this example, the cluster to the left overlapped two
    nodes and the right cluster overlapped only one node, therefore the new
    zone keeps the designation of the left cluster. However, this is merely for
    creating new zone objects.

    Note in the example above that zones are created throughout the grid,
    including boundaries, wherever nodes meet the conditions set within the
    ``zone_function``. Boundaries can be excluded in zones by evaluating if
    nodes are boundaries in ``zone_function``.

    Examples
    --------
    Import modules used in the following examples.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.species_evolution import ZoneController

    The first example uses the default parameters of ZoneController.

    Create a model grid and an elevation field for this grid.

    >>> mg = RasterModelGrid((3, 7))
    >>> z = mg.add_zeros('topographic__elevation', at='node')

    Set elevation to 1 for some nodes.

    >>> z[[9, 10, 11, 12]] = 1

    Define a zone function that returns a boolean array where `True` values
    indicate the nodes where zones can be created.

    >>> def zone_func(grid):
    ...     z = grid.at_node['topographic__elevation']
    ...     return z == 1

    Instantiate ZoneController. Only one zone exists because the nodes that
    were set to one are adjacent to each other in the grid.

    >>> zc = ZoneController(mg, zone_func)
    >>> zc.record_data_frame[['time', 'zone_count']]
       time  zone_count
    0     0           1

    Populate each zone with one species.

    >>> species = zc.populate_zones_uniformly(1)
    >>> len(species)
    1

    A change in elevation is forced to demonstrate a zone fragmentation, and
    then the zones are updated by advancing the record time by 1000.

    >>> z[10] = 0
    >>> zc.run_one_step(1000)

    Two zones now exist because the zone in time 0 fragmented into two zones.

    >>> zc.record_data_frame[['time', 'zone_count', 'fragmentation_count']]
       time  zone_count  fragmentation_count
    0     0           1                  NaN
    1  1000           2                  2.0

    A change in elevation is forced again, this time to demonstrate zone
    capture where multiple zones are overlapped by a zone in the later time
    step. Statistics of the capture can be attained with ``record_data_frame``.

    >>> z[10] = 1
    >>> zc.run_one_step(1000)
    >>> zc.record_data_frame[['time', 'zone_count', 'capture_count',
    ...     'area_captured_sum', 'area_captured_max']]
       time  zone_count  capture_count  area_captured_sum  area_captured_max
    0     0           1            NaN                NaN                NaN
    1  1000           2            0.0                0.0                0.0
    2  2000           1            1.0                2.0                2.0

    The follow example demonstrates non-default ZoneController parameters.

    >>> mg = RasterModelGrid((3, 7))
    >>> z = mg.add_zeros('topographic__elevation', at='node')

    Similar to the prior example, define a zone function that returns a boolean
    array where `True` values indicate the nodes where zones can be created.

    >>> def zone_func(grid):
    ...     z = grid.at_node['topographic__elevation']
    ...     return z == 1

    Set elevation to 1 for nodes so that two clusters of nodes within the zone
    mask exist.
    >>> z[[9, 10, 12]] = 1

    Instantiate ZoneController with options.

    >>> zc = ZoneController(mg, zone_func, minimum_area=2, initial_time=100)

    Only one zone exist, despite two clusters of nodes meeting the zone
    definition, because the ``minimum_area`` was set to 2. Also, the first
    time in the record was set by the ``initial_time`` parameter.

    zc.record_data_frame[['time', 'zone_count']]
       time  zone_count
    0   100           1
    """
    def __init__(
        self, grid, zone_function, minimum_area=0, neighborhood_structure='D8',
        initial_time=0, **kwargs
    ):
        """Initialize a species.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        zone_function : function
            A function that return a mask of the total zone extent. The first
            input parameter of this function must be `grid`.
        minimum_area : float, optional
            The minimum area of the zones that will be created.
        neighborhood_structure : {'D8', 'D4'}, optional
            The structure describes how zones are identified. The default,
            'D8' evaluates the eight neighboring nodes. The diagonal
            neighboring nodes are excluded when 'D4' is selected.
        initial_time : float, int, optional
            The initial time. The unit of time is unspecified within the
            controller. The default is 0.
        kwargs
            Keyword arguments for ``zone_function``. Do not include ``grid``
            in kwargs because ``grid``, the first parameter of this method, is
            automatically added to ``kwargs``.
        """
        # Set parameters.

        self._grid = grid
        self._zone_func = zone_function
        self._zone_params = kwargs
        self._min_area = minimum_area
        self._record = Record(initial_time)
        if neighborhood_structure in ['D8', 'D4']:
            self._neighborhood_struct = neighborhood_structure
        else:
            raise ValueError("`neighborhood_structure` must be 'D8' or 'D4'")

        # Set record initial values.

        self._record.set_value('zone_count', np.nan)
        self._record.set_value('fragmentation_count', np.nan)
        self._record.set_value('capture_count', np.nan)
        self._record.set_value('area_captured_sum', np.nan)
        self._record.set_value('area_captured_max', np.nan)

        # Include `grid` in the zone params dictionary.

        self._zone_params['grid'] = self._grid

        # Set initial zones.

        initial_zone_extent = zone_function(**kwargs)
        self._zones = self._get_zones_with_mask(initial_zone_extent)

        self._record.set_value('zone_count', len(self._zones))

    @property
    def zones(self):
        """The zones of the SpeciesController."""
        return self._zones

    @property
    def record_data_frame(self):
        """A DataFrame of SpeciesEvolver variables over time."""
        return self._record.dataframe

    def populate_zones_uniformly(
        self, count, species_type=ZoneSpecies, **kwargs
    ):
        """Populate each zone with the same species count.

        Parameters
        ----------
        count : int
            The count of species to populate to each zone.
        species_type : type of Species
            A Species type that takes a list of Zones as its first parameter.
        kwargs : dictionary
            Keyword arguments of ``species_type``.
        """
        species = []

        for z in self._zones:
            species.extend([species_type([z], **kwargs) for _ in range(count)])

        return species

    def run_one_step(self, dt):
        """Update the zones for a single timestep.

        This method advances time in the record and determines the connectivity
        of zones between the current and prior time steps.

        Parameters
        ----------
        dt : float
            The model time step duration.
        """
        self._record.advance_time(dt)

        # Resolve the spatiotemporal connectivity of the prior time step zones
        # to the new zones.

        prior_zones = self._zones

        zone_mask = self._zone_func(**self._zone_params)
        new_zones = self._get_zones_with_mask(zone_mask)

        self._zones = _update_zones(
            self._grid, prior_zones, new_zones, self._record
        )

        self._record.set_value('zone_count', len(self._zones))

    def _get_zones_with_mask(self, mask):
        """Get zones using a mask.

        Parameters
        ----------
        mask : ndarray
            A boolean array with the grid number of nodes where `True` values
            are nodes within the extent of all the zones to be created.

        Returns
        -------
        list of Zones
            The discrete zones identified in the mask.
        """
        grid = self._grid

        # Label clusters of `True` values in `mask`.

        if self._neighborhood_struct == 'D8':
            s = 3 * [[1, 1, 1]]
        elif self._neighborhood_struct == 'D4':
            s = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]

        cluster_arr, cluster_ct = label(mask.reshape(grid.shape), structure=s)

        # Create zones for clusters.

        zones = []

        for i in range(1, cluster_ct + 1):
            mask = cluster_arr == i

            if mask.sum() * grid.cellarea >= self._min_area:
                # Create a zone if cluster area exceeds the minimum area.
                zones.append(Zone(mask))

        return zones
