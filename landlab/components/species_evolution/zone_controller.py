#!/usr/bin/env python
"""Zone-based SpeciesController of SpeciesEvolver."""
from scipy.ndimage.measurements import label

from landlab.components.species_evolution import ZoneSpecies, Zone
from .record import Record
from .zone import _update_zones


class ZoneController(object):
    """Controls zones and populates them with species.

    This object manages 'zones' that are used to evaluate the spatial aspect
    of species. A zone represents a portion of a model grid. It is made up of
    spatially continuous grid nodes.

    This controller creates zones using the parameter, ``function`` set at
    initialization. This function identifies all of the grid nodes where zones
    are to be created. A zone is created for each cluster of spatially
    continuous nodes. The function is called and zones are updated when this
    controller is initialized and when the ``run_one_step`` method is called.

    The structure of an example model grid is diagrammed below to demonstrate
    how zones are created. The grid contains six columns and six rows. In this
    example, the function returns True where node values are greater than 0.
    Nodes marked with an ``*`` are nodes where a zone exists and the mask is
    True. All other nodes are marked with a ``·``. A zone is created for
    each cluster of continuous nodes where the mask is True.

    values         function
    evaluated      result
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

    Typically, ``ZoneSpecies`` are used with this controller, and the following
    paragraphs make that assumption. Speciation occurs when the species exists
    in more than one zone once the allopatric wait time has been exceeded in
    that zone (see ZoneSpecies documention for more about allopatric wait
    time).

    Here a different example grid demonstrates the temporal connectivity of
    zones. The grid represents the time, ``T0`` with the nodes of a zone
    marked with ``x``. In the rest of the document, examples will use D8
    neighborhoods and a minimum zone area of 0.

        T0
    · · · · · ·
    · · · · · ·
    · x x x x ·
    · x x x x ·
    · · · · · ·
    · · · · · ·

    Below are variations of the grid at a later time, ``T1`` in three
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
    Multiple zones overlapping a one zone in the prior time step can be
    interpreted as a zone that fragmented, which may affect resident species.
    The number of zone fragmentions can be viewed in the ``record`` attribute.

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

    The grid below continues from T1e. The continuous nodes overlapped two
    zones in T1e. When multiple zones overlap, one zone is assumed to be the
    prior zone and the others are considered captured zones. The number of zone
    captures can be viewed in the ``record`` attribute.

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
    overlap the a zone in the current time step. The zone in the prior time
    step that overlaps the current time step zone the most becomes the zone in
    the earlier time step. In this example, the cluster to the left overlapped
    two nodes and the right cluster overlapped only one node, therefore the new
    zone keeps the designation of the left cluster. However, this is merely for
    creating new zones objects.
    """
    def __init__(self, grid, zone_function, minimum_area=0,
                 neighborhood_structure='D8', initial_time=0, **kwargs):
        """Initialize a species.

        Parameters
        ----------
        species_evolver : SpeciesEvolver
            A Landlab SpeciesEvolver instance.
        zone_function : function
            A function that return a mask of the species range.
        minimum_area : float, optional
            The minimum area of zones.
        neighborhood_structure : {'D8', 'D4'}, optional
            The structure describes how clusters are identified. The default,
            'D8' evaluates the eight neighboring nodes. The diagonal
            neighboring nodes are excluded when 'D4' is selected.
        initial_time : float, int, optional
            The initial time. The unit of time is unspecified within the
            controller. The default is 0.
        kwargs
            Keyword arguments for ``zone_function``. Note that the grid of
            `SpeciesEvolver` is automatically added to `kwargs`.
        """
        # Set parameters.

        self._grid = grid
        self._zone_func = zone_function
        self._zone_params = kwargs
        self._min_area = minimum_area
        self._neighborhood_struct = neighborhood_structure
        self._record = Record(initial_time)

        # Include `grid` in the zone params dictionary.

        self._zone_params['grid'] = self._grid

        # Set initial zones.

        initial_range = zone_function(**kwargs)
        self._zones = self._get_zones_with_mask(initial_range)

    @property
    def zones(self):
        """The zones of the SpeciesController."""
        return self._zones

    @property
    def record(self):
        """A DataFrame of SpeciesEvolver variables over time."""
        return self._record.dataframe

    def populate_zones_uniformly(self, count, species_type=ZoneSpecies,
                                 **kwargs):
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

        determines the connectivity of zones between `time` and the
                prior time step

        Parameters
        ----------
        dt : float
            The model time step duration. The first time step begins at 0.
            Following time steps are advanced by ``dt``.
        """
        self._record.iterate_time(dt)

        # Create an add on to insert into `record`.

        add_on = self._record.get_empty_add_on()

        # Resolve the spatiotemporal connectivity of the prior time step zones
        # to the new zones.

        prior_zones = self._zones

        zone_mask = self._zone_func(**self._zone_params)
        new_zones = self._get_zones_with_mask(zone_mask)

        self._zones = _update_zones(self._grid, prior_zones, new_zones, add_on)

        add_on['zone_count'] = len(self._zones)
        self._record.insert_add_on(add_on)

    def _get_zones_with_mask(self, range_mask):
        """Get zones using a range mask.

        Parameters
        ----------
        range_mask : ndarray
            Grid number of nodes array of booleans where `True` values
            indicates a node within the species range.

        Returns
        -------
        list of Zones
            The discrete zones identified in the range.
        """
        grid = self._grid

        # Label clusters of 'True' values in `range_mask`.

        if self._neighborhood_struct == 'D4':
            s = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]
        elif self._neighborhood_struct == 'D8':
            s = 3 * [[1, 1, 1]]
        cluster_array, cluster_count = label(range_mask.reshape(grid.shape),
                                             structure=s)

        # Create zones for clusters.

        zones = []

        for i in range(1, cluster_count + 1):
            mask = cluster_array == i

            if mask.sum() * grid.cellarea > self._min_area:
                # Create a zone if cluster area exceeds the minimum.
                zones.append(Zone(mask))

        return zones
