#!/usr/bin/env python
"""Zone-based SpeciesController of SpeciesEvolver."""
from scipy.ndimage.measurements import label

from landlab.components.species_evolution import (_SpeciesController,
                                                  ZoneSpecies, Zone)
from .zone import _update_zones


class ZoneSpeciesController(_SpeciesController):
    """A controller for zone-based species.

    This implementation of SpeciesController manages a collection of species
    using 'zones' that characterize the geographic aspect of species. A zone
    represents a portion of a model grid and it is made up of spatially
    continous nodes. The nodes of all zones of a certain type are identified
    using a function set with ``zone_function``. This function is called by
    SpeciesEvolver at each time step to identify all the zone nodes and then
    delineate individual zones.

    The structure of a model grid is diagramed below. The grid contains six
    columns and six rows. In this example, zones are designated where ``field``
    is greater than 0. Nodes marked with an ``*`` are nodes that meet the zone
    designation. All other nodes are marked with a ``·``. A zone is created for
    each cluster of designated nodes.

                   zone
    field          designation
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
    are not included. A minimum zone area can be enforced

    Individual zones are
    identified using the ``neighborhood_structure`` keyword.


    D8             D4             D8, min=2      D4, min=2
    · · · · · ·    · · · · · ·    · · · · · ·    · · · · · ·
    · · + + + ·    · · + + + ·    · · + + + ·    · · + + + ·
    · · + + · ·    · · + + · ·    · · + + · ·    · · + + · ·
    · · · · · ·    · · · · · ·    · · · · · ·    · · · · · ·
    · x x · o ·    · @ @ · o ·    · x x · · ·    · x x · · ·
    x · · · · ·    x · · · · ·    x · · · · ·    · · · · · ·

    Speciation occurs when the species exists in more than one zone, although
    speciation is delayed until the allopatric wait time is exceeded in a zone.


    The grid represents the time, ``T0``.

        T0
    · · · · · ·
    · · · · · ·
    · · x x · ·
    · · x x · ·
    · · · · · ·
    · · · · · ·

    Below are variations of the grid at time 1, ``T1`` in three variations where
    each contains one zone. In ``T1a``, ``T1b``, and ``T1c`` the zone stayed
    the stayed, moved and changed size, respectively. Species migrated with the
    zone. However, in ``T1c``, no nodes overlaps, therefore it is assumed species
    did not disperse from the zone in T0 to the zone in T1d.

    T1a            T1b            T1c            T1d
    · · · · · ·    · · · · · ·    · · · · · ·    · · · · · ·
    · · · · · ·    · · · · · ·    · · · * · *    · · · · · ·
    · · * * · ·    · · · · · ·    · · * * * *    · · · · · ·
    · · * * · ·    · * * · · ·    · · · · · ·    · · · · · ·
    · · · · · ·    · * * · · ·    · · · · · ·    * * · · · ·
    · · · · · ·    · · · · · ·    · · · · · ·    * * · · · ·

    Below are more ``T1`` variations, now demonstrating two zones.

    T1e
    · · · · · ·
    · · · · · ·
    · x · · x x
    x x · · x x
    x x · · x ·
    x · · · · ·


    """
    def __init__(self, species_evolver, zone_function, **kwargs):
        """Initialize a species.

        Parameters
        ----------
        species_evolver : SpeciesEvolver
            A Landlab SpeciesEvolver instance.
        range_function : function
            A function that return a mask of the species range.
        kwargs
            Keyword arguments for ``zone_function``.
        """
        super(ZoneSpeciesController, self).__init__(species_evolver)

        # Set parameters.

        self._zone_func = zone_function
        self._zone_params = kwargs
        self._species = []
        self._extinct_species = []

        # Set initial zones.

        initial_range = zone_function(**kwargs)
        self._zones = self._get_zones_with_mask(initial_range, **kwargs)

    @property
    def extant_species(self):
        """The species that exist at the current model time."""
        return self._species

    @property
    def extinct_species(self):
        """The species that no longer exist at the current model time."""
        return self._extinct_species

    @property
    def zones(self):
        """The zones of the SpeciesController."""
        return self._zones

    def _get_zones_with_mask(self, range_mask, minimum_area=0,
                             neighborhood_structure='D8', **kwargs):
        """Get zones using a range mask.

        Parameters
        ----------
        range_mask : ndarray
            Grid number of nodes array of booleans where `True` values
            indicates a node within the species range.
        minimum_area : float, optional
            The minimum area of zones to be returned.
        neighborhood_structure : {'D8', 'D4'}, optional
            The structure describes how clusters are identified. The default,
            'D8' evaluates the eight neighboring nodes. The diagonal
            neighboring nodes are excluded when 'D4' is selected.

        Returns
        -------
        list of Zones
            The discrete zones identified in the range.
        """
        grid = self._grid

        # Label clusters of 'True' values in `range_mask`.

        if neighborhood_structure == 'D4':
            s = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]
        elif neighborhood_structure == 'D8':
            s = 3 * [[1, 1, 1]]
        cluster_array, cluster_count = label(range_mask.reshape(grid.shape),
                                             structure=s)

        # Create zones for clusters.

        zones = []

        for i in range(1, cluster_count + 1):
            mask = cluster_array == i

            if mask.sum() * grid.cellarea > minimum_area:
                # Create a zone if cluster area exceeds the minimum.
                zones.append(Zone(mask))

        return zones

    def populate_zones_uniformly(self, species_count, species_type=ZoneSpecies,
                                 **kwargs):
        """Populate each zone with the same species count.

        Parameters
        ----------
        species_count : int
            The count of species to populate to each zone.
        species_type : type of Species
            A Species type that takes a list of Zones as its first parameter.
        kwargs : dictionary
            Keyword arguments of ``species_type``.
        """
        species = []

        for z in self._zones:
            species.extend([species_type([z], **kwargs) for _ in range(species_count)])

        self._delegate._introduce_species(species)

        self._species.extend(species)

    def _get_surviving_species(self, dt, time, record_add_on):
        """Run the evolutionary processes of the species.

        This method:
            -   determines the connectivity of zones between `time` and the
                prior time step
            -   disperses the species into the zones that are
                connected
            -   speciates following the speciation rules
            -   updates the record add on
            -   provides SpeciesEvolver with the child species

        Extinction effectively occurs when the species attribute, `zones`
        returns an empty list.

        Parameters
        ----------
        time : float
            The time in the simulation.
        dt : float
            The model time step duration.
        record_add_on : defaultdict
            A dictionary to pass values to the SpeciesEvolver record.

        Returns
        -------
        list of Species
            A list of SpeciesEvolver species objects that are the child species
            that result from the evolutionary processes run. An empty list
            indicates no child species.
        """
        # Resolve the spatiotemporal connectivity of the prior time step zones
        # to the new zones. Connectivity is described by `paths`.

        prior_zones = self._zones

        zone_mask = self._zone_func(**self._zone_params)
        new_zones = self._get_zones_with_mask(zone_mask, **self._zone_params)

        self._zones = _update_zones(self._grid, time, prior_zones, new_zones,
                                    record_add_on)

        # Evolve species.

        surviving_species = []

        for es in self.extant_species:
            species_persists, child_species = es._evolve(time, dt,
                                                         record_add_on)

            if species_persists:
                surviving_species.append(es)
            else:
                self._extinct_species.append(es)

            if len(child_species) > 0:
                self._delegate._introduce_species(child_species)
                surviving_species.extend(child_species)

        self._species = surviving_species

        return surviving_species
