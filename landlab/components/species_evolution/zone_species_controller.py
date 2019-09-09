#!/usr/bin/env python
"""Generic Species type of SpeciesEvolver."""
from scipy.ndimage.measurements import label

from landlab.components.species_evolution import (_SpeciesController,
                                                  ZoneSpecies, Zone)
import zone


class ZoneSpeciesController(_SpeciesController):
    """A controller for zone-based species.

    This class manages a collection of species objects introduced to
    SpeciesEvolver. It features 'zones' that are used to characterize the
    geographic aspect of species. The range of the species is defined by a
    function. The range is updated using this function at each time step. A
    zone is created for each discrete range portion.

    Speciation occurs when the species exists in more than one zone, although
    speciation is delayed until the allopatric wait time is exceeded in a zone.
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
            Keyword arguments for `zone_function`.
        """
        _SpeciesController.__init__(self, species_evolver)

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
        return self._species

    @property
    def extinct_species(self):
        return self._extinct_species

    @property
    def zones(self):
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

        range_mask[grid.node_is_boundary(grid.nodes.flatten())] = False

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

    def populate_each_zone(self, species_count, species_type=ZoneSpecies):
        """Populate each zone with a number of species.

        Parameters
        ----------
        species_count : int
            The count of species to populate to each zone.
        species_type :

        """
        for z in self._zones:
            species = [species_type([z]) for _ in range(species_count)]

            for s in species:
                self._delegate.introduce_species(s)

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
        boolean
            The value informs SpeciesEvolver if the species persists at `time`.
            True indicates that the species persists. False indicates the
            species has become extinct.
        Species list
            A list of SpeciesEvolver species objects that are the child species
            that result from the evolutionary processes run. An empty list
            indicates no child species.
        """
        # Resolve the spatiotemporal connectivity of the prior time step zones
        # to the new zones. Connectivity is described by `paths`.

        prior_zones = self._zones

        zone_mask = self._zone_func(**self._zone_params)
        new_zones = self._get_zones_with_mask(zone_mask, **self._zone_params)

        self._zones = zone._resolve_paths(self._grid, time, prior_zones,
                                          new_zones, record_add_on)

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
                surviving_species.extend(child_species)

            # Set id for child species.

            for cs in child_species:
                clade = cs.parent_species.clade
                cs._identifier = self._delegate._get_unused_species_id(clade)

        self._species = surviving_species

        return surviving_species
