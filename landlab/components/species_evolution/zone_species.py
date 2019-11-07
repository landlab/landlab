#!/usr/bin/env python
"""ZoneSpecies object of SpeciesEvolver."""
import numpy as np

from .base_species import _Species


class Population(object):
    """A population of ZoneSpecies.

    Instances of this class are intended to be created and used directly by
    ZoneSpecies.
    """
    def __init__(self, species, zone, time_to_allopatric_speciation=None):
        """Initialize a population.

        Parameters
        ----------
        species : ZoneSpecies
            The population belongs to ``species``.
        zone : Zone
            The population inhabits ``zone``.
        time_to_allopatric_speciation : float, optional
            The time until speciation due to allopatry occurs. The default
            value of `None` indicates the population is not tending towards
            speciation.
        """
        self._species = species
        self.zone = zone
        self._time_to_allopatric_speciation = time_to_allopatric_speciation

    @property
    def zone(self):
        """The zone of the population."""
        return self._zone

    @zone.setter
    def zone(self, new_zone):
        """Set the population zone and update the species of the zone."""
        if new_zone == None:
            self._zone.species.discard(self._species)
        elif new_zone != None:
            new_zone.species.add(self._species)

        self._zone = new_zone


class ZoneSpecies(_Species):
    """A species based in zones.

    This implementation of SpeciesEvolver species is associated with ``Zone``
    objects. A species can reside in one or multiple zones. Internally, a
    ZoneSpecies is made up of Population objects where each population
    represents the species in a zone. The total extent of all populations of
    the species is depicted by the ``range_mask`` attribute.

    The evolution of this species type is carried out by the following methods,
    in the following order: ``_evaluate_dispersal``, ``_evaluate_speciation``,
    and ``_evaluate_extinction``. Dispersal and speciation are evaluated for
    each population. Extinction is evaluated for the entire species. These
    methods are intended to be expanded or overwritten.

    Dispersal is dictated by the spatial connectivity of zones over two
    successive time steps resolved by zones. Population zones are updated prior
    to ``evaluate_dispersal`` in the ``_update_populations`` method. However,
    the zones thus dispersal can be modified in ``evaluate_dispersal`` to
    modify the species spatial distribution prior to speciation and extinction.

    Speciation occurs when a species exists in multiple zones following a delay
    set by the `allopatric_wait_time` initialization parameter. A clock counts
    down for each zone of the species where the species dispersed beyond one
    zone. The clock starts when a zone in the earlier time step overlaps
    multiple zones in the later time step. The clock is tracked by
    `_zones['time_to_allopatric_speciation'].

    Extinction occurs when the species exists in no zones. Extinction also
    occurs if ``pseudoextinction`` is True and the species speciates.
    """
    def __init__(self, initial_zones, parent_species=None,
                 allopatric_wait_time=0, pseudoextinction=True):
        """Initialize a species.

        Parameters
        ----------
        initial_zones : list of Zones
            A list of SpeciesEvolver Zones where the species is located.
        parent_species : ZoneSpecies, optional
            A SpeciesEvolver species that is the parent species. The default
            value, 'None' indicates no parent species.
        allopatric_wait_time : float, optional
            The delay in model time between geographic seperation and
            speciation. Speciation occurs at the first time step when the delay
            is exceeded. The default value of 0 indicates speciation occurs at
            the same time step when geographic serperation occurs.
        pseudoextinction : boolean, optional
            When 'True', species become extinct when it speciates into child
            species.
        """
        super(ZoneSpecies, self).__init__()

        self._parent_species = parent_species

        self._allopatric_wait_time = allopatric_wait_time
        self._pseudoextinction = pseudoextinction

        # Set initial populations.

        self._populations = []

        for zone in initial_zones:
            self._populations.append(Population(self, zone))

        self._mask_len = len(initial_zones[0].mask)

    @property
    def range_mask(self):
        """A mask representing the geographic extent of the species.

        The mask is an array with a length of grid number of nodes. The species
        exists at nodes where mask elements are ``True``.
        """
        template = np.zeros(self._mask_len, bool)
        if self._populations:
            masks = [pop.zone.mask for pop in self._populations]
            if len(masks) == 1:
                masks.append(template)
            return np.any(masks, 0)
        else:
            return template

    def _evolve_stage_1(self, dt, record_add_on):
        """Run evolutionary processes in preperation of stage 2.

        The species populations are updated in this stage so that the zones are
        updated for all zones at the current time prior to stage 2.

        Parameters
        ----------
        dt : float
            The model time step duration.
        record_add_on : defaultdict
            A dictionary to pass values to the SpeciesEvolver record.
        """
        updated_pops = []

        for pop in self._populations:
            updated_pops.extend(self._disperse_population(pop))

        self._populations = updated_pops

    def _evolve_stage_2(self, dt, record_add_on):
        """Complete evolutionary processes for the time.

        In this stage, dispersal resolved during stage 1 can be modified by
        extending or overriding the dispersal evaluation method, as can
        speciation and extinction that are also evaluated in this stage.

        Parameters
        ----------
        dt : float
            The model time step duration.
        record_add_on : defaultdict
            A dictionary to pass values to the SpeciesEvolver record.

        Returns
        -------
        tuple of boolean and list of Species
            The first element of the returned tuple is boolean indicating if
            the species persists at the current time step. The second element
            of the returned tuple is a list of the child species. An empty list
            indicates no child species.
        """
        # Evaluate dispersal by population.

        for pop in self._populations:
            self._evaluate_dispersal(pop)

        # Handle speciation.

        child_species = []

        for pop in self._populations:
            speciates = self._evaluate_speciation(dt, pop)
            if speciates:
                child_species.append(self._speciate([pop._zone]))

        child_count = len(child_species)
        speciated = child_count > 0

        # Handle extinction.

        extinct = self._evaluate_extinction(speciated)
        pseudoextinct = extinct and speciated and len(self._populations) > 0

        if extinct and self._populations:
            for pop in self._populations:
                pop.zone = None
            self._populations = []

        # Update record add on.

        record_add_on['speciation_count'] += child_count
        record_add_on['extinction_count'] += extinct
        record_add_on['pseudoextinction_count'] += int(pseudoextinct)

        return not extinct, child_species

    def _disperse_population(self, population):
        """Update the population zone and create new zones if necessary.

        Handle species by the number of successors.

        Parameters
        ----------
        population : Population
            A population of the species.

        Returns
        -------
        list of Populations
            The populations that remain extant.
        """
        successors = population.zone.successors
        n_succ = len(successors)

        updated_populations = []

        if n_succ == 0:
            # Remove `population` from its zone because it now exists nowhere.
            population.zone = None

        elif n_succ == 1 and population.zone == successors[0]:
            # Population remains in the same zone.
            updated_populations = [population]

        elif n_succ == 1:
            # Population moved into a different zone.
            successor = successors[0]
            if self not in successor._species:
                # Add population to the successor zone only if the species does
                # not already inhabit this zone.
                population.zone = successor
                updated_populations = [population]

        else:
            # Population moves into multiple zones.
            for zone in successors:
                if population.zone == zone:
                    # The zone and population remain unchanged.
                    updated_populations.append(population)

                    if self._pseudoextinction:
                        population._time_to_allopatric_speciation = self._allopatric_wait_time

                elif self not in zone._species:
                    zone_pop = Population(self, zone,
                                          self._allopatric_wait_time)
                    updated_populations.append(zone_pop)

                else:
                    Exception('Dispersal condition not determined.')

        return updated_populations

    def _speciate(self, zones):
        """Get the species resulting from speciation in zones.

        The default implementation returns a species of the same type as
        `self`. This child species is initialized with the values of `self` for
        the default initialization parameters. At minimum in derived
        implementations, the parent species of the child should be set to
        `self` to correctly record the lineage.

        Parameters
        ----------
        zones : list of Zones
            The zones where the new species will exist.

        Returns
        -------
        Species
            The child species.
        """
        species_type = type(self)

        return species_type(zones, parent_species=self,
                            allopatric_wait_time=self._allopatric_wait_time,
                            pseudoextinction=self._pseudoextinction)

    def _evaluate_dispersal(self, population):
        """Set the range of the species as it results from dispersal.

        The default implementation does not modify population dispersal in
        stage 1.
        """
        ...  # pragma: no cover

    def _evaluate_speciation(self, dt, population):
        """Determine if speciation occurs.

        One condition is evaluated in the default implementation. `True` is
        returned if the time to allopatric species of `zone` has reached or
        fallen below 0.

        Parameters
        ----------
        zone : Zone
            The zone to evaluate.

        Returns
        -------
        boolean
            `True` indicates the species speciates. `False` indicates no
            speciation.
        """
        speciate = False

        speciation_triggered = population._time_to_allopatric_speciation is not None

        if speciation_triggered:
            speciate = population._time_to_allopatric_speciation <= 0

            if not speciate:
                population._time_to_allopatric_speciation -= dt

        return speciate

    def _evaluate_extinction(self, speciation_occurred):
        """Determine if extinction occurs.

        Parameters
        ----------
        speciation_occurred : boolean
            `True` indicates child species were produced.

        Returns
        -------
        boolean
            `True` indicates the species has become extinct. `False` indicates
            the species persists.
        """
        no_populations = len(self._populations) == 0

        pseudoextinct = self._pseudoextinction and speciation_occurred

        return pseudoextinct or no_populations
