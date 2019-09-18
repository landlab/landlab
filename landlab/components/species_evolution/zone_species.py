#!/usr/bin/env python
"""ZoneSpecies object of SpeciesEvolver."""

import numpy as np

from landlab.components.species_evolution import _Species


class _Population(object):
    """A ZoneSpecies population of a zone."""

    def __init__(self, species, zone, time_to_allopatric_speciation=None):
        self._zone = zone
        self._zone._species.append(species)
        self._time_to_allopatric_speciation = time_to_allopatric_speciation

    def _remove_species(self, species):
        """Remove a species from a zone."""
        self._zone._species.remove(species)
        self._zone = None


class ZoneSpecies(_Species):
    """A species based in zones.

    A species can reside in one or multiple zones. The total extent of the
    species is depicted by the ``range_mask`` attribute. The evolution of this
    species type is carried out in the following three processes: dispersal,
    speciation, and extinction.

    Dispersal is dictated by the spatial connectivity of zones over two
    successive time steps.



    Speciation occurs when a species exists in multiple
    zones following a delay set by the `allopatric_wait_time` initialization
    parameter. A clock counts down for each zone of the species where the
    species dispersed beyond one zone. The clock starts when a zone in the
    earlier time step overlaps multiple zones in the later time step. The clock
    is tracked by `_zones['time_to_allopatric_speciation'].

    Extinction occurs when the species exists in no zones. Extinction also
    occurs if ``pseudoextinction`` is True and the
    species speciates.
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

        pops = []

        for zone in initial_zones:
            pops.append(_Population(self, zone))
            zone._species.append(self)

        self._populations = pops

        self._mask_len = len(initial_zones[0].mask)

    @property
    def range_mask(self):
        """A mask representing the geographic extent of the species.

        Returns
        -------
        ndarray
            An array with a length of grid number of nodes. The species exists
            at nodes where mask elements are ``True``.
        """
        template = np.zeros(self._mask_len, bool)
        if self._populations:
            masks = [pop._zone.mask for pop in self._populations]
            if len(masks) == 1:
                masks.append(template)
            return np.any(masks, 0)
        else:
            return template

    def _evolve(self, time, dt, record_add_on):
        """Run the evolutionary processes of the species.

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
        tuple of boolean and Species list
            Tuple of  A list of SpeciesEvolver species objects that are the child species
            that result from the macroevolutionary processes run. An empty list
            indicates no child species.
        """
        # Handle dispersal and update zone data structure.

        updated_pops = []

        for pop in self._populations:
            new_pops = self._evaluate_dispersal(pop, time)
            updated_pops.extend(new_pops)

        self._populations = updated_pops

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
                pop._remove_species(self)
            self._populations = []

        # Update record add on.

        record_add_on['speciation_count'] += child_count
        record_add_on['extinction_count'] += extinct
        record_add_on['pseudoextinction_count'] += int(pseudoextinct)

        return not extinct, child_species

    def _evaluate_dispersal(self, population, time):
        """Set the range of the species as it results from dispersal.

        The default implementation returns all of the potential zones.

        Parameters
        ----------
        prior_zone : Zone
            A zone of the species in the prior time step.
        potential_zones : list of Zones
            The  zones where the species can potentially disperse.

        Returns
        -------
        list of Zones
            The zones where the species will disperse.
        """
        destination_zones = population._zone.path[time][1]

        # Handle `population`.

        updated_populations = []

        if len(destination_zones) == 0:
            # Remove `population` from its zone because it now exists nowhere.
            population._remove_species(self)

        # Handle new populations.

        for zone in destination_zones:
            if population._zone == zone:
                # The zone and population remain unchanged.
                updated_populations.append(population)

                if self._pseudoextinction and len(destination_zones) > 1:
                    population._time_to_allopatric_speciation = self._allopatric_wait_time

            elif self not in zone._species:
                zone_pop = _Population(self, zone, self._allopatric_wait_time)
                updated_populations.append(zone_pop)

            else:
                Exception('Dispersal condition not determined.')

        return updated_populations

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
