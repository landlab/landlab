#!/usr/bin/env python
"""Generic ZoneSpecies object of SpeciesEvolver."""
class ZoneSpecies(object):
    """A generic zone-based species.

    The default implementation evolves species by methods

    `_evaluate_dispersal`, `_evaluate_speciation`, and `evaluate_extinction`.

    Dispersal is dictated by the spatial connectivity of zones over two
    successive time steps. Speciation occurs when a species exists in multiple
    zones following a delay set by the `allopatric_wait_time` initialization
    parameter. A clock counts down for each zone of the species where the
    species dispersed beyond one zone. The clock starts when a zone in the
    earlier time step overlaps multiple zones in the later time step. The clock
    is tracked by `_zones['time_to_allopatric_speciation`].


    Extinction occurs when the species exists in no
    zones. Extinction also occurs if *pseudoextinction* is `True` and the
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
        self._zones = dict.fromkeys(initial_zones, {})
        self._parent_species = parent_species
        self._allopatric_wait_time = allopatric_wait_time
        self._pseudoextinction = pseudoextinction
        self._identifier = None

    @property
    def identifier(self):
        """Get the species identifier.

        Returns
        -------
        tuple
            The identifier of the species. The first element is the clade of a
            species represented by a string. The second element is the specie
            number represented by an integer.
        """
        return self._identifier

    @property
    def parent_species(self):
        """The parent species.

        Returns
        -------
        species
            The species object that produced the species. A value of `None`
            indicates no parent species.
        """

        return self._parent_species

    @property
    def clade(self):
        """Get the species clade identifier.

        Returns
        -------
        string
            The clade identifier of the species.
        """
        return self._identifier[0]

    @property
    def zones(self):
        """The zones of the species.

        Returns
        -------
        list of Zones
            A zone object list of the species.
        """
        return list(self._zones.keys())

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

        updated_dict = {}

        for zone in self._zones:
            possible_destinations = zone.path[time][1]
            destinations = self._evaluate_dispersal(zone,
                                                    possible_destinations)

            zone_dict = self._create_data_structure_for_zones(zone,
                                                              destinations, dt)
            updated_dict = {**updated_dict, **zone_dict}

        self._zones = updated_dict

        # Handle speciation.

        child_species = []

        for zone in self.zones:
            speciates = self._evaluate_speciation(zone)
            if speciates:
                child_species.append(self._speciate([zone]))

        child_count = len(child_species)
        speciated = child_count > 0

        # Handle extinction.

        extinct = self._evaluate_extinction(speciated)
        pseudoextinct = extinct and speciated and len(self._zones) > 0

        # Update record add on.

        record_add_on['speciation_count'] += child_count
        record_add_on['extinction_count'] += extinct
        record_add_on['pseudoextinction_count'] += pseudoextinct

        if extinct and self._zones:
            self._zones = {}

        return not extinct, child_species

    def _create_data_structure_for_zones(self, prior_zone, destinations, dt):
        """Create a zone data structure for species destination zones.

        Parameters
        ----------
        new_zones : zone list
            The zones to insert into the species zone data structure.
        """
        new_dict = dict.fromkeys(destinations, {})

        key = 'time_to_allopatric_speciation'
        disperses_multiple_zones = len(destinations) > 1

        for zone in destinations:
            if zone in self._zones.keys():
                new_dict[zone] = {**new_dict[zone], **self._zones[zone]}

            # Update remaining allopatric wait time.

            splinter_zone = zone != prior_zone or self._pseudoextinction

            if disperses_multiple_zones and splinter_zone:
                new_dict[zone][key] = self._allopatric_wait_time
            elif key in new_dict.keys():
                new_dict[zone][key] -= dt

        return new_dict

    def _evaluate_dispersal(self, prior_zone, potential_zones):
        """Determine the zones where the species will disperse.

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
        destination_zones = potential_zones
        return destination_zones

    def _evaluate_speciation(self, zone):
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
        key = 'time_to_allopatric_speciation'
        speciation_triggered = key in self._zones[zone].keys()

        if speciation_triggered:
            remaining_time = self._zones[zone][key]
            time_is_reached = remaining_time - self._allopatric_wait_time <= 0
        else:
            time_is_reached = False

        return time_is_reached

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
        extant_in_no_zones = not self.zones

        pseudoextinct = self._pseudoextinction and speciation_occurred

        return pseudoextinct or extant_in_no_zones
