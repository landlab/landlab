#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ZoneTaxon object of SpeciesEvolver."""
import numpy as np

from .base_taxon import Taxon
from .zone import Connection


class Population(object):
    """A population of ZoneTaxon.

    Instances of this class are intended to be created and used directly by
    ZoneTaxon.
    """

    def __init__(self, taxon, zone, time_to_allopatric_speciation=None):
        """Initialize a population.

        Parameters
        ----------
        taxon : ZoneTaxon
            The population belongs to ``taxon``.
        zone : Zone
            The population inhabits ``zone``.
        time_to_allopatric_speciation : float, optional
            The time until speciation due to allopatry. The default value of
            `None` indicates the population is not tending towards speciation.
        """
        self._taxon = taxon
        self.zone = zone
        self._time_to_allopatric_speciation = time_to_allopatric_speciation

    @property
    def zone(self):
        """The zone of the population."""
        return self._zone

    @zone.setter
    def zone(self, new_zone):
        """Set the population zone and update the taxon of the zone."""
        if new_zone is None:
            self._zone.taxa.remove(self._taxon)
        elif new_zone is not None:
            new_zone.taxa.append(self._taxon)

        self._zone = new_zone


class ZoneTaxon(Taxon):
    """A taxon based in a zone.

    A ``ZoneTaxon`` is composed of members of a lower taxonomic level. These
    members are stored in the ``populations`` attribute as ``Population``
    objects. Taxonomic rank is not considered in the class despite the use of
    the term, 'speciation', which is used generally here to indicate creation
    of a child taxon.

    Each population is associated with a ``Zone`` object that characterizes the
    geographic aspect of the population. The zone of a population, as well as
    all zones of a model, are created and updated using a ``ZoneController``.
    At model time steps, the connectivity of zones over time is obtained using
    the ``Zone`` object. The total geographic extent of all populations of the
    taxon is depicted by the ``range_mask`` attribute.

    The evolution of this taxon type is carried out in two stages during a
    model time step. Macroevolutionary processes are called during these
    stages. Dispersal and speciation are evaluated for each population.
    Extinction is evaluated for the taxon. These methods, further described
    below, are intended to be expanded or overridden when needed. Additional
    processes can be called by expanding or overridding the evolve method.

    The method, ``_evaluate_dispersal`` is called in the first stage. Dispersal
    is dictated by the spatial connectivity of the taxon zone at model time
    steps over the prior and current times. Dispersal is called in a stage
    prior to other processes so that all taxa are positioned in their landscape
    determined locations prior to other processes.

    The methods, ``_evaluate_speciation`` and ``_evaluate_extinction`` are
    called in the second stage. Speciation occurs when members of a broader
    taxon exist in multiple zones following a delay set by the
    ``allopatric_wait_time`` initialization parameter. A clock counts down for
    each zone of the taxon where the taxon dispersed beyond one zone. The clock
    starts when a zone in the earlier time step overlaps multiple zones in the
    later time step.

    Extinction occurs when the taxon no longer is associated with a zone. This
    can occur when the zone in the prior time step overlaps no zones in the
    current time step, signifying the habitat of the taxon disappeared and no
    suitable habitat was connected over time. The persistance of the taxon can
    also come to an end when the taxon speciates and ``pseudoextinction`` is
    True.
    """

    def __init__(
        self, zones, parent=None, allopatric_wait_time=0, pseudoextinction=True
    ):
        """Initialize a taxon.

        Parameters
        ----------
        zones : list of Zones
            The SpeciesEvolver Zones where the taxon is located.
        parent : ZoneTaxon, optional
            A SpeciesEvolver taxon that is the parent taxon. The default value,
            'None' indicates no parent.
        allopatric_wait_time : float, optional
            The delay in model time between geographic seperation and
            speciation. Speciation occurs at the first time step when the delay
            is reached or exceeded. The default value of 0 indicates speciation
            occurs at the same time step when geographic serperation occurs.
        pseudoextinction : boolean, optional
            When 'True', taxon becomes extinct when it produces child taxa.
        """
        super(ZoneTaxon, self).__init__()

        self.parent = parent
        self._allopatric_wait_time = allopatric_wait_time
        self._pseudoextinction = pseudoextinction

        # Set initial populations.

        self._populations = []

        for zone in zones:
            self._populations.append(Population(self, zone))

        self._mask_len = len(zones[0].mask)

    @Taxon.extant.setter
    def extant(self, value):
        """ """
        self._extant = value

        if not value:
            # Ensure the taxon is not associated with zones when not extant.
            for pop in self._populations:
                pop.zone = None
            self._populations = []

    @property
    def range_mask(self):
        """A mask representing the geographic extent of the taxon.

        The mask is an array with a length of grid number of nodes. The taxon
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

    def _evolve(self, dt, stage, record):
        """Run evolutionary processes for the time.

        In this stage, dispersal resolved during stage 1 can be modified by
        extending or overriding the dispersal evaluation method, as can
        speciation and extinction that are also evaluated in this stage.

        The attribute, ``extant`` is updated by this method.

        Parameters
        ----------
        dt : float
            The model time step duration.
        stage : int
            The evolution stage of the time step.
        record : defaultdict
            The SpeciesEvolver record.

        Returns
        -------
        boolean
           Indicates if the taxon is still evolving. When `False` is returned,
           this method will not be called for the taxon in subsequent stages in
           the current model time step.
        """
        if stage == 0:
            updated_pops = []

            for pop in self._populations:
                updated_pops.extend(self._disperse(pop))

            self._populations = updated_pops

        elif stage == 1:
            # Evaluate dispersal.

            for pop in self._populations:
                self._evaluate_dispersal(dt, pop)

            # Evaluate speciation.

            children = []

            for pop in self._populations:
                speciates = self._evaluate_speciation(dt, pop)
                if speciates:
                    children.append(self._produce_child_taxon([pop.zone]))

            child_count = len(children)
            speciated = child_count > 0

            # Evaluate extinction.

            pseudoextinct = self._pseudoextinction and speciated
            extinct = self._evaluate_extinction(dt) or pseudoextinct
            self.extant = not extinct

            # Update the record.

            record.increment_value("speciations", child_count)
            record.increment_value("extinctions", int(extinct))
            record.increment_value("pseudoextinctions", int(pseudoextinct))

        return stage < 1

    def _disperse(self, population):
        """Update the population zone and create new populations if warranted.

        Dispersal is handled by the connection type of the population zone.

        Parameters
        ----------
        population : Population
            A population of the taxon.
        """
        conn_type = population.zone._conn_type

        if conn_type == Connection.ONE_TO_NONE:
            # Remove the zone from the taxon. The taxon exist no where.
            population.zone = None
            updated_populations = []

        elif conn_type == Connection.ONE_TO_ONE:
            # Population remains in the same zone.
            updated_populations = [population]

        elif conn_type == Connection.MANY_TO_ONE:
            # Taxon moves into one zone.
            if population.zone == population.zone.successors[0]:
                # Population remains in the same zone.
                population._time_to_allopatric_speciation = None
                updated_populations = [population]

            elif self not in population.zone.successors[0].taxa:
                # Population moves into a zone not yet inhabited by taxon.
                population.zone = population.zone.successors[0]
                updated_populations = [population]

            else:
                # Population moves into a zone already inhabited by taxon. The
                # population conceptually merges by forgetting this population.
                updated_populations = []

        elif conn_type in [Connection.ONE_TO_MANY, Connection.MANY_TO_MANY]:
            # Taxon moves into multiple zones.

            updated_populations = []
            awt = self._allopatric_wait_time

            for zone in population.zone.successors:
                if population.zone == zone:
                    # Population remains in the same zone.
                    updated_populations.append(population)

                    if self._pseudoextinction:
                        population._time_to_allopatric_speciation = awt

                elif self not in zone.taxa:
                    # Create population for this zone new to the taxon.
                    zone_pop = Population(self, zone, awt)
                    updated_populations.append(zone_pop)

        return updated_populations

    def _produce_child_taxon(self, zone):
        """Get the taxon resulting from speciation in zones.

        The default implementation returns a taxon of the same type as
        `self`. This child taxon is initialized with the values of `self` for
        the default initialization parameters. At minimum in derived
        implementations, the parent taxon of the child should be set to
        `self` to correctly construct the lineage.

        Parameters
        ----------
        zones : list of Zones
            The zones where the new taxon will exist.

        Returns
        -------
        Taxon
            The child taxon.
        """
        taxon_type = type(self)

        child_taxon = taxon_type(
            zone,
            parent=self,
            allopatric_wait_time=self._allopatric_wait_time,
            pseudoextinction=self._pseudoextinction,
        )

        return child_taxon

    def _evaluate_allopatric_speciation(self, dt, population):
        """Determine if speciation occurs.

        Speciation is triggered during stage 1 dispersal. This speciation
        method decrements the time to allopatric speciation of the populations
        that speciation was triggered. Once the time to allopatric speciation
        is reached or exceeded, this method will return `True` signaling the
        population to speciate.

        Parameters
        ----------
        dt : float, int
            The time step duration to increment time to allopatric speciation,
            if speciation was previously triggered.
        population : Population
            A population of the taxon.

        Returns
        -------
        boolean
            `True` indicates the taxon speciates. `False` indicates no
            speciation.
        """
        triggered = population._time_to_allopatric_speciation is not None

        if triggered:
            # Speciation was previously triggered. Determine if allopatric
            # wait time was reached or exceeded.
            speciate = population._time_to_allopatric_speciation <= 0

            if not speciate:
                # Decrement the time to speciation.
                population._time_to_allopatric_speciation -= dt
        else:
            # Speciation was not previously triggered.
            speciate = False

        return speciate

    def _evaluate_dispersal(self, dt, population):
        """Set the range of the taxon as it results from dispersal.

        Population dispersal is principally determined in stage 1 by the
        method, ``_disperse``. This evaluation method is called by the taxon
        evolve method in stage 2 and allows modification of the stage 1
        dispersal. The default implementation of this method does not modify
        stage 1 dispersal. It is intended to be overridden when needed.

        Parameters
        ----------
        dt : float, int
            The model time step duration.
        population : Population
            A population of the taxon.
        """
        # pragma: no cover

    def _evaluate_speciation(self, dt, population):
        """Determine if speciation occurs.

        This method is called by the taxon stage 2 evolve method.

        Parameters
        ----------
        dt : float, int
            The model time step duration.
        population : Population
            A population of the taxon.

        Returns
        -------
        boolean
            `True` indicates the taxon produces a child taxon. `False`
            indicates no child taxon is produced.
        """
        allo_speciate = self._evaluate_allopatric_speciation(dt, population)

        return allo_speciate

    def _evaluate_extinction(self, dt):
        """Determine if extinction occurs.

        Extinction occurs if no populations exist.

        Parameters
        ----------
        dt : float, int
            The model time step duration.

        Returns
        -------
        boolean
            `True` indicates the taxon has become extinct. `False` indicates
            the taxon persists.
        """
        no_populations = len(self._populations) == 0

        return no_populations
