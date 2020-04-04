#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ZoneTaxon object of SpeciesEvolver."""
import numpy as np
from pandas import Series

from .base_taxon import Taxon


class ZoneTaxon(Taxon):
    """A taxon based in zones.

    A ``ZoneTaxon`` is composed of members of a lower taxonomic level that each
    exists within a ``Zone`` object. Taxonomic rank is not considered by this
    class despite the use of the term, 'speciation', which is used herein to
    generally describe creation of a child taxon object.

    All zones of the taxon can be obtained with the attribute, ``zones`` that
    are the objects that manage the geographic aspect of taxon member
    populations. The total geographic extent of all populations is depicted by
    the ``range_mask`` attribute. The zones of a ZoneTaxon instance are created
    and updated using a ``ZoneController``. At model time steps, the
    connectivity of zones over time is obtained using attributes of the
    ``Zone`` object.

    The evolution of this taxon type is carried out in two stages during a
    model time step. In the first stage, the zones of the taxon are updated
    as the result of zone connectivity between the prior and current step in
    the method, ``_update_zones``. This method is the primary implementation of
    taxon dispersal and it is called in a stage prior to other evolutionary
    processes so that all taxa are positioned in their landscape locations
    prior to the other processes.

    In the second stage, processes are carried out in methods that are readily
    expanded or overridden when needed. The primary methods of second stage
    macroevolution are ``_evaluate_dispersal``, ``_evaluate_speciation``, and
    ``_evaluate_extinction``. The evaluate dispersal method is intended to
    modify dispersal conducted in the first stage and it has no effect unless
    it is expanded or overridden to have an effect. Processes other than those
    listed above can be called by expanding or overridding the ``_evolve``
    method.

    The taxon is allopatric when it is associated with/exists within multiple
    zones (signifying multiple member populations). A timer is started when a
    taxon becomes allopatric. Allopatric speciation occurs once the timer
    reaches or exceeds the ``time_to_allopatric_speciation`` initialization
    parameter. If the initialization parameter, ``persists_post_speciation``
    is True (default), a child taxon is created in each zone except one zone
    (the largest by area) that becomes the sole zone of the taxon. If
    ``persists_post_speciation`` is set to False, a child taxon is created in
    each and every zone, and the parent no longer occupies any zones, and
    therefore the parent taxon is no longer extant.

    Extinction occurs when the taxon is no longer associated with any zones.
    This occurs when zones in the prior time step do not overlap zones in the
    current time step, signifying the geographic range of the taxon is no more.
    A taxon can become no longer extant also when the taxon speciates and
    ``persists_post_speciation`` is False signifying that the parent taxon
    has evolved into multiple taxon distinct from the original taxon.

    The following columns will be added to the ``record_data_frame`` of the
    SpeciesEvolver instance that tracks objects of this Taxon: 'speciations'
    and 'extinctions', which are the counts of these variables at time steps.
    Another column, 'pseudoextinctions' will be included when
    ``persists_post_speciation`` is False. This variable is the count of
    occurrences when a parent taxon became non-extant due to speciation and not
    because of an absence of zones.
    """

    def __init__(
        self,
        zones,
        parent=None,
        time_to_allopatric_speciation=0,
        persists_post_speciation=True,
    ):
        """Initialize a taxon.

        Parameters
        ----------
        zones : list of Zones
            The initial SpeciesEvolver Zones where the taxon is located.
        parent : Taxon, optional
            A SpeciesEvolver taxon that is the parent taxon. The default value,
            'None' indicates no parent.
        time_to_allopatric_speciation : float, int, optional
            The delay in model time to speciate following taxon geographic
            fragmentation, indicated by multiple objects in the attribute,
            ``zones``. Speciation occurs at the time step when the delay is
            reached or exceeded. The default value of 0 indicates speciation
            occurs at the same time step as geographic fragmentation.
        persists_post_speciation : boolean, optional
            When 'True', the default, taxon persists despite speciation. When
            'False' and following speciation, the taxon is no longer extant.
        """
        super().__init__()

        self.parent = parent
        self._tas = time_to_allopatric_speciation
        self._pps = persists_post_speciation

        # Store zones that each represent an instance of a narrower taxonomic
        # level, e.g. a population.

        self._zones = zones

        # Set taxon time in allopatry.

        self._time_in_allopatry = None
        self._update_allopatry_state()

    @Taxon.extant.setter
    def extant(self, state):
        """Set the living state of the taxon."""
        self._extant = state

        if not state:
            # Ensure the taxon is not associated with zones when it is no
            # longer extant.
            self._zones = []

    @property
    def range_mask(self):
        """A mask representing the geographic extent of the taxon.

        The mask is an array with a length of grid number of nodes. The taxon
        exists at nodes where mask elements are ``True``. The mask of a
        ZoneTaxon object is the union of all of its zone masks.
        """
        masks = [zone.mask for zone in self.zones]
        mask = np.any(masks, 0)

        return mask

    @property
    def zones(self):
        """The zones of the taxon."""
        return self._zones

    def _evolve(self, dt, stage, record):
        """Run a step of evolutionary processes.

        Dispersal resolved during stage 1 can be modified by extending or
        overriding the dispersal evaluation method, as can speciation and
        extinction that are also evaluated in this stage.

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
        list of Taxon
            The children produced by the taxon at a given stage. The ``evolve``
            method of child taxon will be called in stages following the stage
            the child taxon was produced. An empty list indicates no child
            taxon.
        """
        if stage == 0:
            self._update_zones()
            child_taxa = []

        elif stage == 1:
            # Evaluate macroevolutionary processes now that zones of all taxon
            # objects were updated in stage 1.

            self._evaluate_dispersal(dt)
            self._update_allopatry_state(dt)

            child_taxa = self._evaluate_speciation(dt)
            child_count = len(child_taxa)

            extinct = self._evaluate_extinction(dt)
            pseudoextinct = child_count > 1 and not self._pps
            self.extant = not extinct and not pseudoextinct

            # Update the record.

            record.increment_value("speciations", child_count)
            record.increment_value("extinctions", int(extinct and not pseudoextinct))

            if not self._pps:
                record.increment_value("pseudoextinctions", int(pseudoextinct))

        return stage < 1, child_taxa

    def _update_zones(self):
        """Update the zones of the taxon.

        Dispersal is represented by setting taxon zones to the zones of the
        current time step that overlap the taxon zones of the prior time step
        (`successors of a zone`).
        """
        successors = []

        for zone in self._zones:
            successors.extend(zone.successors)

        self._zones = Series(successors).drop_duplicates().tolist()

    def _update_allopatry_state(self, dt=None):
        """Update taxon time in allopatry.

        Parameter, ``dt`` can optionally be set to increment time in allopatry
        given that the taxon is already allopatric.

        Parameters
        ----------
        dt : float, int, optional
            The model time step duration.
        """
        if len(self.zones) < 2:
            self._time_in_allopatry = None
        elif self._time_in_allopatry is None:
            self._time_in_allopatry = 0
        elif dt is not None:
            self._time_in_allopatry += dt

    def _produce_child_taxon(self, zones):
        """Get the taxon resulting from speciation.

        This method returns a taxon of the same type as the object that
        speciates. This child taxon is initialized with the initialization
        parameter values of the parent object. At minimum in derived
        implementations of this method, the parent taxon of the child should be
        set to `self` to correctly construct the lineage.

        Parameters
        ----------
        zones : list of Zones
            The zones where the child taxon will exist.

        Returns
        -------
        Taxon
            The child taxon.
        """
        taxon_type = type(self)

        child_taxon = taxon_type(
            zones,
            parent=self,
            time_to_allopatric_speciation=self._tas,
            persists_post_speciation=self._pps,
        )

        return child_taxon

    def _evaluate_allopatric_speciation(self, dt):
        """Return child taxa if the taxon is allopatric.

        A child taxon is returned for each zone except the largest zone if
        ``persists_post_speciation`` is True. A child taxon is returned for all
        zones if ``persists_post_speciation`` is False. The ``zones`` attribute
        is set to a list containing the largest zone or no zone when
        ``persists_post_speciation`` is True or False, respectively. This
        method is called by the ``_evaluate_speciation`` method in the second
        stage of taxon evolution.

        Parameter ``dt`` is not used in the ZoneTaxon implementation of this
        method. It is included to follow the pattern of including this
        parameter in ZoneTaxon evolution methods.

        Parameters
        ----------
        dt : float, int
            The model time step duration.

        Returns
        -------
        list of taxon objects
            The taxon objects produced by allopatric speciation. An empty list
            indicates no child objects and no allopatric speciation.
        """
        zones = self.zones
        allopatric = self._time_in_allopatry is not None
        children = []

        if allopatric and self._time_in_allopatry >= self._tas:
            if self._pps:
                # The zone/member with the greatest zone area remains
                # associated with the object.
                idx = np.argmax([zone.area for zone in zones])
                largest_zone = zones.pop(idx)

            for zone in zones:
                child = self._produce_child_taxon([zone])
                children.append(child)

            if self._pps:
                self._zones = [largest_zone]
            else:
                self._zones = []

            self._update_allopatry_state()

        return children

    def _evaluate_dispersal(self, dt):
        """Modify taxon dispersal.

        Population dispersal is principally determined in stage 1 by the
        method, ``_update_zones``. This evaluation method is called by the
        taxon evolve method in stage 2 and allows modification of the stage 1
        dispersal. This method implemented in ZoneTaxon does not modify stage 1
        dispersal. It is intended to be overridden when needed.

        Parameters
        ----------
        dt : float, int
            The model time step duration.
        """
        # pragma: no cover

    def _evaluate_speciation(self, dt):
        """Return child taxa if speciation occurs.

        This method is called by the taxon stage 2 evolve method. The default
        implementation of this method solely gets any taxon objects resulting
        from the method, ``_evaluate_allopatric_speciation``. Other modes of
        speciation, including sympatric, can be evaluated here by expanded this
        ``_evaluate_speciation`` method.

        Parameters
        ----------
        dt : float, int
            The model time step duration.

        Returns
        -------
        list of taxon objects
            The taxon objects produced by allopatric speciation. An empty list
            indicates no child objects and no allopatric speciation.
        """
        child_taxa = self._evaluate_allopatric_speciation(dt)

        return child_taxa

    def _evaluate_extinction(self, dt):
        """Determine if extinction occurs.

        Extinction occurs if no zone/member populations exist. Other conditions
        of extinction can be included by expanding or overridding this method.

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
        taxon_occupies_no_zones = len(self.zones) == 0

        return taxon_occupies_no_zones
