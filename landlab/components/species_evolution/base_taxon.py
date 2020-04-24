#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Base Taxon of SpeciesEvolver."""
from abc import ABC, abstractmethod


class Taxon(ABC):
    """Base Taxon of SpeciesEvolver.

    A SpeciesEvolver Taxon represents a group of organisms. Examples of groups
    represented by this object include an analogue group (e.g., seed plants in
    general, a specific trout species), a taxonomic level (e.g., phylum,
    species, population). More generally, Taxon subclasses act as approaches to
    simulate evolution of biologic groups. Taxon was designed to represent
    biologic groups, although they could be individual organisms.

    This class is intended to be subclassed. Subclasses must implement the
    attributes and methods of this base class that are designated as abstract.
    The methods must take the same parameters and both the attributes and
    methods must return the values described in their docstrings.
    """

    def __init__(self):
        """Instantiate a Taxon object.

        A base taxon object cannot be instantiated. This initialization method
        can be called by `super` in the initialization method of a subclass to
        set initial values of required private properties.
        """
        self._extant = True
        self._tid = None
        self._parent = None

    def __repr__(self):
        return "<{}, tid={}>".format(self.__class__.__name__, self._tid)

    @property
    def tid(self):
        """The identifier of the taxon.

        The identifier is an integer automatically and uniquely assigned by
        SpeciesEvolver once the component begins tracking the taxon. It is
        read-only as it should not be changed once this parameter is set.
        """
        return self._tid

    @property
    def extant(self):
        """The living state of the taxon.

        The taxon lives at the current model time if ``True``. The taxon is
        extinct as of the current model time if ``False``.
        """
        return self._extant

    @extant.setter
    def extant(self, state):
        """Set the living state of the taxon."""
        self._extant = state

    @property
    def parent(self):
        """The parent taxon.

        The parent is the taxon object that produced this object. A value of
        ``None`` indicates no parent taxon.
        """
        return self._parent

    @parent.setter
    def parent(self, taxon):
        """Set the parent taxon."""
        self._parent = taxon

    @property
    @abstractmethod
    def range_mask(self):
        """A mask of the taxon geographic extent.

        The range mask is a boolean numpy array where True values indicate
        where the taxon is located in the model grid associated with a
        SpeciesEvolver instance.

        This property must be implemented in a subclass.
        """
        # pragma: no cover

    @abstractmethod
    def _evolve(self, dt, stage, record):
        """Run the evolutionary processes of the taxon.

        SpeciesEvolver loops through the evolution processes of extant taxa in
        stages during the ``run_one_step`` method of the component. Therefore
        if a taxon type requires all other taxa to undergo some processing
        before an evolution process, then the taxon can evolve at a later
        stage using the ``stage`` parameter. Taxon subclasses should be
        designed with as few stages in this method as possible.

        This method at each stage must return both a boolean indicating if the
        taxon has additional stages to run and a list of child taxa produced
        during that evolution stage. The ``evolve`` method of child taxon will
        be called in stages subsequent to the stage the child taxon was
        produced. An empty list indicates no child taxon.

        See this method in ``ZoneTaxon`` for an example implementation.

        This method must be implemented in a subclass.

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
        # pragma: no cover
