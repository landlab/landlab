#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Base Taxon of SpeciesEvolver."""
from abc import ABC, abstractmethod


class Taxon(ABC):
    """Base Taxon of SpeciesEvolver.

    A SpeciesEvolver Taxon represents a group of organisms. Examples of the
    kind of group represented by a subclass include an analog group (e.g., a
    specific trout species, seed plants), a taxonomic level (e.g., phylum,
    species, population). More generally, Taxon subclasses act as approaches to
    simulate evolution of biologic groups. Taxon can be made of individual
    organisms, although SpeciesEvolver currently has no built in functionality
    for individuals.

    This class is intended to be subclassed. Subclasses must implement the
    properties and methods of this base class that are designated as abstract.
    The methods must take the same parameters and both the properties and
    methods must return the values described in their docstrings.
    """

    def __init__(self):
        self._extant = True
        self._uid = None
        self._parent = None
        self._children = []

    def __repr__(self):
        return "<{}, uid={}>".format(self.__class__.__name__, self.uid)

    @property
    def uid(self):
        """The unique identifier of the taxon.

        The identifier is a unique integer automatically assigned by
        SpeciesEvolver once the component begins tracking the taxon. It is
        read-only as it should not be changed once it is assigned.
        """
        return self._uid

    @property
    def extant(self):
        """The living state of the taxon.

        The taxon lives at the current model time if ``True``. The taxon is
        extinct as of the current model time if ``False``.
        """
        return self._extant

    @extant.setter
    def extant(self, value):
        """Set the living state of the taxon."""
        self._extant = value

    @property
    def parent(self):
        """The parent taxon.

        The parent is the taxon object that produced this object. A value of
        ``None`` indicates no parent taxon.
        """
        return self._parent

    @parent.setter
    def parent(self, value):
        """Set the parent taxon.

        This method also appends this taxon to the list of children belonging
        to the parent taxon.
        """
        self._parent = value

        if value is not None:
            value.children.append(self)

    @property
    def children(self):
        """The immediate descendents of the taxon.

        The children are the objects produced by this taxon. An empty list
        indicates no child taxon.
        """
        return self._children

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
    def _evolve(self, dt, stage, record, id):
        """Run the evolutionary processes of the taxon.

        SpeciesEvolver loops through the evolution processes of extant taxa in
        stages during the ``run_one_step`` method of the component. Therefore
        if a taxon type requires all other taxa to undergo some processing
        before an evolution process, then the taxon can evolve at a later
        stage using the ``stage`` parameter. This method must return a boolean
        indicating if the taxon has completed evolution of all stages. See this
        method implemented in ``ZoneTaxon`` for an example.

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
        """
        # pragma: no cover
