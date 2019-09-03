#!/usr/bin/env python

from abc import ABC, abstractmethod

from landlab.components.species_evolution.species_evolution import SpeciesDataMixIn

class _SpeciesController(SpeciesDataMixIn, ABC):
    """The base Species controller.

    All species classes must include these class attributes: `identifier`,
    `clade`, `parent_species`, and `extant`; and these methods: `_evolve`.
    """

    def __init__(self, species_evolver):
        # Set parameters.

        self._delegate = species_evolver._delegate
        species_evolver._species_controllers.append(self)
        self._grid = species_evolver._grid

#    @abstractmethod
#    def _evolve(self, time, dt, record_add_on):
#        """Abstract method."""
#        ...  # pragma: no cover
