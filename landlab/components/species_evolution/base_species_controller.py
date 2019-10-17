#!/usr/bin/env python
from abc import ABC, abstractmethod


class _SpeciesController(ABC):
    """The base Species controller.

    A SpeciesController manages a collection of species introduced to
    SpeciesEvolver.
    """

    def __init__(self, species_evolver):
        """Instantiate a SpeciesController.

        Parameters
        ----------
        species_evolver: SpeciesEvolver
            An instance of the SpeciesEvolver component.
        """
        self._se = species_evolver
        species_evolver._species_controllers.append(self)

    @abstractmethod
    def _get_surviving_species(self):
        """Abstract method."""
        ...  # pragma: no cover
