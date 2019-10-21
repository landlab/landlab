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
        self._grid = species_evolver.grid
        species_evolver._add_species_controller(self)

        # Reference SpeciesEvolver methods needed by this controller.
        self._introduce_species = species_evolver._introduce_species

    @abstractmethod
    def _get_surviving_species(self):
        """Abstract method."""
        ...  # pragma: no cover
