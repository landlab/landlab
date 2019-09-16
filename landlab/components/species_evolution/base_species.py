#!/usr/bin/env python
"""Base Species of SpeciesEvolver."""
from abc import ABC, abstractmethod


class _Species(ABC):
    """Base Species."""

    def __init__(self, parent_species=None):
        self._parent_species = parent_species
        self._identifier = None

    @property
    @abstractmethod
    def range_mask(self):
        """Abstract method."""
        ...  # pragma: no cover
