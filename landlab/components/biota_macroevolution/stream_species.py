"""BiotaEvolver StreamSpecies object.
"""

from landlab import FieldError
import numpy as np
from .species import Species
from watershed import get_watershed_mask, get_watershed_outlet


class StreamSpecies(Species):

    def __init__(self, timestep, habitat_patches, parent_species_id=-1,
                 min_area_for_extinction=0):
        """Initialize a stream species.

        Parameters
        ----------
        timestep : float
            Initial timestep of species.
        habitat_patches : HabitatPatch list
            A list of BiotaEvolver HabitatPatch objects.
        parent_species_id : UUID
            The identifier of the parent species. An id of -1 indicates no
            parent species.
        extinction_area_threshold : float
            The species will become extinct at and below this area threshold.
        """

        super().__init__(timestep, habitat_patches, parent_species_id)

        self.min_area_for_extinction = min_area_for_extinction

    def disperse(self, grid, destination_patches):
        self.habitat_patches = destination_patches

    def speciate(self, destination_patches, timestep):
        new_species = []
        for d in destination_patches not in self.habitat_patches:
            # Bifurcate species.
            new = StreamSpecies(timestep, d, parent_species_id=self.identifier)
            new_species.append(new)

        return new_species
