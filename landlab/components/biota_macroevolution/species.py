"""BiotaEvolver Species object.
"""

import numpy as np
from uuid import uuid4


class Species(object):

    def __init__(self, timestep, habitat_patches, parent_species_id=-1):
        """Initialize a species.

        Parameters
        ----------
        timestep : float
            Initial timestep of species.
        habitat_patches : HabitatPatch list
            A list of BiotaEvolver HabitatPatch objects.
        parent_species_id : UUID
            The identifier of the parent species. An id of -1 indicates no
            parent species.
        """

        self.timesteps_existed = [timestep]
        self.range_mask = []
        if isinstance(habitat_patches, list):
            self.habitat_patches = habitat_patches
        else:
            self.habitat_patches = [habitat_patches]

        self.parent_species_id = parent_species_id
        self.identifier = uuid4()

    def disperse(self, destination_patches, grid):
        self.habitat_patches = destination_patches

    def speciate(self):
        pass