"""BiotaEvolver Species object.
"""

from landlab.components.biota_macroevolution import HabitatPatchVector
from uuid import uuid4


class Species(object):

    def __init__(self, initial_time, habitat_patches, parent_species_id=-1):
        """Initialize a species.

        Parameters
        ----------
        initial_time : float
            Initial time of the species.
        habitat_patches : HabitatPatch list
            A list of BiotaEvolver HabitatPatch objects.
        parent_species_id : UUID
            The identifier of the parent species. An id of -1 indicates no
            parent species.
        """

        self.times_existed = [initial_time]
        self.range_mask = []
        if isinstance(habitat_patches, list):
            self.habitat_patches = habitat_patches
        else:
            self.habitat_patches = [habitat_patches]

        self.parent_species_id = parent_species_id
        self.identifier = uuid4()

    def evolve(self, time, patch_vectors):

        child_species = []
        print(555,len(patch_vectors))
        for v in patch_vectors:
            if v.cardinality in [HabitatPatchVector.ONE_TO_ONE,
                                 HabitatPatchVector.MANY_TO_ONE]:
                self.habitat_patches = v.destinations[0]
                self.times_existed.append(time)
            elif v.cardinality in [HabitatPatchVector.ONE_TO_MANY,
                                   HabitatPatchVector.MANY_TO_MANY]:
                child_species.append(Species(time, v.origin,parent_species_id=
                                             self.identifier))
                for d in v.destinations:
                    child_species.append(Species(time, v.origin,
                                                 parent_species_id=
                                                 self.identifier))

        return child_species

    def disperse(self, destination_patches, grid):
        self.habitat_patches = destination_patches

    def speciate(self):
        pass
