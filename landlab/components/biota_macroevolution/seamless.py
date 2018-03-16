"""TODO: Description.

Component written by Nathan Lyons beginning August 20, 2017.
"""

from landlab import Component
from landlab.components.biota_macroevolution import plot
import numpy as np
from pprint import pprint
from uuid import UUID


class BiotaEvolver(Component):
    """TODO: Description.

    This component is adapted from SEAMLESS (Spatially Explicit Area Model of
    Landscape Evolution by SimulationS). See Albert et al., Systematic Biology 66, 2017.

    The primary method of this class is :func:`run_one_step`.

    Construction:

        BiotaEvolver(grid)

    Parameters
    ----------
    grid : ModelGrid
        A Landlab ModelGrid.
    """

    _name = 'BiotaEvolver'

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}

    def __init__(self, grid):
        """Initialize BiotaEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        """

        # Store grid and set parameters.
        self._grid = grid
        self.at_timestep = {'time': [], 'habitat_patches': []}
        self._timestep = -1
        self.species = []

    @property
    def species_ids(self):
        """Get the species ids."""
        ids = []
        for species in self.species:
            ids.append(species.identifier)
        return ids

    @property
    def time(self):
        return self.at_timestep['time'][self._timestep]

#    @property
#    def timesteps(self):
#        number_of_timesteps = len(self.at_timestep['time'])
#        return range(number_of_timesteps)

    # Update methods

    def run_one_step(self, time, new_habitat_patches, **kwargs):
        """Run the macroevolution processes for a single timestep.

        Parameters
        ----------
        time : float
            The current time in the simulation.
        new_habitat_patches : HabitatPatch list
            A list of BiotaEvolver HabitatPatch objects.
        """
        self._timestep = len(self.at_timestep['time']) + 1
        self.evolve(time, new_habitat_patches, **kwargs)

    def evolve(self, time, new_habitat_patches, **kwargs):
        """Run the macroevolution processes.

        Parameters
        ----------
        time : float
            The time in the simulation.
        new_habitat_patches : HabitatPatch list
            A list of BiotaEvolver HabitatPatch objects.
        """
        if self._timestep == 0:
            # Set non-default keys of at_timestep to null value.
            default_keys = ['time', 'habitat_patches']
            nondefault_keys = set(self.at_timestep.keys()) - set(default_keys)
            for key in nondefault_keys:
                self.at_timestep[key].append(-1)

            updated_habitat_patches = new_habitat_patches

        elif self._timestep > 0:
            vectors = self._get_habitat_patch_vectors(time,
                                                      new_habitat_patches,
                                                      **kwargs)

            self._update_species(vectors)

            # Flatten and get unique habitat patch vector destinations.
            destinations = [v.destinations for v in vectors]
            updated_habitat_patches = list(set(sum(destinations, [])))

        else:
            raise ValueError('BiotaEvolver timestep must be a positive '
                             'number.')

        self.at_timestep['habitat_patches'].append(updated_habitat_patches)
        self.at_timestep['time'].append(time)

    def _get_habitat_patch_vectors(self, time, new_habitat_patches, **kwargs):
        prior_patches = self.at_timestep['habitat_patches'][-1]
        patch_types = set([type(p) for p in new_habitat_patches])
        vectors = []

        for pt in patch_types:
            priors_with_type = list(filter(lambda p: isinstance(p, pt),
                                           prior_patches))
            new_with_type = list(filter(lambda p: isinstance(p, pt),
                                        new_habitat_patches))
            pt_vectors = pt._get_habitat_patch_vectors(priors_with_type,
                                                       new_with_type, time,
                                                       **kwargs)
            vectors.extend(pt_vectors)

        return vectors

    def _update_species(self, patch_vectors):
        # Update only the species extant in the prior time.
        prior_time = self.at_timestep['time'][-1]
        prior_species = self.species_at_time(prior_time)

        # Run macroevolution processes for each extant species.
        origins = list(filter(None, [v.origin for v in patch_vectors]))

        for species in prior_species:
            # Get vectors that include the origin of this species.
            indices = np.where(np.isin(origins, species.habitat_patches))[0]

            if indices:
                child_species = species.evolve(self.time, [patch_vectors[i] for i in indices])
                self.species.extend(child_species)

    # Species methods

    def add_new_species(self, species):
        self.species.append(species)

    def species_at_time(self, time):
        species = []
        for s in self.species:
            if time in s.times_existed:
                species.append(s)
        return species

    def species_at_node(self, node):
        species_list = []
        for species in self.species:
            if species.range_mask[node]:
                species_list.append(species)
        return species_list

    def species_with_id(self, identifier):
        _identifier = UUID(identifier)
        return self._object_with_identifier(self.species, _identifier)

    def get_tree(self, selected_species=None):

        if selected_species == None:
            selected_species = self.species

        tree = {}

        time_reversed = self.at_timestep['time'][::-1]

        for t in time_reversed:

            # Get species that existed in step.
            tree[t] = {}
            for s in selected_species:
                if t in s.times_existed:
                    pid = s.parent_species_id

                    if pid not in tree[t].keys():
                        tree[t][pid] = []
                    tree[t][pid].append(s)

        return tree

    # Print methods.

    def print_phylogeny(self):
        for species in self.species:
            print('\n')
            print(species.identifier)
            pprint(species.phylogeny)

    def print_tree(self, selected_species=None):
        tree = self.get_tree(selected_species)
        for time, parent_id in tree.items():
            print(time)
            for parent_id, species_list in tree[time].items():
                print('    ', parent_id)
                for species in species_list:
                    print('        ', species.identifier,
                          species.times_existed)

    # Convenience methods.

    def _unique(self, values):
        unique = np.unique(values)
        unique_no_nulls = self._remove_value(-1, unique)
        return unique_no_nulls

    def _remove_value(self, value, array):
        if any(value == array):
            indices_to_delete = np.where(array == value)
            array = np.delete(array, indices_to_delete)

        return array

    def _object_with_identifier(self, item_list, identifier):
        for item in item_list:
            if item.identifier == identifier:
                return item
                break
