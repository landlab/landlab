"""TODO: Description.

Component written by Nathan Lyons beginning August 20, 2017.
"""

from landlab import Component
from landlab.components.biota_macroevolution import BiotaEvolverObject
from landlab.core.messages import warning_message
import numpy as np
from pprint import pprint
from string import ascii_uppercase
from uuid import UUID


class BiotaEvolver(Component, BiotaEvolverObject):
    """TODO: Description.

    This component is adapted from SEAMLESS (Spatially Explicit Area Model of
    Landscape Evolution by SimulationS).See Albert et al., Systematic Biology 66, 2017.

    The primary method of this class is :func:`run_one_step`.

    Construction:

        BiotaEvolver(grid)
    """

    _name = 'BiotaEvolver'

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}

    def __init__(self, grid, map_vars=None, **kwds):
        """Initialize BiotaEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        """
        Component.__init__(self, grid, map_vars=None, **kwds)
        BiotaEvolverObject.__init__(self)

        self._clades = {}

    @property
    def species_ids(self):
        """Get the species ids."""
        ids = []
        for species in self.species:
            ids.append(species.identifier)
        return ids

    # Update methods

    def run_one_step(self, time, new_habitat_patches, **kwargs):
        """Run the macroevolution processes for a single timestep.

        Parameters
        ----------
        time : float
            The time in the simulation.
        new_habitat_patches : HabitatPatch list
            A list of BiotaEvolver HabitatPatch objects.
        """
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
        # Time tuple: first element is prior time, second is the next time.
        times = (self.record.time__last, time)

        self.record.setdefault(times[1], {})

        if self.record.count == 0:
            warning_message('Species must be introduced.')
        else:
            vectors = self._get_habitat_patch_vectors(times,
                                                      new_habitat_patches,
                                                      **kwargs)

            time_species = self._get_advancing_species(times,
                                                       vectors)

            # Flatten and get unique habitat patch vector destinations.
            destinations = [v.destinations for v in vectors]
            updated_habitat_patches = list(set(sum(destinations, [])))

        self.record[times[1]].update({'habitat_patches':
            updated_habitat_patches, 'species': time_species})

    def _get_habitat_patch_vectors(self, times, new_habitat_patches, **kwargs):
        prior_patches = self.record[times[0]]['habitat_patches']
        patch_types = set([type(p) for p in new_habitat_patches])
        vectors = []

        for pt in patch_types:
            priors_with_type = list(filter(lambda p: isinstance(p, pt),
                                           prior_patches))
            new_with_type = list(filter(lambda p: isinstance(p, pt),
                                        new_habitat_patches))
            output = pt._get_habitat_patch_vectors(priors_with_type,
                                                   new_with_type, times[1],
                                                   self._grid, **kwargs)

            # Parse habitat patch vector output.
            if(isinstance(output, tuple)):
                pt_vectors = output[0]
                self.record[times[1]].update(output[1])
            else:
                pt_vectors = output

            vectors.extend(pt_vectors)

        return vectors

    def _get_advancing_species(self, times, patch_vectors):
        # Update only the extant species.
        extant_species = self.record[times[0]]['species']

        # Run macroevolution processes for each extant species.
        origins = list(filter(None, [v.origin for v in patch_vectors]))

        advancing_species = []

        for es in extant_species:
            # Get vectors that include the origin of this species.
            species_patches = es.record[es.record.time__last]['habitat_patches']
            indices = np.where(np.isin(origins, species_patches))[0]

            if len(indices) > 0:
                pv = [patch_vectors[i] for i in indices]
                child_species = es.run_macroevolution_processes(times[1], pv)

                # Set child species id.
                for cs in child_species:
                    sid = self._get_unused_species_id(es.clade)
                    cs._identifier = sid

                advancing_species.extend(child_species)

        return advancing_species

    # Species methods

    def introduce_species(self, species, time):
        """Add a species to BiotaEvolver.

        Parameters
        ----------
        species : BiotaEvolver Species
            The species to introduce.
        time : float
            The time in the simulation.
        """
        cid = self._get_unused_clade_id()
        sid = self._get_unused_species_id(cid)
        species._identifier = sid

        self.record.setdefault(time, {})

        self.record[time].setdefault('species', [])
        self.record[time]['species'].append(species)

        self.record[time].setdefault('habitat_patches', [])
        species_patches = species.record[time]['habitat_patches']
        self.record[time]['habitat_patches'].extend(species_patches)

    def species_at_node(self, node):
        species_list = []
        for species in self.species:
            if species.range_mask[node]:
                species_list.append(species)
        return species_list

    def species_with_id(self, identifier):
        _identifier = UUID(identifier)
        return self._object_with_identifier(self.species, _identifier)

    def get_tree(self):
        """Get phylogenetic tree.
        """
        tree = {}

        for t in self.record.times:
            tree[t] = {}
            for s in self.record[t]['species']:
                p = s.parent_species
                tree[t].setdefault(p, [])
                tree[t][p].append(s)

        return tree

    def print_tree(self):
        tree = self.get_tree()

        readable_tree = {}
        for time in tree.keys():
            for p in tree[time].keys():
                if p == -1:
                    pid = p
                else:
                    pid = p.identifier
                species_ids = [s.identifier for s in tree[time][p]]
                readable_tree.setdefault(time, {})
                readable_tree[time][pid] = species_ids

        pprint(readable_tree)

    # Convenience methods.

    def _object_with_identifier(self, item_list, identifier):
        for item in item_list:
            if item.identifier == identifier:
                return item
                break

    def _get_unused_clade_id(self):
        used_ids = list(self._clades.keys())

        alphabet = list(ascii_uppercase)
        clade_id = np.setdiff1d(alphabet, used_ids)

        duplicator = alphabet
        while len(clade_id) == 0:
           duplicator = np.core.defchararray.add(alphabet, duplicator)
           clade_id = np.setdiff1d(duplicator, used_ids)

        self._clades[clade_id[0]] = 0

        return clade_id[0]

    def _get_unused_species_id(self, clade_id):
        self._clades[clade_id] += 1
        return (clade_id, self._clades[clade_id])
