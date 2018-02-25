"""TODO: Description.
"""

from landlab import FieldError
import numpy as np
from .species import Species
from watershed import get_watershed_mask, get_watershed_outlet


class StreamSpecies(Species):
    
    def __init__(self, timestep, range_mask, parent_species_id=-1):
        """Initialize a stream species.
        """
        
        super().__init__(timestep, range_mask, parent_species_id)
    
    def disperse(self, grid, destination_patches):        
        stream_nodes = grid.at_node['stream']
        
        for patch in destination_patches:
            
        
        nodes = np.all([self.range_mask, stream_nodes], 0)
        
        self.habitat_patches = patches
        
        
    
    def update_range(self, grid, timestep, extinction_area_threshold):
        
        if 'drainage_area' not in grid.at_node.keys():
            raise FieldError('A grid "drainage_area" field is required to '
                             'update a stream species range.')
        if 'stream' not in grid.at_node.keys():
            raise FieldError('A grid "stream" field is required to '
                             'update a stream species range.')
        
        A = grid.at_node['drainage_area']
        stream_nodes = grid.at_node['stream']
        
        # The while loops continues until all nodes are moved from
        # unresolved to updated. 
        unresolved_nodes = self.range_mask
        updated_nodes = np.zeros(grid.number_of_nodes, dtype=bool)
        
        child_species = []
        
        barrier_created = False
        
        captured_area = 0
        
        while any(unresolved_nodes):
        
            # Each loop processes the largest watershed of unresolved nodes.
            max_A = max(A[unresolved_nodes])
            max_A_array = np.all([A == max_A, unresolved_nodes], 0)
            
            # Find the outlet of the watershed that contains the species.
            prior_outlet = np.where(max_A_array)[0][0]
            updated_outlet = get_watershed_outlet(grid, prior_outlet)

            new_range = get_watershed_mask(grid, updated_outlet)
            new_range_streams = np.all([new_range, stream_nodes], 0)
            updated_nodes = np.any([updated_nodes, new_range_streams], 0)
            
            updated_max_A = max(A[new_range_streams])
            
            captured_area += updated_max_A - max_A
            
            # Identify species nodes not yet updated.
            inverted_updated_nodes = np.invert(updated_nodes)
            unresolved_nodes = np.all([inverted_updated_nodes, self.range_mask,
                                       stream_nodes], 0)

            # Evolve species
            # TODO: Barriers could be over identified for species with
            # disconnected regions.
            streams_elsewhere = any(stream_nodes[unresolved_nodes])
            
            if barrier_created or streams_elsewhere:
                # Barrier created: create new child species, this species does
                # not persist.
                new_species = Species(timestep, new_range_streams,
                                      self.identifier)
                child_species.append(new_species)
                barrier_created = True
            elif max_A > extinction_area_threshold:
                # No changes: species persists. Otherwise, species becomes
                # extinct.
                self.timesteps_existed.append(timestep)
            
        captured_mask = np.all([np.invert(self.range_mask), updated_nodes], 0)
        captured_nodes = np.where(captured_mask)[0]
        
        self.range_mask = updated_nodes

        return child_species, captured_area
