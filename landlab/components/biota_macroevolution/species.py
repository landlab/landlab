"""TODO: Description.
"""

from copy import deepcopy
from landlab import FieldError
import numpy as np
from uuid import uuid4
from watershed import get_watershed_array

class Species(object):
    
    def __init__(self, timestep, nodes, parent_species_id=-1):
        """Initialize Species.
        """
        
        if timestep == 0:
            self.timesteps_existed = []
        else:
            self.timesteps_existed = [timestep]
        self.nodes = nodes
        self.parent_species_id = parent_species_id
        self.identifier = uuid4()
    
    def evolve(self, grid, timestep, min_stream_drainage_area, extinction_area_threshold):
        
        if 'drainage_area' not in grid.at_node.keys():
            raise FieldError('A grid "drainage_area" field is required to '
                             'update a species range.')
        
        if 'flow__receiver_node' not in grid.at_node.keys():
            raise FieldError('A grid "flow__receiver_node" field is required '
                             'to update a species range.')
        
        A = grid.at_node['drainage_area']
        stream_nodes = A >= min_stream_drainage_area
        receiver_at_node = grid.at_node['flow__receiver_node']
        
        # The while loops continues until all nodes are moved from
        # unresolved to updated. 
        unresolved_nodes = self.nodes
        updated_nodes = np.zeros(grid.number_of_nodes, dtype=bool)
        
        child_species = []
        
        barrier_created = False
        
        while any(unresolved_nodes):
        
            # Each loop processes the largest watershed of unresolved nodes.
            max_A = max(A[unresolved_nodes])
            max_A_array = np.all([A == max_A, unresolved_nodes], 0)
            
            # Find the outlet of the watershed that contains the species.
            prior_outlet = np.where(max_A_array)[0][0]
            giver_node = prior_outlet
            receiver_node = receiver_at_node[prior_outlet]

            if prior_outlet == receiver_node:
                updated_outlet = receiver_node
            else:
                outlet_not_found = True
                while outlet_not_found:
                    if np.any([grid.node_is_boundary(receiver_node),
                               receiver_node == receiver_at_node[receiver_node]]):
                        outlet_not_found = False
                        updated_outlet = giver_node
                    else:
                        giver_node = deepcopy(receiver_node)
                        receiver_node = receiver_at_node[receiver_node]

            new_range = get_watershed_array(grid, updated_outlet)
            new_range_streams = np.all([new_range, stream_nodes], 0)
            updated_nodes = np.any([updated_nodes, new_range_streams], 0)
            
            # Identify species nodes not yet updated.
            inverted_updated_nodes = np.invert(updated_nodes)
            unresolved_nodes = np.all([inverted_updated_nodes, self.nodes, stream_nodes], 0)
            
            # Evolve species
            # TODO: Barriers could be over identified for species with
            # disconnected regions.
            streams_elsewhere = any(A[unresolved_nodes] >= min_stream_drainage_area)
            
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
            
        captured_mask = np.all([np.invert(self.nodes), updated_nodes], 0)
        captured_nodes = np.where(captured_mask)[0]
        
        self.nodes = updated_nodes
        print(999,len(child_species))
        return child_species, captured_nodes
