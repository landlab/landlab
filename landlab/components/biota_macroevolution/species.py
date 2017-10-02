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
    
    def update_range(self, grid, min_stream_drainage_area):
        
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
            
        self.nodes = updated_nodes