"""TODO: Description.
"""

from uuid import uuid4


class Species(object):
    
    def __init__(self, timestep, range_mask, parent_species_id=-1):
        """Initialize a species.
        """
        
        if timestep == 0:
            self.timesteps_existed = []
        else:
            self.timesteps_existed = [timestep]
        self.range_mask = range_mask
        self.parent_species_id = parent_species_id
        self.identifier = uuid4()
    
    def update_range(self, grid, timestep, extinction_area_threshold):
        
        child_species = []
        captured_nodes = []
        
        return child_species, captured_nodes

    def disperse(self, grid, stream_min_area, destination_patches):
        stream_nodes = grid.at_node['stream']
        
        for patch in destination_patches:
            
        
        nodes = np.all([self.range_mask, stream_nodes], 0)
        
        self.range_mask = nodes