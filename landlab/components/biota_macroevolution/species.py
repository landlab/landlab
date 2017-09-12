"""TODO: Description.
"""

from uuid import uuid4

class Species(object):
    
    def __init__(self, timestep, regions, parent_species_id=-1):
        """Initialize Species.
        """
        
        self.identifier = uuid4()
        
        self.phylogeny = {}  
        self.record_timestep(timestep, regions)
        
        self.parent_species_id = parent_species_id
    
    @property
    def phylogeny_timesteps(self):
        """Get the timesteps when the species existed."""
        return self.phylogeny.keys()
    
    def record_timestep(self, timestep, regions):
        self.phylogeny[timestep] = {'regions': regions}
        
    def exists_at_timestep(self, timestep):
        return timestep in self.phylogeny.keys()
    