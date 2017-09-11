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
        
    def record_timestep(self, timestep, regions):
        self.phylogeny[timestep] = {'regions': regions}
    
    def increment_timestep(self, timestep):
        if timestep < 1:
            previous_timestep = 0
        else:
            previous_timestep = timestep - 1
            
        if previous_timestep in self.phylogeny.keys():
            self.record_timestep(timestep,
                                 self.phylogeny[previous_timestep]['regions'])
        
    def exists_at_timestep(self, timestep):
        return timestep in self.phylogeny.keys()