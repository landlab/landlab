"""TODO: Description.
"""

from uuid import uuid4


class Habitat(object):
    
    def __init__(self, time, nodes):
        """
        """
        
        self.at_time = {time: [nodes]}
        self.identifier = uuid4()
        self.species = []
        
        self.patches = {'id': uuid4(), 'nodes_at_time': {}}
        