from .leaf_retention import *

#Duration classes and selection method
class Duration(object):
    def __init__(self, **kwargs):
        pass

class Annual(Duration):
    def __init__(self):
        super().__init__()

    def enter_dormancy(self):
        print('I die at dormancy')
    
    def emerge(self,emerge_size=[0.01,0.1,0.5]):
        print('I start as a seedling between 100 and 200% of the minimum size')
        #emerge size should be the min size for annuals
        return self.initialize_biomass(emerge_size)
    
    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5]}):
        #This function provides the range of mass an initial annual plant can have
        #Since initalization of annuals is the same as emergence, use the same function
        init_size_min=sum(grow_params['plant_part_min'])
        init_size_max=2*init_size_min
        return init_size_min, init_size_max

class Perennial(Duration):
    def __init__(self, retention_val):
        self.retention=self.select_retention_class(retention_val)
    
    def select_retention_class(self, retention_val):
        retention={
            'evergreen':Evergreen,
            'deciduous':Deciduous
        }
        return retention[retention_val]

    def enter_dormancy(self):
        print('I move carbs among live parts around during dormancy')
    
    def emerge(self, emerge_min=[0.01,0.1,0.5]):
        print('I will transfer some stored carbohydrate to photosynthetic parts')

    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        #This function provides the range of mass an initial perennial plant can have
        init_min_mass=sum(grow_params['plant_part_min'])
        init_max_mass=sum(grow_params['plant_part_max'])
        return [init_min_mass, init_max_mass]