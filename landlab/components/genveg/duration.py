from .leaf_retention import *
import numpy as np
#from landlab.components.genveg import PlantGrowth 

#Duration classes and selection method
class Duration(object):
    def __init__(self, **kwargs):
        pass

class Annual(Duration):
    def __init__(self):
        super().__init__()

    def enter_dormancy(self, grow_params = {'senescence_start':270,'growing_season_end':285}, current_day=0, plants=(np.recarray((0,),
            dtype=[('species','U10'),('pid',int),('cell_index',int),('root_biomass',float),('leaf_biomass',float),('stem_biomass',float)]))):
        #on or after senescence_day, the plant needs to lost 2% of its daily biomass after calculating new total biomass
        if current_day >= grow_params['senescence_start'] and current_day < grow_params['growing_season_end']:
            plants['root_biomass'] = plants['root_biomass'] - (plants['root_biomass']*0.02)
            plants['leaf_biomass'] = plants['leaf_biomass'] - (plants['leaf_biomass'] * 0.02)
            plants['stem_biomass'] = plants['stem_biomass'] - (plants['stem_biomass'] * 0.02)
        #on growing season end, the total biomass needs to be set to 0
        if current_day >= grow_params['growing_season_end']:
            plants['root_biomass'] = 0
            plants['leaf_biomass'] = 0
            plants['stem_biomass'] = 0 
        return plants
        print('I die at dormancy and my biomass is: {}'.format(total_biomass))
    
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