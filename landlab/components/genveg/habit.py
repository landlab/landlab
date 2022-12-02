from .duration import *
#Growth habit classes and selection method
#Growth habit uses duration properties to assign dormancy and emergence methods
class Habit(object):
    def __init__(self,duration_val, retention_val):
        self.duration=self.select_duration_class(duration_val, retention_val)

    def select_duration_class(self, duration_val, retention_val):
        duration={
            'annual':Annual,
            'perennial':Perennial(retention_val)
        }
        return duration[duration_val]   

class Forbherb(Habit):
    def __init__(self, duration_val, retention_val='None'):
        super().__init__(duration_val, retention_val)
        self.green_parts=('leaf', 'stem')

    def enter_dormancy(self):
        self.duration.enter_dormancy()

    def emerge(self):
        # Use this to move carbohydrate to aboveground biomass from storage organs and roots
        self.duration.emerge()

    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        self.duration.initialize_biomass(grow_params)

class Graminoid(Habit):
    def __init__(self, duration_val, retention_val='None'):
        super().__init__(duration_val,retention_val)

        self.green_parts=('leaf','stem')

    def enter_dormancy(self, grow_params = {'senescence_start':270,'growing_season_end':285}, current_day=0, plants=(np.recarray((0,),
            dtype=[('species','U10'),('pid',int),('cell_index',int),('root_biomass',float),('leaf_biomass',float),('stem_biomass',float)]))):
        self.duration.enter_dormancy(grow_params, current_day,plants)

    def emerge(self):
        self.duration.emerge()

    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        init_mass_min, init_mass_max=self.duration.initialize_biomass(grow_params)
        return init_mass_min, init_mass_max

class Shrub(Habit):
    def __init__(self, duration_val, retention_val='None'):
        super().__init__(duration_val,retention_val)
        
        self.green_parts=('leaf')

class Tree(Habit):
    def __init__(self, duration_val, retention_val):
        super().__init__(duration_val, retention_val='None')

        self.green_parts=('leaf')

class Vine(Habit):
    pass
