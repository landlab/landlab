from .duration import *
#Growth habit classes and selection method
#Growth habit uses duration properties to assign dormancy and emergence methods
class Habit(object):
    def __init__(self, species_grow_params, green_parts, duration_val, retention_val):
        self.duration=self.select_duration_class(species_grow_params, green_parts, duration_val, retention_val)
        
    def select_duration_class(self, species_grow_params, green_parts, duration_val, retention_val):
        if duration_val=='perennial':
            duration_val=retention_val
        duration={
            'annual':Annual(species_grow_params),
            'deciduous':Deciduous(species_grow_params, green_parts),
            'evergreen':Evergreen(species_grow_params)
        }
        return duration[duration_val]

    def calc_lateral_width(self, volume, plants):
        plants['root_sys_width']=0.08+0.24*volume
        return plants
    
    def emerge(self, plants):
        plants=self.duration.emerge(plants)
        return plants

    def senesce(self, plants):
        plants=self.duration.senesce(plants)
        return plants
    
    def set_initial_biomass(self, plants, in_growing_season):
        plants=self.duration.set_initial_biomass(plants, in_growing_season)
        return plants
    
    def enter_dormancy(self, plants):
        plants=self.duration.enter_dormancy(plants)
        return plants
    
class Forbherb(Habit):
    def __init__(self, species_grow_params, duration_val):
        green_parts=('leaf', 'stem')
        retention_val='deciduous'
        super().__init__(species_grow_params, green_parts, duration_val, retention_val)
        
class Graminoid(Habit):
    def __init__(self, species_grow_params,duration_val):
        green_parts=('leaf', 'stem')
        retention_val='deciduous'
        super().__init__(species_grow_params, green_parts, duration_val, retention_val)

class Shrub(Habit):
    def __init__(self, species_grow_params, duration_val, retention_val):
        green_parts=('leaf')
        super().__init__(species_grow_params, green_parts, duration_val, retention_val)
    
    def calc_lateral_width(self, volume, plants):
        plants['lateral_width']=0.35+0.31*volume
        return plants
        
class Tree(Habit):
    def __init__(self, species_grow_params, duration_val, retention_val):
        green_parts=('leaf')
        super().__init__(species_grow_params, green_parts, duration_val, retention_val)

    def calc_lateral_width(self, volume, plants):
        plants['lateral_width']=0.35+0.31*volume
        return plants

class Vine(Habit):
    def __init__(self, species_grow_params, duration_val, retention_val):
        green_parts=('leaf')
        super().__init__(species_grow_params, green_parts, duration_val, retention_val)