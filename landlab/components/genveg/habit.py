from .duration import *
#Growth habit classes and selection method
#Growth habit uses duration properties to assign dormancy and emergence methods
class Habit(object):
    def __init__(self, duration_val, retention_val='deciduous'):
        self.duration=self.select_duration_class(duration_val, retention_val)

    def select_duration_class(self, duration_val, retention_val):
        duration={
            'annual':Annual(),
            'perennial':Perennial(retention_val)
        }
        return duration[duration_val]  

    def senesce(self, plants):
        plants=self.duration.senesce(plants)
        return plants
    
    def enter_dormancy(self, plants):
        plants=self.duration.enter_dormancy(plants)
        return plants
    
    def set_init_biomass_range(self, grow_params):
        init_mass_min, init_mass_max=self.duration.set_init_biomass_range(grow_params)
        return init_mass_min, init_mass_max


class Forbherb(Habit):
    def __init__(self, duration_val):
        super().__init__(duration_val)
        self.green_parts=('leaf', 'stem')

    def emerge(self):
        # Use this to move carbohydrate to aboveground biomass from storage organs and roots
        self.duration.emerge()


class Graminoid(Habit):
    def __init__(self, duration_val):
        super().__init__(duration_val)

        self.green_parts=('leaf','stem')

    def emerge(self):
        self.duration.emerge()

class Shrub(Habit):
    def __init__(self, duration_val, retention_val):
        super().__init__(duration_val, retention_val)
        
        self.green_parts=('leaf')

class Tree(Habit):
    def __init__(self, duration_val, retention_val):
        super().__init__(duration_val, retention_val)

        self.green_parts=('leaf')

class Vine(Habit):
    pass
