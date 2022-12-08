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

class Forbherb(Habit):
    def __init__(self, duration_val):
        super().__init__(duration_val)
        self.green_parts=('leaf', 'stem')

    def enter_dormancy(
            self,
            duration_params,
            current_jday,
            plants
        ):
        plants=self.duration.enter_dormancy(duration_params, current_jday, plants)
        return plants

    def emerge(self):
        # Use this to move carbohydrate to aboveground biomass from storage organs and roots
        self.duration.emerge()

    def initialize_biomass(self, grow_params):
        self.duration.initialize_biomass(grow_params)

class Graminoid(Habit):
    def __init__(self, duration_val):
        super().__init__(duration_val)

        self.green_parts=('leaf','stem')

    def enter_dormancy(
            self,
            duration_params,
            current_jday,
            plants
        ):
        plants=self.duration.enter_dormancy(duration_params, current_jday, plants)
        return plants

    def emerge(self):
        self.duration.emerge()

    def initialize_biomass(self, grow_params):
        init_mass_min, init_mass_max=self.duration.initialize_biomass(grow_params)
        return init_mass_min, init_mass_max

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
