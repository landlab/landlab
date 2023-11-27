from .duration import *


# Growth habit classes and selection method
# Growth habit uses duration properties to assign dormancy and emergence methods
class Habit(object):
    def __init__(self, species_grow_params, durationdict, duration_val, green_parts):
        self.duration = self.select_duration_class(
            species_grow_params,
            duration_val,
            green_parts,
            durationdict["senesce_rate"],
        )

    def select_duration_class(
        self,
        species_grow_params,
        duration_val,
        green_parts=("root", "leaf", "stem", "reproductive"),
        senesce_rate=0.0,
    ):
        duration = {
            "annual": Annual(species_grow_params, senesce_rate),
            "perennial deciduous": Deciduous(
                species_grow_params, green_parts, senesce_rate
            ),
            "perennial evergreen": Evergreen(species_grow_params),
        }
        return duration[duration_val]

    def emerge(self, plants):
        plants = self.duration.emerge(plants)
        return plants

    def senesce(self, plants, ns_green_mass, persistent_mass):
        plants = self.duration.senesce(plants, ns_green_mass, persistent_mass)
        return plants

    def set_initial_biomass(self, plants, in_growing_season):
        plants = self.duration.set_initial_biomass(plants, in_growing_season)
        return plants

    def enter_dormancy(self, plants):
        plants = self.duration.enter_dormancy(plants)
        return plants


class Forbherb(Habit):
    def __init__(self, species_grow_params, durationdict, duration_val):
        green_parts = ("leaf", "stem")
        super().__init__(species_grow_params, durationdict, duration_val, green_parts)


class Graminoid(Habit):
    def __init__(self, species_grow_params, durationdict, duration_val):
        green_parts = ("leaf", "stem")
        super().__init__(species_grow_params, durationdict, duration_val, green_parts)

    def set_initial_height(self, max_height, min_height, arr_size):
        height = min_height + rng.rayleigh(scale=0.5, size=arr_size) * max_height
        return height


class Shrub(Habit):
    def __init__(self, species_grow_params, durationdict, duration_val):
        green_parts = "leaf"
        super().__init__(species_grow_params, durationdict, duration_val, green_parts)


class Tree(Habit):
    def __init__(self, species_grow_params, durationdict, duration_val):
        green_parts = "leaf"
        super().__init__(species_grow_params, durationdict, duration_val, green_parts)


class Vine(Habit):
    def __init__(self, species_grow_params, durationdict, duration_val):
        green_parts = "leaf"
        super().__init__(species_grow_params, durationdict, duration_val, green_parts)
