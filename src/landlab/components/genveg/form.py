import numpy as np

from .dispersal import Clonal
from .dispersal import Random
from .dispersal import Seed

rng = np.random.default_rng()


class Abg_repro:
    def __init__(self, params):
        self.abg_parts = ("leaf", "stem", "reproductive")
        self.dead_abg_parts = ("dead_leaf", "dead_stem", "dead_reproductive")
        sums = {}
        sum_vars = [

            ["max_abg_biomass", "plant_part_max", self.abg_parts],
            ["min_abg_biomass", "plant_part_min", self.abg_parts],
        ]
        for sum_var in sum_vars:
            sums[sum_var[0]] = 0
            for part in sum_var[2]:
                sums[sum_var[0]] += params[
                    "grow_params"
                ][sum_var[1]][part]
        self.max_abg_biomass = sums["max_abg_biomass"]
        self.min_abg_biomass = sums["min_abg_biomass"]


class Blg_repro:
    def __init__(self, params):
        self.abg_parts = ("leaf", "stem")
        self.dead_abg_parts = ("dead_leaf", "dead_stem")
        sums = {}
        sum_vars = [

            ["max_abg_biomass", "plant_part_max", self.abg_parts],
            ["min_abg_biomass", "plant_part_min", self.abg_parts],
        ]
        for sum_var in sum_vars:
            sums[sum_var[0]] = 0
            for part in sum_var[2]:
                sums[sum_var[0]] += params[
                    "grow_params"
                ][sum_var[1]][part]
        self.max_abg_biomass = sums["max_abg_biomass"]
        self.min_abg_biomass = sums["min_abg_biomass"]


# Growth form classes and selection method
class Bunch(Seed):
    def __init__(self, params):
        super().__init__(params)

    def branch(self):
        print("Limited lateral branching due to clumping")

    def disperse(self, plants):
        return plants


class Colonizing(Random):
    def __init__(self, params):
        super().__init__(params)

    def branch(self):
        print("No branching annual")

    def disperse(self, plants):
        return plants


class Multiplestems(Seed):
    def __init__(self, params):
        super().__init__(params)

    def branch(self):
        print("Create two or more main stems at or near soil surface")

    def disperse(self, plants):
        return plants


class Rhizomatous(Clonal):
    def __init__(self, params):
        super().__init__(params)

    def set_initial_branches(self, max_branches, arr_size):
        n_branches = np.ceil(rng.rayleigh(scale=0.26, size=arr_size) * max_branches)
        return n_branches

    def branch(self):
        print("Tiller via rhizomes")


class Singlecrown(Seed):
    def __init__(self, params):
        super().__init__(params)

    def branch(self):
        print("Herbaceous plant with one persistent base")

    def disperse(self, plants):
        return plants


class Singlestem(Seed):
    def __init__(self, params):
        super().__init__(params)

    def branch(self):
        print("Plant develops one stem like a tree or a corn plant")

    def disperse(self, plants):
        return plants


class Stoloniferous(Clonal):
    def __init__(self, params):
        super().__init__(params)

    def branch(self):
        print("Branches via stolons")


class Thicketforming(Seed):
    def __init__(self, params):
        super().__init__(params)

    def branch(self):
        print("Limited lateral branching due to dense thickets")

    def disperse(self, plants):
        return plants
