import numpy as np

rng = np.random.default_rng()


# Dispersal classes and selection method
class Repro:
    def __init__(self, params):
        grow_params = params["grow_params"]
        self.min_size = grow_params["growth_biomass"]["min"]


class Clonal(Repro):
    def __init__(self, params):
        disperse_params = params["dispersal_params"]
        super().__init__(params)
        self.unit_cost = disperse_params["unit_cost_dispersal"]
        self.max_dist_dispersal = disperse_params["max_dist_dispersal"]

    def disperse(self, plants):
        """
        This method determines the potential location and cost
        of clonal dispersal. A separate method in integrator.py
        determines if dispersal is successful based if the potential location
        is occupied.
        """
        max_runner_length = np.zeros_like(plants["root"])
        available_carb = plants["reproductive"] - 2 * self.min_size
        max_runner_length[available_carb > 0] = (
            available_carb[available_carb > 0] / self.unit_cost
        )
        runner_length = rng.uniform(
            low=0.05,
            high=self.max_dist_dispersal,
            size=plants.size,
        )
        pup_azimuth = np.deg2rad(rng.uniform(low=0.01, high=360, size=plants.size))

        filter = np.nonzero(runner_length <= max_runner_length)

        plants["pup_x_loc"][filter] = (
            runner_length[filter] * np.cos(pup_azimuth[filter])
            + plants["x_loc"][filter]
        )
        plants["pup_y_loc"][filter] = (
            runner_length[filter] * np.sin(pup_azimuth)[filter]
            + plants["y_loc"][filter]
        )
        plants["pup_cost"][filter] = runner_length[filter] * self.unit_cost

        return plants


class Seed(Repro):
    def __init__(self, params):
        pass

    def disperse(self):
        print("New plant emerges from seed some distance within parent plant")


class Random(Repro):
    def __init__(self, params):
        pass

    def disperse(self):
        print("New plant randomly appears")
