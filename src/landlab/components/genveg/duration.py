import numpy as np
from scipy.optimize import fsolve

rng = np.random.default_rng()


# Duration class and child classes Annual and Perennial
# (Child classes are Deciduous and Evergreen)
class Duration:
    def __init__(
        self,
        species_grow_params,
        duration_params,
        dt,
        green_parts,
        persistent_parts,
    ):
        self.dt = dt
        self.growdict = species_grow_params
        self.allocation_coeffs = {
            "root_to_leaf": species_grow_params["root_to_leaf"],
            "root_to_stem": species_grow_params["root_to_stem"],
        }
        self.green_parts = green_parts
        self.persistent_parts = persistent_parts
        self.death_rate = duration_params["death_rate"]

    def enter_dormancy(self, plants):
        print("I kill green parts at end of growing season")
        for part in self.green_parts:
            plants["dead_" + str(part)] += plants[part]
            plants[part] = np.zeros_like(plants[part], dtype=np.float64)
        return plants

    def senesce(self, plants, ns_green_mass=None, persistent_mass=None):
        # use senesce rate to move portion of nsc to persistent parts, then remove
        filter = np.nonzero((ns_green_mass > 0.001) | (persistent_mass > 0.001))
        for part in self.persistent_parts:
            plants[part][filter] += (
                self.death_rate[part]
                * ns_green_mass[filter]
                * plants[part][filter]
                / persistent_mass[filter]
            )
        for part in self.green_parts:
            plants[part][filter] -= self.death_rate[part] * plants[part][filter]
        return plants

    def set_new_biomass(self, plants):
        print("I create new plants")
        total_biomass_ideal = rng.uniform(
            low=self.growdict["growth_min_biomass"],
            high=2 * self.growdict["growth_min_biomass"],
            size=plants.size,
        )
        (
            plants["root_biomass"],
            plants["leaf_biomass"],
            plants["stem_biomass"],
        ) = self._solve_biomass_allocation(total_biomass_ideal)
        plants["repro_biomass"] = np.zeros_like(plants["root_biomass"])
        plants["plant_age"] = np.zeros_like(plants["root_biomass"])
        return plants

    def _solve_biomass_allocation(self, total_biomass):
        # Initialize arrays to calculate root, leaf and stem biomass from total
        root = []
        leaf = []
        stem = []

        # Loop through grid array
        for total_biomass_in_cell in total_biomass:
            solver_guess = np.full(3, np.log10(total_biomass_in_cell / 3))
            part_biomass_log10 = fsolve(
                self._solverFuncs,
                solver_guess,
                (self.allocation_coeffs, total_biomass_in_cell),
            )
            part_biomass = 10**part_biomass_log10

            root.append(part_biomass[0])
            leaf.append(part_biomass[1])
            stem.append(part_biomass[2])

        # Convert to numpy array
        root = np.array(root)
        leaf = np.array(leaf)
        stem = np.array(stem)
        return root, leaf, stem

    def _solverFuncs(self, solver_guess, solver_coeffs, total_biomass):
        root_part_log10 = solver_guess[0]
        leaf_part_log10 = solver_guess[1]
        stem_part_log10 = solver_guess[2]
        plant_part_biomass_log10 = np.empty([(3)])

        plant_part_biomass_log10[0] = (
            10**root_part_log10
            + 10**leaf_part_log10
            + 10**stem_part_log10
            - total_biomass
        )
        plant_part_biomass_log10[1] = (
            solver_coeffs["root_to_leaf"]["a"]
            + solver_coeffs["root_to_leaf"]["b1"] * root_part_log10
            + solver_coeffs["root_to_leaf"]["b2"] * root_part_log10**2
            - leaf_part_log10
        )
        plant_part_biomass_log10[2] = (
            solver_coeffs["root_to_stem"]["a"]
            + solver_coeffs["root_to_stem"]["b1"] * root_part_log10
            + solver_coeffs["root_to_stem"]["b2"] * root_part_log10**2
            - stem_part_log10
        )
        return plant_part_biomass_log10


class Annual(Duration):
    def __init__(self, species_grow_params, duration_params, dt):
        green_parts = ("root", "leaf", "stem", "reproductive")
        persistent_parts = ()
        super().__init__(
            species_grow_params, duration_params, dt, green_parts, persistent_parts
        )

    def senesce(self, plants, mass_of_green=None, mass_of_persistent=None):
        print("I start to lose biomass during senescence periood")
        for part in self.green_parts:
            plants[part] -= plants[part] * self.death_rate[part] * self.dt
        return plants

    def emerge(self, plants, available_mass, persistent_total_mass):
        print("I emerge from dormancy and I am an annual")
        plants = self.set_new_biomass(plants)
        return plants

    def set_initial_biomass(self, plants, in_growing_season):
        if in_growing_season:
            plants = self.set_new_biomass(plants)
        else:
            plants = self.enter_dormancy(plants)
        return plants


class Perennial(Duration):
    def __init__(
        self,
        species_grow_params,
        duration_params,
        dt,
        green_parts=(),
        persistent_parts=("root", "leaf", "stem"),
    ):
        self.persistent_parts = persistent_parts
        super().__init__(
            species_grow_params, duration_params, dt, green_parts, persistent_parts
        )

    def senesce(self, plants, ns_green_mass=None, persistent_mass=None):
        # use senesce rate to move portion of nsc to persistent parts, then remove
        filter = np.nonzero((ns_green_mass > 0.001) | (persistent_mass > 0.001))
        for part in self.persistent_parts:
            plants[part][filter] += (
                ns_green_mass[filter]
                * plants[part][filter]
                / persistent_mass[filter]
            )
        for part in self.green_parts:
            plants[part][filter] -= self.death_rate[part] * plants[part][filter]
        return plants

    def set_new_biomass(self, plants):
        return super().set_new_biomass(plants)

    def set_initial_biomass(self, plants, in_growing_season):
        plants = self.set_new_biomass(plants)
        return plants


class Evergreen(Perennial):
    def __init__(self, species_grow_params, duration_params, dt):
        self.keep_green_parts = True
        super().__init__(species_grow_params, duration_params, dt)

    def set_initial_biomass(self, plants, in_growing_season):
        return super().set_initial_biomass(plants, in_growing_season)

    def set_new_biomass(self, plants):
        return super().set_new_biomass(plants)

    def emerge(self, plants, available_mass, persistent_total_mass):
        return plants

    def senesce(self, plants, mass_of_green=None, mass_of_persistent=None):
        return plants


class Deciduous(Perennial):
    def __init__(self, species_grow_params, duration_params, dt, green_parts):
        self.keep_green_parts = False
        all_veg_sources = ("root", "leaf", "stem", "reproductive")
        persistent_parts = []
        for part in all_veg_sources:
            if part not in green_parts:
                persistent_parts.append(part)
        persistent_parts = tuple(persistent_parts)

        super().__init__(
            species_grow_params, duration_params, dt, green_parts, persistent_parts
        )

    def emerge(self, plants, available_mass, persistent_total_mass):
        print("I emerge from dormancy and I am a deciduous perennial")
        # next steps are to clean this up using same approach as dormancy

        total_mass_new_green = np.zeros_like(plants["root_biomass"])
        new_green_biomass = {}
        for part in self.green_parts:
            new_green_biomass[part] = rng.uniform(
                low=self.growdict["plant_part_min"][part],
                high=self.growdict["plant_part_min"][part] * 2,
                size=plants.size,
            )
            total_mass_new_green += new_green_biomass[part]
        adjusted_total_new_green = np.minimum(available_mass, total_mass_new_green)

        for part in self.green_parts:
            plants[part] = (
                adjusted_total_new_green / total_mass_new_green
            ) * new_green_biomass[part]
            new_green_biomass[part] = plants[part]

        total_mass_emerged = sum(new_green_biomass.values())
        for part in self.persistent_parts:
            deleted_mass = total_mass_emerged * plants[part] / persistent_total_mass
            plants[part] -= deleted_mass
        return plants

    def enter_dormancy(self, plants):
        return super().enter_dormancy(plants)

    def set_initial_biomass(self, plants, in_growing_season):
        if not in_growing_season:
            for part in self.green_parts:
                plants[part] = np.zeros_like(plants[part])
            for part in self.persistent_parts:
                plants[part] = rng.uniform(
                    low=self.growdict["plant_part_min"][part],
                    high=self.growdict["plant_part_min"][part] * 2,
                    size=plants.size,
                )
        else:
            for part in (self.green_parts + self.persistent_parts):
                plants[part] = rng.uniform(
                    low=self.growdict["plant_part_min"][part],
                    high=self.growdict["plant_part_min"][part] * 2,
                    size=plants.size,
                )
            plants["repro_biomass"] = (
                self.growdict["plant_part_min"]["reproductive"]
                + rng.rayleigh(scale=0.2, size=plants.size)
                * self.growdict["plant_part_max"]["reproductive"]
            )
        return plants
