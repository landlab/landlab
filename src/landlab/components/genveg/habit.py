import numpy as np

from .check_objects import UnitTestChecks
from .duration import Annual
from .duration import Deciduous
from .duration import Evergreen

rng = np.random.default_rng()


# Growth habit classes and selection method
# Growth habit uses duration properties to assign dormancy and emergence methods
class Habit:
    def __init__(self, params, dt=1, green_parts=(None)):
        self.grow_params = params["grow_params"]
        duration_params = params["duration_params"]
        self.duration = self._select_duration_class(
            self.grow_params,
            duration_params,
            dt,
            params["plant_factors"]["duration"],
            green_parts,
        )
        self.morph_params = self._calc_derived_morph_params(params)

    def _calc_canopy_area_from_shoot_width(self, shoot_sys_width):
        UnitTestChecks().is_negative_present(shoot_sys_width, "shoot_sys_width")
        canopy_area = 0.25 * np.pi * shoot_sys_width**2
        return canopy_area

    def _calc_crown_volume(self, shoot_width, basal_width, shoot_height):
        shoot_rad = shoot_width / 2
        basal_rad = basal_width / 2
        volume = (
            1
            / 3
            * np.pi
            * shoot_height
            * (shoot_rad**2 + shoot_rad * basal_rad + basal_rad**2)
        )
        return volume

    def calc_root_sys_width(self, shoot_sys_width, basal_width, shoot_sys_height=1):
        volume = self._calc_crown_volume(shoot_sys_width, basal_width, shoot_sys_height)
        root_sys_width = 0.08 + 0.24 * volume
        root_sys_width[root_sys_width > self.morph_params["max_root_sys_width"]] = (
            self.morph_params["max_root_sys_width"]
        )
        root_sys_width[root_sys_width < self.morph_params["min_root_sys_width"]] = (
            self.morph_params["max_root_sys_width"]
        )
        return root_sys_width

    def _calc_shoot_width_from_canopy_area(self, canopy_area):
        shoot_width_cm = np.sqrt(4 * canopy_area / np.pi)
        shoot_width = shoot_width_cm / 100
        return shoot_width

    def _select_duration_class(
        self,
        species_grow_params,
        duration_params,
        dt,
        duration_val,
        green_parts=(None),
    ):
        duration = {
            "annual": Annual(species_grow_params, duration_params, dt),
            "perennial deciduous": Deciduous(
                species_grow_params, duration_params, dt, green_parts
            ),
            "perennial evergreen": Evergreen(species_grow_params, duration_params, dt),
        }
        return duration[duration_val]

    def _calc_derived_morph_params(self, params):
        # Calculate growth habit specific parameters depending on morphology
        # allometry method and save them as a dictionary to the habit class
        params["morph_params"]["max_canopy_area"] = (
            np.pi / 4 * params["morph_params"]["max_shoot_sys_width"] ** 2
        )
        params["morph_params"]["min_canopy_area"] = (
            np.pi / 4 * params["morph_params"]["min_shoot_sys_width"] ** 2
        )
        params["morph_params"]["max_root_area"] = (
            np.pi / 4 * params["morph_params"]["max_root_sys_width"] ** 2
        )
        params["morph_params"]["min_root_area"] = (
            np.pi / 4 * params["morph_params"]["min_root_sys_width"] ** 2
        )
        return params["morph_params"]

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # These dimensions are empirically derived allometric relationships.
        # Using ones for placeholder right now
        log_basal_width = basal_width = log_height = height = log_canopy_area = (
            canopy_area
        ) = np.ones_like(abg_biomass)
        basal_width = np.exp(log_basal_width)
        height = np.exp(log_height)
        canopy_area = np.exp(log_canopy_area)
        shoot_sys_width = self._calc_shoot_width_from_canopy_area(canopy_area)
        return basal_width, shoot_sys_width, height

    def emerge(self, plants):
        plants = self.duration.emerge(plants)
        return plants

    def estimate_abg_biomass_from_morphology(self, plants):
        # Edit this to a common form
        log_basal_width_cm = np.log(
            (plants["basal_width"] * 100),
            out=np.zeros_like(plants["basal_width"], dtype=np.float64),
            where=(plants["basal_width"] > 0.0),
        )
        log_abg_biomass = (
            self.morph_params["basal_coeffs"]["a"]
            + self.morph_params["basal_coeffs"]["b"] * log_basal_width_cm
        )
        est_abg_biomass = np.exp(
            log_abg_biomass,
            out=np.zeros_like(log_abg_biomass, dtype=np.float64),
            where=(plants["basal_width"] > 0.0),
        )
        return est_abg_biomass

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
    def __init__(self, params, dt):
        green_parts = ("leaf", "stem")
        super().__init__(params, dt, green_parts)


class Graminoid(Habit):
    def __init__(
        self,
        params,
        dt,
        empirical_coeffs={
            "perennial": {
                "C3": {
                    "basal_coeffs": {"a": 0.20, "b": 1.17},
                    "height_coeffs": {"a": np.log(0.232995), "b": 0.619077},
                    "canopy_coeffs": {"a": np.log(0.0597478), "b": 1.31244},
                },
                "C4": {
                    "basal_coeffs": {"a": 0.20, "b": 1.17},
                    "height_coeffs": {"a": np.log(0.2776634), "b": 0.4176197},
                    "canopy_coeffs": {"a": np.log(0.0516016), "b": 1.38799},
                },
            },
            "annual": {
                "C3": {
                    "basal_coeffs": {"a": 0.20, "b": 1.17},
                    "height_coeffs": {"a": np.log(0.8548639), "b": 0.9187837},
                    "canopy_coeffs": {"a": np.log(0.030101), "b": 1.24346249},
                },
                "C4": {
                    "basal_coeffs": {"a": 0.20, "b": 1.17},
                    "height_coeffs": {"a": np.log(0.2776634), "b": 0.4176197},
                    "canopy_coeffs": {"a": np.log(0.157143), "b": 1.3828},
                },
            },
        },
    ):
        green_parts = ("leaf", "stem")
        duration_val = params["plant_factors"]["duration"]
        photo_val = params["plant_factors"]["p_type"]

        if params["morph_params"]["allometry_method"] == "min-max":
            max_basal_dia_cm = params["morph_params"]["max_basal_dia"] * 100
            min_basal_dia_cm = params["morph_params"]["min_basal_dia"] * 100
            (
                params["morph_params"]["basal_coeffs"]["a"],
                params["morph_params"]["basal_coeffs"]["b"],
            ) = self._calc2_allometry_coeffs(
                min_basal_dia_cm,
                max_basal_dia_cm,
                params["grow_params"]["min_abg_biomass"],
                params["grow_params"]["max_abg_biomass"],
            )
            (
                params["morph_params"]["height_coeffs"]["a"],
                params["morph_params"]["height_coeffs"]["b"],
            ) = self._calc2_allometry_coeffs(
                min_basal_dia_cm,
                max_basal_dia_cm,
                params["morph_params"]["min_height"],
                params["morph_params"]["max_height"],
            )
            min_canopy_area = self._calc_canopy_area_from_shoot_width(
                params["morph_params"]["min_shoot_sys_width"]
            )
            max_canopy_area = self._calc_canopy_area_from_shoot_width(
                params["morph_params"]["max_shoot_sys_width"]
            )
            (
                params["morph_params"]["canopy_coeffs"]["a"],
                params["morph_params"]["canopy_coeffs"]["b"],
            ) = self._calc2_allometry_coeffs(
                min_basal_dia_cm,
                max_basal_dia_cm,
                min_canopy_area,
                max_canopy_area,
            )

        else:
            if params["morph_params"]["allometry_method"] == "default":
                params["morph_params"]["basal_coeffs"] = empirical_coeffs[duration_val][
                    photo_val
                ]["basal_coeffs"]
                params["morph_params"]["height_coeffs"] = empirical_coeffs[
                    duration_val
                ][photo_val]["height_coeffs"]
                params["morph_params"]["canopy_coeffs"] = empirical_coeffs[
                    duration_val
                ][photo_val]["canopy_coeffs"]

            (
                params["morph_params"]["min_basal_dia"],
                params["morph_params"]["min_shoot_sys_width"],
                params["morph_params"]["min_shoot_sys_height"],
            ) = self.calc_abg_dims_from_biomass(
                params["grow_params"]["min_abg_biomass"]
            )
            (
                params["morph_params"]["max_basal_dia"],
                params["morph_params"]["max_shoot_sys_width"],
                params["morph_params"]["max_shoot_sys_height"],
            ) = self._calc_abg_dims_from_biomass(
                params["grow_params"]["max_abg_biomass"]
            )
            params["morph_params"]["min_basal_dia"]

        super().__init__(params, dt, green_parts)

    def _calc2_allometry_coeffs(self, x_min, x_max, y_min, y_max):
        ln_x_min = np.log(x_min)
        ln_x_max = np.log(x_max)
        ln_y_min = np.log(y_min)
        ln_y_max = np.log(y_max)
        b = (ln_y_max - ln_y_min) / (ln_x_max - ln_x_min)
        a = ln_y_max - b * ln_x_max
        return (a, b)

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # These dimensions are empirically derived allometric relationships for grasses.
        log_basal_width_cm = log_height = height = canopy_area = np.zeros_like(
            abg_biomass
        )
        filter = np.nonzero(abg_biomass > self.grow_params["min_abg_biomass"])
        log_basal_width_cm[filter] = (
            np.log(abg_biomass[filter]) - self.morph_params["basal_coeffs"]["a"]
        ) / self.morph_params["basal_coeffs"]["b"]
        basal_width_cm = np.exp(log_basal_width_cm)
        basal_width = basal_width_cm / 100
        canopy_area = self._calc_canopy_area(basal_width)
        log_height[filter] = (
            self.morph_params["height_coeffs"]["a"] - log_basal_width_cm[filter]
        ) / self.morph_params["height_coeffs"]["b"]
        height = np.exp(log_height)
        shoot_sys_width = self._calc_shoot_width_from_canopy_area(canopy_area)
        return basal_width, shoot_sys_width, height

    def _calc_canopy_area(self, basal_width):
        filter = np.nonzero(basal_width >= 0.0)
        log_canopy_area = log_basal_width_cm = np.zeros_like(basal_width)
        log_basal_width_cm[filter] = np.log(basal_width[filter] * 100)
        log_canopy_area[filter] = (
            log_basal_width_cm[filter] - self.morph_params["canopy_coeffs"]["a"]
        ) / self.morph_params["canopy_coeffs"]["b"]
        canopy_area = np.exp(log_canopy_area)
        return canopy_area

    def estimate_abg_biomass_from_morphology(self, plants):
        log_basal_width_cm = np.log(
            (plants["basal_width"] * 100),
            out=np.zeros_like(plants["basal_width"], dtype=np.float64),
            where=(plants["basal_width"] > 0.0),
        )
        log_abg_biomass = (
            self.morph_params["basal_coeffs"]["a"]
            + self.morph_params["basal_coeffs"]["b"] * log_basal_width_cm
        )
        est_abg_biomass = np.exp(
            log_abg_biomass,
            out=np.zeros_like(log_abg_biomass, dtype=np.float64),
            where=(plants["basal_width"] > 0.0),
        )
        return est_abg_biomass


class Shrub(Habit):
    def __init__(self, params, dt):
        green_parts = ("leaf")
        super().__init__(params, dt, green_parts)


class Tree(Habit):
    def __init__(self, params, dt):
        green_parts = ("leaf")
        super().__init__(params, dt, green_parts)


class Vine(Habit):
    def __init__(self, params, dt):
        green_parts = ("leaf")
        super().__init__(params, dt, green_parts)
