import numpy as np
from scipy import linalg

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
        self.duration_params = params["duration_params"]
        self.morph_params = params["morph_params"]
        self.duration = self._select_duration_class(
            self.grow_params,
            self.duration_params,
            dt,
            params["plant_factors"]["duration"],
            green_parts,
        )

    def _calc2_allometry_coeffs(self, x_min, x_max, y_min, y_max):
        """
        Calculates coefficients for natural log relationships
        log(y) = a + b(log(x))
        """
        ln_x_min = np.log(x_min)
        ln_x_max = np.log(x_max)
        ln_y_min = np.log(y_min)
        ln_y_max = np.log(y_max)
        b = (ln_y_max - ln_y_min) / (ln_x_max - ln_x_min)
        a = ln_y_max - b * ln_x_max
        return (a, b)

    def _apply2_allometry_eq_for_xs(self, ys, predict_var):
        # Assuming you want to solve for ys rather than xs
        ln_ys = np.log(ys)
        a = self.morph_params[predict_var]["a"]
        b = self.morph_params[predict_var]["b"]
        ln_xs = (ln_ys - a) / b
        return np.exp(ln_xs)

    def _apply2_allometry_eq_for_ys(self, xs, predict_var):
        # Assuming you want to solve for ys rather than xs
        ln_xs = np.log(xs)
        a = self.morph_params[predict_var]["a"]
        b = self.morph_params[predict_var]["b"]
        ln_ys = a + (b * ln_xs)
        return np.exp(ln_ys)

    def _calc3_allometry_coeffs(self, x_min, x_mean, x_max, y_min, y_mean, y_max, z_min, z_mean, z_max):
        """
        Calculates coefficients for natural log relationships
        log(z) = a + b(log(x)) + c(log(y))
        """
        ln_x_min = np.log(x_min)
        ln_x_mean = np.log(x_mean)
        ln_x_max = np.log(x_max)
        ln_y_min = np.log(y_min)
        ln_y_mean = np.log(y_mean)
        ln_y_max = np.log(y_max)
        ln_z_min = np.log(z_min)
        ln_z_mean = np.log(z_mean)
        ln_z_max = np.log(z_max)
        eqs = np.array([
            [1, ln_x_min.item(), ln_y_min.item()],
            [1, ln_x_mean.item(), ln_y_mean.item()],
            [1, ln_x_max.item(), ln_y_max.item()]
        ])
        z_est = np.array([ln_z_min.item(), ln_z_mean.item(), ln_z_max.item()])
        out = linalg.solve(eqs, z_est)
        return (out.item(0), out.item(1), out.item(2))

    def _apply3_allometry_eq_for_zs(self, xs, ys, predict_var):
        ln_xs = np.log(xs)
        ln_ys = np.log(ys)
        a = self.morph_params[predict_var]["a"]
        b = self.morph_params[predict_var]["b"]
        c = self.morph_params[predict_var]["c"]
        ln_zs = a + b * ln_xs + c * ln_ys
        return np.exp(ln_zs)

    def _apply3_allometry_eq_for_ys(self, xs, zs, predict_var):
        ln_xs = np.log(xs)
        ln_zs = np.log(zs)
        a = self.morph_params[predict_var]["a"]
        b = self.morph_params[predict_var]["b"]
        c = self.morph_params[predict_var]["c"]
        ln_ys = (ln_zs - a - b * ln_xs) / c
        return np.exp(ln_ys)

    def _apply3_allometry_eq_for_xs(self, ys, zs, predict_var):
        ln_ys = np.log(ys)
        ln_zs = np.log(zs)
        a = self.morph_params[predict_var]["a"]
        b = self.morph_params[predict_var]["b"]
        c = self.morph_params[predict_var]["c"]
        ln_xs = (ln_zs - a - c * ln_ys) / b
        return np.exp(ln_xs)

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
        root_sys_width[root_sys_width > self.morph_params["root_sys_width"]["max"]] = (
            self.morph_params["root_sys_width"]["max"]
        )
        root_sys_width[root_sys_width < self.morph_params["root_sys_width"]["min"]] = (
            self.morph_params["root_sys_width"]["min"]
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
        # Assume log-log scaling for all params with biomass
        params["morph_params"]["basal_dia_cm"] = {}
        params["morph_params"]["basal_dia_cm"]["min"] = params["morph_params"]["basal_dia"]["min"] * 100
        params["morph_params"]["basal_dia_cm"]["max"] = params["morph_params"]["basal_dia"]["max"] * 100
        (
            params["morph_params"]["basal_coeffs"]["a"],
            params["morph_params"]["basal_coeffs"]["b"],
        ) = self._calc2_allometry_coeffs(
            params["morph_params"]["basal_dia_cm"]["min"],
            params["morph_params"]["basal_dia_cm"]["max"],
            params["grow_params"]["abg_biomass"]["min"],
            params["grow_params"]["abg_biomass"]["max"],
        )
        (
            params["morph_params"]["height_coeffs"]["a"],
            params["morph_params"]["height_coeffs"]["b"],
        ) = self._calc2_allometry_coeffs(
            params["morph_params"]["height"]["min"],
            params["morph_params"]["height"]["max"],
            params["grow_params"]["abg_biomass"]["min"],
            params["grow_params"]["abg_biomass"]["max"],
        )
        params["morph_params"]["canopy_area"] = {}
        params["morph_params"]["canopy_area"]["min"] = self._calc_canopy_area_from_shoot_width(
            params["morph_params"]["shoot_sys_width"]["min"]
        )
        params["morph_params"]["canopy_area"]["max"] = self._calc_canopy_area_from_shoot_width(
            params["morph_params"]["shoot_sys_width"]["max"]
        )
        (
            params["morph_params"]["canopy_coeffs"]["a"],
            params["morph_params"]["canopy_coeffs"]["b"],
        ) = self._calc2_allometry_coeffs(
            params["morph_params"]["canopy_area"]["min"],
            params["morph_params"]["canopy_area"]["max"],
            params["grow_params"]["abg_biomass"]["min"],
            params["grow_params"]["abg_biomass"]["max"],
        )
        return params

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # These dimensions are empirically derived allometric relationships.
        # Using ones for placeholder right now
        basal_dia_cm = height = canopy_area = np.ones_like(abg_biomass)
        filter = np.nonzero(abg_biomass > self.grow_params["abg_biomass"]["min"])
        basal_dia_cm[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], "basal_coeffs")
        basal_width = basal_dia_cm / 100
        height[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], "height_coeffs")
        canopy_area[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], "canopy_coeffs")
        shoot_sys_width = self._calc_shoot_width_from_canopy_area(canopy_area)
        return basal_width, shoot_sys_width, height

    def emerge(self, plants):
        plants = self.duration.emerge(plants)
        return plants

    def estimate_abg_biomass_from_cover(self, plants):
        # redo this to back calculate percent cover and height to basal diameter and abg biomass
        log_canopy_area = np.log(
            0.25 * np.pi * (plants["shoot_sys_width"])**2,
            out=np.zeros_like(plants["shoot_sys_width"], dtype=np.float64),
            where=(plants["shoot_sys_width"] > 0.0),
        )
        log_abg_biomass = (
            self.morph_params["canopy_coeffs"]["a"]
            + self.morph_params["canopy_coeffs"]["b"] * log_canopy_area
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
    def __init__(
        self,
        params,
        dt,
        empirical_coeffs={
            "height_coeffs": {
                "a": (np.log(1.403) - 0.370),
                "b": 1.903,
            },
            "canopy_coeffs": {
                "a": -2.057,
                "b": 1.741,
                "c": 0.945,
            },
            "basal_coeffs": {
                "a": (np.log(1.0787) - 2.757),
                "b": 2.474,
                "c": 1,
            },
        }
    ):
        params = self._calc_derived_morph_params(params)
        green_parts = ("leaf", "stem")
        super().__init__(params, dt, green_parts)

    def _calc_derived_morph_params(self, params):
        # Forbs and herbs morphology is not well defined
        # From Miao et al 2008 and Patzig et al 2020
        # - height is a function of aboveground biomass
        # - cover (canopy area) is a function of biomass and height
        # - basal area is a function of aboveground biomass and height
        # Need to better define defaults for monocot and dicot from literature.
        # Using max-min scaling only here
        params["morph_params"]["basal_area_cm2"] = {}
        params["morph_params"]["basal_area_cm2"]["min"] = self._calc_canopy_area_from_shoot_width(params["morph_params"]["basal_dia"]["min"] * 100)
        params["morph_params"]["basal_area_cm2"]["mean"] = self._calc_canopy_area_from_shoot_width(params["morph_params"]["basal_dia"]["mean"] * 100)
        params["morph_params"]["basal_area_cm2"]["max"] = self._calc_canopy_area_from_shoot_width(params["morph_params"]["basal_dia"]["max"] * 100)
        params["morph_params"]["canopy_area"] = {}
        params["morph_params"]["canopy_area"]["min"] = self._calc_canopy_area_from_shoot_width(params["morph_params"]["shoot_sys_width"]["min"])
        params["morph_params"]["canopy_area"]["mean"] = self._calc_canopy_area_from_shoot_width(params["morph_params"]["shoot_sys_width"]["mean"])
        params["morph_params"]["canopy_area"]["max"] = self._calc_canopy_area_from_shoot_width(params["morph_params"]["shoot_sys_width"]["max"])
        (
            params["morph_params"]["height_coeffs"]["a"],
            params["morph_params"]["height_coeffs"]["b"],
        ) = self._calc2_allometry_coeffs(
            params["morph_params"]["height"]["min"],
            params["morph_params"]["height"]["max"],
            params["grow_params"]["abg_biomass"]["min"],
            params["grow_params"]["abg_biomass"]["max"],
        )
        (
            params["morph_params"]["basal_coeffs"]["a"],
            params["morph_params"]["basal_coeffs"]["b"],
            params["morph_params"]["basal_coeffs"]["c"]
        ) = self._calc3_allometry_coeffs(
            params["morph_params"]["basal_area_cm2"]["min"],
            params["morph_params"]["basal_area_cm2"]["mean"],
            params["morph_params"]["basal_area_cm2"]["max"],
            params["morph_params"]["height"]["min"],
            params["morph_params"]["height"]["mean"],
            params["morph_params"]["height"]["max"],
            params["grow_params"]["abg_biomass"]["min"],
            params["grow_params"]["abg_biomass"]["mean"],
            params["grow_params"]["abg_biomass"]["max"],
        )
        (
            params["morph_params"]["canopy_coeffs"]["a"],
            params["morph_params"]["canopy_coeffs"]["b"],
            params["morph_params"]["canopy_coeffs"]["c"]
        ) = self._calc3_allometry_coeffs(
            params["morph_params"]["canopy_area"]["min"],
            params["morph_params"]["canopy_area"]["mean"],
            params["morph_params"]["canopy_area"]["max"],
            params["morph_params"]["height"]["min"],
            params["morph_params"]["height"]["mean"],
            params["morph_params"]["height"]["max"],
            params["grow_params"]["abg_biomass"]["min"],
            params["grow_params"]["abg_biomass"]["mean"],
            params["grow_params"]["abg_biomass"]["max"],
        )
        return params

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # These dimensions are empirically derived allometric relationships for grasses.
        basal_area_cm2 = height = canopy_area = np.zeros_like(
            abg_biomass
        )
        filter = np.nonzero(abg_biomass > self.grow_params["abg_biomass"]["min"])
        height[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], "height_coeffs")
        basal_area_cm2[filter] = self._apply3_allometry_eq_for_xs(height[filter], abg_biomass[filter], "basal_area_coeffs")
        basal_width = self._calc_shoot_width_from_canopy_area(basal_area_cm2) / 100
        canopy_area[filter] = self._apply3_allometry_eq_for_xs(height[filter], abg_biomass[filter], "canopy_coeffs")
        shoot_sys_width = self._calc_shoot_width_from_canopy_area(canopy_area)
        return basal_width, height, shoot_sys_width


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
                    "canopy_coeffs": {"a": np.log(0.23702483), "b": 0.72682, "c": 0.9459644},
                },
                "C4": {
                    "basal_coeffs": {"a": 0.20, "b": 1.17},
                    "height_coeffs": {"a": np.log(0.2776634), "b": 0.4176197},
                    "canopy_coeffs": {"a": np.log(0.06669907), "b": 1.3043469, "c": 0.2002879},
                },
            },
            "annual": {
                "C3": {
                    "basal_coeffs": {"a": 0.20, "b": 1.17},
                    "height_coeffs": {"a": np.log(0.1476171), "b": 0.6995105},
                    "canopy_coeffs": {"a": np.log(0.12826361), "b": 0.7134629, "c": 0.7576721},
                },
                "C4": {
                    "basal_coeffs": {"a": 0.20, "b": 1.17},
                    "height_coeffs": {"a": np.log(0.4204882), "b": 0.5194908},
                    "canopy_coeffs": {"a": np.log(0.25749493), "b": 1.0866763, "c": 0.5700335},
                },
            },
        },
    ):
        # Graminoid morphology
        # - basal diameter is a function of aboveground biomass
        # - height is a function of basal diameter
        # - canopy area is a function of height and basal diameter
        # - stem diameter and number of stems is a function of basal diameter
        # Look at how we can use this to loop calculation
        # allometry_relationships = {
        #    {"basal_dia": ["abg_biomass"]},
        #    {"height": ["basal_dia"]},
        #    {"canopy_area": ["basal_dia", "height"]}
        # }
        green_parts = ("leaf", "stem")
        params = self._calc_derived_morph_params(params, empirical_coeffs)
        super().__init__(params, dt, green_parts)

    def _calc_derived_morph_params(self, params, empirical_coeffs):
        # Gao et al 2024 scales grasses basal diameter to biomass, height to basal diameter
        # and canopy area to basal area
        duration_val = params["plant_factors"]["duration"]
        photo_val = params["plant_factors"]["p_type"]
        params["morph_params"]["basal_dia_cm"] = {}
        if params["morph_params"]["allometry_method"] == "min-max":
            params["morph_params"]["basal_dia_cm"]["min"] = params["morph_params"]["basal_dia"]["min"] * 100
            params["morph_params"]["basal_dia_cm"]["mean"] = params["morph_params"]["basal_dia"]["mean"] * 100
            params["morph_params"]["basal_dia_cm"]["max"] = params["morph_params"]["basal_dia"]["max"] * 100
            (
                params["morph_params"]["basal_coeffs"]["a"],
                params["morph_params"]["basal_coeffs"]["b"],
            ) = self._calc2_allometry_coeffs(
                params["morph_params"]["basal_dia_cm"]["min"],
                params["morph_params"]["basal_dia_cm"]["max"],
                params["grow_params"]["abg_biomass"]["min"],
                params["grow_params"]["abg_biomass"]["max"],
            )
            (
                params["morph_params"]["height_coeffs"]["a"],
                params["morph_params"]["height_coeffs"]["b"],
            ) = self._calc2_allometry_coeffs(
                params["morph_params"]["basal_dia_cm"]["min"],
                params["morph_params"]["basal_dia_cm"]["max"],
                params["morph_params"]["height"]["min"],
                params["morph_params"]["height"]["max"],
            )
            params["morph_params"]["canopy_area"] = {}
            params["morph_params"]["canopy_area"]["min"] = self._calc_canopy_area_from_shoot_width(
                params["morph_params"]["shoot_sys_width"]["min"]
            )
            params["morph_params"]["canopy_area"]["mean"] = self._calc_canopy_area_from_shoot_width(
                params["morph_params"]["shoot_sys_width"]["mean"]
            )
            params["morph_params"]["canopy_area"]["max"] = self._calc_canopy_area_from_shoot_width(
                params["morph_params"]["shoot_sys_width"]["max"]
            )
            (
                params["morph_params"]["canopy_coeffs"]["a"],
                params["morph_params"]["canopy_coeffs"]["b"],
                params["morph_params"]["canopy_coeffs"]["c"]
            ) = self._calc3_allometry_coeffs(
                params["morph_params"]["basal_dia_cm"]["min"],
                params["morph_params"]["basal_dia_cm"]["mean"],
                params["morph_params"]["basal_dia_cm"]["max"],
                params["morph_params"]["height"]["min"],
                params["morph_params"]["height"]["mean"],
                params["morph_params"]["height"]["max"],
                params["morph_params"]["canopy_area"]["min"],
                params["morph_params"]["canopy_area"]["mean"],
                params["morph_params"]["canopy_area"]["max"],
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
                params["morph_params"]["basal_dia"]["min"],
                params["morph_params"]["height"]["min"],
                params["morph_params"]["shoot_sys_width"]["min"],
            ) = self.calc_abg_dims_from_biomass(
                params["grow_params"]["abg_biomass"]["min"]
            )
            (
                params["morph_params"]["basal_dia"]["max"],
                params["morph_params"]["height"]["max"],
                params["morph_params"]["shoot_sys_width"]["max"],
            ) = self.calc_abg_dims_from_biomass(
                params["grow_params"]["abg_biomass"]["max"]
            )
            (
                params["morph_params"]["basal_dia"]["mean"],
                params["morph_params"]["height"]["mean"],
                params["morph_params"]["shoot_sys_width"]["mean"],
            ) = self.calc_abg_dims_from_biomass(
                params["grow_params"]["abg_biomass"]["mean"]
            )
        return params

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # These dimensions are empirically derived allometric relationships for grasses.
        basal_dia_cm = height = canopy_area = np.zeros_like(
            abg_biomass
        )
        filter = np.nonzero(abg_biomass > self.grow_params["abg_biomass"]["min"])
        basal_dia_cm[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], "basal_coeffs")
        basal_width = basal_dia_cm / 100
        height[filter] = self._apply2_allometry_eq_for_ys(basal_dia_cm[filter], "height_coeffs")
        canopy_area[filter] = self._apply3_allometry_eq_for_zs(basal_dia_cm[filter], height[filter], "canopy_coeffs")
        shoot_sys_width = self._calc_shoot_width_from_canopy_area(canopy_area)
        return basal_width, height, shoot_sys_width

    def estimate_abg_biomass_from_cover(self, plants):
        # Edit to derive back from percent cover assuming cover for grasses is more reflective of basal diameter than canopy area
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
    def __init__(self, params, dt, empirical_coeffs={
            "basal_coeffs": {
                "a": (np.log(1.0787) - 2.757),
                "b": 2.474,
            },
            "shoot_sys_width_coeffs": {
                "a": -2.057,
                "b": 1.741,
                "c": 0.945,
            },
            "height_coeffs": {
                "a": (np.log(1.403) - 0.370),
                "b": 1.903,
                "c": 0.652,
            }}):
        green_parts = ("leaf")
        params = self._calc_derived_morph_params(params, empirical_coeffs)
        super().__init__(params, dt, green_parts)
        # Shrub morhology - note these relationships will compound uncertainty
        # - basal diameter is a function of aboveground biomass - in cm
        # - shoot system width is a function of aboveground biomass and basal diameter in m
        # - height is a function of aboveground biomass, basal diameter, amd crown diameter in m

    def _calc_derived_morph_params(self, params, empirical_coeffs):
        # Conti et al 2018 scales shrubs basal diameter to biomass, canopy diameter
        # (shoot system width) to biomass and basal diameter,
        # and height to biomass and canopy diameter
        # Note these relationships are designed to predict aboveground biomass from
        # field measured allometric variables and their combination in this form will
        # result in compounded uncertainty and should not be considered accurate for all
        # species and individuals
        if params["morph_params"]["allometry_method"] == "min-max":
            params["morph_params"]["basal_dia_cm"]["min"] = params["morph_params"]["basal_dia"]["min"] * 100
            params["morph_params"]["basal_dia_cm"]["mean"] = params["morph_params"]["basal_dia"]["mean"] * 100
            params["morph_params"]["basal_dia_cm"]["max"] = params["morph_params"]["basal_dia"]["max"] * 100
            (
                params["morph_params"]["basal_coeffs"]["a"],
                params["morph_params"]["basal_coeffs"]["b"],
            ) = self._calc2_allometry_coeffs(
                params["morph_params"]["basal_dia_cm"]["min"],
                params["morph_params"]["basal_dia_cm"]["max"],
                params["grow_params"]["abg_biomass"]["min"],
                params["grow_params"]["abg_biomass"]["max"],
            )
            (
                params["morph_params"]["canopy_coeffs"]["a"],
                params["morph_params"]["canopy_coeffs"]["b"],
                params["morph_params"]["canopy_coeffs"]["c"]
            ) = self._calc3_allometry_coeffs(
                params["morph_params"]["basal_dia_cm"]["min"],
                params["morph_params"]["basal_dia_cm"]["mean"],
                params["morph_params"]["basal_dia_cm"]["max"],
                params["morph_params"]["shoot_sys_width"]["min"],
                params["morph_params"]["shoot_sys_width"]["mean"],
                params["morph_params"]["shoot_sys_width"]["max"],
                params["grow_params"]["abg_biomass"]["min"],
                params["grow_params"]["abg_biomass"]["mean"],
                params["grow_params"]["abg_biomass"]["max"],
            )
            (
                params["morph_params"]["height_coeffs"]["a"],
                params["morph_params"]["height_coeffs"]["b"],
                params["morph_params"]["height_coeffs"]["c"]
            ) = self._calc3_allometry_coeffs(
                params["morph_params"]["shoot_sys_width"]["min"],
                params["morph_params"]["shoot_sys_width"]["mean"],
                params["morph_params"]["shoot_sys_width"]["max"],
                params["morph_params"]["height"]["min"],
                params["morph_params"]["height"]["mean"],
                params["morph_params"]["height"]["max"],
                params["grow_params"]["abg_biomass"]["min"],
                params["grow_params"]["abg_biomass"]["mean"],
                params["grow_params"]["abg_biomass"]["max"],
            )
        else:
            if params["morph_params"]["allometry_method"] == "default":
                params["morph_params"]["basal_coeffs"] = empirical_coeffs["basal_coeffs"]
                params["morph_params"]["height_coeffs"] = empirical_coeffs["height_coeffs"]
                params["morph_params"]["canopy_coeffs"] = empirical_coeffs["canopy_coeffs"]
            (
                params["morph_params"]["basal_dia"]["min"],
                params["morph_params"]["height"]["min"],
                params["morph_params"]["shoot_sys_width"]["min"],
            ) = self.calc_abg_dims_from_biomass(
                params["grow_params"]["abg_biomass"]["min"]
            )
            (
                params["morph_params"]["basal_dia"]["max"],
                params["morph_params"]["height"]["max"],
                params["morph_params"]["shoot_sys_width"]["max"],
            ) = self.calc_abg_dims_from_biomass(
                params["grow_params"]["abg_biomass"]["max"]
            )
            (
                params["morph_params"]["basal_dia"]["mean"],
                params["morph_params"]["height"]["mean"],
                params["morph_params"]["shoot_sys_width"]["mean"],
            ) = self.calc_abg_dims_from_biomass(
                params["grow_params"]["abg_biomass"]["mean"]
            )
        return params

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # These dimensions are empirically derived allometric relationships for grasses.
        basal_dia_cm = height = shoot_sys_width = np.zeros_like(
            abg_biomass
        )
        filter = np.nonzero(abg_biomass > self.grow_params["abg_biomass"]["min"])
        basal_dia_cm[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], "basal_coeffs")
        basal_width = basal_dia_cm / 100
        shoot_sys_width[filter] = self._apply3_allometry_eq_for_ys(abg_biomass[filter], basal_dia_cm[filter], abg_biomass[filter], "canopy_coeffs")
        height[filter] = self._apply3_allometry_eq_for_ys(shoot_sys_width[filter], abg_biomass[filter], "height_coeffs")
        return basal_width, height, shoot_sys_width


class Tree(Habit):
    def __init__(self, params, dt, empirical_coeffs={}):
        green_parts = ("leaf")
        params = self._calc_derived_morph_params(params)
        super().__init__(params, dt, green_parts)


class Vine(Habit):
    def __init__(self, params, dt, empirical_coeffs={}):
        green_parts = ("leaf")
        params = self._calc_derived_morph_params(params)
        super().__init__(params, dt, green_parts)
