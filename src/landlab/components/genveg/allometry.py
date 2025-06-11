import numpy as np


class Biomass:
    def __init__(self, params, empirical_coeffs={}, cm=False):
        self.abg_biomass = params["grow_params"]["abg_biomass"]
        params["morph_params"]["empirical_coeffs"] = empirical_coeffs
        params["morph_params"]["empirical_coeffs"] = self._set_allometric_coeffs(params)
        self.cm = cm
        self.morph_params = params["morph_params"]

    def _set_allometric_coeffs(self, params):
        if params["morph_params"]["allometry_method"] != "min-max":
            return params["morph_params"]["empirical_coeffs"]
        else:
            xs_coeffs_list = [
                ["basal_dia", "basal_dia_coeffs"],
                ["height", "height_coeffs"],
                ["canopy_area", "canopy_area_coeffs"],
            ]
            for x in xs_coeffs_list:
                params["morph_params"]["empirical_coeffs"][x[1]] = {}
                x_min = params["morph_params"][x[0]]["min"]
                x_max = params["morph_params"][x[0]]["max"]
                a, b = self._calc2_allometry_coeffs(
                    x_min, x_max, self.abg_biomass["min"], self.abg_biomass["max"]
                )
                params["morph_params"]["empirical_coeffs"][x[1]]["a"] = a
                params["morph_params"]["empirical_coeffs"][x[1]]["b"] = b
            return params["morph_params"]["empirical_coeffs"]

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
        a = self.morph_params["empirical_coeffs"][predict_var]["a"]
        b = self.morph_params["empirical_coeffs"][predict_var]["b"]
        ln_xs = (ln_ys - a) / b
        return np.exp(ln_xs)

    def _apply2_allometry_eq_for_ys(self, xs, predict_var):
        # Assuming you want to solve for ys rather than xs
        ln_xs = np.log(xs)
        a = self.morph_params["empirical_coeffs"][predict_var]["a"]
        b = self.morph_params["empirical_coeffs"][predict_var]["b"]
        ln_ys = a + (b * ln_xs)
        return np.exp(ln_ys)

    def calc_abg_dims(self, abg_biomass, cm=False):
        # These dimensions are empirically derived allometric relationships.
        filter = np.nonzero(abg_biomass > self.abg_biomass["min"])
        basal_dia = np.zeros_like(abg_biomass)
        height = np.zeros_like(abg_biomass)
        canopy_area = np.zeros_like(abg_biomass)
        basal_dia[filter] = self._apply2_allometry_eq_for_xs(
            abg_biomass[filter], "basal_dia_coeffs"
        )
        if cm is True:
            basal_dia /= 100
        height[filter] = self._apply2_allometry_eq_for_xs(
            abg_biomass[filter], "height_coeffs"
        )
        canopy_area[filter] = self._apply2_allometry_eq_for_xs(
            abg_biomass[filter], "canopy_area_coeffs"
        )
        shoot_width = np.sqrt(4 * canopy_area / np.pi)
        return (basal_dia, shoot_width, height)

    def calc_abg_biomass_from_dim(self, dim, var_name, cm=False):
        abg_biomass = np.zeros_like(dim)
        filter = np.nonzero(dim > self.morph_params[var_name]["min"])
        if cm is True:
            dim *= 100
        abg_biomass[filter] = self._apply2_allometry_eq_for_ys(
            dim[filter], var_name + "_coeffs"
        )
        return abg_biomass


class Dimensional(Biomass):
    def __init__(self, params, empirical_coeffs, cm=True):
        # Formulation from LPJ-GUESS and BIOM-E are power formulations
        # This generally follows meta-analysis showing biomass is a function of abg volume
        # canopy area = a1 * basal diameter^b1 which is approx equiv to
        # log(canopy area) = a1 + b1 log(basal diameter)
        # height = a2 * basal diameter^b2
        # log(height) = a2 + b2 log(basal diameter)
        # cm paramater indicates allometric relationship assumes basal diameter is in cm
        # literature typically uses cm for basal diameter or diameter breast height
        super().__init__(params, empirical_coeffs, cm)

    def calc_abg_dims(self, abg_biomass, cm=True):
        basal_dia = np.zeros_like(abg_biomass)
        shoot_sys_width = np.zeros_like(abg_biomass)
        height = np.zeros_like(abg_biomass)
        filter = np.nonzero(abg_biomass > self.abg_biomass["min"])
        basal_dia[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], "basal_dia_coeffs")
        if cm is True:
            basal_dia /= 100
        shoot_sys_width[filter] = self._calc_shoot_width_from_basal_dia(
            basal_dia[filter], cm=cm
        )
        height[filter] = self._calc_height_from_basal_dia(basal_dia[filter], cm=cm)
        return (basal_dia, shoot_sys_width, height)

    def _calc_shoot_width_from_basal_dia(self, basal_dia, cm=False):
        canopy_area = np.zeros_like(basal_dia)
        filter = np.nonzero(basal_dia > self.morph_params["basal_dia"]["min"])
        if cm is True:
            basal_dia *= 100
        canopy_area[filter] = self._apply2_allometry_eq_for_ys(
            basal_dia[filter], "canopy_area_coeffs"
        )
        shoot_width = np.sqrt(4 * canopy_area / np.pi)
        return shoot_width

    def _calc_height_from_basal_dia(self, basal_dia, cm=False):
        height = np.zeros_like(basal_dia)
        filter = np.nonzero(basal_dia > self.morph_params["basal_dia"]["min"])
        if cm is True:
            basal_dia *= 100
        height[filter] = self._apply2_allometry_eq_for_ys(
            basal_dia[filter], "height_coeffs"
        )
        return height

    def _calc_basal_dia_from_shoot_width(self, shoot_width, cm=False):
        basal_dia = np.zeros_like(shoot_width)
        filter = np.nonzero(shoot_width > self.morph_params["shoot_sys_width"]["min"])
        canopy_area = 0.25 * np.pi * shoot_width**2
        basal_dia[filter] = self._apply2_allometry_eq_for_xs(
            canopy_area[filter], "canopy_area_coeffs"
        )
        if cm is True:
            basal_dia /= 100
        return basal_dia


class Multi_Dimensional(Biomass):
    def __init__(self, params):
        # Not implemented yet but saved here for future versions allowing
        # for more complex user-defined relationships
        super().__init__(params)

    def _apply3_allometry_eq_for_zs(self, xs, ys, predict_var):
        ln_xs = np.log(xs)
        ln_ys = np.log(ys)
        a = self.morph_params["empirical_coeffs"][predict_var]["a"]
        b = self.morph_params["empirical_coeffs"][predict_var]["b"]
        c = self.morph_params["empirical_coeffs"][predict_var]["c"]
        ln_zs = a + b * ln_xs + c * ln_ys
        return np.exp(ln_zs)

    def _apply3_allometry_eq_for_ys(self, xs, zs, predict_var):
        ln_xs = np.log(xs)
        ln_zs = np.log(zs)
        a = self.morph_params["empirical_coeffs"][predict_var]["a"]
        b = self.morph_params["empirical_coeffs"][predict_var]["b"]
        c = self.morph_params["empirical_coeffs"][predict_var]["c"]
        ln_ys = (ln_zs - a - b * ln_xs) / c
        return np.exp(ln_ys)

    def _apply3_allometry_eq_for_xs(self, ys, zs, predict_var):
        ln_ys = np.log(ys)
        ln_zs = np.log(zs)
        a = self.morph_params["empirical_coeffs"][predict_var]["a"]
        b = self.morph_params["empirical_coeffs"][predict_var]["b"]
        c = self.morph_params["empirical_coeffs"][predict_var]["c"]
        ln_xs = (ln_zs - a - c * ln_ys) / b
        return np.exp(ln_xs)
