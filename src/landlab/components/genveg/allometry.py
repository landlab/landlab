import numpy as np


class Biomass:
    def __init__(self, params):
        self.abg_biomass = params["grow_params"]["abg_biomass"]
        xs_coeffs_list = [
            ["basal_dia","basal_coeffs"],
            ["height","height_coeffs"],
            ["shoot_sys_width","canopy_coeffs"],
        ]
        for x in xs_coeffs_list:
            x_min = params["morph_params"][x[0]]["min"]
            x_max = params["morph_params"][x[0]]["max"]
            a, b = self._calc2_allometry_coeffs(x_min, x_max, self.abg_biomass["min"], self.abg_biomass["max"])
            params["morph_params"][x[1]]["a"] = a
            params["morph_params"][x[1]]["b"] = b
        self.morph_params = params["morph_params"]


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

    def _calc_abg_dim_from_biomass(self, abg_biomass, variable, cm=False):
        # These dimensions are empirically derived allometric relationships.
        dimension = np.zeros_like(abg_biomass)
        filter = np.nonzero(abg_biomass > self.abg_biomass["min"])
        dimension[filter] = self._apply2_allometry_eq_for_xs(abg_biomass[filter], variable)
        if cm is True:
            return dimension / 100
        else:
            return dimension

    def _calc_abg_biomass_from_dim(self, dim, variable, cm=False):
        abg_biomass = np.zeros_like(dim)
        filter = np.nonzero(dim > self.morph_params[variable]["min"])
        if cm is True:
            dim *= 100
        abg_biomass[filter] = self._apply2_allometry_eq_for_ys(dim[filter], variable)
        return abg_biomass


class Dimensional(Biomass):
    def __init__(self, params):
        # Formulation from LPJ-GUESS and BIOM-E are power formulations
        # This generally follows meta-analysis showing biomass is a function of abg volume
        # crown area = a1 * basal diameter^b1 which is approx equiv to
        # log(crown area) = a1 + b1 log(basal diameter)
        # height = a2 * basal diameter^b2
        # log(height) = a2 + b2 log(basal diameter)
        # cm paramater indicates allometric relationship assumes basal diameter is in cm
        super.__init__(params)

    def _calc_shoot_width_from_basal_dia(self, basal_dia, cm=False):
        canopy_area = np.zeros_like(basal_dia)
        filter = np.nonzero(basal_dia > self.morph_params["basal_dia"]["min"])
        if cm is True:
            basal_dia *= 100
        canopy_area[filter] = self._apply2_allometry_eq_for_ys(basal_dia[filter], "shoot_sys_width_coeffs")
        shoot_width = np.sqrt(4 * canopy_area / np.pi)
        return shoot_width

    def _calc_height_from_basal_dia(self, basal_dia, cm=False):
        height = np.zeros_like(basal_dia)
        filter = np.nonzero(basal_dia > self.morph_params["basal_dia"]["min"])
        if cm is True:
            basal_dia *= 100
        height[filter] = self._apply2_allometry_eq_for_ys(basal_dia[filter], "height_coeffs")
        return height

    def _calc_basal_dia_from_shoot_width(self, shoot_width, cm=False):
        basal_dia = np.zeros_like(shoot_width)
        filter = np.nonzero(shoot_width > self.morph_params["shoot_sys_width"]["min"])
        basal_dia[filter] = self._apply2_allometry_eq_for_xs(shoot_width[filter], "shoot_sys_width_coeffs")
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