import numpy as np

# Growth form classes and selection method


class PlantShape(object):
    def __init__(self, morph_params, grow_params):
        self.morph_params = morph_params
        self.grow_params = grow_params
        self.min_crown_area = self.calc_crown_area_from_shoot_width(
            self.morph_params["min_shoot_sys_width"]
        )
        self.max_crown_area = self.calc_crown_area_from_shoot_width(
            self.morph_params["max_shoot_sys_width"]
        )

    def calc_root_sys_width(self, shoot_sys_width, shoot_sys_height=1):
        volume = self.calc_crown_volume(shoot_sys_width, shoot_sys_height)
        root_sys_width = 0.08 + 0.24 * volume
        root_sys_width[root_sys_width > self.morph_params["max_root_sys_width"]] = (
            self.morph_params["max_root_sys_width"]
        )
        root_sys_width[root_sys_width < self.morph_params["min_root_sys_width"]] = (
            self.morph_params["max_root_sys_width"]
        )
        return root_sys_width

    def calc_crown_area_from_shoot_width(self, shoot_sys_width):
        crown_area = 0.25 * np.pi * shoot_sys_width**2
        return crown_area


    def calc_crown_volume(self, shoot_sys_width, shoot_sys_height):
        volume = np.pi / 4 * shoot_sys_width**2 * shoot_sys_height
        return volume

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # These dimensions are empirically derived allometric relationships for grasses. This should be moved into habit and the shape classes deleted.
        log_basal_width = crown_area = basal_width = plant_height = np.zeros_like(
            abg_biomass
        )
        filter = np.nonzero(abg_biomass > self.grow_params["min_abg_biomass"])
        log_basal_width[filter] = (
            np.log(abg_biomass[filter]) - self.basal_width_coeffs["a"]
        ) / self.basal_width_coeffs["b"]
        basal_width = np.exp(log_basal_width)
        height = self.height_coeffs["a"] * basal_width ** self.height_coeffs["b"]
        crown_area[filter] = (
            self.crown_coeffs["a"]
            * basal_width ** self.crown_coeffs["b"]
            * height ** self.crown_coeffs["c"]
        )
        shoot_sys_width = (4 * crown_area / np.pi) ** 0.5
        return basal_width, shoot_sys_width, plant_height


class Climbing(PlantShape):
    def __init__(self, morph_params, grow_params):
        pass


class Conical(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)

    def calc_crown_volume(self, shoot_sys_width, shoot_sys_height):
        volume = np.pi / 12 * shoot_sys_width**2 * shoot_sys_height
        return volume


class Decumbent(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)

    def calc_crown_volume(self, shoot_sys_width, shoot_sys_height):
        volume = np.pi / 3 * shoot_sys_width**2 * shoot_sys_height
        return volume


class Erect(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)

    def abg_biomass_transform(self, abg_biomass):
        return np.log10(abg_biomass / 1000)

    def calc_root_sys_width(self, shoot_sys_width, root_sys_width=np.nan):
        return shoot_sys_width


class Irregular(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)


class Oval(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)

    def calc_crown_volume(self, shoot_sys_width, shoot_sys_height):
        volume = np.pi / 6 * shoot_sys_width**2 * shoot_sys_height
        return volume


class Prostrate(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)


class Rounded(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)

    def calc_crown_volume(self, shoot_sys_width, shoot_sys_height):
        volume = np.pi / 6 * shoot_sys_width**3
        return volume


class Semierect(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)


class Vase(PlantShape):
    def __init__(self, morph_params, grow_params):
        super().__init__(morph_params, grow_params)
