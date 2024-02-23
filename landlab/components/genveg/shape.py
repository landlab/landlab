import numpy as np
from scipy import interpolate

# Growth form classes and selection method

# how do these composition classes need to relate to each other? we need to use properties from one composition class to in methods of another.
# Need to distinguish live plant size from dead. How should we handle this? Make a dead flag? Make a dead array? We have plants with live and dead parts.


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
        # Calculate log10 linear equation to calculate relationship between aboveground biomass and aspect ratio
        # This is not working as expected
        width_x0 = self.abg_biomass_transform(self.grow_params["min_abg_biomass"])
        width_x1 = self.abg_biomass_transform(self.grow_params["max_abg_biomass"])
        width_y0 = np.log10(self.min_crown_area)
        width_y1 = np.log10(self.max_crown_area)
        width_m = (width_y1 - width_y0) / (width_x1 - width_x0)
        width_b = width_y0 - (width_m * width_x0)
        self.crown_area_coeffs = {"m": width_m, "b": width_b}
        height_x0 = width_x0
        height_x1 = width_x1
        height_y0 = np.log10(self.morph_params["min_height"])
        height_y1 = np.log10(self.morph_params["max_height"])
        height_m = (height_y1 - height_y0) / (height_x1 - height_x0)
        height_b = height_y0 - height_m * height_x0
        self.height_coeffs = {"m": height_m, "b": height_b}

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

    def calc_vital_volume_from_biomass(self, abg_biomass):
        log_vital_volume = (
            np.log10(abg_biomass / 1000)
            / np.log10(self.grow_params["max_abg_biomass"] / 1000)
        ) * np.log10(self.morph_params["max_vital_volume"])
        return 10**log_vital_volume

    def abg_biomass_transform(self, abg_biomass):
        return np.log10(abg_biomass / 1000)

    def calc_crown_volume(self, shoot_sys_width, shoot_sys_height):
        volume = np.pi / 4 * shoot_sys_width**2 * shoot_sys_height
        return volume

    def calc_abg_dims_from_biomass(self, abg_biomass):
        # interpolation function not working as expected
        # shoot sys width(t) should be dependent on shoot_sys_width (t-1)
        log_crown_area = vital_volume = log_plant_height = crown_area = plant_height = (
            np.zeros_like(abg_biomass)
        )
        filter = np.nonzero(abg_biomass > self.grow_params["min_abg_biomass"])
        # log_aspect_ratio[filter]=self.aspect_ratio_interp_func(np.log10(abg_biomass[filter]/1000))
        log_crown_area[filter] = self.crown_area_coeffs["b"] + (
            self.crown_area_coeffs["m"]
            * self.abg_biomass_transform(abg_biomass[filter])
        )
        crown_area[filter] = 10 ** log_crown_area[filter]
        shoot_sys_width = (4 * crown_area / np.pi) ** 0.5
        vital_volume[filter] = self.calc_vital_volume_from_biomass(abg_biomass[filter])
        # shoot_sys_width = ((4 * vital_volume) / (np.pi * aspect_ratio)) ** (1 / 3)
        log_plant_height[filter] = self.height_coeffs["b"] + self.height_coeffs[
            "m"
        ] * self.abg_biomass_transform(abg_biomass[filter])
        plant_height[filter] = 10 ** log_plant_height[filter]
        return shoot_sys_width, plant_height


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
