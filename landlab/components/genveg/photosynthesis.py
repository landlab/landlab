import numpy as np


class Photosynthesis(object):
    def __init__(self):
        pass

    def photosynthesize(
        self, grid_par_W_per_sqm, growdict, last_biomass, _glu_req, daylength
    ):
        print("I am photosynthesizing during the growing season")
        par_micromol_per_sqm_s = (grid_par_W_per_sqm) * 4.6
        # from Charisma instructions: tells how much of the light a plant is going to get as PAR in microeinsteins based on how many leaves are on the plant
        intSolarRad = par_micromol_per_sqm_s * np.exp(
            -(growdict["k_light_extinct"]) * last_biomass["leaf_biomass"]
        )
        # amount of light absorbed, per half saturaion constants from Charisma eq. 3. the monod or michaelis/menten function is adequate for describing the photosynthetic response to light
        intLightpH = intSolarRad / (intSolarRad + growdict["light_half_sat"])
        # pMax is the maximum rate of photosynthesis, species specific
        photosynthesis = (growdict["p_max"]) * intLightpH
        # calculates total biomass gained across plant (twlvg is amount of leaver/green matter): you feed the model total biomass and then from that we determine how much leaf mass there is and so then basically an average of how much that average leaf will produce multiplied by the number of leaves, this is assuming that all leaves are mature
        dtgaCollapsed = photosynthesis * last_biomass["leaf_biomass"]
        # total biomass for day length
        assimilatedCH2O = dtgaCollapsed * daylength
        # converts carbohydrates to glucose where photosynthesis unit is glucose and then we later convert that glucose to biomass in another section
        gphot = assimilatedCH2O * (30 / 44)
        delta_tot = np.zeros_like(_glu_req)
        delta_tot[_glu_req != 0] = gphot[_glu_req != 0] / _glu_req[_glu_req != 0]
        return delta_tot

        def estimate_PAR(self, jday, ktp=0.78, rel_humidity=0.5, albedo=0.25):
            log_mod_coeffs_low = (2.0196, -5.6485, 1.3469, 0.7309, 0.3045)
            log_mod_coeffs_high = (1.2438, -2.3335, 0.7046, 0.4107, -1.9484)
            pass


class C3(Photosynthesis):
    def __init__(self):
        super().__init__()


class C4(Photosynthesis):
    def __init__(self):
        super().__init__()


class Cam(Photosynthesis):
    def __init__(self):
        super().__init__()
