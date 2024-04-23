"""
Species class definition, composition classes, and factory methods
to generate species classes.
These are used by PlantGrowth to differentiate plant properties
and processes for species.
"""

from .habit import Forbherb, Graminoid, Shrub, Tree, Vine
from .form import (
    Bunch,
    Colonizing,
    Multiplestems,
    Rhizomatous,
    Singlecrown,
    Singlestem,
    Stoloniferous,
    Thicketforming,
)
from .shape import (
    Climbing,
    Conical,
    Decumbent,
    Erect,
    Irregular,
    Oval,
    Prostrate,
    Rounded,
    Semierect,
    Vase,
)
from .photosynthesis import C3, C4, Cam
import numpy as np
from sympy import symbols, diff, lambdify, log

rng = np.random.default_rng()


# Define species class that inherits composite class methods
class Species(object):
    def __init__(self, species_params, latitude):
        self.all_parts = list(
            species_params["grow_params"]["glucose_requirement"].keys()
        )
        self.growth_parts = self.all_parts.copy()
        self.growth_parts.remove("reproductive")
        self.abg_parts = self.growth_parts.copy()
        self.abg_parts.remove("root")
        self.dead_parts = [
            "dead_root",
            "dead_leaf",
            "dead_stem",
            "dead_reproductive",
        ]
        self.dead_abg_parts = self.dead_parts.copy()
        self.dead_abg_parts.remove("dead_root")
        self.dead_abg_parts.remove("dead_reproductive")

        self.validate_plant_factors(species_params["plant_factors"])
        self.validate_duration_params(species_params["duration_params"])

        species_params = self.calculate_derived_params(species_params)

        self.species_plant_factors = species_params["plant_factors"]
        self.species_duration_params = species_params["duration_params"]
        self.species_grow_params = species_params["grow_params"]
        self.species_photo_params = species_params["photo_params"]
        self.species_dispersal_params = species_params["dispersal_params"]
        self.species_mort_params = species_params["mortality_params"]
        self.species_morph_params = species_params["morph_params"]

        self.populate_biomass_allocation_array()

        self.habit = self.select_habit_class()
        self.form = self.select_form_class()
        self.shape = self.select_shape_class()
        self.photosynthesis = self.select_photosythesis_type(latitude)

    def validate_plant_factors(self, plant_factors):
        plant_factor_options = {
            "species": [],
            "growth_habit": ["forb_herb", "graminoid", "shrub", "tree", "vine"],
            "monocot_dicot": ["monocot", "dicot"],
            "angio_gymno": ["angiosperm", "gymnosperm"],
            "duration": ["annual", "perennial deciduous", "perennial evergreen"],
            "growth_form": [
                "bunch",
                "colonizing",
                "multiple_stems",
                "rhizomatous",
                "single_crown",
                "single_stem",
                "stoloniferous",
                "thicket_forming",
            ],
            "shape": [
                "climbing",
                "columnar",
                "conical",
                "decumbent",
                "erect",
                "irregular",
                "oval",
                "prostrate",
                "rounded",
                "semi_erect",
                "vase",
            ],
            "p_type": ["C3", "C4"],
        }

        for key in plant_factors:
            try:
                opt_list = plant_factor_options[key]
                if opt_list:
                    if plant_factors[key] not in opt_list:
                        msg = "Invalid " + str(key) + " option"
                        raise ValueError(msg)
            except ValueError:
                print(
                    "Unexpected variable name in species parameter dictionary."
                    "Please check input parameter file"
                )

    def validate_duration_params(self, duration_params):
        if (duration_params["growing_season_start"] < 0) | (
            duration_params["growing_season_start"] > 366
        ):
            msg = "Growing season beginning must be integer values between 1-365"
            raise ValueError(msg)
        elif (
            duration_params["growing_season_end"]
            < duration_params["growing_season_start"]
        ) | (duration_params["growing_season_end"] > 366):
            msg = (
                "Growing season end must be between 1-365"
                "and greater than the growing season beginning"
            )
            raise ValueError(msg)
        elif (
            duration_params["senescence_start"]
            < duration_params["growing_season_start"]
        ) | (
            duration_params["senescence_start"] > duration_params["growing_season_end"]
        ):
            msg = "Start of senescence must be within the growing season"
            raise ValueError(msg)

    def calculate_derived_params(self, species_params):
        species_params["morph_params"]["max_crown_area"] = (
            np.pi / 4 * species_params["morph_params"]["max_shoot_sys_width"] ** 2
        )
        species_params["morph_params"]["min_crown_area"] = (
            np.pi / 4 * species_params["morph_params"]["min_shoot_sys_width"] ** 2
        )
        species_params["morph_params"]["max_root_area"] = (
            np.pi / 4 * species_params["morph_params"]["max_root_sys_width"] ** 2
        )
        species_params["morph_params"]["min_root_area"] = (
            np.pi / 4 * species_params["morph_params"]["min_root_sys_width"] ** 2
        )
        species_params["morph_params"]["max_vital_volume"] = (
            species_params["morph_params"]["max_crown_area"]
            * species_params["morph_params"]["max_height"]
        )
        species_params["morph_params"]["area_per_stem"] = (
            species_params["morph_params"]["max_crown_area"]
            / species_params["morph_params"]["max_n_stems"]
        )
        species_params["morph_params"]["min_abg_aspect_ratio"] = (
            species_params["morph_params"]["max_height"]
            / species_params["morph_params"]["min_shoot_sys_width"]
        )
        species_params["morph_params"]["max_abg_aspect_ratio"] = (
            species_params["morph_params"]["max_height"]
            / species_params["morph_params"]["max_shoot_sys_width"]
        )

        sum_vars = [
            ["max_total_biomass", "plant_part_max", self.all_parts],
            ["max_growth_biomass", "plant_part_max", self.growth_parts],
            ["max_abg_biomass", "plant_part_max", self.abg_parts],
            ["min_total_biomass", "plant_part_min", self.all_parts],
            ["min_growth_biomass", "plant_part_min", self.growth_parts],
            ["min_abg_biomass", "plant_part_min", self.abg_parts],
            [
                "min_nsc_biomass",
                "min_nsc_content",
                self.growth_parts,
            ],  # this is dynamic and varying with plant size
        ]
        for sum_var in sum_vars:
            species_params["grow_params"][sum_var[0]] = 0
            for part in sum_var[2]:
                species_params["grow_params"][sum_var[0]] += species_params[
                    "grow_params"
                ][sum_var[1]][part]

        species_params["morph_params"]["biomass_packing"] = (
            species_params["grow_params"]["max_growth_biomass"]
            / species_params["morph_params"]["max_vital_volume"]
        )

        species_params["duration_params"]["senesce_rate"] = 0.9 / (
            species_params["duration_params"]["growing_season_end"]
            - species_params["duration_params"]["senescence_start"]
        )

        seasonal_nsc_assim_rates = [
            "winter_nsc_rate",
            "spring_nsc_rate",
            "summer_nsc_rate",
            "fall_nsc_rate",
        ]
        seasonal_days = [
            ["growing_season_end", "growing_season_start"],
            ["growing_season_start", "peak_biomass"],
            ["peak_biomass", "senescence_start"],
            ["senescence_start", "growing_season_end"],
        ]
        species_params["duration_params"]["nsc_rate_change"] = {}
        for idx, season in enumerate(seasonal_nsc_assim_rates):
            end_date = species_params["duration_params"][seasonal_days[idx][1]]
            start_date = species_params["duration_params"][seasonal_days[idx][0]]
            species_params["duration_params"]["nsc_rate_change"][season] = {}
            if start_date > end_date:
                season_length = (365 - start_date) + end_date
            else:
                season_length = end_date - start_date
            for part in self.all_parts:
                rate_change_nsc = (
                    species_params["grow_params"]["incremental_nsc"][part][idx]
                    - species_params["grow_params"]["incremental_nsc"][part][idx - 1]
                ) / season_length  # figure out how to loop this
                species_params["duration_params"]["nsc_rate_change"][season][
                    part
                ] = rate_change_nsc

        return species_params

    def calculate_lai(self, leaf_area, shoot_sys_width):
        crown_area = self.shape.calc_crown_area_from_shoot_width(shoot_sys_width)
        lai = np.divide(
            leaf_area,
            crown_area,
            out=np.zeros_like(leaf_area),
            where=~np.isclose(crown_area, np.zeros_like(shoot_sys_width)),
        )
        return lai

    def select_photosythesis_type(self, latitude):
        photosynthesis_options = {"C3": C3, "C4": C4, "cam": Cam}
        return photosynthesis_options[self.species_plant_factors["p_type"]](
            latitude, photo_params=self.species_photo_params
        )

    def select_habit_class(self):
        habit = {
            "forb_herb": Forbherb,
            "graminoid": Graminoid,
            "shrub": Shrub,
            "tree": Tree,
            "vine": Vine,
        }
        return habit[self.species_plant_factors["growth_habit"]](
            self.species_grow_params,
            self.species_duration_params,
            self.species_plant_factors["duration"],
        )

    def select_form_class(self):
        form = {
            "bunch": Bunch(self.species_dispersal_params, self.species_grow_params),
            "colonizing": Colonizing(
                self.species_dispersal_params, self.species_grow_params
            ),
            "multiple_stems": Multiplestems(
                self.species_dispersal_params, self.species_grow_params
            ),
            "rhizomatous": Rhizomatous(
                self.species_dispersal_params, self.species_grow_params
            ),
            "single_crown": Singlecrown(
                self.species_dispersal_params, self.species_grow_params
            ),
            "single_stem": Singlestem(
                self.species_dispersal_params, self.species_grow_params
            ),
            "stoloniferous": Stoloniferous(
                self.species_dispersal_params, self.species_grow_params
            ),
            "thicket_forming": Thicketforming(
                self.species_dispersal_params, self.species_grow_params
            ),
        }
        return form[self.species_plant_factors["growth_form"]]

    def select_shape_class(self):
        shape = {
            "climbing": Climbing(self.species_morph_params, self.species_grow_params),
            "conical": Conical(self.species_morph_params, self.species_grow_params),
            "decumbent": Decumbent(self.species_morph_params, self.species_grow_params),
            "erect": Erect(self.species_morph_params, self.species_grow_params),
            "irregular": Irregular(self.species_morph_params, self.species_grow_params),
            "oval": Oval(self.species_morph_params, self.species_grow_params),
            "prostrate": Prostrate(self.species_morph_params, self.species_grow_params),
            "rounded": Rounded(self.species_morph_params, self.species_grow_params),
            "semi_erect": Semierect(
                self.species_morph_params, self.species_grow_params
            ),
            "vase": Vase(self.species_morph_params, self.species_grow_params),
        }
        return shape[self.species_plant_factors["shape"]]

    def populate_biomass_allocation_array(self):
        root2leaf = self.species_grow_params["root_to_leaf"]
        root2stem = self.species_grow_params["root_to_stem"]
        prior_root_biomass = np.arange(
            start=self.species_grow_params["plant_part_min"]["root"],
            stop=self.species_grow_params["plant_part_max"]["root"] + 0.1,
            step=0.1,
        )
        length_of_array = len(prior_root_biomass)
        place_zeros = np.zeros(length_of_array)
        biomass_allocation_map = np.column_stack(
            (
                prior_root_biomass,
                place_zeros,
                place_zeros,
                place_zeros,
                place_zeros,
                place_zeros,
                place_zeros,
            )
        )
        biomass_allocation_map = list(map(tuple, biomass_allocation_map))
        self.biomass_allocation_array = np.array(
            biomass_allocation_map,
            dtype=[
                ("prior_root_biomass", float),
                ("total_biomass", float),
                ("delta_leaf_unit_root", float),
                ("delta_stem_unit_root", float),
                ("leaf_mass_frac", float),
                ("stem_mass_frac", float),
                ("abg_biomass", float),
            ],
        )

        # set up sympy equations
        rootsym = symbols("rootsym")
        dleaf = diff(
            10
            ** (
                root2leaf["a"]
                + root2leaf["b1"] * log(rootsym, 10)
                + root2leaf["b2"] * (log(rootsym, 10)) ** 2
            ),
            rootsym,
        )
        dstem = diff(
            10
            ** (
                root2stem["a"]
                + root2stem["b1"] * log(rootsym, 10)
                + root2stem["b2"] * (log(rootsym, 10)) ** 2
            ),
            rootsym,
        )
        # Generate numpy expressions and solve for rate change
        # in leaf and stem biomass per unit mass of root
        fleaf = lambdify(rootsym, dleaf, "numpy")
        fstem = lambdify(rootsym, dstem, "numpy")
        self.biomass_allocation_array["delta_leaf_unit_root"] = fleaf(
            self.biomass_allocation_array["prior_root_biomass"]
        )
        self.biomass_allocation_array["delta_stem_unit_root"] = fstem(
            self.biomass_allocation_array["prior_root_biomass"]
        )
        _leaf_biomasss = 10 ** (
            root2leaf["a"]
            + root2leaf["b1"] * np.log10(prior_root_biomass)
            + root2leaf["b2"] * (np.log10(prior_root_biomass)) ** 2
        )
        _stem_biomass = 10 ** (
            root2stem["a"]
            + root2stem["b1"] * np.log10(prior_root_biomass)
            + root2stem["b2"] * (np.log10(prior_root_biomass)) ** 2
        )
        self.biomass_allocation_array["total_biomass"] = (
            self.biomass_allocation_array["prior_root_biomass"]
            + _leaf_biomasss
            + _stem_biomass
        )
        self.biomass_allocation_array["leaf_mass_frac"] = (
            _leaf_biomasss / self.biomass_allocation_array["total_biomass"]
        )
        self.biomass_allocation_array["stem_mass_frac"] = (
            _stem_biomass / self.biomass_allocation_array["total_biomass"]
        )
        self.biomass_allocation_array["abg_biomass"] = _leaf_biomasss + _stem_biomass

    def allocate_biomass_dynamically(self, _live_biomass, delta_tot):
        ###
        # This method allocates new biomass according to the size-dependent
        # biomass allocation array calculated upon initiation of the PlantGrowth class.
        # The array is only valid for actively growing plants so this method is only
        # used during the growing season.
        # After initial allocation, storage redistribution, reallocation, and the
        # minimum size check methods are called to adjust the biomass in each part
        # before saving.
        # Required parameters are the new net biomass generated for the day
        ###

        _new_biomass = _live_biomass
        growth_biomass = self.sum_plant_parts(_live_biomass, parts="growth")

        # Interpolate values from biomass allocation array
        delta_leaf_unit_root = np.interp(
            _live_biomass["root_biomass"],
            self.biomass_allocation_array["prior_root_biomass"],
            self.biomass_allocation_array["delta_leaf_unit_root"],
        )
        delta_stem_unit_root = np.interp(
            _live_biomass["root_biomass"],
            self.biomass_allocation_array["prior_root_biomass"],
            self.biomass_allocation_array["delta_stem_unit_root"],
        )
        filter = np.nonzero(delta_tot > 0)
        frac_to_growth = np.ones_like(_live_biomass["root"])
        frac_to_repro = np.zeros_like(_live_biomass["root"])
        mass_ratio = growth_biomass / self.species_grow_params["max_growth_biomass"]
        frac_to_growth[filter] = 1 / (1 + 0.002 * np.exp(9.50 * mass_ratio[filter]))
        frac_to_repro[filter] = 1 - frac_to_growth[filter]
        _new_biomass["reproductive"] += (
            delta_tot
            * frac_to_repro
            / self.species_grow_params["glucose_requirement"]["reproductive"]
        )

        # Calculate allocation
        _glu_req_sum = np.zeros_like(_new_biomass["root"])
        for part in self.growth_parts:
            _glu_req_sum = (
                self.species_grow_params["glucose_requirement"][part]
                * _new_biomass[part]
                / growth_biomass
            )
        filter = np.nonzero(_glu_req_sum > 0)
        delta_tot_growth = np.zeros_like(_new_biomass["root"])
        delta_tot_growth[filter] = (
            delta_tot[filter] * frac_to_growth[filter] / _glu_req_sum[filter]
        )

        root = delta_tot_growth / (1 + delta_leaf_unit_root + delta_stem_unit_root)
        delta = {
            "root": root,
            "leaf": delta_leaf_unit_root * root,
            "stem": delta_stem_unit_root * root,
        }
        for part in self.growth_parts:
            # update biomass in plant array
            _new_biomass[part] = _live_biomass[part] + delta[part]
        # Adjust biomass allocation among storage and growth parts
        _new_biomass = self.adjust_biomass_allocation_towards_ideal(_new_biomass)
        return _new_biomass

    def branch(self):
        self.form.branch()

    def disperse(self, plants, jday):
        # decide how to parameterize reproductive schedule, make repro event
        plants = self.form.disperse(plants)
        return plants

    def enter_dormancy(
        self, plants, jday
    ):  # calculate sum of green parts and sum of persistant parts
        end_plants = plants.copy()
        plants = self.habit.enter_dormancy(plants)
        plants = self.update_dead_biomass(plants, end_plants)
        return plants

    def emerge(self, plants, jday):
        ns_conc = self.get_daily_nsc_concentration(jday)
        available_stored_biomass = np.zeros_like(plants["root"])
        total_persistent_biomass = np.zeros_like(plants["root"])
        for part in self.habit.duration.persistent_parts:
            avail_nsc_content = (
                ns_conc[part] - self.species_grow_params["min_nsc_content"][part]
            ) * np.ones_like(plants[part])
            avail_nsc_content[avail_nsc_content < 0] = 0.0

            available_stored_biomass += plants[part] * avail_nsc_content
            total_persistent_biomass += plants[part]
        plants = self.habit.duration.emerge(
            plants, available_stored_biomass, total_persistent_biomass
        )
        plants = self.update_morphology(plants)

        return plants

    def get_daily_nsc_concentration(self, _current_jday):
        d = _current_jday
        days = self.species_duration_params
        rate = self.species_duration_params["nsc_rate_change"]
        nsc_content = {}
        day_conditions = [
            d < days["growing_season_start"],
            (d >= days["growing_season_start"]) & (d < days["peak_biomass"]),
            (d >= days["peak_biomass"]) & (d < days["senescence_start"]),
            (d >= days["senescence_start"]) & (d < days["growing_season_end"]),
            d >= days["growing_season_end"],
        ]

        for part in self.all_parts:
            nsc_content_opts_b = [
                (
                    self.species_grow_params["incremental_nsc"][part][3]
                    + (self.species_grow_params["nsc_content"][part] * 1000) ** 0.5
                )
                ** 2,
                (
                    self.species_grow_params["incremental_nsc"][part][0]
                    + (self.species_grow_params["nsc_content"][part] * 1000) ** 0.5
                )
                ** 2,
                (
                    self.species_grow_params["incremental_nsc"][part][1]
                    + (self.species_grow_params["nsc_content"][part] * 1000) ** 0.5
                )
                ** 2,
                (
                    self.species_grow_params["incremental_nsc"][part][2]
                    + (self.species_grow_params["nsc_content"][part] * 1000) ** 0.5
                )
                ** 2,
                (
                    self.species_grow_params["incremental_nsc"][part][3]
                    + (self.species_grow_params["nsc_content"][part] * 1000) ** 0.5
                )
                ** 2,
            ]

            nsc_content_opts_mx = [
                rate["winter_nsc_rate"][part] * (d + 365 - days["growing_season_end"]),
                rate["spring_nsc_rate"][part] * (d - days["growing_season_start"]),
                rate["summer_nsc_rate"][part] * (d - days["peak_biomass"]),
                rate["fall_nsc_rate"][part] * (d - days["senescence_start"]),
                rate["winter_nsc_rate"][part] * (d - days["growing_season_end"]),
            ]
            nsc_content[part] = (
                (np.select(day_conditions, nsc_content_opts_b)) ** 0.5
                + np.select(day_conditions, nsc_content_opts_mx)
            ) ** 2 / 1000
        return nsc_content

    def litter_decomp(self, _new_biomass):
        decay_rate = self.species_morph_params["biomass_decay_rate"]
        for part in self.dead_parts:
            filter = np.nonzero(_new_biomass[part] > 0.0)
            part_init_mass = np.zeros_like(_new_biomass[part])
            part_init_mass[filter] = _new_biomass[part][filter] / np.exp(
                -decay_rate * _new_biomass[str(part) + "_age"][filter]
            )
            _new_biomass[part] = part_init_mass * np.exp(
                -decay_rate * (_new_biomass[str(part) + "_age"] + self.dt.astype(float))
            )
            _new_biomass[str(part) + "_age"] += self.dt.astype(float) * np.ones_like(
                _new_biomass[part]
            )
        for part in self.dead_parts:
            filter = np.nonzero(
                np.isnan(_new_biomass[part])
                | (_new_biomass[part] < 0)
                | np.isinf(_new_biomass[part])
            )
            _new_biomass[part][filter] = np.zeros_like(_new_biomass[part][filter])
        return _new_biomass

    def mortality(self, plants, _in_growing_season):
        old_biomass = plants.copy()
        plants = self.calculate_whole_plant_mortality(plants, _in_growing_season)
        plants = self.calculate_shaded_leaf_mortality(plants)
        plants = self.update_dead_biomass(plants, old_biomass)
        return plants

    def calculate_whole_plant_mortality(self, plants, _in_growing_season):
        mortdict = self.species_mort_params
        # set flags for three types of mortality periods
        mort_period_bool = {
            "during growing season": _in_growing_season is True,
            "during dormant season": _in_growing_season is False,
            "year-round": True,
        }
        factors = mortdict["mort_variable_name"]
        for fact in factors:
            # Determine if mortality factor is applied
            run_mort = mort_period_bool[mortdict["period"][fact]]
            if not run_mort:
                continue
            else:
                try:
                    # Assign mortality predictor from grid to plant
                    pred = self._grid["cell"][factors[fact]][plants["cell_index"]]
                    coeffs = mortdict["coeffs"][fact]
                    # Calculate the probability of survival and cap from 0-1
                    prob_survival = 1 / (1 + coeffs[0] * np.exp(-coeffs[1] * pred))
                    prob_survival[np.isnan(prob_survival)] = 1.0
                    prob_survival[prob_survival < 0] = 0
                    prob_survival_daily = prob_survival ** (
                        1 / (mortdict["duration"][fact] / self.dt.astype(int))
                    )
                    daily_survival = prob_survival_daily > rng.random(pred.shape)
                    for part in self.all_parts:
                        plants["dead_" + str(part)] = plants["dead_" + str(part)] + (
                            plants[part] * (np.invert(daily_survival).astype(int))
                        )
                        plants[part] = plants[part] * daily_survival.astype(int)

                except KeyError:
                    msg = f"No data available for mortality factor {factors[fact]}"
                    raise ValueError(msg)
        return plants

    def calculate_shaded_leaf_mortality(self, plants):
        # Based on Teh code equation 7.18 and 7.20 (pg. 154)
        lai = self.calculate_lai(plants["total_leaf_area"], plants["shoot_sys_width"])

        excess_lai = (
            lai - self.species_morph_params["lai_cr"]
        ) / self.species_morph_params["lai_cr"]
        shaded_leaf = np.nonzero((excess_lai > 0) & (plants["leaf"] > 0))
        D_shade = np.zeros_like(plants["total_leaf_area"])
        D_shade[shaded_leaf] = 0.03 * excess_lai[shaded_leaf]
        D_shade[D_shade > 0.03] = 0.03
        leaf_loss = plants["leaf_biomass"] * D_shade
        plants["dead_leaf"][shaded_leaf] += leaf_loss[shaded_leaf]
        plants["leaf"][shaded_leaf] -= leaf_loss[shaded_leaf]
        return plants

    def photosynthesize(
        self,
        _par,
        _min_temperature,
        _max_temperature,
        cell_lai,
        _last_biomass,
        _current_day,
    ):
        ind_lai = self.calculate_lai(
            _last_biomass["total_leaf_area"], _last_biomass["shoot_sys_width"]
        )
        lai = np.maximum(cell_lai, ind_lai)
        carb_generated = self.photosynthesis.photosynthesize(
            _par,
            _min_temperature,
            _max_temperature,
            lai,
            _last_biomass,
            _current_day,
        )
        random_water_stress = rng.normal(loc=0.65, scale=0.12, size=carb_generated.size)
        random_water_stress[random_water_stress < 0] = 0.0
        random_water_stress[random_water_stress > 1] = 1.0
        adj_carb_generated = random_water_stress * carb_generated
        return adj_carb_generated

    def respire(self, _min_temperature, _max_temperature, _last_biomass):
        _temperature = (_min_temperature + _max_temperature) / 2
        growdict = self.species_grow_params
        _new_biomass = _last_biomass.copy()
        # respiration coefficient temp dependence from Teh 2006
        temp_adj = 2 ** ((_temperature - 25) / 10)
        for part in self.all_parts:
            delta_respire = np.zeros_like(_last_biomass["root"])
            filter = np.nonzero(_last_biomass[part] > 0)
            delta_respire[filter] = (
                temp_adj[filter]
                * growdict["respiration_coefficient"][part]
                * _last_biomass[part][filter]
            ) / growdict["glucose_requirement"][part]
            _new_biomass[part][filter] -= delta_respire[filter]

        return _new_biomass

    def senesce(self, plants, jday):
        ns_conc = self.get_daily_nsc_concentration(jday)
        ns_green_mass = 0.0
        for part in self.habit.duration.green_parts:
            ns_green_mass += plants[part] * ns_conc[part]
        persistent_mass = self.sum_plant_parts(plants, parts="persistent")
        plants = self.habit.senesce(
            plants,
            ns_green_mass=ns_green_mass,
            persistent_mass=persistent_mass,
        )
        return plants

    def set_initial_biomass(self, plants, in_growing_season):
        growdict = self.species_grow_params
        morphdict = self.species_morph_params
        plants["repro_biomass"] = (
            growdict["plant_part_min"]["reproductive"]
            + rng.rayleigh(scale=0.2, size=plants.size)
            * growdict["plant_part_max"]["reproductive"]
        )
        crown_area = self.shape.calc_crown_area_from_shoot_width(
            plants["shoot_sys_width"]
        )
        plants["shoot_sys_height"] = self.habit.set_initial_height(
            morphdict["min_height"], morphdict["max_height"], crown_area.size
        )
        log_vital_volume = np.log10(crown_area * plants["shoot_sys_height"])
        log_abg_biomass_ideal = (
            log_vital_volume / np.log10(morphdict["max_vital_volume"])
        ) * np.log10(growdict["max_abg_biomass"] / 1000)
        total_biomass = np.interp(
            ((10**log_abg_biomass_ideal) * 1000),
            self.biomass_allocation_array["abg_biomass"],
            self.biomass_allocation_array["total_biomass"],
        )
        (
            plants["root_biomass"],
            plants["leaf_biomass"],
            plants["stem_biomass"],
        ) = self.habit.duration._solve_biomass_allocation(total_biomass)
        plants["root_sys_width"] = self.shape.calc_root_sys_width(
            plants["shoot_sys_width"], plants["shoot_sys_height"]
        )
        plants["n_stems"] = self.form.set_initial_branches(
            morphdict["max_n_stems"], crown_area.size
        )
        plants = self.habit.duration.set_initial_biomass(plants, in_growing_season)
        return plants

    def set_new_biomass(self, plants):
        plants = self.habit.duration.set_new_biomass(plants)
        return plants

    def update_morphology(self, plants):
        abg_biomass = self.sum_plant_parts(plants, parts="aboveground")
        dead_abg_biomass = self.sum_plant_parts(plants, parts="dead_aboveground")
        total_abg_biomass = abg_biomass + dead_abg_biomass
        dims = self.shape.calc_abg_dims_from_biomass(total_abg_biomass)
        plants["shoot_sys_width"] = dims[0]
        plants["shoot_sys_height"] = dims[1]
        plants["root_sys_width"] = self.shape.calc_root_sys_width(
            plants["shoot_sys_width"]
        )
        dead_leaf_area = plants["total_leaf_area"] - plants["live_leaf_area"]
        dead_leaf_area[dead_leaf_area < 0] = 0.0
        filter = np.nonzero(dead_leaf_area > 0)
        plants["live_leaf_area"] = (
            plants["leaf"] * self.species_morph_params["sp_leaf_area"]
        )

        cohort_dead_leaf_area = np.zeros_like(dead_leaf_area)
        cohort_dead_leaf_area[filter] = dead_leaf_area[filter] / np.exp(
            -self.species_morph_params["biomass_decay_rate"]
            * plants["dead_leaf_age"][filter]
        )
        dead_leaf_area_ratio = np.ones_like(plants["live_leaf_area"])
        dead_leaf_area_ratio[filter] = (
            dead_leaf_area[filter] / cohort_dead_leaf_area[filter]
        )
        plants["total_leaf_area"] = plants["live_leaf_area"] + (
            plants["dead_leaf"]
            * self.species_morph_params["sp_leaf_area"]
            * dead_leaf_area_ratio
            / 3
        )
        return plants

    def update_dead_biomass(self, _new_biomass, old_biomass):
        for part in self.all_parts:
            part_biomass_change = _new_biomass[part] - old_biomass[part]
            filter = np.nonzero(part_biomass_change < 0.0)
            _new_biomass["dead_" + str(part)][filter] -= part_biomass_change[filter]
            _new_biomass["dead_" + str(part) + "_age"][filter] = 0.0
        return _new_biomass

    def sum_plant_parts(self, _new_biomass, parts="total"):
        parts_choices = {
            "total": self.all_parts,
            "growth": self.growth_parts,
            "aboveground": self.abg_parts,
            "persistent": self.habit.duration.persistent_parts,
            "green": self.habit.duration.green_parts,
            "dead": self.dead_parts,
            "dead_aboveground": self.dead_abg_parts,
        }

        parts_dict = parts_choices[parts]
        _new_tot = np.zeros_like(_new_biomass["root_biomass"])
        for part in parts_dict:
            _new_tot += _new_biomass[part]
        return _new_tot

    def calculate_dead_age(self, age_t1, mass_t1, mass_t2):
        age_t2 = np.zeros_like(age_t1)
        filter = np.nonzero(mass_t2 > 0)
        age_t2[filter] = (
            (age_t1[filter] * mass_t1[filter])
            + ((mass_t2[filter] - mass_t1[filter]) * np.zeros_like(age_t1[filter]))
        ) / (mass_t2[filter])
        return age_t2
