"""
Species class definition, composition classes, and factory methods to generate species classes. 
These are used by PlantGrowth to differentiate plant properties and processes for species.
"""
from .habit import *
from .form import *
from .shape import *
from .photosynthesis import *
import numpy as np
from sympy import symbols, diff, lambdify, log


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
            except:
                msg = "Unexpected variable name in species parameter dictionary. Please check input parameter file."
                raise ValueError(msg)

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
            msg = "Growing season end must be between 1-365 and greater than the growing season beginning"
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

    def calculate_lai(self, leaf_biomass, shoot_sys_width):
        return (
            self.species_morph_params["sp_leaf_area"]
            * leaf_biomass
            / (100 * 0.25 * np.pi * shoot_sys_width**2)
        )

    def select_photosythesis_type(self, latitude):
        photosynthesis_options = {"C3": C3, "C4": C4, "cam": Cam}
        return photosynthesis_options[self.species_plant_factors["p_type"]](latitude)

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
        # Generate numpy expressions and solve for rate change in leaf and stem biomass per unit mass of root
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

    def branch(self):
        self.form.branch()

    def disperse(self, plants, jday):
        # decide how to parameterize reproductive schedule, make repro event
        # right now we are just taking 20% of available storage and moving to
        # growth_biomass = self.sum_plant_parts(plants, parts="growth")
        # filter = np.nonzero(
        #    growth_biomass
        #    >= (
        #        self.species_dispersal_params["min_size_dispersal"]
        #        * self.species_grow_params["max_growth_biomass"]
        #    )
        # )
        # ns_conc = self.get_daily_nsc_concentration(jday)
        # available_stored_biomass = np.zeros_like(plants["root"])
        # for part in self.growth_parts:
        #    avail_nsc_content = np.ones_like(plants[part]) * (
        #        ns_conc[part] - self.species_grow_params["min_nsc_content"][part]
        #    )
        #    avail_nsc_content[avail_nsc_content < 0] = 0.0
        #    available_stored_biomass += plants[part] * avail_nsc_content
        # plants["repro_biomass"][filter] += 0.05 * available_stored_biomass[filter]
        # for part in self.growth_parts:
        #    plants[part][filter] -= (0.05 * (available_stored_biomass[filter])) * (
        #        plants[part][filter] / growth_biomass[filter]
        #    )
        plants = self.form.disperse(plants)
        return plants

    def enter_dormancy(
        self, plants, jday
    ):  # calculate sum of green parts and sum of persistant parts
        end_dead_age = plants["dead_age"]
        end_dead_bio = self.sum_plant_parts(plants, parts="dead")
        plants = self.habit.enter_dormancy(plants)
        new_dead_bio = self.sum_plant_parts(plants, parts="dead")
        plants["dead_age"] = self.calculate_dead_age(
            end_dead_age, end_dead_bio, new_dead_bio
        )
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
            (d >= days["growing_season_start"]) & (d < days["reproduction_start"]),
            (d >= days["reproduction_start"]) & (d < days["senescence_start"]),
            (d >= days["senescence_start"]) & (d < days["growing_season_end"]),
            d >= days["growing_season_end"],
        ]

        for part in self.all_parts:
            nsc_content_opts_b = [
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
                (
                    self.species_grow_params["incremental_nsc"][part][0]
                    + (self.species_grow_params["nsc_content"][part] * 1000) ** 0.5
                )
                ** 2,
            ]

            nsc_content_opts_mx = [
                rate["winter_nsc_rate"][part] * (d + 365 - days["growing_season_end"]),
                rate["spring_nsc_rate"][part] * (days["growing_season_start"] - d),
                rate["summer_nsc_rate"][part] * (days["reproduction_start"] - d),
                rate["fall_nsc_rate"][part] * (days["senescence_start"] - d),
                rate["winter_nsc_rate"][part] * (days["growing_season_end"] - d),
            ]
            nsc_content[part] = (
                (np.select(day_conditions, nsc_content_opts_b)) ** 0.5
                + np.select(day_conditions, nsc_content_opts_mx)
            ) ** 2 / 1000
        return nsc_content

    def litter_decomp(self, _new_biomass):
        decay_rate = self.species_morph_params["biomass_decay_rate"]
        sum_dead_mass = self.sum_plant_parts(_new_biomass, parts="dead")
        cohort_init_mass = sum_dead_mass / np.exp(
            -decay_rate * _new_biomass["dead_age"]
        )
        filter = np.nonzero(sum_dead_mass > 0.0)
        for part in self.dead_parts:
            part_init_mass = np.zeros_like(_new_biomass["dead_age"])
            part_init_mass[filter] = (
                cohort_init_mass[filter]
                * _new_biomass[part][filter]
                / sum_dead_mass[filter]
            )
            _new_biomass[part] = part_init_mass * np.exp(
                -decay_rate * (_new_biomass["dead_age"] + self.dt.astype(float))
            )
        _new_biomass["dead_age"] += self.dt.astype(float) * np.ones_like(
            _new_biomass["dead_age"]
        )
        for part in self.dead_parts:
            filter = np.nonzero(_new_biomass[part] < 0)
            _new_biomass[part][filter] = np.zeros_like(_new_biomass[part][filter])
        return _new_biomass

    def mortality(self, plants, _in_growing_season):
        ###EMILY - the cleanest way to handle this may be to move the code here for whole plant mortality to a
        # separate method (function) and call it then call your leaf mortality method
        mortdict = self.species_mort_params
        #
        # This is just used for whole plant mortality. You can leave it here so you
        # don't have to pass the entire mortdict but you will need to pass the mort_period_bool
        # and factors to the whole plant mortality function
        # set flags for three types of mortality periods
        mort_period_bool = {
            "during growing season": _in_growing_season == True,
            "during dormant season": _in_growing_season == False,
            "year-round": True,
        }
        factors = mortdict["mort_variable_name"]
        # Leave this here since it is used for the dead plant part mass balance
        old_dead_bio = self.sum_plant_parts(plants, parts="dead")
        old_dead_age = plants["dead_age"]

        ####move to whole plant mortality
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
        # leave this part here since we can use the same routine to move the dead leaves
        # into the dead biomass pool and calculate the weighted dead age - which is used for decomp
        new_dead_bio = self.sum_plant_parts(plants, parts="dead")
        plants["dead_age"] = self.calculate_dead_age(
            old_dead_age, old_dead_bio, new_dead_bio
        )
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
        lai = self.calculate_lai(
            _last_biomass["leaf_biomass"], _last_biomass["shoot_sys_width"]
        )
        delta_tot = self.photosynthesis.photosynthesize(
            _par,
            _min_temperature,
            _max_temperature,
            cell_lai,
            _last_biomass,
            _current_day,
        )
        return delta_tot

    def respire(self, _min_temperature, _max_temperature, _last_biomass):
        _temperature = (_min_temperature + _max_temperature) / 2
        growdict = self.species_grow_params
        _new_biomass = _last_biomass.copy()
        # respiration coefficient temp dependence from Teh 2006
        temp_adj = 2 ** ((_temperature - 25) / 10)
        total_delta_respire = np.zeros_like(_last_biomass["root"])
        for part in self.all_parts:
            delta_respire = np.zeros_like(_last_biomass["root"])
            filter = np.nonzero(_last_biomass[part] > 0)
            delta_respire = (
                temp_adj[filter]
                * growdict["respiration_coefficient"][part]
                * _last_biomass[part][filter]
            ) / growdict["glucose_requirement"][part]
            total_delta_respire[filter] += delta_respire[filter]
            _new_biomass[part][filter] -= delta_respire[filter]

        # print("Respiration")
        # print(total_delta_respire)

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
        # dead_abg_biomass = self.sum_plant_parts(plants, parts="dead_aboveground")
        # total_abg_biomass = abg_biomass + dead_abg_biomass
        total_abg_biomass = abg_biomass
        dims = self.shape.calc_abg_dims_from_biomass(total_abg_biomass)
        plants["shoot_sys_width"] = dims[0]
        plants["shoot_sys_height"] = dims[1]
        plants["root_sys_width"] = self.shape.calc_root_sys_width(
            plants["shoot_sys_width"]
        )
        plants["leaf_area"] = plants["leaf"] * self.species_morph_params["sp_leaf_area"]
        return plants

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
        filter = np.where(mass_t2 > 0)
        age_t2[filter] = (
            (age_t1[filter] * mass_t1[filter])
            + ((mass_t2[filter] - mass_t1[filter]) * np.zeros_like(age_t1[filter]))
        ) / (mass_t2[filter])
        return age_t2
