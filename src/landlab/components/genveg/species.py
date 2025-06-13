"""
Species class definition, composition classes, and factory methods
to generate species classes.
These are used by PlantGrowth to differentiate plant properties
and processes for species.
"""

import numpy as np
from sympy import diff
from sympy import lambdify
from sympy import log
from sympy import symbols

from .check_objects import UnitTestChecks
from .form import Bunch
from .form import Colonizing
from .form import Multiplestems
from .form import Rhizomatous
from .form import Singlecrown
from .form import Singlestem
from .form import Stoloniferous
from .form import Thicketforming
from .habit import Forbherb
from .habit import Graminoid
from .habit import Shrub
from .habit import Tree
from .habit import Vine
from .photosynthesis import C3
from .photosynthesis import C4
from .photosynthesis import Cam

rng = np.random.default_rng()


# Define species class that inherits composite class methods
class Species:
    def __init__(self, species_params, latitude, dt=1):
        self.dt = dt
        self.validate_plant_factors(species_params["plant_factors"])
        self.validate_duration_params(species_params["duration_params"])
        self.define_plant_parts(species_params)
        species_params = self.calculate_derived_params(species_params)
        self.form = self.select_form_class(species_params)
        self.habit = self.select_habit_class(species_params)
        self.photosynthesis = self.select_photosythesis_type(species_params, latitude)
        # check these below to see if we need to save the composition dictionary not the original
        self.species_plant_factors = species_params["plant_factors"]
        self.species_duration_params = species_params["duration_params"]
        self.species_grow_params = species_params["grow_params"]
        self.species_photo_params = species_params["photo_params"]
        self.species_dispersal_params = species_params["dispersal_params"]
        self.species_mort_params = species_params["mortality_params"]
        self.species_morph_params = self.habit.morph_params
        self.populate_biomass_allocation_array()

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
            "storage": ["aboveground", "belowground"],
            "p_type": ["C3", "C4", "Cam"],
        }

        for key in plant_factors:
            try:
                opt_list = plant_factor_options[key]
                if opt_list:
                    if plant_factors[key] not in opt_list:
                        msg = "Invalid " + str(key) + " option"
                        raise ValueError(msg)
            except KeyError:
                raise KeyError(
                    "Unexpected variable name in species parameter dictionary. Please check input parameter file"
                )

    def validate_duration_params(self, duration_params):
        if (duration_params["growing_season_start"] < 1) | (
            duration_params["growing_season_start"] > 365
        ):
            msg = "Growing season beginning must be integer values between 1-365"
            raise ValueError(msg)

        if (
            duration_params["growing_season_end"]
            < duration_params["growing_season_start"]
        ) | (duration_params["growing_season_end"] > 365) | (duration_params["growing_season_end"] < 1):
            msg = (
                "Growing season end must be between 1-365 "
                "and greater than the growing season beginning"
            )
            raise ValueError(msg)

        if (
            duration_params["senescence_start"]
            < duration_params["growing_season_start"]
        ) | (
            duration_params["senescence_start"] > duration_params["growing_season_end"]
        ):
            msg = "Start of senescence must be within the growing season"
            raise ValueError(msg)

    def define_plant_parts(self, species_params):
        self.all_parts = list(
            species_params["grow_params"]["glucose_requirement"].keys()
        )
        self.growth_parts = self.all_parts.copy()
        self.growth_parts.remove("reproductive")
        self.dead_parts = [
            "dead_root",
            "dead_leaf",
            "dead_stem",
            "dead_reproductive",
        ]

        if species_params["plant_factors"]["storage"] == "aboveground":
            self.abg_parts = ("leaf", "stem", "reproductive")
            self.dead_abg_parts = ("dead_leaf", "dead_stem", "dead_reproductive")
        else:
            self.abg_parts = ("leaf", "stem")
            self.dead_abg_parts = ("dead_leaf", "dead_stem")

    def calc_area_of_circle(self, diameter):
        return np.pi / 4 * diameter**2

    def calc_volume_cylinder(self, area, height):
        return area * height

    def calc_param_ratio(self, numerator, denominator):
        return numerator / denominator

    def calculate_derived_params(self, species_params):
        morph_params = species_params["morph_params"]
        # Area of circle calcuations
        # check for negative values
        for m_params in [
            "shoot_sys_width",
            "root_sys_width",
        ]:
            for val in [
                "max",
                "mean",
                "min",
            ]:
                UnitTestChecks().is_negative_present(morph_params[m_params][val], val)
        species_params["morph_params"]["canopy_area"] = {}
        species_params["morph_params"]["canopy_area"]["max"] = self.calc_area_of_circle(
            diameter=morph_params["shoot_sys_width"]["max"]
        )
        species_params["morph_params"]["canopy_area"]["mean"] = self.calc_area_of_circle(
            diameter=morph_params["shoot_sys_width"]["mean"]
        )
        species_params["morph_params"]["canopy_area"]["min"] = self.calc_area_of_circle(
            diameter=morph_params["shoot_sys_width"]["min"]
        )
        species_params["morph_params"]["root_area"] = {}
        species_params["morph_params"]["root_area"]["max"] = self.calc_area_of_circle(
            diameter=morph_params["root_sys_width"]["max"]
        )
        species_params["morph_params"]["root_area"]["mean"] = self.calc_area_of_circle(
            diameter=morph_params["root_sys_width"]["mean"]
        )
        species_params["morph_params"]["root_area"]["min"] = self.calc_area_of_circle(
            diameter=morph_params["root_sys_width"]["min"]
        )

        # volume of a cylinder
        # check for negative values
        UnitTestChecks().is_negative_present(morph_params["height"]["max"], "max")
        UnitTestChecks().is_negative_present(
            species_params["morph_params"]["canopy_area"]["max"], "max_canopy_area"
        )

        # ratio calculations
        # check if zero
        for denominator in [
            "n_stems",
            "shoot_sys_width",
            "basal_dia",
        ]:
            for val in [
                "max",
                "mean",
                "min",
            ]:
                UnitTestChecks().is_zero(morph_params[denominator][val], val)
        species_params["morph_params"]["area_per_stem"] = self.calc_param_ratio(
            morph_params["canopy_area"]["max"], morph_params["n_stems"]["max"]
        )

        sum_vars = [
            ["total_biomass", "max", "plant_part_max", self.all_parts],
            ["growth_biomass", "max", "plant_part_max", self.growth_parts],
            ["abg_biomass", "max", "plant_part_max", self.abg_parts],
            ["abg_biomass", "mean", "plant_part_mean", self.abg_parts],
            ["total_biomass", "min", "plant_part_min", self.all_parts],
            ["growth_biomass", "min", "plant_part_min", self.growth_parts],
            ["abg_biomass", "min", "plant_part_min", self.abg_parts],
            [
                "nsc_biomass",
                "min",
                "min_nsc_content",
                self.growth_parts,
            ],  # this is dynamic and varying with plant size
        ]
        for sum_var in sum_vars:
            species_params["grow_params"][sum_var[0]] = {}
        for sum_var in sum_vars:
            species_params["grow_params"][sum_var[0]][sum_var[1]] = 0
            for part in sum_var[3]:
                species_params["grow_params"][sum_var[0]][sum_var[1]] += species_params[
                    "grow_params"
                ][sum_var[2]][part]
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
                UnitTestChecks().is_zero(season_length, f"season length for {season}")
                rate_change_nsc = (
                    species_params["grow_params"]["incremental_nsc"][part][idx]
                    - species_params["grow_params"]["incremental_nsc"][part][idx - 1]
                ) / season_length  # figure out how to loop this
                species_params["duration_params"]["nsc_rate_change"][season][
                    part
                ] = rate_change_nsc
        return species_params

    def calculate_lai(self, leaf_area, shoot_sys_width):
        canopy_area = self.habit._calc_canopy_area_from_shoot_width(shoot_sys_width)
        # value check for leaf_area
        UnitTestChecks().is_negative_present(leaf_area, "leaf_area")

        lai = np.divide(
            leaf_area,
            canopy_area,
            out=np.zeros_like(leaf_area),
            where=~np.isclose(canopy_area, np.zeros_like(shoot_sys_width)),
        )
        return lai

    def select_photosythesis_type(self, species_params, latitude):
        p_type = species_params["plant_factors"]["p_type"]
        photosynthesis_options = {"C3": C3, "C4": C4, "cam": Cam}
        return photosynthesis_options[p_type](
            latitude, photo_params=species_params["photo_params"]
        )

    def select_habit_class(self, species_params):
        habit_type = species_params["plant_factors"]["growth_habit"]
        habit = {
            "forb_herb": Forbherb,
            "graminoid": Graminoid,
            "shrub": Shrub,
            "tree": Tree,
            "vine": Vine,
        }
        return habit[habit_type](species_params, self.dt)

    def select_form_class(self, species_params):
        form_type = species_params["plant_factors"]["growth_form"]
        form = {
            "bunch": Bunch,
            "colonizing": Colonizing,
            "multiple_stems": Multiplestems,
            "rhizomatous": Rhizomatous,
            "single_crown": Singlecrown,
            "single_stem": Singlestem,
            "stoloniferous": Stoloniferous,
            "thicket_forming": Thicketforming,
        }
        return form[form_type](species_params)

    def populate_biomass_allocation_array(self):
        # This method precalculates the biomass allocation array based on plant
        # type (angiosperm/gymnosperm, monocot/dicot) or based on user-defined
        # coefficients
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
                -decay_rate[part] * _new_biomass[str(part) + "_age"][filter]
            )
            _new_biomass[part] = part_init_mass * np.exp(
                -decay_rate[part]
                * (_new_biomass[str(part) + "_age"] + self.dt.astype(float))
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

    def mortality(self, plants, grid, _in_growing_season):
        old_biomass = plants.copy()
        mortdict = self.species_mort_params
        # set flags for three types of mortality periods
        mort_period_bool = {
            "during growing season": _in_growing_season is True,
            "during dormant season": _in_growing_season is False,
            "year-round": True,
        }
        factors = mortdict["mort_variable_name"]
        for key, factor in factors.items():
            run_mort = mort_period_bool[mortdict["period"][key]]
            if not run_mort:
                continue
            else:
                try:
                    pred = grid["cell"][factor][plants["cell_index"]]
                    plants = self.calculate_whole_plant_mortality(plants, pred, key)
                except KeyError:
                    msg = f"No data available for mortality factor {factor}"
                    raise ValueError(msg)
        if _in_growing_season is True:
            plants = self.calculate_shaded_leaf_mortality(plants)
        plants = self.update_dead_biomass(plants, old_biomass)
        return plants

    def calculate_whole_plant_mortality(self, plants, grid_value, key):
        mortdict = self.species_mort_params
        factor_coeffs = mortdict["coeffs"][key]
        # Calculate the probability of survival and cap from 0-1
        prob_survival = 1 / (1 + np.exp(-factor_coeffs[0] * (grid_value - factor_coeffs[1])))
        prob_survival[np.isnan(prob_survival)] = 1.0
        prob_survival[prob_survival < 0] = 0
        prob_survival_daily = prob_survival ** (
            1 / (mortdict["duration"][key] / self.dt.astype(int))
        )
        daily_survival = prob_survival_daily > rng.random(grid_value.shape)
        for part in self.all_parts:
            plants["dead_" + str(part)] = plants["dead_" + str(part)] + (
                plants[part] * (np.invert(daily_survival).astype(int))
            )
            plants[part] = plants[part] * daily_survival.astype(int)
        return plants

    def calculate_shaded_leaf_mortality(self, plants):
        # Teh, Christopher BS. Introduction to mathematical modeling of
        # crop growth: How the equations are derived and assembled into
        # a computer model. Dissertation. com, 2006 based on equation
        # 7.18 and 7.20 (pg. 154)
        leaf_death_rate = self.species_duration_params["death_rate"]["leaf"]
        lai = self.calculate_lai(plants["total_leaf_area"], plants["shoot_sys_width"])
        excess_lai = (
            lai - self.species_morph_params["lai_cr"]
        ) / self.species_morph_params["lai_cr"]
        shaded_leaf = np.nonzero((excess_lai > 0) & (plants["leaf"] > 0))
        D_shade = np.zeros_like(plants["total_leaf_area"])
        D_shade[shaded_leaf] = leaf_death_rate * excess_lai[shaded_leaf]
        D_shade[D_shade > leaf_death_rate] = leaf_death_rate
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
        _relative_water_content,
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
        carb_generated_photo_adj = (
            carb_generated
            * _relative_water_content
            / self.photosynthesis.crit_water_content
        )
        carb_generated_photo_adj[
            carb_generated_photo_adj > carb_generated
        ] = carb_generated[carb_generated_photo_adj > carb_generated]
        return carb_generated_photo_adj

    def respire(self, _min_temperature, _max_temperature, _rel_sat, _last_biomass):
        """
        This function calculates maintenance respiration
        """
        _temperature = (_min_temperature + _max_temperature) / 2
        growdict = self.species_grow_params
        _new_biomass = _last_biomass.copy()
        # respiration coefficient temp dependence from Teh 2006
        temp_adj = 2 ** ((_temperature - 25) / 10)
        for part in self.all_parts:
            delta_respire = np.zeros_like(_last_biomass["root"])
            hypoxic_ratio = growdict["hypoxic_ratio"][part]
            if part in self.abg_parts:
                hypoxic_ratio = 1
            filter = np.nonzero(_last_biomass[part] > 0)
            hypoxia_adjustment = _rel_sat * (hypoxic_ratio - 1) + 1
            delta_respire[filter] = (
                temp_adj[filter]
                * hypoxia_adjustment[filter]
                * growdict["respiration_coefficient"][part]
                * _last_biomass[part][filter]
                / growdict["glucose_requirement"][part]
            )
            _new_biomass[part][filter] -= delta_respire[filter]
        return _new_biomass

    def senesce(self, plants, jday):
        old_biomass = plants.copy()
        death_rate = self.species_duration_params["death_rate"]
        ns_conc = self.get_daily_nsc_concentration(jday)
        ns_green_mass_lost = 0.0
        for part in self.habit.duration.green_parts:
            ns_green_mass_lost += plants[part] * ns_conc[part] * death_rate[part]
        persistent_mass = self.sum_plant_parts(plants, parts="persistent")
        plants = self.habit.senesce(
            plants,
            ns_green_mass=ns_green_mass_lost,
            persistent_mass=persistent_mass,
        )
        plants = self.update_dead_biomass(plants, old_biomass)
        return plants

    def set_initial_biomass(self, plants, in_growing_season):
        est_abg_biomass = self.habit.estimate_abg_biomass_from_cover(plants)
        log_abg_biomass = np.log(
            est_abg_biomass,
            out=np.zeros_like(est_abg_biomass, dtype=np.float64),
            where=(est_abg_biomass > 0.0),
        )
        total_biomass = np.interp(
            (10 ** (log_abg_biomass)),
            self.biomass_allocation_array["abg_biomass"],
            self.biomass_allocation_array["total_biomass"],
        )
        (
            plants["root_biomass"],
            plants["leaf_biomass"],
            plants["stem_biomass"],
        ) = self.habit.duration._solve_biomass_allocation(total_biomass)
        plants["root_sys_width"] = self.habit.calc_root_sys_width(
            plants["shoot_sys_width"], plants["shoot_sys_height"]
        )
        plants = self.habit.duration.set_initial_biomass(
            plants, in_growing_season
        )
        return plants

    def set_new_biomass(self, plants):
        plants = self.habit.duration.set_new_biomass(plants)
        return plants

    def update_morphology(self, plants):
        # Right now this only tracks live biomass - should we track all?
        # Assuming dead leaf area is 95% of live leaf area
        abg_biomass = self.sum_plant_parts(plants, parts="aboveground")
        plants["basal_dia"], plants["shoot_sys_width"], plants["shoot_sys_height"] = (
            self.habit.calc_abg_dims_from_biomass(abg_biomass)
        )
        plants["root_sys_width"] = self.habit.calc_root_sys_width(
            plants["shoot_sys_width"], plants["basal_dia"], plants["shoot_sys_height"]
        )
        dead_leaf_area = np.zeros_like(plants["total_leaf_area"])
        filter = np.nonzero(plants["dead_leaf"] > 0)
        plants["live_leaf_area"] = (
            plants["leaf"] * self.species_morph_params["sp_leaf_area"]
        )
        dead_leaf_area[filter] = (
            0.95 * plants["dead_leaf"][filter]
            * self.species_morph_params["sp_leaf_area"]
        )
        plants["total_leaf_area"] = plants["live_leaf_area"] + dead_leaf_area
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
        age_t2 = age_t1.copy()
        filter = np.nonzero(mass_t2 > 0)
        age_t2[filter] = (
            (age_t1[filter] * mass_t1[filter])
        ) / (mass_t2[filter])
        return age_t2
