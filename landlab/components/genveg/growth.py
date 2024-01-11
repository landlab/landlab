"""
Growth component of GenVeg - this is the main driver of vegetation growth and
is driven by a photosynthesis model. Vegetation growth depends on the availability
of carbohydrate produced by photosynthetically active plant parts. 
"""
from landlab.components.genveg.species import Species
from landlab.data_record import DataRecord
import numpy as np

rng = np.random.default_rng()


class PlantGrowth(Species):
    """
    Add Intro Stuff here
    _name = "PlantGrowth"
    _unit_agnostic = False
    _cite_as =
    @article{piercygv,
        author = {Piercy, C.D.; Swannack, T.M.; Carrillo, C.C.; Russ, E.R.; Charbonneau, B. M.]
    }
    #Add all variables to be saved or chosen as outputs here
    _info = {
        "vegetation__total_biomass": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units":"g",
            "mapping":"cell",
            "doc": "Total plant biomass for the plant class at the end of the time step"
        },
    }
    """

    def __init__(
        self,
        grid,
        dt,
        rel_time,
        _current_jday,
        species_params={
            "col_params": {},
            "disp_params": {},
            "duration_params": {
                "growing_season_start": 91,
                "growing_season_end": 290,
                "senescence_start": 228,
            },
            "grow_params": {
                "respiration_coefficient": [0.015, 0.015, 0.03],
                "glucose_requirement": [1.444, 1.513, 1.463],
                "k_light_extinct": 0.02,
                "light_half_sat": 9,
                "p_max": 0.055,
                "root_to_leaf_coeffs": [0.031, 0.951, 0],
                "root_to_stem_coeffs": [-0.107, 1.098, 0.0216],
                "plant_part_min": [0.01, 0.1, 0.5],
            },
            "mort_params": {
                "s1_days": 365,
                "s1_name": "Mortality factor",
                "s1_pred": [1, 2, 3, 4],
                "s1_rate": [0, 0.1, 0.9, 1],
                "s1_weight": [1000, 1, 1, 1000],
            },
            "plant_factors": {
                "species": "Corn",
                "growth_form": 1,
                "monocot_dicot": "monocot",
                "angio_gymno": "angiosperm",
                "annual_perennial": "annual",
                "p_type": "C3",
            },
            "size_params": {
                "max_height_stem": 2.5,
                "max_mass_stem": 72,
                "max_n_stems": 3,
                "max_plant_density": 1,
            },
            "stor_params": {"r_wint_die": 0.25, "r_wint_stor": 0.25},
        },
        **kwargs,
    ):
        """Instantiate PlantGrowth
        Parameters
        ----------
        grid: RasterModelGrid
            A Landlab ModelGrid

        dt: NumPy time delta, required,
            time step interval

        rel_time: int, required,
            number of time steps elapsed

        _current_jday: int, required
            day of the year assuming Jan 1 is 1

        **kwargs to send to init
            plants: Numpy structured array of individual plants, optional
                with columns
                species: string, plant species names
                pid: int, plant ID
                cell_index: int, index of cell location on grid
                root_biomass: float, title='root', plant live root biomass in g
                leaf_biomass: float, title='stem', plant live leaf biomass in g
                stem_biomass: float, title='stem', plant live stem biomass in g
                storage_biomass: float, title='storage', plant live stem biomass in g
                repro_biomass: float, title='reproductive', plant live reproductive biomass in g
                plant_age: int, plant age in days

            species_params: dict, optional,
                a nested dictionary of named vegetation parameters for the
                species or community and process of interest with below sub-dictionaries.
                plant_factors: dict, required,
                    dictionary of plant characteristics describing the
                    species or community of interest with below keys
                    species: string, required,
                        name of species or community used to identify plant
                    growth_form: string, required,
                        USDA plant growth habit,
                        graminoid, forb/herb, shrub, tree, vine
                    monocot_dicot: string, required,
                        should be monocot or dicot
                    angio_gymno: string, required,
                        should be angiosperm or gymnosperm
                    annual_perennial: string, required,
                        plant growth duration, annual (1 year) or
                        perennial (multiple years)
                    p_type: string, required,
                        photosythesis type, either 'C3', 'C4', or 'CAM'
                    leaf_retention: string, required,
                        evergreen or deciduous (annuals are deciduous)
                duration_params: dict, required,
                    dictionary of parameters defining the growing season,
                    growing_season_start: int, required,
                        growing season start day of year,
                        must be between 1-365
                    growing_season_end: int, required,
                        growing season end day of year,
                        must be between 1-365
                    senesecence_start: int, required,
                        start of senescence period after plant reaches peak biomass,
                        must be between gs_start and gs_end
                grow_params: dict, required,
                    dictionary of paramaters required to simulate plant growth

                    respiration_coefficient: dict, required,
                        respiration coefficient with keys
                            'root': float
                            'leaf': float
                            'stem': float
                            'reproductive': float
                    glucose_requirements: dict, required,
                        glucose requirement
                    le_k: float, required,
                        light extinction coefficient
                    hi: float, required,
                        something
                    p_max: float, required,
                        maximum photosyntehtic output
        """
        # Initialize species object to get correct species parameter list
        self._grid = grid
        (_, _latitude) = self._grid.xy_of_reference
        self._lat_rad = np.radians(_latitude)
        super().__init__(species_params, self._lat_rad)
        self.species_name = self.species_plant_factors["species"]

        self.dt = dt
        self.time_ind = 1

        event_flags = self.set_event_flags(_current_jday)
        _in_growing_season = event_flags.pop("_in_growing_season")
        self.plants = kwargs.get(
            "plant_array",
            self._init_plants_from_grid(_in_growing_season, kwargs["species_cover"]),
        )
        self.call = []

        # Create empty Datarecord to store plant data
        # Instantiate data record
        self.record_plants = DataRecord(
            self._grid,
            time=[rel_time],
            items={
                "grid_element": np.repeat(["cell"], self.plants["pid"].size).reshape(
                    self.plants["pid"].size, 1
                ),
                "element_id": np.reshape(
                    self.plants["cell_index"], (self.plants["pid"].size, 1)
                ),
            },
            data_vars={
                "vegetation__species": (
                    ["item_id", "time"],
                    np.reshape(self.plants["species"], (self.plants["pid"].size, 1)),
                ),
                "vegetation__root_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["root_biomass"], (self.plants["pid"].size, 1)
                    ),
                ),
                "vegetation__leaf_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["leaf_biomass"], (self.plants["pid"].size, 1)
                    ),
                ),
                "vegetation__stem_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["stem_biomass"], (self.plants["pid"].size, 1)
                    ),
                ),
                "vegetation__repro_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["repro_biomass"], (self.plants["pid"].size, 1)
                    ),
                ),
                "vegetation__dead_root_biomass": (
                    ["item_id", "time"],
                    np.reshape(self.plants["dead_root"], (self.plants["pid"].size, 1)),
                ),
                "vegetation__dead_leaf_biomass": (
                    ["item_id", "time"],
                    np.reshape(self.plants["dead_leaf"], (self.plants["pid"].size, 1)),
                ),
                "vegetation__dead_stem_biomass": (
                    ["item_id", "time"],
                    np.reshape(self.plants["dead_stem"], (self.plants["pid"].size, 1)),
                ),
                "vegetation__dead_repro_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_reproductive"], (self.plants["pid"].size, 1)
                    ),
                ),
                "vegetation__shoot_sys_width": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["shoot_sys_width"], (self.plants["pid"].size, 1)
                    ),
                ),
                "vegetation__leaf_area": (
                    ["item_id", "time"],
                    np.reshape(self.plants["leaf_area"], (self.plants["pid"].size, 1)),
                ),
                "vegetation__plant_age": (
                    ["item_id", "time"],
                    np.reshape(self.plants["plant_age"], (self.plants["pid"].size, 1)),
                ),
            },
            attrs={
                "vegetation__species": "species name, string",
                "vegetation__root_biomass": "g",
                "vegetation__leaf_biomass": "g",
                "vegetation__stem_biomass": "g",
                "vegetation__repro_biomass": "g",
                "vegetation__dead_root_biomass": "g",
                "vegetation__dead_leaf_biomass": "g",
                "vegetation__dead_stem_biomass": "g",
                "vegetation__dead_repro_biomass": "g",
                "vegetation__leaf_area": "sq m",
                "vegetation__shoot_sys_width": "m",
                "vegetation__plant_age": "days",
            },
        )
        self.plants["item_id"] = self.record_plants.item_coordinates
        # Set constants for PAR formula
        self._wgaus = [0.2778, 0.4444, 0.2778]
        self._xgaus = [0.1127, 0.5, 0.8873]
        self.delta_tot = []

    def species_plants(self):
        return self.plants

    def species_get_variable(self, var_name):
        return self.species_grow_params

    def update_plants(self, var_names, pids, var_vals):
        updated_plants = self.plants
        for idx, var_name in enumerate(var_names):
            updated_plants[var_name][np.isin(self.plants["pid"], pids)] = var_vals[idx]
        self.plants = updated_plants
        return updated_plants

    def add_new_plants(self, new_plants_list):
        old_plants = self.plants
        last_pid = self.plants["pid"][-1]
        pids = np.arange(last_pid + 1, last_pid + new_plants_list.size + 1)
        new_plants_list["pid"] = pids
        new_plants_list["item_id"] = pids
        np.concatenate((old_plants, new_plants_list), axis=0)
        self.plants = old_plants
        return self.plants

    def _grow(self, _current_jday):
        ###This is the primary method in PlantGrowth and is run within each
        # GenVeg run_one_step at each timestep. This method applies new environmental
        # conditions daily from the grid, determines what processes run, and implements
        # them in order to update the plant array.

        # set up shorthand aliases and reset
        growdict = self.species_grow_params
        _last_biomass = self.plants
        _last_dead_biomass = self.sum_plant_parts(_last_biomass, parts="dead")
        _new_biomass = _last_biomass.copy()

        # Decide what processes happen today
        event_flags = self.set_event_flags(_current_jday)
        processes = {
            "_in_growing_season": self.photosynthesize,
            "_in_senescence_period": self.senesce,
            "_in_reproductive_period": self.disperse,
            "_is_emergence_day": self.emerge,
            "_is_dormant_day": self.enter_dormancy,
        }

        # Run mortality and decompose litter each day
        _new_biomass = self.mortality(_new_biomass, event_flags["_in_growing_season"])
        _new_biomass = self.litter_decomp(_new_biomass)

        # Limit growth processes only to live plants
        _total_biomass = self.sum_plant_parts(_new_biomass, parts="total")
        filter = np.nonzero(_total_biomass > 0.0)
        _last_live_biomass = _last_biomass[filter]
        _live_biomass = _new_biomass[filter]
        _last_total_biomass = _total_biomass[filter]

        # calculate variables needed to run plant processes
        _par = self._grid["cell"]["radiation__par_tot"][_last_biomass["cell_index"]][
            filter
        ]
        _min_temperature = self._grid["cell"]["air__min_temperature_C"][
            _last_biomass["cell_index"]
        ][filter]
        _max_temperature = self._grid["cell"]["air__max_temperature_C"][
            _last_biomass["cell_index"]
        ][filter]
        _cell_lai = self._grid["cell"]["vegetation__cell_lai"][
            _last_biomass["cell_index"]
        ][filter]
        _new_live_biomass = self.respire(
            _min_temperature, _max_temperature, _last_live_biomass
        )

        # Change this so for positive delta_tot we allocate by size and
        if event_flags["_in_growing_season"]:
            delta_photo = processes["_in_growing_season"](
                _par,
                _min_temperature,
                _max_temperature,
                _cell_lai,
                _new_live_biomass,
                _current_jday,
            )
            # we will need to add a turnover rate estiamte. Not sure where to include it though.
            # delta_tot=delta_tot-self.species_grow_params['biomass_turnover_rate']*self.dt
            _new_live_biomass = self.allocate_biomass_dynamically(
                _live_biomass, delta_photo
            )
        # else:
        #    _new_live_biomass = self.allocate_biomass_proportionately(
        #        _live_biomass, _last_total_biomass, delta_biomass_respire
        #    )
        _new_live_biomass = self.kill_small_plants(_new_live_biomass)
        event_flags.pop("_in_growing_season")
        # Run all other processes that need to occur
        for process in event_flags.items():
            if process[1]:
                _new_live_biomass = processes[process[0]](
                    _new_live_biomass, _current_jday
                )
        _new_live_biomass["plant_age"] += self.dt.astype(float) * np.ones_like(
            _new_live_biomass["plant_age"]
        )

        _new_biomass[filter] = _new_live_biomass
        _new_biomass = self.update_morphology(_new_biomass)
        _new_biomass = self.update_dead_biomass(
            _new_biomass, _last_biomass, _last_dead_biomass
        )
        _new_biomass = self.remove_plants(_new_biomass)

        # print('Completed run_one_step for '+self.species_name+' on day '+str(_current_jday))
        self.plants = _new_biomass

    def _init_plants_from_grid(self, in_growing_season, species_cover):
        ###This method initializes the plants in the PlantGrowth class
        # from the vegetation fields stored on the grid. This method
        # is only called if no initial plant array is parameterized
        # as part of the PlantGrowth initialization.
        # Required parameters are a boolean inidicating if the plants are
        # in the active growing season.
        ###
        dtypes = [
            ("species", "U10"),
            ("pid", int),
            ("cell_index", int),
            ("x_loc", float),
            ("y_loc", float),
            (("root", "root_biomass"), float),
            (("leaf", "leaf_biomass"), float),
            (("stem", "stem_biomass"), float),
            (("reproductive", "repro_biomass"), float),
            ("dead_root", float),
            ("dead_stem", float),
            ("dead_leaf", float),
            ("dead_reproductive", float),
            ("dead_age", float),
            ("shoot_sys_width", float),
            ("root_sys_width", float),
            ("shoot_sys_height", float),
            ("root_sys_depth", float),
            ("leaf_area", float),
            ("plant_age", float),
            ("n_stems", int),
            ("pup_x_loc", float),
            ("pup_y_loc", float),
            ("pup_cost", float),
            ("item_id", int),
        ]
        pidval = 0
        plantlist = []
        # Loop through grid cells
        for cell_index in range(self._grid.number_of_cells):
            cell_plants = self._grid["cell"]["vegetation__plant_species"][cell_index]
            cell_cover = species_cover[cell_index]
            # Loop through list of plants stored on grid cell
            for plant in cell_plants:
                if plant == self.species_plant_factors["species"]:
                    plant_cover = cell_cover[plant]
                    cover_area = (
                        plant_cover * self._grid.area_of_cell[cell_index] * 0.907
                    )
                    plant_shoot_widths = []
                    while cover_area > (
                        1.2 * self.species_morph_params["min_crown_area"]
                    ):
                        plant_width = rng.uniform(
                            low=self.species_morph_params["min_shoot_sys_width"],
                            high=self.species_morph_params["max_shoot_sys_width"],
                            size=1,
                        )
                        cover_area -= np.pi / 4 * plant_width**2
                        if cover_area > 0:
                            plant_shoot_widths.append(plant_width)
                        else:
                            breakpoint
                    for new_plant_width in plant_shoot_widths:
                        plantlist.append(
                            (
                                plant,
                                pidval,
                                cell_index,
                                np.nan,
                                np.nan,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                new_plant_width,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0,
                                np.nan,
                                np.nan,
                                np.nan,
                                0,
                            )
                        )
                        pidval += 1
        plant_array = np.array(plantlist, dtype=dtypes)
        plant_array = self.set_initial_biomass(plant_array, in_growing_season)
        return plant_array

    def allocate_biomass_dynamically(self, _live_biomass, delta_tot):
        ###This method allocates new biomass according to the size-dependent
        # biomass allocation array calculated upon initiation of the PlantGrowth class.
        # The array is only valid for actively growing plants so this method is only
        # used during the growing season.
        # After initial allocation, storage redistribution, reallocation, and the
        # minimum size check methods are called to adjust the biomass in each part before
        # saving.
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
        # filter = np.nonzero(
        #    growth_biomass
        #    >= (
        #        self.species_dispersal_params["min_size_dispersal"]
        #        * self.species_grow_params["max_growth_biomass"]
        #    )
        # )
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
        # _new_biomass = self.redistribute_storage_biomass(_new_biomass)
        _new_biomass = self.adjust_biomass_allocation_towards_ideal(_new_biomass)
        return _new_biomass

    def allocate_biomass_proportionately(
        self, _last_biomass, _total_biomass, delta_tot
    ):
        ###This method allocates new net biomass amongst growth parts
        # proportionately based on the relative size of the part. This
        # method is used outside of the growing season since some plant parts
        # may not be present while the plant is dormant. The storage redistribution
        # method is called after initial biomass allocation to redistribute storage
        # biomass to dormant growth parts as needed.
        # Required parameter is a numpy array of net biomass change to be applied
        # to the plant and it returns the structured array _new_biomass.
        _new_biomass = _last_biomass
        _total_biomass = self.sum_plant_parts(_last_biomass, parts="total")
        filter = np.nonzero(_total_biomass != 0)
        for part in self.all_parts:
            _new_biomass[part][filter] = (
                _last_biomass[part][filter] / _total_biomass[filter]
            ) * delta_tot[filter] + _last_biomass[part][filter]
        return _new_biomass

    def adjust_biomass_allocation_towards_ideal(self, _new_biomass):
        ###This method adjusts biomass allocation towards the ideal allocation
        # proportions based on the plant size. If parts of the plant are
        # removed via herbivory or damage, this allows the plant to utilize
        # other stored resources to regrow the damaged parts.
        _total_biomass = self.sum_plant_parts(_new_biomass, parts="growth")
        _min_leaf_mass_frac = (
            self.species_grow_params["plant_part_min"]["leaf"] / _total_biomass
        )
        _min_stem_mass_frac = (
            self.species_grow_params["plant_part_min"]["stem"] / _total_biomass
        )

        current_leaf_mass_frac = np.divide(
            _new_biomass["leaf_biomass"],
            _total_biomass,
            out=np.zeros_like(_total_biomass),
            where=~np.isclose(_total_biomass, np.zeros_like(_total_biomass)),
        )
        current_stem_mass_frac = np.divide(
            _new_biomass["stem_biomass"],
            _total_biomass,
            out=np.zeros_like(_total_biomass),
            where=~np.isclose(_total_biomass, np.zeros_like(_total_biomass)),
        )

        ideal_leaf_mass_frac = np.interp(
            _total_biomass,
            self.biomass_allocation_array["total_biomass"],
            self.biomass_allocation_array["leaf_mass_frac"],
        )
        ideal_stem_mass_frac = np.interp(
            _total_biomass,
            self.biomass_allocation_array["total_biomass"],
            self.biomass_allocation_array["stem_mass_frac"],
        )

        current_diff_leaf = ideal_leaf_mass_frac - current_leaf_mass_frac
        current_diff_stem = ideal_stem_mass_frac - current_stem_mass_frac

        _new_leaf_mass_frac = ideal_leaf_mass_frac - current_diff_leaf * (
            1
            - np.exp(
                -self.species_duration_params["senesce_rate"] * self.dt.astype(int)
            )
        )
        _new_stem_mass_frac = ideal_stem_mass_frac - current_diff_stem * (
            1
            - np.exp(
                -self.species_duration_params["senesce_rate"] * self.dt.astype(int)
            )
        )

        _new_leaf_mass_frac[
            _new_leaf_mass_frac < _min_leaf_mass_frac
        ] = _min_leaf_mass_frac[_new_leaf_mass_frac < _min_leaf_mass_frac]
        _new_stem_mass_frac[
            _new_stem_mass_frac < _min_stem_mass_frac
        ] = _min_stem_mass_frac[_new_stem_mass_frac < _min_stem_mass_frac]

        _new_biomass["root_biomass"] = (
            1 - _new_leaf_mass_frac - _new_stem_mass_frac
        ) * _total_biomass
        _new_biomass["leaf_biomass"] = _new_leaf_mass_frac * _total_biomass
        _new_biomass["stem_biomass"] = _new_stem_mass_frac * _total_biomass
        return _new_biomass

    def set_event_flags(self, _current_jday):
        ###This method sets event flags so required processes are run based
        # on the day of year.
        durationdict = self.species_duration_params
        flags_to_test = {
            "_in_growing_season": bool(
                (_current_jday > durationdict["growing_season_start"])
                & (_current_jday < durationdict["growing_season_end"])
            ),
            "_is_emergence_day": bool(
                _current_jday == durationdict["growing_season_start"]
            ),
            "_in_reproductive_period": bool(
                (_current_jday >= durationdict["reproduction_start"])
                & (_current_jday < durationdict["reproduction_end"])
            ),
            "_in_senescence_period": bool(
                (_current_jday >= durationdict["senescence_start"])
                & (_current_jday < durationdict["growing_season_end"])
            ),
            "_is_dormant_day": bool(
                _current_jday == durationdict["growing_season_end"]
            ),
        }
        return flags_to_test

    def kill_small_plants(self, _new_biomass):
        ###This method moved live biomass to dead biomass is the plant
        # is too small to grow.

        # Edit this and figure out when to kill plants
        for part in self.all_parts:
            _new_biomass[part][
                _new_biomass[part] < self.species_grow_params["plant_part_min"][part]
            ]
        min_size = self.species_grow_params["min_growth_biomass"]
        total_biomass = self.sum_plant_parts(_new_biomass, parts="growth")
        dead_plants = np.nonzero(total_biomass < min_size)
        if dead_plants[0].size > 0:
            print(str(dead_plants[0].size) + " were too small to survive")

        for part in self.all_parts:
            _new_biomass["dead_" + str(part)][dead_plants] += _new_biomass[part][
                dead_plants
            ]
            _new_biomass[part][dead_plants] = np.zeros_like(
                _new_biomass[part][dead_plants]
            )
            _new_biomass[part][_new_biomass[part] < 0] = np.zeros_like(
                _new_biomass[part][_new_biomass[part] < 0]
            )
        return _new_biomass

    def update_dead_biomass(self, _new_biomass, old_biomass, old_dead_biomass):
        for part in self.all_parts:
            part_biomass_change = _new_biomass[part] - old_biomass[part]
            filter = np.nonzero(part_biomass_change < 0.0)
            _new_biomass["dead_" + str(part)][filter] += -part_biomass_change[filter]

        _new_dead_biomass = self.sum_plant_parts(_new_biomass, parts="dead")
        _new_biomass["dead_age"] = self.calculate_dead_age(
            _new_biomass["dead_age"], old_dead_biomass, _new_dead_biomass
        )
        return _new_biomass

    def remove_plants(self, _new_biomass):
        ###Plants that have too little dead biomass remaining to track
        # are removed from the plant array and no longer tracked.
        min_size_dead = 0.1
        min_size_live = self.species_grow_params["min_growth_biomass"]
        total_live_biomass = self.sum_plant_parts(_new_biomass, parts="growth")
        total_dead_biomass = self.sum_plant_parts(_new_biomass, parts="dead")
        remove_plants = np.nonzero(
            (total_live_biomass < min_size_live) & (total_dead_biomass < min_size_dead)
        )
        # print("There were " + str(_new_biomass.size) + " plants")
        _new_biomass = np.delete(_new_biomass, remove_plants, axis=None)
        # print("After removal, there were " + str(_new_biomass.size) + (" plants"))
        return _new_biomass

    def save_plant_output(self, rel_time, save_params):
        ###This method saves plant properties at the required time step
        # future work versions will save additional variables based on user input.
        prev_time = self.record_plants.latest_time
        # new_ids=np.where(np.isin(self.plants['pid'],self.record_plants.dataset['item_id'].values[:,prev_time], invert=True))
        # for i in new_ids:
        #    self.record_plants.add_item(time=[rel_time],
        #        new_item={
        #            'grid_element'=['cell'],
        #            'element_id'=self.plants['cell_index'][i],
        #        },

        #    )
        self.record_plants.add_record(time=np.array([rel_time]))
        self.record_plants.ffill_grid_element_and_id()

        item_ids = self.plants["item_id"]

        self.record_plants.dataset["vegetation__species"].values[
            item_ids, self.time_ind
        ] = self.plants["species"]
        self.record_plants.dataset["vegetation__root_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["root_biomass"]
        self.record_plants.dataset["vegetation__leaf_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["leaf_biomass"]
        self.record_plants.dataset["vegetation__stem_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["stem_biomass"]
        self.record_plants.dataset["vegetation__repro_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["repro_biomass"]
        self.record_plants.dataset["vegetation__dead_root_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_root"]
        self.record_plants.dataset["vegetation__dead_leaf_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_leaf"]
        self.record_plants.dataset["vegetation__dead_stem_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_stem"]
        self.record_plants.dataset["vegetation__dead_repro_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_reproductive"]
        self.record_plants.dataset["vegetation__leaf_area"].values[
            item_ids, self.time_ind
        ] = self.plants["leaf_area"]
        self.record_plants.dataset["vegetation__shoot_sys_width"].values[
            item_ids, self.time_ind
        ] = self.plants["shoot_sys_width"]
        self.record_plants.dataset["vegetation__plant_age"].values[
            item_ids, self.time_ind
        ] = self.plants["plant_age"]
        #    item_id= np.reshape(self.plants['item_id'],(self.plants['item_id'].size,1)),
        #    new_record={
        #        'vegetation__species':(['item_id','time'],self.plants['species']),
        #        #'vegetation__root_biomass':(['item_id','time'],np.reshape(self.plants['root_biomass'],(self.plants['cell_index'].size,1))),
        #        #'vegetation__leaf_biomass':(['item_id','time'],np.reshape(self.plants['leaf_biomass'],(self.plants['cell_index'].size,1))),
        #        #'vegetation__stem_biomass':(['item_id','time'],np.reshape(self.plants['stem_biomass'],(self.plants['cell_index'].size,1)))
        #    }
        # )
        self.time_ind += 1
