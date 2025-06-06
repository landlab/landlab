"""
Growth component of GenVeg - this is the main driver of vegetation growth and
is driven by a photosynthesis model. Vegetation growth depends on the availability
of carbohydrate produced by photosynthetically active plant parts.
"""

import numpy as np

from landlab.components.genveg.species import Species
from landlab.data_record import DataRecord

rng = np.random.default_rng()


class PlantGrowth(Species):
    """
    Add Intro Stuff here
    _name = "PlantGrowth"
    _unit_agnostic = False
    _cite_as =
    @article{piercygv,
        author = {
            Piercy, C.D.;
            Swannack, T.M.;
            Carrillo, C.C.;
            Russ, E.R.;
            Charbonneau, B. M.;
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
                repro_biomass: float, title='reproductive',
                                plant live reproductive biomass in g
                plant_age: int, plant age in days

            species_params: dict, optional,
                a nested dictionary of named vegetation parameters for the
                species or community and process of interest with below sub-dictionaries
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
        self.dt = dt
        super().__init__(species_params, self._lat_rad, self.dt)
        self.species_name = self.species_plant_factors["species"]
        self.time_ind = 1
        event_flags = self.set_event_flags(_current_jday)
        _in_growing_season = event_flags.pop("_in_growing_season")
        max_plants = np.round(
            self._grid.number_of_cells
            * self.species_morph_params["max_plant_density"]
            * self._grid.area_of_cell
        ).astype(int)
        self.no_data_scalar = (
            "N/A",
            999999,
            999999,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            999999,
            np.nan,
            np.nan,
            np.nan,
            999999,
        )
        self.dtypes = [
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
            ("dead_leaf", float),
            ("dead_stem", float),
            ("dead_reproductive", float),
            ("dead_root_age", float),
            ("dead_leaf_age", float),
            ("dead_stem_age", float),
            ("dead_reproductive_age", float),
            ("shoot_sys_width", float),
            ("basal_dia", float),
            ("root_sys_width", float),
            ("shoot_sys_height", float),
            ("root_sys_depth", float),
            ("total_leaf_area", float),
            ("live_leaf_area", float),
            ("plant_age", float),
            ("n_stems", int),
            ("pup_x_loc", float),
            ("pup_y_loc", float),
            ("pup_cost", float),
            ("item_id", int),
        ]
        mask_scalar = 1
        empty_list = []
        mask = []
        for i in range(max_plants[0]):
            empty_list.append(self.no_data_scalar)
            mask.append(mask_scalar)
        self.plants = np.ma.array(empty_list, mask=mask, dtype=self.dtypes)
        self.plants.fill_value = self.no_data_scalar
        try:
            (init_plants, self.n_plants) = kwargs.get(
                ("plant_array", "n_plants"),
                self._init_plants_from_grid(
                    _in_growing_season, kwargs["species_cover"]
                ),
            )
        except KeyError:
            msg = "GenVeg requires a pre-populated plant array or a species cover."
            raise ValueError(msg)

        self.plants[: self.n_plants] = init_plants

        self.call = []
        # Create empty Datarecord to store plant data
        # Instantiate data record
        self.record_plants = DataRecord(
            self._grid,
            time=[rel_time],
            items={
                "grid_element": np.repeat(["cell"], self.n_plants).reshape(
                    self.n_plants, 1
                ),
                "element_id": np.reshape(
                    self.plants["cell_index"][: self.n_plants],
                    (self.n_plants, 1),
                ),
            },
            data_vars={
                "vegetation__species": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["species"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__root_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["root_biomass"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__leaf_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["leaf_biomass"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__stem_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["stem_biomass"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__repro_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["repro_biomass"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__dead_root_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_root"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__dead_leaf_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_leaf"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__dead_stem_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_stem"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__dead_repro_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_reproductive"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__shoot_sys_width": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["shoot_sys_width"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__total_leaf_area": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["total_leaf_area"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
                ),
                "vegetation__plant_age": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["plant_age"][: self.n_plants],
                        (self.n_plants, 1),
                    ),
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
                "vegetation__total_leaf_area": "sq m",
                "vegetation__shoot_sys_width": "m",
                "vegetation__plant_age": "days",
            },
        )
        self.plants["item_id"][: self.n_plants] = self.record_plants.item_coordinates
        # Set constants for PAR formula
        self._wgaus = [0.2778, 0.4444, 0.2778]
        self._xgaus = [0.1127, 0.5, 0.8873]
        self.delta_tot = []

    def species_plants(self):
        unmasked_rows = np.nonzero(self.plants["pid"] >= 0)
        return self.plants[unmasked_rows].filled()

    def species_get_variable(self, var_name):
        return self.species_grow_params

    def update_plants(self, var_names, pids, var_vals):
        for idx, var_name in enumerate(var_names):
            self.plants[var_name][np.isin(self.plants["pid"], pids)] = var_vals[idx]
        return self.plants

    def add_new_plants(self, new_plants_list, _rel_time):

        # Reassess this. We need the INDEX of the last nanmax PID
        last_pid = np.ma.max(self.plants["pid"])
        pids = np.arange(last_pid + 1, last_pid + 1 + new_plants_list.size)
        new_plants_list["pid"] = pids
        new_plants_list["item_id"] = pids
        (n_new_plants,) = new_plants_list.shape
        start_index = np.flatnonzero(self.plants["pid"] == last_pid).astype(int) + 1
        end_index = n_new_plants + start_index[0]
        self.plants[start_index[0] : end_index] = new_plants_list
        self.n_plants += n_new_plants
        self.record_plants.add_item(
            time=np.array([_rel_time]),
            new_item={
                "grid_element": np.repeat(["cell"], n_new_plants).reshape(
                    n_new_plants, 1
                ),
                "element_id": np.reshape(
                    self.plants["cell_index"][start_index[0] : end_index],
                    (n_new_plants, 1),
                ),
            },
            new_item_spec={
                "vegetation__species": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["species"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__root_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["root_biomass"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__leaf_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["leaf_biomass"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__stem_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["stem_biomass"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__repro_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["repro_biomass"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__dead_root_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_root"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__dead_leaf_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_leaf"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__dead_stem_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_stem"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__dead_repro_biomass": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["dead_reproductive"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__shoot_sys_width": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["shoot_sys_width"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__total_leaf_area": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["total_leaf_area"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
                "vegetation__plant_age": (
                    ["item_id", "time"],
                    np.reshape(
                        self.plants["plant_age"][start_index[0] : end_index],
                        (n_new_plants, 1),
                    ),
                ),
            },
        )
        return self.plants

    def _grow(self, _current_jday, _grid_par, _grid_relative_water_content, _grid_relative_saturation):
        # This is the primary method in PlantGrowth and is run within each
        # GenVeg run_one_step at each timestep. This method applies new environmental
        # conditions daily from the grid, determines what processes run, and implements
        # them in order to update the plant array.

        # set up shorthand aliases and reset
        # look at checking to see if we can use recordmask here
        _last_biomass = self.plants[~self.plants["pid"].mask].copy()
        _new_biomass = self.plants[~self.plants["pid"].mask]
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
        _new_biomass = self.mortality(_new_biomass)
        _new_biomass = self.litter_decomp(_new_biomass)
        # Limit growth processes only to live plants
        _total_biomass = self.sum_plant_parts(_new_biomass, parts="total")
        filter = np.nonzero(_total_biomass > 0.0)
        _new_live_biomass = _new_biomass[filter]

        # calculate variables needed to run plant processes
        _par = _grid_par[_last_biomass["cell_index"]][filter]
        _relative_water_content = _grid_relative_water_content[
            _last_biomass["cell_index"]
        ][filter]
        _rel_sat = _grid_relative_saturation[_last_biomass["cell_index"]][filter]
        _min_temperature = self._grid["cell"]["air__min_temperature_C"][
            _last_biomass["cell_index"]
        ][filter]
        _max_temperature = self._grid["cell"]["air__max_temperature_C"][
            _last_biomass["cell_index"]
        ][filter]
        _cell_lai = self._grid["cell"]["vegetation__leaf_area_index"][
            _last_biomass["cell_index"]
        ][filter]
        _new_live_biomass = self.respire(
            _min_temperature, _max_temperature, _rel_sat, _new_live_biomass
        )

        # Change this so for positive delta_tot we allocate by size and
        if event_flags["_in_growing_season"]:
            carb_generated_photo = processes["_in_growing_season"](
                _par,
                _min_temperature,
                _max_temperature,
                _cell_lai,
                _relative_water_content,
                _new_live_biomass,
                _current_jday,
            )

            # Future add turnover rate
            _new_live_biomass = self.allocate_biomass_dynamically(
                _new_live_biomass, carb_generated_photo
            )
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
        _new_live_biomass = self.update_morphology(_new_live_biomass)
        _new_biomass[filter] = _new_live_biomass
        _new_biomass = self.update_dead_biomass(_new_biomass, _last_biomass)
        self.plants[~self.plants["pid"].mask] = _new_biomass
        self.plants, self.n_plants = self.remove_plants()

    def _init_plants_from_grid(self, in_growing_season, species_cover):
        """
        This method initializes the plants in the PlantGrowth class
        # from the vegetation fields stored on the grid. This method
        # is only called if no initial plant array is parameterized
        # as part of the PlantGrowth initialization.
        # Required parameters are a boolean inidicating if the plants are
        # in the active growing season.
        """
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
                    plant_basal_dias = []
                    while cover_area > (
                        1.2 * self.species_morph_params["min_canopy_area"]
                    ):
                        plant_basal_dia = rng.uniform(
                            low=self.species_morph_params["min_basal_dia"],
                            high=self.species_morph_params["max_basal_dia"],
                            size=1,
                        )
                        plant_canopy_area = self.habit._calc_canopy_area(
                            plant_basal_dia
                        )
                        plant_shoot_width = (
                            self.habit._calc_shoot_width_from_canopy_area(
                                plant_canopy_area
                            )
                        )
                        cover_area -= (
                            np.pi
                            / 4
                            * np.sqrt(plant_shoot_width * plant_basal_dia) ** 2
                        )
                        if cover_area > 0:
                            plant_basal_dias.append(plant_basal_dia)
                            plant_shoot_widths.append(plant_shoot_width)
                        else:
                            breakpoint
                    for index, new_plant_width in enumerate(plant_shoot_widths):
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
                                0.0,
                                0.0,
                                0.0,
                                new_plant_width,
                                plant_basal_dias[index],
                                0.0,
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
        plant_array = np.array(plantlist, dtype=self.dtypes)
        plant_array = self.set_initial_biomass(plant_array, in_growing_season)
        return (plant_array, pidval)

    def allocate_biomass_proportionately(
        self, _last_biomass, _total_biomass, delta_tot
    ):
        """
        This method allocates new net biomass amongst growth parts
        proportionately based on the relative size of the part. This
        method is used outside of the growing season since some plant parts
        may not be present while the plant is dormant. The storage redistribution
        method is called after initial biomass allocation to redistribute storage
        biomass to dormant growth parts as needed.
        Required parameter is a numpy array of net biomass change to be applied
        to the plant and it returns the structured array _new_biomass.
        """
        _new_biomass = _last_biomass
        _total_biomass = self.sum_plant_parts(_last_biomass, parts="total")
        filter = np.nonzero(_total_biomass != 0)
        for part in self.all_parts:
            _new_biomass[part][filter] = (
                _last_biomass[part][filter] / _total_biomass[filter]
            ) * delta_tot[filter] + _last_biomass[part][filter]
        return _new_biomass

    def adjust_biomass_allocation_towards_ideal(self, _new_biomass):
        """
        This method adjusts biomass allocation towards the ideal allocation
        proportions based on the plant size. If parts of the plant are
        removed via herbivory or damage, this allows the plant to utilize
        other stored resources to regrow the damaged parts.
        """
        _total_biomass = self.sum_plant_parts(_new_biomass, parts="growth")
        _min_leaf_mass_frac = (
            self.species_grow_params["plant_part_min"]["leaf"] / _total_biomass
        )
        _min_stem_mass_frac = (
            self.species_grow_params["plant_part_min"]["stem"] / _total_biomass
        )

        _min_root_mass_frac = (
            self.species_grow_params["plant_part_min"]["root"] / _total_biomass
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

        _new_root_mass_frac = 1 - _new_leaf_mass_frac - _new_stem_mass_frac
        _new_leaf_mass_frac[_new_leaf_mass_frac < _min_leaf_mass_frac] = (
            _min_leaf_mass_frac[_new_leaf_mass_frac < _min_leaf_mass_frac]
        )
        _new_stem_mass_frac[_new_stem_mass_frac < _min_stem_mass_frac] = (
            _min_stem_mass_frac[_new_stem_mass_frac < _min_stem_mass_frac]
        )

        root_diff = _new_root_mass_frac - _min_root_mass_frac
        filter = np.nonzero(root_diff < 0)
        _new_root_mass_frac[filter] = _min_root_mass_frac[filter]
        _leaf_allocation = _new_leaf_mass_frac / (
            _new_leaf_mass_frac + _new_stem_mass_frac
        )
        _stem_allocation = _new_stem_mass_frac / (
            _new_leaf_mass_frac + _new_stem_mass_frac
        )
        _new_leaf_mass_frac[filter] = _new_leaf_mass_frac[filter] + (
            root_diff[filter] * _leaf_allocation[filter]
        )
        _new_stem_mass_frac[filter] = _new_stem_mass_frac[filter] + (
            root_diff[filter] * _stem_allocation[filter]
        )
        _new_biomass["root_biomass"] = (_new_root_mass_frac) * _total_biomass
        _new_biomass["leaf_biomass"] = _new_leaf_mass_frac * _total_biomass
        _new_biomass["stem_biomass"] = _new_stem_mass_frac * _total_biomass
        return _new_biomass

    def set_event_flags(self, _current_jday):
        """
        This method sets event flags so required processes are run based
        on the day of year.
        """
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
        # This method moved live biomass to dead biomass is the plant
        # is too small to grow.
        min_size = self.species_grow_params["min_growth_biomass"]
        total_biomass = self.sum_plant_parts(_new_biomass, parts="growth")
        dead_plants = np.nonzero(total_biomass < min_size)
        if dead_plants[0].size > 0:
            print(str(dead_plants[0].size) + " were too small to survive")

        for part in self.all_parts:
            _new_biomass[part][dead_plants][
                np.isnan(_new_biomass[part][dead_plants])
                | (_new_biomass[part][dead_plants] < 0)
            ] = 0.0
            _new_biomass["dead_" + str(part)][dead_plants] += _new_biomass[part][
                dead_plants
            ]
            _new_biomass[part][dead_plants] = 0.0
            _new_biomass[part][_new_biomass[part] < 0] = 0.0
        return _new_biomass

    def remove_plants(self):
        # Plants that have too little dead biomass remaining to track
        # are removed from the plant array and no longer tracked.
        min_size_dead = 0.1
        min_size_live = self.species_grow_params["min_growth_biomass"]
        total_live_biomass = self.sum_plant_parts(self.plants, parts="growth")
        total_dead_biomass = self.sum_plant_parts(self.plants, parts="dead")
        remove_plants = np.flatnonzero(
            (total_dead_biomass < min_size_dead) & (total_live_biomass < min_size_live)
        )
        self.plants[remove_plants] = self.no_data_scalar
        self.plants[remove_plants] = np.ma.masked
        remove_array_length = np.count_nonzero(remove_plants)
        self.n_plants -= remove_array_length
        return self.plants, self.n_plants

    def save_plant_output(self, rel_time, save_params):
        # This method saves plant properties at the required time step
        # future work versions will save additional variables based on user input.
        self.record_plants.add_record(time=np.array([rel_time]))
        self.record_plants.ffill_grid_element_and_id()

        item_ids = self.plants["item_id"][~self.plants["item_id"].mask]

        self.record_plants.dataset["vegetation__species"].values[
            item_ids, self.time_ind
        ] = self.plants["species"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__root_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["root_biomass"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__leaf_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["leaf_biomass"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__stem_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["stem_biomass"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__repro_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["repro_biomass"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__dead_root_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_root"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__dead_leaf_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_leaf"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__dead_stem_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_stem"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__dead_repro_biomass"].values[
            item_ids, self.time_ind
        ] = self.plants["dead_reproductive"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__total_leaf_area"].values[
            item_ids, self.time_ind
        ] = self.plants["total_leaf_area"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__shoot_sys_width"].values[
            item_ids, self.time_ind
        ] = self.plants["shoot_sys_width"][~self.plants["item_id"].mask]
        self.record_plants.dataset["vegetation__plant_age"].values[
            item_ids, self.time_ind
        ] = self.plants["plant_age"][~self.plants["item_id"].mask]
        self.time_ind += 1
