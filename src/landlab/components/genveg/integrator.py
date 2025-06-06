"""
Plant integration component of GenVeg - this is the part of GenVeg
that handles interactions between plants and plants and the physical grid
"""

import matplotlib as plt
import numpy as np
import pandas as pd

# from landlab.components import Radiation
from landlab import Component

from .growth import PlantGrowth

rng = np.random.default_rng()


class GenVeg(Component, PlantGrowth):
    """
    Add Intro Stuff here
    """

    _name = "GenVeg"

    _unit_agnostic = False

    _cite_as = """
    @article{piercygv,
        author = {
            Piercy, C.D.;
            Swannack, T.M.;
            Russ, E.R.;
            Catlett, A.R.
    }
    """
    # Add all variables to be saved or chosen as outputs here
    _info = {
        "plant__total_biomass": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "g",
            "mapping": "cell",
            "doc": "Total plant biomass at the end of the time step",
        },
    }

    def __init__(
        self, grid, dt, current_day, vegparams, plant_array=np.empty((0, 30), dtype=[])
    ):
        # save grid object to class
        super().__init__(grid)
        # Check to see is the grid has the right data and assign defaults
        self.current_day = current_day
        self.check_for_grid_fields()
        self._max_water_availability = self._field_capacity - self._wilt_pt
        (_, _latitude) = self._grid.xy_of_reference
        self._lat_rad = np.radians(_latitude)
        # Set initial time variables
        self.dt = dt
        self.start_date = current_day
        self.time_ind = 0
        # self.neighbors=self._grid.looped_neighbors_at_cell()
        self.nodes = self._grid.node_at_cell
        _current_jday = self._calc_current_jday()
        rel_time = self._calc_rel_time()
        # Create empty array to store PlantGrowth objects
        plantspecies = []
        _ = self._grid.add_zeros("vegetation__total_biomass", at="cell", clobber=True)
        _ = self._grid.add_zeros("vegetation__n_plants", at="cell", clobber=True)
        _ = self._grid.add_zeros("vegetation__plant_height", at="cell", clobber=True)
        _ = self._grid.add_zeros("vegetation__lai", at="cell", clobber=True)

        # Instantiate a PlantGrowth object and
        # summarize number of plants and biomass per cell
        if plant_array.size == 0:
            cover_allocation = []
            for cell_index in range(self._grid.number_of_cells):
                species_list = self._grid.at_cell["vegetation__plant_species"][
                    cell_index
                ]
                species_list = species_list[~np.isin(species_list, "null")]
                cell_cover = self._grid.at_cell["vegetation__cover_fraction"][
                    cell_index
                ]
                number_of_species = len(species_list)
                cover_species = rng.uniform(low=0.5, high=1.0, size=number_of_species)
                cover_sum = sum(cover_species)
                species_cover_allocation = cover_species / cover_sum
                cover_allocation.append(
                    dict(zip(species_list, (species_cover_allocation * cell_cover)))
                )

            # for each species in the parameters file
            for species in vegparams:
                if species == "null":
                    continue
                species_dict = vegparams[species]
                species_obj = PlantGrowth(
                    self._grid,
                    self.dt,
                    _current_jday,
                    rel_time,
                    species_dict,
                    species_cover=cover_allocation,
                )
                plantspecies.append(species_obj)
        else:
            for species in vegparams:
                plant_array = plant_array[plant_array["species"] == species]
                species_dict = vegparams[species]
                species_canopy_area = np.pi / 4 * plant_array["shoot_sys_width"] ** 2
                species_basal_area = np.pi / 4 * plant_array["basal_dia"] ** 2
                species_abg_area = np.sqrt(species_basal_area * species_canopy_area)
                species_percent_cover = (
                    self.calculate_grid_vars(
                        plant_array["cell_index"], species_abg_area
                    )
                    / self._grid.area_of_cell
                )
                species_cover = [
                    dict(
                        zip(
                            [self._grid.at_cell["vegetation__plant_species"][i]],
                            [species_percent_cover[i]],
                        )
                    )
                    for i in range(self._grid.number_of_cells)
                ]

                species_obj = PlantGrowth(
                    self._grid,
                    self.dt,
                    _current_jday,
                    rel_time,
                    species_dict,
                    species_cover=species_cover,
                    plant_array=plant_array,
                )
                plantspecies.append(species_obj)
        self.plant_species = plantspecies

        # Summarize biomass and number of plants across grid
        all_plants = self.combine_plant_arrays()
        tot_bio_species = (
            all_plants["root_biomass"]
            + all_plants["leaf_biomass"]
            + all_plants["stem_biomass"]
        )
        abg_area = np.pi / 4 * all_plants["shoot_sys_width"] ** 2
        cell_biomass = self.calculate_grid_vars(
            all_plants["cell_index"], tot_bio_species
        )
        cell_plant_count = self.calculate_grid_vars(all_plants["cell_index"])

        frac_cover = (
            self.calculate_grid_vars(
                all_plants["cell_index"],
                abg_area,
            )
            / self._grid.area_of_cell
        )
        cell_leaf_area = self.calculate_grid_vars(
            all_plants["cell_index"],
            all_plants["total_leaf_area"],
        )
        plant_height = np.zeros_like(frac_cover)
        n_of_plants = cell_plant_count.astype(np.float64)
        cells_with_plants = np.where(n_of_plants > 0.0)
        sum_plant_height = self.calculate_grid_vars(
            all_plants["cell_index"],
            all_plants["shoot_sys_height"],
        )
        plant_height[cells_with_plants] = (
            sum_plant_height[cells_with_plants] / n_of_plants[cells_with_plants]
        )

        self._grid.at_cell["vegetation__live_biomass"] = cell_biomass
        self._grid.at_cell["vegetation__plant_count"] = cell_plant_count
        self._grid.at_cell["vegetation__cover_fraction"] = frac_cover
        self._grid.at_cell["vegetation__plant_height"] = plant_height
        self._grid.at_cell["vegetation__leaf_area_index"] = (
            cell_leaf_area / self._grid.area_of_cell
        )

        # add location information for each plant
        for cell_index in range(self._grid.number_of_cells):
            cell_corners = self._grid.corners_at_cell[cell_index]
            x_vertices = self._grid.x_of_corner[cell_corners]
            y_vertices = self._grid.y_of_corner[cell_corners]
            min_x = np.min(x_vertices)
            max_x = np.max(x_vertices)
            corner_vertices = np.array(list(zip(x_vertices, y_vertices)))
            cell_plants = all_plants[all_plants["cell_index"] == cell_index]
            cell_poly = plt.pyplot.Polygon(corner_vertices)

            # Check if point falls in cell
            for idx, plant in enumerate(cell_plants):
                unoccupied_center = False
                radius = plant["basal_dia"] / 2
                while unoccupied_center is False:
                    x = rng.uniform(low=min_x + radius, high=max_x - radius, size=1)
                    y_lims = self.get_cell_boundary_points(corner_vertices, x)
                    y = rng.uniform(
                        low=y_lims[0] + radius, high=y_lims[1] - radius, size=1
                    )
                    point = (*x, *y)
                    if cell_poly.contains_point(point):
                        [unoccupied_center] = self.check_if_loc_unocc(
                            [point], [radius], cell_plants, "above"
                        )
                        if (unoccupied_center is True) or (idx == 0):
                            cell_plants[idx]["x_loc"] = x
                            cell_plants[idx]["y_loc"] = y
                            break

            for species_obj in self.plant_species:
                species = species_obj.species_plant_factors["species"]
                update_plants = cell_plants[cell_plants["species"] == species]
                update_plants = species_obj.update_morphology(update_plants)
                species_obj.update_plants(
                    ["x_loc", "y_loc", "shoot_sys_width", "shoot_sys_height"],
                    update_plants["pid"],
                    np.vstack(
                        (
                            update_plants["x_loc"],
                            update_plants["y_loc"],
                            update_plants["shoot_sys_width"],
                            update_plants["shoot_sys_height"],
                        )
                    ),
                )

    def calculate_grid_vars(self, indices, grid_var=None):
        """
        Uses the bincount method to calculate plant counts when no
        value for grid_var is provided and sum variables by grid cell
        when a grid variable is provided.
        """
        obs = np.nonzero(indices >= 0.0)
        if grid_var is None:
            weight_var = grid_var
        else:
            weight_var = grid_var[obs]
        var_out = np.bincount(
            indices[obs],
            weights=weight_var,
            minlength=self._grid.number_of_cells,
        )
        return var_out.astype(np.float64)

    def check_for_grid_fields(
        self,
        soil_texture_defaults={
            "porosity": {
                "sand": 0.43,
                "sandy loam": 0.39,
                "sandy clay loam": 0.41,
                "loam": 0.42,
                "silt loam": 0.43,
                "silt": 0.40,
                "silty clay loam": 0.47,
                "clay loam": 0.44,
                "clay": 0.49,
            },
            "field_capacity": {
                "sand": 0.09,
                "sandy loam": 0.16,
                "sandy clay loam": 0.26,
                "loam": 0.25,
                "silt loam": 0.29,
                "silt": 0.29,
                "silty clay loam": 0.37,
                "clay loam": 0.34,
                "clay": 0.42,
            },
            "wilting_point": {
                "sand": 0.02,
                "sandy loam": 0.07,
                "sandy clay loam": 0.16,
                "loam": 0.12,
                "silt loam": 0.11,
                "silt": 0.06,
                "silty clay loam": 0.20,
                "clay loam": 0.20,
                "clay": 0.28,
            },
        },
    ):
        try:
            self.plants_on_grid = self._grid["cell"]["vegetation__plant_species"]
        except KeyError:
            msg = "GenVeg requires initial distribution of plant species at-cell field."
            raise ValueError(msg)
        # Check to see if grid contains required environmental fields
        try:
            self.min_air_temp = self._grid["cell"]["air__min_temperature_C"][:]
            self.max_air_temp = self._grid["cell"]["air__max_temperature_C"][:]
        except KeyError:
            msg = (
                "GenVeg requires min and max air temperatures "
                "in Celcius for each time step."
            )
            raise KeyError(msg)

        try:
            self._par = self._grid["cell"]["radiation__total_par"][:]
        except KeyError:
            msg = (
                "GenVeg requires incoming PAR for each timestep. "
                "Empiricial estimation will be used for the run."
            )
            print(msg)
            # add radiation here
        else:
            self._par_method = "direct_input"
        try:
            self._soil_water = self._grid["cell"]["soil_water__volume_fraction"][
                :
            ]
        except KeyError:
            msg = (
                "Soil moisture field is not found. "
                "GenVeg will use random soil moisture values"
            )
            print(msg)
            self._soil_water = rng.uniform(
                low=0.3,
                high=1.0,
                size=self._grid["cell"]["radiation__total_par"][:].size,
            )

        try:
            self._wilt_pt = self._grid["cell"]["soil__wilting_point"][:].copy()
            self._field_capacity = self._grid["cell"]["soil__field_capacity"][:].copy()
            self._porosity = self._grid["cell"]["soil__porosity"]
        except KeyError:
            msg = "Default soil physical properties will be used based on soil texture"
            print(msg)
            try:
                self._wilt_pt = np.vectorize(
                    soil_texture_defaults["wilting_point"].get
                )(self._grid["cell"]["surface__soil_texture"])
                self._field_capacity = np.vectorize(
                    soil_texture_defaults["field_capacity"].get
                )(self._grid["cell"]["surface__soil_texture"])
            except KeyError:
                msg = "No soil texture provided so assuming values for silt loam"
                print(msg)
                self._wilt_pt = np.ones_like(self._par) * (
                    (soil_texture_defaults["wilting_point"].get)("silt loam")
                )
                self._field_capacity = np.ones_like(self._par) * (
                    soil_texture_defaults["field_capacity"].get
                )("silt loam")

    def get_int_output(self):
        print(self.species_cover_allocation)

    def run_one_step(self):
        _current_jday = self._calc_current_jday()
        cell_biomass = np.zeros_like(self._grid.at_cell["vegetation__live_biomass"])
        cell_plant_count = np.zeros_like(self._grid.at_cell["vegetation__plant_count"])
        _available_water_cell = (
            np.minimum(self._soil_water, self._field_capacity) - self._wilt_pt
        )
        _max_water_available = self._field_capacity - self._wilt_pt
        _available_water_frac = _available_water_cell / _max_water_available
        _frac_above_fc = np.subtract(
            self._soil_water,
            self._field_capacity,
            out=np.zeros_like(self._soil_water),
            where=(self._soil_water > self._field_capacity)
        )
        _rel_saturation = _frac_above_fc / (1 - self._field_capacity)
        all_plants = []
        for species_obj in self.plant_species:
            species_obj._grow(_current_jday, self._par, _available_water_frac, _rel_saturation)

        all_plants = self.combine_plant_arrays()
        all_plants = self.check_for_dispersal_success(all_plants)

        tot_bio_species = (
            all_plants["root_biomass"]
            + all_plants["leaf_biomass"]
            + all_plants["stem_biomass"]
        )
        abg_area = (
            np.pi
            / 4
            * (np.sqrt(all_plants["shoot_sys_width"] * all_plants["basal_dia"])) ** 2
        )
        cell_biomass = self.calculate_grid_vars(
            all_plants["cell_index"], tot_bio_species
        )
        cell_plant_count = self.calculate_grid_vars(all_plants["cell_index"])
        cell_percent_cover = (
            self.calculate_grid_vars(
                all_plants["cell_index"],
                abg_area,
            )
            / self._grid.area_of_cell
        )
        cell_leaf_area = self.calculate_grid_vars(
            all_plants["cell_index"],
            all_plants["total_leaf_area"],
        )
        cell_leaf_area[cell_leaf_area < 0] = 0.0
        cell_leaf_area[np.isnan(cell_leaf_area)] = 0.0
        plant_height = np.zeros_like(cell_percent_cover)
        n_of_plants = cell_plant_count.astype(np.float64)
        cells_with_plants = np.nonzero(n_of_plants > 0.0)
        sum_plant_height = self.calculate_grid_vars(
            all_plants["cell_index"],
            all_plants["shoot_sys_height"],
        )
        plant_height[cells_with_plants] = (
            sum_plant_height[cells_with_plants] / n_of_plants[cells_with_plants]
        )
        self._grid.at_cell["vegetation__live_biomass"] = cell_biomass
        self._grid.at_cell["vegetation__plant_count"] = cell_plant_count
        self._grid.at_cell["vegetation__cover_fraction"] = cell_percent_cover
        self._grid.at_cell["vegetation__plant_height"] = plant_height
        self._grid.at_cell["vegetation__leaf_area_index"] = np.divide(
            cell_leaf_area,
            self._grid.area_of_cell,
            np.zeros_like(self._grid.at_cell["vegetation__live_biomass"]),
            where=~np.isclose(
                cell_leaf_area,
                np.zeros_like(self._grid.at_cell["vegetation__live_biomass"]),
            ),
        )
        self.current_day += 1

    def _calc_current_jday(self):
        jday_td = self.current_day - np.datetime64(
            str(self.current_day.astype("datetime64[Y]")) + "-01-01"
        )
        _current_jday = jday_td.astype(int)
        return _current_jday

    def _calc_available_water(self):
        water_ratio = (self._soil_water - self._wilt_pt) / self._max_water_availability
        water_ratio[water_ratio > 1] = 1
        return water_ratio

    def _calc_rel_time(self):
        return (self.current_day - self.start_date).astype(float)

    def get_cell_boundary_points(self, vertices, x_value):
        # Create a Path object from the polygon
        # Initialize a list to store the intersection points
        intersection_points = []

        # Loop through all the vertices of the polygon
        for i in range(len(vertices)):
            x1, y1 = vertices[i]
            x2, y2 = vertices[(i + 1) % len(vertices)]
            if x1 <= x_value and x2 >= x_value or x1 >= x_value and x2 <= x_value:
                # Calculate the intersection point between the polygon edge
                # and the vertical line at x_value
                y_intersection = (x_value - x1) * (y2 - y1) / (x2 - x1) + y1
                intersection_points.append(y_intersection)

        if len(intersection_points) == 0:
            return None
        else:
            return min(intersection_points), max(intersection_points)

    def check_if_loc_unocc(self, plant_loc, plant_width, all_plants, check_type):
        plants_with_locations = all_plants[~np.isnan(all_plants["x_loc"])]

        # This code looks for the cells around the plant cell_index
        area = {"above": "shoot_sys_width", "below": "root_sys_width"}
        is_center_unocc = []
        for idx, loc in enumerate(plant_loc):
            distance = (
                np.sqrt(
                    (loc[0] - plants_with_locations["x_loc"]) ** 2
                    + (loc[1] - plants_with_locations["y_loc"]) ** 2
                )
                - plant_width[idx] / 2
            )
            no_conflict = distance > plants_with_locations[area[check_type]] / 2
            is_center_unocc.append(np.all(no_conflict))
        return is_center_unocc

    def check_for_dispersal_success(self, all_plants):
        new_pups = all_plants[~np.isnan(all_plants["pup_x_loc"])]

        if new_pups.size != 0:
            pup_locs = tuple(zip(new_pups["pup_x_loc"], new_pups["pup_y_loc"]))
            pup_widths = np.zeros_like(new_pups["pup_x_loc"])
            loc_unoccupied = self.check_if_loc_unocc(
                pup_locs, pup_widths, all_plants, "below"
            )
            new_pups = new_pups[np.nonzero(loc_unoccupied)]
            new_pups["x_loc"] = new_pups["pup_x_loc"]
            new_pups["y_loc"] = new_pups["pup_y_loc"]

        for species_obj in self.plant_species:
            species = species_obj.species_name
            species_new_pups = new_pups[
                new_pups["species"] == species
            ]  # This is a slice
            species_plants = all_plants[all_plants["species"] == species]
            if species_new_pups.size != 0:
                species_parents = species_new_pups.copy()
                species_new_pups = species_obj.habit.duration.set_new_biomass(
                    species_new_pups
                )
                species_new_pups = species_obj.update_morphology(species_new_pups)
                species_new_pups["plant_age"] = np.zeros_like(species_new_pups["root"])
                species_new_pups["cell_index"] = self._grid.cell_at_node[
                    self._grid.find_nearest_node(
                        (species_new_pups["x_loc"], species_new_pups["y_loc"]),
                        mode="clip",
                    )
                ]
                species_parents["reproductive"] = species_parents["reproductive"] - (
                    species_new_pups["root"]
                    + species_new_pups["leaf"]
                    + species_new_pups["stem"]
                    + species_parents["pup_cost"]
                )

                species_obj.update_plants(
                    ["reproductive"],
                    species_parents["pid"],
                    species_parents["reproductive"],
                )
                _rel_time = self._calc_rel_time
                species_obj.add_new_plants(species_new_pups, _rel_time)

            species_obj.update_plants(
                ["pup_x_loc", "pup_y_loc", "pup_cost"],
                species_plants["pid"],
                np.vstack(
                    (
                        np.full_like(species_plants["root"], np.nan),
                        np.full_like(species_plants["root"], np.nan),
                        np.full_like(species_plants["root"], np.nan),
                    )
                ),
            )

        return self.combine_plant_arrays()

    def combine_plant_arrays(self):
        all_plants = []
        for species_obj in self.plant_species:
            array_out = species_obj.species_plants()
            plant_entries = array_out[: species_obj.n_plants]
            all_plants.append(np.ravel(plant_entries))

        all_plants_array = all_plants[0]
        for i in range(1, len(all_plants)):
            all_plants_array = np.concatenate((all_plants_array, all_plants[i]))
        return all_plants_array

    def view_record_grid(
        self,
    ):
        view = self.record_grid.dataset.to_dataframe()
        return view

    def print_test_output(self):
        pass
        # return self.test_output

    def save_output(self, save_params=["root_biomass", "leaf_biomass", "stem_biomass"]):
        rel_time = self._calc_rel_time()
        for species_obj in self.plant_species:
            species_obj.species_plants()
            species_obj.save_plant_output(rel_time, save_params)
        self.time_ind += 1

    def get_plant_output(self, species="all"):
        if species == "all":
            out_df = pd.DataFrame()
            for species_obj in self.plant_species:
                species_df = species_obj.record_plants.dataset.to_dataframe()
                species_df.reset_index(inplace=True)
                species_df.set_index(
                    ["time", "vegetation__species", "item_id"], inplace=True
                )
                out_df = pd.concat([out_df, species_df])
        else:
            for species_obj in self.plant_species:
                if species_obj.species_name == species:
                    out_df = species_obj.record_plants.dataset.to_dataframe()
        return out_df
