"""
Plant integration component of GenVeg - this is the part of GenVeg 
that handles interactions between plants and plants and the physical grid
"""

# from landlab.components import Radiation
from landlab import Component
import numpy as np
import pandas as pd
from .growth import PlantGrowth
import matplotlib as plt

rng = np.random.default_rng()


class GenVeg(Component, PlantGrowth):
    """
    Add Intro Stuff here
    """

    _name = "GenVeg"

    _unit_agnostic = False

    _cite_as = """
    @article{piercygv,
        author = {Piercy, C.D.; Swannack, T.M.; Carrillo, C.C.; Russ, E.R.; Charbonneau, B. M.]
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
            "doc": "Total plant biomass for the plant class at the end of the time step",
        },
    }

    def __init__(self, grid, dt, current_day, vegparams):
        # save grid object to class
        super().__init__(grid)
        # Check to see if there are plants on the grid

        try:
            self.plants_on_grid = self._grid["cell"]["vegetation__plant_species"]
        except KeyError:
            msg = "GenVeg requires initial distribution of plant species at-cell field."
            raise ValueError(msg)
        # Check to see if grid contains required environmental fields
        try:
            self.min_air_temp = self._grid["cell"]["air__min_temperature_C"][:].copy()
            self.max_air_temp = self._grid["cell"]["air__max_temperature_C"][:].copy()
        except KeyError:
            msg = "GenVeg requires min and max air temperatures in Celcius for each time step."
            raise ValueError(msg)

        try:
            self._par = self._grid["cell"]["radiation__par_tot"][:].copy()
        except RuntimeWarning:
            msg = "GenVeg requires incoming PAR for each timestep. Empiricial estimation will be used for the run."
            print(msg)
            # try:
            #    self.__albedo_bare=self._grid['cell']['bare_ground_albedo'][:].copy()
            # except KeyError:
            #   msg=('Empirical estimation of PAR requires bare ground albedo at-cell field.')
            #   raise ValueError(msg)
            #            # From chapter 2 of Teh 2006 (pg. 31; equation 2.13)
            # try:

            # except KeyError:
            #    msg=('Empirical estimation of PAR requires lat/long of grid xy reference.')
            # try:
            #    self._clear_sky_index=self._grid['cell']['air__clear_sky_index'][:].copy()
            # except KeyError:
            #    msg=('Empirical estimation of PAR requires a clear sky index value at-cell field')
            # except RuntimeError:

            #    raise RuntimeError(msg)
            # else:
            #    self._par_method='empirical_estimation'
        else:
            self._par_method = "direct_input"

            (_, _latitude) = self._grid.xy_of_reference
            self._lat_rad = np.radians(_latitude)

        # Set initial time variables
        self.dt = dt
        self.current_day = current_day
        self.start_date = current_day
        self.time_ind = 0
        # self.neighbors=self._grid.looped_neighbors_at_cell()
        self.nodes = self._grid.node_at_cell
        _current_jday = self._calc_current_jday()
        rel_time = self._calc_rel_time()

        ##Need to modify to allow user to input plant array for hotstarting
        ##Instantiate a plantgrowth object and summarize number of plants and biomass per cell
        # Create empty array to store PlantGrowth objects
        plantspecies = []
        _ = self._grid.add_zeros("vegetation__total_biomass", at="cell", clobber=True)
        _ = self._grid.add_zeros("vegetation__n_plants", at="cell", clobber=True)
        _ = self._grid.add_zeros("vegetation__plant_height", at="cell", clobber=True)
        _ = self._grid.add_zeros("vegetation__cell_lai", at="cell", clobber=True)
        self.species_cover_allocation = []

        for cell_index in range(self._grid.number_of_cells):
            species_list = self._grid.at_cell["vegetation__plant_species"][cell_index]
            cell_cover = self._grid.at_cell["vegetation__percent_cover"][cell_index]
            number_of_species = len(species_list)
            cover_species = rng.uniform(low=0.5, high=1.0, size=number_of_species)
            cover_species[species_list == "null"] = 0.0
            cover_sum = sum(cover_species)
            species_cover_allocation = cover_species / cover_sum

            self.species_cover_allocation.append(
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
                species_cover=self.species_cover_allocation,
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
        cell_biomass = np.bincount(
            all_plants["cell_index"],
            weights=tot_bio_species,
            minlength=self._grid.number_of_cells,
        )
        cell_plant_count = np.bincount(
            all_plants["cell_index"], minlength=self._grid.number_of_cells
        )
        frac_cover = (
            np.bincount(
                all_plants["cell_index"],
                weights=abg_area,
                minlength=self._grid.number_of_cells,
            )
            / self._grid.area_of_cell
        )
        cell_leaf_area = np.bincount(
            all_plants["cell_index"],
            weights=all_plants["total_leaf_area"],
            minlength=self._grid.number_of_cells,
        )
        plant_height = np.zeros_like(frac_cover)
        n_of_plants = cell_plant_count.astype(np.float64)
        cells_with_plants = np.where(n_of_plants > 0.0)
        sum_plant_height = np.bincount(
            all_plants["cell_index"],
            weights=all_plants["shoot_sys_height"],
            minlength=self._grid.number_of_cells,
        )
        plant_height[cells_with_plants] = (
            sum_plant_height[cells_with_plants] / n_of_plants[cells_with_plants]
        )

        self._grid.at_cell["vegetation__total_biomass"] = cell_biomass
        self._grid.at_cell["vegetation__n_plants"] = cell_plant_count
        self._grid.at_cell["vegetation__percent_cover"] = frac_cover
        self._grid.at_cell["vegetation__plant_height"] = plant_height
        self._grid.at_cell["vegetation__cell_lai"] = (
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
                radius = plant["shoot_sys_width"] / 2
                while unoccupied_center == False:
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
                        if (unoccupied_center == True) or (idx == 0):
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

    def get_int_output(self):
        print(self.species_cover_allocation)

    def run_one_step(self):
        _current_jday = self._calc_current_jday()
        cell_biomass = np.zeros_like(self._grid.at_cell["vegetation__total_biomass"])
        cell_plant_count = np.zeros_like(self._grid.at_cell["vegetation__n_plants"])

        all_plants = []
        for species_obj in self.plant_species:
            species_obj._grow(_current_jday)

        all_plants = self.combine_plant_arrays()
        all_plants = self.check_for_dispersal_success(all_plants)

        tot_bio_species = (
            all_plants["root_biomass"]
            + all_plants["leaf_biomass"]
            + all_plants["stem_biomass"]
        )
        abg_area = np.pi / 4 * all_plants["shoot_sys_width"] ** 2
        cell_biomass = np.bincount(
            all_plants["cell_index"],
            weights=tot_bio_species,
            minlength=self._grid.number_of_cells,
        )
        cell_plant_count = np.bincount(
            all_plants["cell_index"], minlength=self._grid.number_of_cells
        )
        cell_percent_cover = (
            np.bincount(
                all_plants["cell_index"],
                weights=abg_area,
                minlength=self._grid.number_of_cells,
            )
            / self._grid.area_of_cell
        )
        cell_leaf_area = np.bincount(
            all_plants["cell_index"],
            weights=all_plants["total_leaf_area"],
            minlength=self._grid.number_of_cells,
        )
        cell_leaf_area[cell_leaf_area < 0] = 0
        cell_leaf_area[np.isnan(cell_leaf_area)] = 0
        plant_height = np.zeros_like(cell_percent_cover)
        n_of_plants = cell_plant_count.astype(np.float64)
        cells_with_plants = np.nonzero(n_of_plants > 0.0)
        sum_plant_height = np.bincount(
            all_plants["cell_index"],
            weights=all_plants["shoot_sys_height"],
            minlength=self._grid.number_of_cells,
        )
        plant_height[cells_with_plants] = (
            sum_plant_height[cells_with_plants] / n_of_plants[cells_with_plants]
        )
        self._grid.at_cell["vegetation__total_biomass"] = cell_biomass
        self._grid.at_cell["vegetation__n_plants"] = cell_plant_count
        self._grid.at_cell["vegetation__percent_cover"] = cell_percent_cover
        self._grid.at_cell["vegetation__plant_height"] = plant_height
        self._grid.at_cell["vegetation__cell_lai"] = np.divide(
            cell_leaf_area,
            self._grid.area_of_cell,
            np.zeros_like(self._grid.at_cell["vegetation__total_biomass"]),
            where=~np.isclose(
                cell_leaf_area,
                np.zeros_like(self._grid.at_cell["vegetation__total_biomass"]),
            ),
        )
        self.current_day += 1

    def _calc_current_jday(self):
        jday_td = self.current_day - np.datetime64(
            str(self.current_day.astype("datetime64[Y]")) + "-01-01"
        )
        _current_jday = jday_td.astype(int)
        return _current_jday

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
                # Calculate the intersection point between the polygon edge and the vertical line at x_value
                y_intersection = (x_value - x1) * (y2 - y1) / (x2 - x1) + y1
                intersection_points.append(y_intersection)

        if len(intersection_points) == 0:
            return None
        else:
            return min(intersection_points), max(intersection_points)

    def check_if_loc_unocc(self, plant_loc, plant_width, all_plants, check_type):
        plants_with_locations = all_plants[~np.isnan(all_plants["x_loc"])]

        # This code looks for the cells around the plant cell_index
        #    parent_cell=new_pups['cell_index']
        #    query_cells=self._grid.looped_neighbors_at_cell[pup['cell_index']]
        #    query_cells.append(parent_cell)
        #    plants_to_check=all_plants[np.isin(new_pups['cell_index'], query_cells)]
        # check to see if we can vectorize this method
        # We can also downselect so we are only checking plants in the surrounding cells
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
            # print("There were " + str(new_pups.size) + " potential new plants")
            new_pups = new_pups[np.nonzero(loc_unoccupied)]
            # print("Only " + str(new_pups.size) + " were successful")
            new_pups["x_loc"] = new_pups["pup_x_loc"]
            new_pups["y_loc"] = new_pups["pup_y_loc"]

        for species_obj in self.plant_species:
            species = species_obj.species_name
            species_new_pups = new_pups[new_pups["species"] == species]
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
                species_obj.add_new_plants(species_new_pups)
                print("Successful dispersal occurred")

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
            all_plants.append(np.ravel(array_out))

        all_plants_array = all_plants[0]
        for i in range(1, len(all_plants)):
            all_plants_array = np.concatenate((all_plants_array, all_plants[i]))
        # all_plants_array=np.ravel(all_plants_array)
        return all_plants_array

        # need list of all plants within neighborhood plus their location and radius
        # for pup in new_pups:
        #    parent_cell=new_pups['cell_index']
        #    query_cells=self._grid.looped_neighbors_at_cell[pup['cell_index']]
        #    query_cells.append(parent_cell)
        #    plants_to_check=all_plants[np.isin(new_pups['cell_index'], query_cells)]
        #    plant_centers=zip(plants_to_check['pup_x_loc'], plants_to_check['pup_y_loc'])
        #    plant_radii=plants_to_check['root_sys_width']/2
        #    parent_loc=zip(pup['x_loc'], pup['y_loc'])
        #    pup_loc=zip(pup['pup_x_loc'], pup['pup_y_loc'])
        #    shortest_distance=self._calc_shortest_path(parent_loc, pup_loc, plant_centers, plant_radii)

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

    def _calc_shortest_path(
        self, parent_loc, pup_loc, obstacle_centers, obstacle_radii
    ):
        min_distance = np.inf
        min_path = None

        for idx, obs_center in enumerate(obstacle_centers):
            obs_radius = obstacle_radii(idx)

            runner_direction = parent_loc - pup_loc
            obs_direction = parent_loc - obs_center

            angle = np.arccos(
                np.dot(runner_direction, obs_direction)
                / (np.linalg.norm(obs_direction))
            )
            distance_to_obstacle = (
                np.linalg.norm(obs_direction) * np.sin(angle) - obs_radius
            )

            # Calculate the vector from the obstacle center to the closest point on the path
            path_direction = np.array([-obs_direction[1], obs_direction[0]])
            if np.dot(path_direction, runner_direction) < 0:
                path_direction = -path_direction

            # Calculate the closest point on the path to the start point
            path_point = (
                parent_loc
                + path_direction * distance_to_obstacle / np.linalg.norm(path_direction)
            )

            # Calculate the distance from the start point to the closest point on the path
            distance_to_path = np.linalg.norm(path_point - parent_loc)

            # Calculate the distance from the closest point on the path to the end point
            distance_from_path_to_end = np.linalg.norm(pup_loc - path_point)

            # Calculate the total distance
            total_distance = distance_to_path + distance_from_path_to_end

            # Update the minimum distance and path if necessary
            if total_distance < min_distance:
                min_distance = total_distance
                min_path = path_point

        # If there are no obstacles, the shortest path is a straight line
        if min_path is None:
            return np.linalg.norm(pup_loc - parent_loc)

        # Calculate the shortest path from the start point to the closest point on the path
        shortest_path_start = self.shortest_path(
            parent_loc, min_path, obstacle_centers, obstacle_radii
        )

        # Calculate the shortest path from the closest point on the path to the end point
        shortest_path_end = self.shortest_path(
            min_path, pup_loc, obstacle_centers, obstacle_radii
        )

        # Return the total shortest path
        return shortest_path_start + shortest_path_end
