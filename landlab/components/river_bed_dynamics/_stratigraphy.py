"""
Implements a series of functions to create, initialize, and run the stratigraphy
tracking algorithm

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager

"""
import csv

import numpy as np


def checks_correct_equation_to_track_stratigraphy(self):
    """If by mistake track_stratigraphy was set as True but a MPM style equation was selected
    this function returns track_stratigraphy to False

    Examples
    --------

    Import the required libraries and create a grid and configure mandatory fields

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import RiverBedDynamics

    >>> grid = RasterModelGrid((5, 5), xy_spacing=100)

    >>> grid.at_node["topographic__elevation"] = np.full(grid.number_of_nodes, 1.0)
    >>> grid.at_node["topographic__elevation"][2] = 0

    >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.40)
    >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.40)
    >>> grid.at_link["surface_water__velocity"] = np.full(grid.number_of_links, 0.40)
    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])

    Case 1: We choose parker1990 equation and set track_stratigraphy=True.
    So, after instantiation if we check rbd._track_stratigraphy it should be True

    >>> rbd = RiverBedDynamics(
    ...     grid,
    ...     bedload_equation="Parker1990",
    ...     track_stratigraphy=True,
    ... )

    >>> rbd._track_stratigraphy
    True

    Case 2: We choose WilcockAndCrowe equation and set track_stratigraphy=True.
    So, after instantiation if we check rbd._track_stratigraphy it should be True

    >>> rbd = RiverBedDynamics(
    ...     grid,
    ...     bedload_equation="WilcockAndCrowe",
    ...     track_stratigraphy=True,
    ... )

    >>> rbd._track_stratigraphy
    True

    Case 3: We choose FLvB equation and set track_stratigraphy=True. So, after instantiation
    if we check rbd._track_stratigraphy it should be False because that equation does not
    allow tracking the stratigraphy

    >>> rbd = RiverBedDynamics(
    ...     grid,
    ...     bedload_equation="FLvB",
    ...     track_stratigraphy=True,
    ... )

    >>> rbd._track_stratigraphy
    False

    """
    if self._track_stratigraphy:
        if not (
            self._bedload_equation == "Parker1990"
            or self._bedload_equation == "WilcockAndCrowe"
        ):
            self._track_stratigraphy = False


def create_links_dictionary(self):
    """Creates a dictonary containing link's information related to
    surface and subsurface gsd"""
    if self._track_stratigraphy:
        z = self._grid["link"]["topographic__elevation"]
        gsd_F = self._bed_surf__gsd_link
        gsd_Fs = self._bed_subsurf__gsd_link
        active_layer_La = self._bed_surf__act_layer_thick_link
        dzl = self._bed_surf__thick_new_layer_link

        self._link_stratigraphy = {}  # Creates dictionay to store data
        self._link_stratigraphy_temp = {}  # Creates dictionay to store data temporarely

        for i in self._grid.active_links:
            self._link_stratigraphy[i] = []
            self._link_stratigraphy_temp[i] = []

        for i in self._grid.active_links:
            data = np.hstack((self._current_t, z[i], gsd_F[i, :], gsd_Fs[i, :]))
            data = np.reshape(data, [1, data.shape[0]])
            self._link_stratigraphy[i].append(data)

            data = np.hstack(
                (
                    self._current_t,
                    z[i],
                    dzl[i],
                    active_layer_La[i],
                    gsd_F[i, :],
                    gsd_Fs[i, :],
                )
            )
            data = np.reshape(data, [1, data.shape[0]])
            self._link_stratigraphy_temp[i].append(data)


def checks_erosion_or_deposition(self):
    if self._track_stratigraphy:
        # Updates the thickness of the new layers
        self._bed_surf__thick_new_layer_link = (
            self._grid["link"]["topographic__elevation"]
            - self._topogr__elev_subsurf_link
        )
        # Checks if deposited material needs to be updated
        (update_deposited_link_id,) = np.where(
            self._bed_surf__thick_new_layer_link > self._bed_surf_new_layer_thick
        )
        update_deposited_link_id = update_deposited_link_id[
            np.in1d(update_deposited_link_id, self._grid.active_links)
        ]
        self._update_deposited_link_id = update_deposited_link_id[
            ~np.in1d(update_deposited_link_id, self._bed_surf__elev_fix_link_id)
        ]

        if self._update_deposited_link_id.shape[0] > 0:
            self._update_stratigraphy = (
                True  # Allows entering into the stratigraphy function
            )
            self._update_subsurface_deposited = (
                True  # Allows entering into the surface gsd update function
            )

        # Checks if eroded material needs to be updated
        (update_eroded_link_id,) = np.where(
            self._bed_surf__thick_new_layer_link < -self._bed_surf_new_layer_thick
        )
        update_eroded_link_id = update_eroded_link_id[
            np.in1d(update_eroded_link_id, self._grid.active_links)
        ]
        self._update_eroded_link_id = update_eroded_link_id[
            ~np.in1d(update_eroded_link_id, self._bed_surf__elev_fix_link_id)
        ]

        if self._update_eroded_link_id.shape[0] > 0:
            self._update_stratigraphy = (
                True  # Allows entering into the stratigraphy function
            )
            self._update_subsurface_eroded = (
                True  # Allows entering into the surface gsd update function
            )


def evolve(self):
    """This function controls how the stratigraphy is beign updated

    An example where stratigraphy is affected by changes caused by deposition and
    erosion is tested below. This is a very slow test. The conditions of this test
    were designed to check that the code works as expected and therefore may be
    unrealistic.

    Examples
    --------

    As per usual, we define import the required libraries and create a grid and
    configure the mandatory fields

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import RiverBedDynamics
    >>> from . import _stratigraphy as stratigraphy
    >>> import os

    >>> grid = RasterModelGrid((8, 3), xy_spacing=100)

    >>> grid.at_node["topographic__elevation"] = [
    ...     [1.12, 1.00, 1.12],
    ...     [1.12, 1.01, 1.12],
    ...     [1.12, 1.01, 1.12],
    ...     [1.12, 1.01, 1.12],
    ...     [1.12, 1.01, 1.12],
    ...     [1.12, 1.01, 1.12],
    ...     [1.12, 1.01, 1.12],
    ...     [1.12, 1.12, 1.12],
    ... ]

    >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.40)
    >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.40)
    >>> grid.at_link["surface_water__velocity"] = np.full(grid.number_of_links, 0.40)
    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])
    >>> gsd = [[8, 100], [4, 90], [2, 0]]

    >>> fixed_nodes = np.zeros(grid.number_of_nodes)
    >>> fixed_nodes[[1, 4]] = 1

    >>> fixed_bed_gsd_nodes = np.zeros(grid.number_of_nodes)
    >>> fixed_bed_gsd_nodes[[1, 4]] = 1

    >>> qb = np.full(grid.number_of_links, 0.0)
    >>> qb[[28, 33]] = -0.002

    We imposed a bed load rate of 0.002 going south, that is why we used the negative sign.
    At this point we can instantiate the component. We use Parker 1990 eq. MPM equations
    can't be use when tracking the evolution of the surface and subsurface.

    >>> rbd = RiverBedDynamics(
    ...     grid,
    ...     gsd=gsd,
    ...     bedload_equation="Parker1990",
    ...     outlet_boundary_condition="fixedValue",
    ...     bed_surf__elev_fix_node=fixed_nodes,
    ...     bed_surf__gsd_fix_node=fixed_bed_gsd_nodes,
    ...     sed_transp__bedload_rate_fix_link=qb,
    ...     track_stratigraphy=True,
    ...     bed_surf_new_layer_thick=0.02,
    ...     num_cycles_to_process_strat=2,
    ... )

    We will run the model for 1299 s. This is exactly the time required for the first
    link to reach a deposition of 2 cm (Notice bed_surf_new_layer_thick=0.02).
    However, first, we will run until 1298 s, which is 1 second before reaching the
    2 cm threshold

    >>> for t in range(1299):
    ...     rbd._current_t = t
    ...     rbd.run_one_step()
    ...

    Let's take a look at the evolution of the surface gsd of node 23 at times
    0, 25, 115, and 1297. The time indices are 0, 13, 58, and 649 because
    num_cycles_to_process_strat=2

    >>> rbd._link_stratigraphy_temp[23][0][0, 4:6]
    array([ 0.1,  0.9])

    >>> rbd._link_stratigraphy_temp[23][13][0, 4:6]
    array([ 0.1052982,  0.8947018])

    >>> rbd._link_stratigraphy_temp[23][58][0, 4:6]
    array([ 0.12428022,  0.87571978])

    >>> rbd._link_stratigraphy_temp[23][649][0, 4:6]
    array([ 0.3631502,  0.6368498])

    The list rbd._link_stratigraphy_temp[23] has 650 elements, this means that
    it's evolution was recorded every 2 seconds

    >>> len(rbd._link_stratigraphy_temp[23])
    650

    And the elevation after 1298 s.

    >>> z = grid.at_node["topographic__elevation"][16]
    >>> np.round(z, decimals=3)
    1.05

    This is larger that 2 cm, but keep in mind that the calculation is done in
    links, so when mapping, from nodes to links and then from links to nodes,
    small differences appear. This is only important in small grids such as the
    one we are using in this example.

    Now, let's run for one more sec.

    >>> t = 1300
    >>> rbd._current_t = t
    >>> rbd.run_one_step()

    The list rbd._link_stratigraphy_temp[23] has 0 elements, this means that
    a new layer was created

    >>> len(rbd._link_stratigraphy_temp[23])
    0

    Let's see the evolution of the layers in node 23. At time 0 the elevation was:

    >>> rbd._link_stratigraphy[23][0][0, 1]
    1.01

    When was created the new layer? The result is in seconds

    >>> rbd._link_stratigraphy[23][1][0, 0]
    1300.0

    What is the elevation when the new layer was created?

    >>> new_elevation = rbd._link_stratigraphy[23][1][0, 1]
    >>> np.around(new_elevation, decimals=2)
    1.03

    Let's keep running the component, this time with a lower bed load rate such
    that erosion is generated. We will run the model until 3902 s. This is exactly
    the time required for the first link to reach an erosion of 2 cm. However,
    first, we will run until 3901 s, which is 1 second before reaching the
    2 cm threshold

    >>> rbd._sed_transp__bedload_rate_fix_link[[28, 33]] = 0.001
    >>> for t in range(1301, 3901):
    ...     rbd._current_t = t
    ...     rbd.run_one_step()
    ...

    Let's take a look at the evolution of the surface gsd of node 23 at times
    1302, 1330, 1420, and 3900. The time indices are 0, 14, 59, and 1299 because
    this list was deleted previously, when deposition occurred

    >>> rbd._link_stratigraphy_temp[23][0][0, 4:6]
    array([ 0.2368653,  0.7631347])

    >>> rbd._link_stratigraphy_temp[23][14][0, 4:6]
    array([ 0.23001637,  0.76998363])

    >>> rbd._link_stratigraphy_temp[23][59][0, 4:6]
    array([ 0.21015594,  0.78984406])

    >>> rbd._link_stratigraphy_temp[23][1299][0, 4:6]
    array([ 0.09472506,  0.90527494])

    The list rbd._link_stratigraphy_temp[23] has 1300 elements, this means that
    it's evolution was recorded every 2 seconds

    >>> len(rbd._link_stratigraphy_temp[23])
    1300

    And the elevation after 3900 s.

    >>> z = grid.at_node["topographic__elevation"][16]
    >>> np.round(z, decimals=3)
    1.01

    Now, let's run for one more sec.

    >>> t = 3902
    >>> rbd._current_t = t
    >>> rbd.run_one_step()

    The list rbd._link_stratigraphy_temp[23] has 0 elements, this means that
    a new layer was created, or in this case we reached the original subsurface

    >>> len(rbd._link_stratigraphy_temp[23])
    0

    Let's see the evolution of the layers in node 23. At time 0 the elevation was:

    >>> rbd._link_stratigraphy[23][0][0, 1]
    1.01

    When was created the 2nd (deposited) new layer? The result is in seconds

    >>> rbd._link_stratigraphy[23][1][0, 0]
    1300.0

    When was created the 3rd (erored) new layer? The result is in seconds

    >>> rbd._link_stratigraphy[23][2][0, 0]
    3902.0

    So, this node come back to it's original conditions after going through a
    process of deposition and erosion.

    If we want to post-process the stratigraphy data we can write it to a file
    that is by deafult called Stratigraphy_evolution.csv

    >>> stratigraphy.write_evolution(rbd)

    In this case we will delete the file to keep Landlab clean

    >>> os.remove("Stratigraphy_evolution.csv")

    """
    if self._track_stratigraphy:
        # Define variables
        z = self._grid["link"]["topographic__elevation"]
        n_grain_sizes = np.size(self._gsd[:, 0])
        active_layer_La = self._bed_surf__act_layer_thick_link
        dzl = self._bed_surf__thick_new_layer_link
        gsd_F = self._bed_surf__gsd_link
        gsd_Fs = self._bed_subsurf__gsd_link
        self._stratigraphy_cycle += 1

        if self._stratigraphy_cycle >= self._num_cycles_to_process_strat:
            # Add data to dictionary containing temporal information
            for i in self._grid.active_links:
                data = np.hstack(
                    (
                        self._current_t,
                        z[i],
                        dzl[i],
                        active_layer_La[i],
                        gsd_F[i, :],
                        gsd_Fs[i, :],
                    )
                )
                data = np.reshape(data, [1, data.shape[0]])
                self._link_stratigraphy_temp[i].append(data)
                self._stratigraphy_cycle = 0

        # Here we update stratigraphy in case of deposition
        # Enters this cycle only when sediment deposited on a link is larger
        # (deeper) than the user-defined layer thickness
        if self._update_subsurface_deposited:
            for i in self._update_deposited_link_id:
                # Gather information from links that reached the new layer depth
                link_data = np.vstack(self._link_stratigraphy_temp[i])
                mean_subsurface_gsd = np.mean(
                    link_data[:, 4 : 4 + n_grain_sizes - 1], axis=0
                )  # 4 is the number of elements before gsd_F GSD
                self._bed_subsurf__gsd_link[i, :] = mean_subsurface_gsd / np.sum(
                    mean_subsurface_gsd
                )
                self._link_stratigraphy_temp[i] = []

                # Now updates the link history
                data = np.hstack(
                    (
                        self._current_t,
                        z[i],
                        gsd_F[i, :],
                        self._bed_subsurf__gsd_link[i, :],
                    )
                )
                data = np.reshape(data, [1, data.shape[0]])
                self._link_stratigraphy[i].append(data)

                # Updates the surface gsd
                self._bed_surf__gsd_link[i, :] = self._bed_subsurf__gsd_link[i, :]

            # A new layer was created, therefore the thickness is returned to ~zero
            self._bed_surf__thick_new_layer_link[self._update_deposited_link_id] = (
                self._bed_surf__thick_new_layer_link[self._update_deposited_link_id]
                - self._bed_surf_new_layer_thick
            )
            self._topogr__elev_subsurf_link[self._update_deposited_link_id] = (
                self._topogr__elev_subsurf_link[self._update_deposited_link_id]
                + self._bed_surf_new_layer_thick
            )

        # Here we update the stratigraphy in case of erosion
        # Enters this cycle only when sediment eroded on a link is lower (scoured)
        # than the user-defined layer thickness
        if self._update_subsurface_eroded:
            # The layer has been eroded, so we delete the temporal information
            for i in self._update_eroded_link_id:
                self._link_stratigraphy_temp[i] = []

            # Looks for the layer's information in the previous (in time) layer
            for i in self._update_eroded_link_id:
                link_data = np.vstack(self._link_stratigraphy[i])
                if link_data.shape[0] > 1:
                    self._bed_subsurf__gsd_link[i, :] = link_data[
                        -2, (2 + self._gsd.shape[0] - 1) : None
                    ]

                # Now updates the link history
                data = np.hstack(
                    (
                        self._current_t,
                        z[i],
                        gsd_F[i, :],
                        self._bed_subsurf__gsd_link[i, :],
                    )
                )
                data = np.reshape(data, [1, data.shape[0]])
                self._link_stratigraphy[i].append(data)

                # Updates the surface gsd
                self._bed_surf__gsd_link[i, :] = self._bed_subsurf__gsd_link[i, :]

            # A new layer was created, therefore the thickness is returned to ~zero
            self._bed_surf__thick_new_layer_link[self._update_eroded_link_id] = (
                self._bed_surf__thick_new_layer_link[self._update_eroded_link_id]
                + self._bed_surf_new_layer_thick
            )
            self._topogr__elev_subsurf_link[self._update_eroded_link_id] = (
                self._topogr__elev_subsurf_link[self._update_eroded_link_id]
                - self._bed_surf_new_layer_thick
            )

        self._update_subsurface_deposited = False
        self._update_subsurface_eroded = False


def write_evolution(self):
    "Writes the stratigraphy time evolution into a csv file"

    resultsDict = self._link_stratigraphy.copy()
    for key in resultsDict:
        for i, arr in enumerate(resultsDict[key]):
            # Calculate the positions where 666 should be inserted
            insert_pos_1 = 1 + self._gsd.shape[0]
            modified_arr = np.insert(arr[0], insert_pos_1, 0)

            insert_pos_2 = modified_arr.size  # After the first insertion
            modified_arr = np.insert(modified_arr, insert_pos_2, 0)

            # Update the array in the dictionary
            resultsDict[key][i] = modified_arr.reshape(1, -1)

    gsd_values = self._gsd[:, 0].tolist()  # Extract the first column values as a list
    header = (
        ["Link #", "Index", "Time [s]", "z [m]"] + gsd_values + gsd_values
    )  # Create the header

    with open("Stratigraphy_evolution.csv", "w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow(header)

        # Iterate over the dictionary and write rows to the CSV file
        for link, arrays in resultsDict.items():
            for i, arr in enumerate(arrays):
                # Convert each array to a list and flatten it for writing
                data = arr.flatten().tolist()
                row = [link, i] + data
                writer.writerow(row)
