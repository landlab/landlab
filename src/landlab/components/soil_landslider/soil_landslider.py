#!/usr/bin/env python3
"""Grid-based simulation of soil landslides.
This component simulates soil landslides by accomplishing the following tasks:
Calculating the stability of cells according to a factor of safety equations
Triggering  stochastic landslides at critical cells according to a temporal probability
Propagating landslide failure erosional areas uphill of critical cells,
entirely within the uppermost soil layer and
Depositing landslide material downhill according to a multi-path flow routing algorithm

This component is built on the architecture of the BedrockLandslider component
which was devoloped by Benjamin Campforts

Paul Morgan
"""

#### comments with four numpads are by me


import numpy as np

from landlab import Component
from landlab.grid.nodestatus import NodeStatus

from ..bedrock_landslider.cfuncs import _landslide_runout
from ..depression_finder.lake_mapper import _FLOODED

MAX_HEIGHT_SLOPE = 100  # in m


class SoilLandsliderGeo(Component):

    _name = "SoilLandsliderGeo"

    _unit_agnostic = True

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        # Note that this field has to be provided in addition to the \
        # flow__receiver_node and will be used to route sediments over the hillslope
        "hill_flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        # Note that this field has to be provided in addition to the \
        # flow__receiver_proportions and will be used to route sediments
        # over the hillslope
        "hill_flow__receiver_proportions": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
        "hill_topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "LS_sediment__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "node",
            "doc": "Sediment flux originating from landslides \
                (volume per unit time of sediment entering each node)",
        },
        "landslide__erosion": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Total erosion caused by landsliding ",
        },
        "landslide__deposition": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Total deposition of derived sediment",
        },
        "landslide_sediment_point_source": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3",
            "mapping": "node",
            "doc": "Landslide derived sediment, as point sources on all the \
                critical nodes where landslides initiate, \
                before landslide runout is calculated ",
        },
        "rel__wetness": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "frac",
            "mapping": "node",
            "doc": "fraction of vertical water column saturated by groundwater ",
        },
        "factor_of_safety": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "none",
            "mapping": "node",
            "doc": "Factor of safety calculation output ",
        },
    }

    def __init__(
        self,
        grid,
        angle_int_frict=1.0,
        threshold_slope=None,
        cohesion_eff=4e4,
        landslides_return_time=1e5,
        rho_s=2000,
        rho_h2o=1000,
        grav=9.81,
        fraction_fines_LS=0,
        phi=0,
        max_pixelsize_landslide=1e9,
        seed=2021,
        verbose_landslides=False,
        landslides_on_boundary_nodes=True,
        critical_sliding_nodes=None,
        min_deposition_slope=0,
        Rw=0.8,
        Rw_mode="constant",
        max_step_height=False,
        stable_stopper_factor=5,
    ):
        """Initialize the SoilLandsliderGeo model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        angle_int_frict: float, optional
            Materials angle of internal friction in [m/m]
        threshold_slope: float, optional
            Threshold slope used in non-linear deposition scheme [m/m]
            Default value is set to angle_int_frict if not specified
        cohesion_eff : float, optional
            Effective cohesion of material [m L^-1 T^-2].
        landslides_return_time : float, optional
            Return time for stochastic landslide events to occur [T]
        rho_s : float, optional
            Bulk density sediment/soil [m L^-3].
        rho_h20 : float, optional
            Bulk density water [m L^-3].
        fraction_fines_LS : float
            Fraction of permanently suspendable fines in bedrock
            Value must be between 0 and 1 [-].
        phi : float, optional
            Sediment porosity, value must be between 0 and 1 [-].
        max_pixelsize_landslide : int , optional
            Maximum size for landslides in number of pixels
        verbose_landslides : bool , optional
            Print output as number of simulated landslides per timestep
        seed : float , optional
            Provide seed to set stochastic model.
            If not provided, seed is set to 2021.
            Provide None to keep current seed.
        landslides_on_boundary_nodes : bool, optional
            Allow landslides to initiate (critical node) and extend over
            boundary nodes.
        critical_sliding_nodes : list, optional
            Provide list with critical nodes where landslides have to initiate
            This cancels the stochastic part of the algorithm and allows the
            user to form landslides at the provided critical nodes.
        Rw_mode : str, optional
                How to assign Rw values:
                "constant": use the provided Rw value for all nodes (default).
                "calculated": read Rw values from the grid field "rel__wetness".
                None or unrecognized: defaults to constant Rw = 0.8
        Rw : float,optional
            Wetness fraction of  height water column over height soil column (vertical)
            Used only if Rw_mode="constant". Value must be between 0 and 1 [-]
        max_step_height : float or False, optional
            Maximum height of the bedrock elevation above the failure plane for
            the slide to erode the soil on top and continue
            In meters Defaults to False
        stable_stopper_factor: float, optional
            The factor of safety that is considered too stable for a slide to
            propagate uphill through. Note FS<1 is still required to initiate failure.
            Defaults to 5
        """
        super().__init__(grid)

        topo = self.grid.at_node["topographic__elevation"]
        soil = self.grid.at_node["soil__depth"]

        if "bedrock__elevation" not in grid.at_node:
            grid.add_field("bedrock__elevation", topo - soil, at="node", dtype=float)

        # Check consistency of bedrock, soil and topographic elevation fields
        if not np.allclose(
            grid.at_node["bedrock__elevation"] + grid.at_node["soil__depth"],
            grid.at_node["topographic__elevation"],
        ):
            raise RuntimeError(
                "The sum of bedrock elevation and topographic elevation should be equal"
            )

        self.initialize_output_fields()

        # Store grid and parameters
        self._angle_int_frict = angle_int_frict
        if threshold_slope is None:
            self._threshold_slope = angle_int_frict
        else:
            self._threshold_slope = threshold_slope
        self._cohesion_eff = cohesion_eff
        self._rho_s = rho_s
        self._rho_h2o = rho_h2o
        self._grav = grav
        self._fraction_fines_LS = fraction_fines_LS
        self._phi = phi
        self._landslides_return_time = landslides_return_time
        self._max_pixelsize_landslide = max_pixelsize_landslide
        self._verbose_landslides = verbose_landslides
        self._landslides_on_boundary_nodes = landslides_on_boundary_nodes
        self._critical_sliding_nodes = critical_sliding_nodes
        self._min_deposition_slope = min_deposition_slope
        self._Rw = Rw  # see below for edits
        self._max_step_height = max_step_height
        self._stable_stopper_factor = stable_stopper_factor

        # Data structures to store properties of simulated landslides.
        self._landslides_size = []
        self._landslides_volume = []
        self._landslides_volume_sed = []
        #### removed landslide voume bedrock

        # Make sure the relative wetness values/grid exist/work and are less than one
        #
        if Rw_mode == "calculated":
            if "rel__wetness" not in grid.at_node:
                raise ValueError(
                    "Rw_mode='calculated' requires a grid field 'rel__wetness' at nodes."
                )
            relwetness = grid.at_node["rel__wetness"]

        elif Rw_mode == "constant":
            relwetness = grid.add_zeros(
                "rel__wetness", at="node", dtype=float, clobber=True
            )
            relwetness[:] = Rw

        else:
            # Default case: nothing provided, assume value set in defaults above
            relwetness = grid.add_zeros(
                "rel__wetness", at="node", dtype=float, clobber=True
            )
            relwetness[:] = Rw

        # Check input values
        if phi >= 1.0 or phi < 0.0:
            raise ValueError(f"Porosity must be between 0 and 1 ({phi})")

        if fraction_fines_LS > 1.0 or fraction_fines_LS < 0.0:
            raise ValueError(
                f"Fraction of fines must be between 0 and 1 ({fraction_fines_LS})"
            )

        if max_step_height is not False and not isinstance(
            max_step_height, (int, float)
        ):
            raise TypeError(
                f"max_step_height must be False or a number, got {type(max_step_height)}"
            )

        if isinstance(max_step_height, bool) and max_step_height is not False:
            # catches the edge case where someone passes True
            raise TypeError("max_step_height must be False or a number, not True")

        # Set seed
        if seed is not None:
            np.random.seed(seed)

    # Getters for properties
    @property
    def landslides_size(self):
        """
        List with the size of simulated landslides.
        The list is reset every time the _landslide_erosion function is called
        """
        return self._landslides_size

    @property
    def landslides_volume(self):
        """
        List with the volume of simulated landslides.
        The list is reset every time the _landslide_erosion function is called
        """
        return self._landslides_volume

    @property
    def landslides_volume_sed(self):
        """
        List with the volume of sediment eroded by landslides.
        The list is reset every time the _landslide_erosion function is called
        """
        return self._landslides_volume_sed

    def _landslide_erosion(self, dt):
        """
        Calculate soil landsliding for a time period 'dt'.
        Based on the infinite slope factor of safety calculations

        Parameters
        ----------
        dt: float
            The imposed timestep.

        Returns
        -------
        suspended_sed : float
            Volume of suspended sediment.

        """
        # troubleshoot tracking
        # print("starting landslide erosion")

        # Pointers
        topo = self.grid.at_node["topographic__elevation"]
        bed = self.grid.at_node["bedrock__elevation"]
        steepest_slope = self.grid.at_node["topographic__steepest_slope"]
        relwetness = self.grid.at_node["rel__wetness"]
        relwetness[:] = np.clip(relwetness, 0, 1.0)  # clip this at 1

        #### calculate slope in degrees for calculations
        slope = np.rad2deg(np.arctan(steepest_slope))
        soil_d = self.grid.at_node["soil__depth"]
        landslide_sed_in = self.grid.at_node["landslide_sediment_point_source"]
        landslide__ero = self.grid.at_node["landslide__erosion"]

        # Reset LS Plains
        landslide__ero.fill(0.0)
        # Reset landslide sediment point source field
        landslide_sed_in.fill(0.0)

        # Reset data structures to store properties of simulated landslides.
        self._landslides_size = []
        self._landslides_volume = []
        self._landslides_volume_sed = []
        #### remove bedrock landslide volume again

        # Identify flooded nodes
        flood_status = self.grid.at_node["flood_status_code"]
        flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]

        # In the following section the location of critical nodes where
        # landsldies are initatated is calcualted, unless these critical nodes
        # are provided as critical_sliding_nodes
        if self._critical_sliding_nodes is None:
            # find the critical sliding nodes
            # start by calculating the quick infinite slope factor of safety
            # at each node
            # equation goes to zero if there's no soil

            # factor of safety equation
            shearstrenth = self._cohesion_eff + self._grav * soil_d * (
                np.cos(np.deg2rad(slope))
            ) ** 2 * (self._rho_s - relwetness * self._rho_h2o) * np.tan(
                np.arctan(self._angle_int_frict)
            )
            shearstress = (
                self._rho_s
                * self._grav
                * soil_d
                * (np.cos(np.deg2rad(slope)))
                * (np.sin(np.deg2rad(slope)))
            )
            FS = np.divide(
                shearstrenth,
                shearstress,
                where=shearstress > 0,
                out=np.zeros_like(shearstress),
            )
            self.grid.at_node["factor_of_safety"] = FS  # save FS for everwhere

            # sliding (nodes where a slide may be triggered) is set by the temporal probability,
            # but only occurs if the cell is unstable
            # Temporal probability
            temporal_prob = 1 - np.exp(-dt / self._landslides_return_time)
            sliding = np.random.rand(FS.size) < temporal_prob
            # stability is checked at the start of the erosion loop

            # Now, find the critical node, which is the receiver of critical_landslide_nodes
            # Critical nodes must be unique (a given node can have more receivers...)
            critical_landslide_nodes = np.unique(
                self.grid.at_node["flow__receiver_node"][np.where(sliding)]
            )

            # Remove boundary nodes
            if not self._landslides_on_boundary_nodes:
                critical_landslide_nodes = critical_landslide_nodes[
                    ~self.grid.node_is_boundary(critical_landslide_nodes)
                ]
        else:
            critical_landslide_nodes = np.array(self._critical_sliding_nodes)
            # still need to calculate FS though cause we're saving it later.
            shearstrenth = self._cohesion_eff + self._grav * soil_d * (
                np.cos(np.deg2rad(slope))
            ) ** 2 * (self._rho_s - relwetness * self._rho_h2o) * np.tan(
                np.arctan(self._angle_int_frict)
            )
            shearstress = (
                self._rho_s
                * self._grav
                * soil_d
                * (np.cos(np.deg2rad(slope)))
                * (np.sin(np.deg2rad(slope)))
            )
            # FS= shearstrenth/ shearstress
            FS = np.divide(
                shearstrenth,
                shearstress,
                where=shearstress > 0,
                out=np.zeros_like(shearstress),
            )
            self.grid.at_node["factor_of_safety"] = FS

        # output variables
        suspended_sed = 0.0
        if self._verbose_landslides:
            print(f"nbSlides = {len(critical_landslide_nodes)}")

        store_cumul_volume = 0.0

        # troubleshoot tracking
        # print("looping through cells to erode")

        # this while loop loops through all the critical nodes
        while critical_landslide_nodes.size > 0:
            # work with the first of the list
            crit_node = critical_landslide_nodes[0]  # start at first critical node
            crit_node_el = topo[crit_node]
            crit_node_bed = bed[crit_node]  # bedrock elevation at the critical node
            crit_node_soil_d = soil_d[crit_node]  # soil depth at crit node
            crit_node_rel_wetness = relwetness[
                crit_node
            ]  # relative wetness at crit node

            # get 8 neighbors and only keep those to active nodes which are upstream
            neighbors = np.concatenate(
                (
                    self.grid.active_adjacent_nodes_at_node[crit_node],
                    self.grid.diagonal_adjacent_nodes_at_node[crit_node],
                )
            )
            neighbors = neighbors[neighbors != -1]
            neighbors_up = neighbors[topo[neighbors] > crit_node_el]

            x_crit_node = self.grid.node_x[crit_node]
            y_crit_node = self.grid.node_y[crit_node]

            # distance to the uphill neighbors (only the first ones)
            dist_to_initial_node = np.sqrt(
                np.add(
                    np.square(x_crit_node - self.grid.node_x[neighbors_up]),
                    np.square(y_crit_node - self.grid.node_y[neighbors_up]),
                )
            )
            # slope to the first uphill neighbors
            slope_neighbors_to_crit_node = (
                topo[neighbors_up] - crit_node_el
            ) / dist_to_initial_node

            #####
            if self._verbose_landslides:
                print(f"slope_neighbors_to_crit_node {slope_neighbors_to_crit_node}")
            #####

            # here check for the failure plane is stable or unstable FS<1
            # slope failure plane is the steepest slope from the crit node
            slope_fail_critnode_frac = np.max(slope_neighbors_to_crit_node)
            slope_fail_critnode = np.rad2deg(np.arctan(slope_fail_critnode_frac))
            z_fail_critnode = crit_node_soil_d
            relwetness

            # troubleshooting
            # breakpoint()

            shearstrenth = self._cohesion_eff + self._grav * z_fail_critnode * (
                np.cos(np.deg2rad(slope_fail_critnode))
            ) ** 2 * (self._rho_s - crit_node_rel_wetness * self._rho_h2o) * np.tan(
                np.arctan(self._angle_int_frict)
            )
            shearstress = (
                self._rho_s
                * self._grav
                * z_fail_critnode
                * (np.cos(np.deg2rad(slope_fail_critnode)))
                * (np.sin(np.deg2rad(slope_fail_critnode)))
            )
            # FS= shearstrenth/ shearstress
            FS_critnode = np.divide(
                shearstrenth,
                shearstress,
                where=shearstress > 0,
                out=np.zeros_like(shearstress),
            )

            # for troubleshooting
            # print(f'FS_critnode {FS_critnode}')

            # if the critical node is stable then remove it, and restart the loop
            # no landslide there.
            if FS_critnode > 1:
                # then delete this crit node (as long as theres more)
                if critical_landslide_nodes.size > 0:
                    critical_landslide_nodes = np.delete(critical_landslide_nodes, 0, 0)
                continue
            # print(f'FS crit node {FS_critnode}')

            # heres where the ground failure happens
            # identify the sliding angle

            if slope_neighbors_to_crit_node.size > 0:
                # slope slide is the angle to the uphill
                slope_slide = max(slope_neighbors_to_crit_node)
                store_volume_bed = 0.0
                store_volume_sed = 0.0
                upstream_count = 0
                upstream_neighbors = neighbors_up
                if self._verbose_landslides:
                    print(f"upstream_neighbors {upstream_neighbors}")
                if not self._landslides_on_boundary_nodes:
                    upstream_neighbors = upstream_neighbors[
                        ~self.grid.node_is_boundary(upstream_neighbors)
                    ]
                # Fix sliding angle of particular LS
                # sliding angle is the steepest slope to the crit node
                sliding_angle = slope_slide
                nb_landslide_cells = 0

                visited_nodes = set()
                # this stops the code from re-running nodes that aren't
                # eroded all the way to new_el_1
                # If landslides become unrealistically big, exit algorithm
                while upstream_neighbors.size > 0 and (
                    upstream_count <= self._max_pixelsize_landslide
                    and nb_landslide_cells < 1e5
                ):
                    # if the slide has been visited already skip it!
                    if upstream_neighbors[0] in visited_nodes:
                        upstream_neighbors = np.delete(upstream_neighbors, 0, 0)
                        continue

                    visited_nodes.add(
                        upstream_neighbors[0]
                    )  # save that we've run this node

                    distance_to_crit_node = np.sqrt(
                        np.add(
                            np.square(
                                x_crit_node - self.grid.node_x[upstream_neighbors[0]]
                            ),
                            np.square(
                                y_crit_node - self.grid.node_y[upstream_neighbors[0]]
                            ),
                        )
                    )
                    #### new_el is sliding angle up from base of soil at crit node
                    #### every cell has a new_el_1

                    new_el_1 = crit_node_bed + distance_to_crit_node * sliding_angle
                    nb_landslide_cells += 1
                    # troubleshooting
                    # print(f'new_el_1 {new_el_1} upstream_neighbors[0] {upstream_neighbors[0]}')
                    # print(f'new_el_1 {new_el_1} topo upstream_neighbors[0] {topo[upstream_neighbors[0]]}')

                    # this if block is for if the cell topo elevation is within
                    # the failure envelope
                    if new_el_1 < topo[upstream_neighbors[0]]:
                        current_node = upstream_neighbors[0]
                        # add in the max step height constraint
                        # (this may be removed in future versions)
                        if self._max_step_height is not False and bed[
                            [upstream_neighbors[0]]
                        ] > (new_el_1 + self._max_step_height):
                            print("maxstepheight triggered")
                            continue

                        # if the slide hits a too stable cell stop propagating uphill that way slide
                        if FS[current_node] > self._stable_stopper_factor:
                            # print("stablestopper induced")
                            continue

                        # Do actual slide
                        upstream_count += 1
                        sed_landslide_ero = np.clip(
                            min(
                                soil_d[upstream_neighbors[0]],
                                topo[upstream_neighbors[0]] - new_el_1,
                            ),
                            a_min=0,  # trying out changing this from 0 to .1mm to fix a bug
                            a_max=None,
                        )
                        # print(f'cell {upstream_neighbors[0]} sed_landslide_ero {sed_landslide_ero}')

                        soil_d[upstream_neighbors[0]] -= sed_landslide_ero
                        # instead of setting the new elevation to the new el,
                        # set it to how much sediment was eroded
                        topo[upstream_neighbors[0]] = (
                            topo[upstream_neighbors[0]] - sed_landslide_ero
                        )
                        # print(f'cell {upstream_neighbors[0]} updated topo {topo[upstream_neighbors[0]]}')

                        # there is no bedrock erosion
                        vol_sed = (
                            sed_landslide_ero * (1 - self._phi) * (self.grid.dx**2)
                        )
                        # Troubleshooting
                        # print(f'vol_sed {vol_sed}')
                        store_volume_sed = store_volume_sed + vol_sed

                        # update the neighbors to include upstream neighbors to current cell
                        neighbors = np.concatenate(
                            (
                                self.grid.active_adjacent_nodes_at_node[
                                    upstream_neighbors[0]
                                ],
                                self.grid.diagonal_adjacent_nodes_at_node[
                                    upstream_neighbors[0]
                                ],
                            )
                        )
                        neighbors = neighbors[neighbors != -1]
                        neighbors_up = neighbors[topo[neighbors] > crit_node_el]
                        upstream_neighbors = [*upstream_neighbors, *neighbors_up]

                        temp, idx = np.unique(upstream_neighbors, return_index=True)
                        upstream_neighbors = np.array(upstream_neighbors)
                        upstream_neighbors = upstream_neighbors[np.sort(idx)]
                        if not self._landslides_on_boundary_nodes:
                            upstream_neighbors = upstream_neighbors[
                                ~self.grid.node_is_boundary(upstream_neighbors)
                            ]
                        # if one of the LS pixels also appears in critical_landslide_nodes list,
                        # remove it there so that no new landslide is initialized
                        critical_landslide_nodes = critical_landslide_nodes[
                            np.where(critical_landslide_nodes != upstream_neighbors[0])
                        ]

                        #### only sediment erosion
                        landslide__ero[current_node] = sed_landslide_ero
                        # Troubleshooting
                        # print(f'cell {current_node} sed_landslide_ero {sed_landslide_ero}')

                    upstream_neighbors = np.delete(upstream_neighbors, 0, 0)

                store_volume = store_volume_sed + store_volume_bed
                store_cumul_volume += store_volume
                if upstream_count > 0:
                    landslide_sed_in[crit_node] += (store_volume / dt) * (
                        1.0 - self._fraction_fines_LS
                    )
                    suspended_sed += (store_volume / dt) * self._fraction_fines_LS

                    self._landslides_size.append(upstream_count)
                    self._landslides_volume.append(store_volume)
                    self._landslides_volume_sed.append(store_volume_sed)
                    #### no bedload

            if critical_landslide_nodes.size > 0:
                critical_landslide_nodes = np.delete(critical_landslide_nodes, 0, 0)

        if self._verbose_landslides:
            print(f"nbSlides = {len(self._landslides_size)}")

        return suspended_sed

    def _landslide_runout(self, dt):
        """
        Calculate landslide runout using a non-local deposition algorithm based on:
        * Carretier S., Martinod P., Reich M., Godderis Y. (2016) Modelling
          sediment clasts transport during landscape evolution.
          Earth Surf Dyn: 4(1):237–51.
        * Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.
          (2020) HyLands 1.0: a hybrid landscape evolution model to simulate
          the impact of landslides and landslide-derived sediment on landscape
          evolution. Geosci Model Dev: 13(9):3863–86.

        Parameters
        ----------
        dt : float
            Timestep.

        Returns
        -------
        dh_hill : float
            Hillslope erosion over the simulated domain.
        volume_leaving : float
            Total volume of sediment leaving the simulated domain.
        flux_core_nodes : float
            Sediment flux over the simulated domain.

        """
        # troubleshoot tracking
        # print("starting landslide runout")

        topo = self.grid.at_node["topographic__elevation"]
        bed = self.grid.at_node["bedrock__elevation"]
        soil_d = self.grid.at_node["soil__depth"]
        sed_flux = self.grid.at_node["LS_sediment__flux"]
        stack_rev = np.flip(self.grid.at_node["flow__upstream_node_order"])
        landslide_depo = self.grid.at_node["landslide__deposition"]
        landslide_sed_in = self.grid.at_node["landslide_sediment_point_source"]
        node_status = self.grid.status_at_node

        # Only process core nodes
        stack_rev_sel = stack_rev[node_status[stack_rev] == NodeStatus.CORE]
        receivers = self.grid.at_node["hill_flow__receiver_node"]
        fract_receivers = self.grid.at_node["hill_flow__receiver_proportions"]

        # keep only steepest slope
        slope = np.max(self.grid.at_node["hill_topographic__steepest_slope"], axis=1)
        slope[slope < 0] = 0.0

        flux_in = landslide_sed_in * dt  # flux_in, in m3 per timestep

        # L following carretier 2016
        transport_length_hill = np.where(
            slope < self._threshold_slope,
            self.grid.dx / (1 - (slope / self._threshold_slope) ** 2),
            1e6,
        )

        flux_out = np.zeros(topo.shape)
        dh_hill = np.zeros(topo.shape)
        topo_copy = np.array(topo)
        max_depo = np.zeros(topo.shape)
        length_adjacent_cells = np.array(
            [
                self.grid.dx,
                self.grid.dx,
                self.grid.dx,
                self.grid.dx,
                self.grid.dx * np.sqrt(2),
                self.grid.dx * np.sqrt(2),
                self.grid.dx * np.sqrt(2),
                self.grid.dx * np.sqrt(2),
            ]
        )

        _landslide_runout(
            self.grid.dx,
            self._phi,
            self._min_deposition_slope,
            stack_rev_sel,
            receivers,
            fract_receivers,
            flux_in,
            transport_length_hill,
            flux_out,
            dh_hill,
            topo_copy,
            max_depo,
            length_adjacent_cells,
        )
        sed_flux[:] = flux_out

        flux_core_nodes = np.sum(flux_in[self.grid.status_at_node == 0])
        volume_leaving = np.sum(flux_in)  # Qs_leaving # in m3 per timestep

        # Change sediment layer
        soil_d[:] += dh_hill
        topo[:] = bed + soil_d

        # Reset Qs
        landslide_sed_in.fill(0.0)
        # Update deposition field
        landslide_depo[:] = dh_hill

        return dh_hill, volume_leaving, flux_core_nodes

    def run_one_step(self, dt):
        """Advance SoilLandsliderGeo component by one time step of size dt.

        Parameters
        ----------
        dt: float
            The imposed timestep.

        Returns
        -------
        vol_suspended_sediment_yield : float
            volume of sediment evacuated as syspended sediment.
        volume_leaving : float
            Volume of sediment leaving the domain.
        """
        dt = float(dt)

        if self.current_time is None:
            self.current_time = dt
        else:
            self.current_time += dt

        # Landslides
        vol_suspended_sediment_yield = self._landslide_erosion(dt)
        dh_hill, volume_leaving, flux_core_nodes = self._landslide_runout(dt)

        return vol_suspended_sediment_yield, volume_leaving
