"""Grid-based simulation of bedrock landslides.

Benjamin Campforts
"""

import numpy as np

from landlab import Component
from landlab.grid.nodestatus import NodeStatus

from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import _landslide_runout

MAX_HEIGHT_SLOPE = 100  # in m


class BedrockLandslider(Component):
    """Calculate the location and magnitude of episodic bedrock landsliding.

    Landlab component that calculates the location and magnitude of episodic
    bedrock landsliding following the Cullman criterion.
    See the publication:

    Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.
    (2020) HyLands 1.0: a hybrid landscape evolution model to simulate the
    impact of landslides and landslide-derived sediment on landscape evolution.
    Geosci Model Dev: 13(9):3863–86.
    `https://dx.doi.org/10.5194/esurf-6-1-2018 <https://dx.doi.org/10.5194/esurf-6-1-2018>`_

    Campforts, B., Shobe, C. M., Overeem, I., & Tucker, G. E. (2022).
    The Art of Landslides: How Stochastic Mass Wasting Shapes Topography and
    Influences Landscape Dynamics.
    Journal of Geophysical Research: Earth Surface,
    127(8), 1–16. https://doi.org/10.1029/2022JF006745


    Examples
    --------

    >>> import numpy as np
    >>> from numpy import testing
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import PriorityFloodFlowRouter, BedrockLandslider

    Make a ``RasterModelGrid`` and create a plateau.

    * 5x5 grid
    * Initial topography is set to plateau value of 10

    >>> mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> s = mg.add_zeros("soil__depth", at="node")
    >>> b = mg.add_zeros("bedrock__elevation", at="node")

    Make plateau at 10 m

    >>> b += 10

    Lower boundary cell to 0

    >>> b[2] = 0
    >>> z[:] = b + s

    Instantiate the :class:`~.priority_flood_flow_router.PriorityFloodFlowRouter`
    for flow accumulation and the ``BedrockLandslider``

    >>> fd = PriorityFloodFlowRouter(
    ...     mg,
    ...     separate_hill_flow=True,
    ...     suppress_out=True,
    ... )
    >>> hy = BedrockLandslider(mg, landslides_return_time=1)

    Run the flow director and ``BedrockLandslider`` for one timestep

    >>> fd.run_one_step()
    >>> vol_suspended_sediment_yield, volume_leaving = hy.run_one_step(dt=1)

    After one timestep, we can predict exactly where the landslide will occur.
    The return time is set to 1 year so that probability for sliding is 100%.
    The angle of internal friction is 1 m/m, the topographical gradient is 10 m/m.
    At cardinal cells, the sliding plane will be at *(1 + 10) / 2 = 5.5* m/m.
    With a *dx* of 1, the cardinal cell next to the critical sliding node must
    be 5.5 m and the diagonal one at *5.5 * sqrt(2) = 7.8* m

    >>> testing.assert_almost_equal(
    ...     [5.5 * np.sqrt(2), 5.5, 5.5 * np.sqrt(2)], z[6:9], decimal=5
    ... )

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.
    (2020) BedrockLandslider 1.0: a hybrid landscape evolution model to simulate the
    impact of landslides and landslide-derived sediment on landscape evolution.
    Geosci Model Dev: 13(9):3863–86.
    `https://dx.doi.org/10.5194/esurf-6-1-2018 <https://dx.doi.org/10.5194/esurf-6-1-2018>`_

    **Additional References**

    None Listed

    """

    _name = "BedrockLandslider"

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
    }

    _cite_as = """
    @Article{gmd-13-3863-2020,
        AUTHOR = {Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.},
        TITLE = {BedrockLandslider 1.0: a hybrid landscape evolution model to
                 simulate the impact of landslides and landslide-derived sediment
                 on landscape evolution.},
        JOURNAL = {Geoscientific Model Development},
        VOLUME = {13},
        YEAR = {2020},
        NUMBER = {9},
        PAGES = {3863--3886},
        URL = {https://doi.org/10.5194/gmd-13-3863-2020},
        DOI = {10.5194/gmd-13-3863-2020}
    }"""

    def __init__(
        self,
        grid,
        angle_int_frict=1.0,
        threshold_slope=None,
        cohesion_eff=1e4,
        landslides_return_time=1e5,
        rho_r=2700,
        grav=9.81,
        fraction_fines_LS=0,
        phi=0,
        max_pixelsize_landslide=1e9,
        seed=2021,
        verbose_landslides=False,
        landslides_on_boundary_nodes=True,
        critical_sliding_nodes=None,
        min_deposition_slope=0,
    ):
        """Initialize the BedrockLandslider model.

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
            Return time for stochastic landslide events to occur
        rho_r : float, optional
            Bulk density rock [m L^-3].
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
        self._rho_r = rho_r
        self._grav = grav
        self._fraction_fines_LS = fraction_fines_LS
        self._phi = phi
        self._landslides_return_time = landslides_return_time
        self._max_pixelsize_landslide = max_pixelsize_landslide
        self._verbose_landslides = verbose_landslides
        self._landslides_on_boundary_nodes = landslides_on_boundary_nodes
        self._critical_sliding_nodes = critical_sliding_nodes
        self._min_deposition_slope = min_deposition_slope

        # Data structures to store properties of simulated landslides.
        self._landslides_size = []
        self._landslides_volume = []
        self._landslides_volume_sed = []
        self._landslides_volume_bed = []

        # Check input values
        if phi >= 1.0 or phi < 0.0:
            raise ValueError(f"Porosity must be between 0 and 1 ({phi})")

        if fraction_fines_LS > 1.0 or fraction_fines_LS < 0.0:
            raise ValueError(
                f"Fraction of fines must be between 0 and 1 ({fraction_fines_LS})"
            )

        # Set seed
        if seed is not None:
            np.random.seed(seed)

    # Getters for properties
    @property
    def fraction_fines(self):
        """
        Fraction of permanently suspendable fines in bedrock.
        Value must be between 0 and 1 [-].
        """
        return self._fraction_fines_LS

    @property
    def phi(self):
        """
        Sediment porosity, value must be between 0 and 1 [-].
        """
        return self._phi

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

    @property
    def landslides_volume_bed(self):
        """
        List with the volume of bedrock eroded by landslides.
        The list is reset every time the _landslide_erosion function is called
        """
        return self._landslides_volume_bed

    def _landslide_erosion(self, dt):
        """
        Calculate bedrock landsliding for a time period 'dt'.

        Parameters
        ----------
        dt: float
            The imposed timestep.

        Returns
        -------
        suspended_sed : float
            Volume of suspended sediment.

        """
        # Pointers
        topo = self.grid.at_node["topographic__elevation"]
        bed = self.grid.at_node["bedrock__elevation"]
        steepest_slope = self.grid.at_node["topographic__steepest_slope"]
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
        self._landslides_volume_bed = []

        # Identify flooded nodes
        flood_status = self.grid.at_node["flood_status_code"]
        flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]

        # In the following section the location of critical nodes where
        # landsldies are initatated is calcualted, unless these critical nodes
        # are provided as critical_sliding_nodes
        if self._critical_sliding_nodes is None:
            # Calculate gradients
            height_cell = topo - topo[self.grid.at_node["flow__receiver_node"]]

            height_cell[flooded_nodes] = 0
            height_cell[height_cell > MAX_HEIGHT_SLOPE] = MAX_HEIGHT_SLOPE

            angle_int_frict_radians = np.arctan(self._angle_int_frict)
            height_critical = np.divide(
                (4 * self._cohesion_eff / (self._grav * self._rho_r))
                * (np.sin(np.arctan(steepest_slope)) * np.cos(angle_int_frict_radians)),
                1 - np.cos(np.arctan(steepest_slope) - angle_int_frict_radians),
                where=(1 - np.cos(np.arctan(steepest_slope) - angle_int_frict_radians))
                > 0,
                out=np.zeros_like(steepest_slope),
            )
            spatial_prob = np.divide(
                height_cell,
                height_critical,
                where=height_critical > 0,
                out=np.zeros_like(height_critical),
            )
            spatial_prob[np.arctan(steepest_slope) <= angle_int_frict_radians] = 0
            spatial_prob[spatial_prob > 1] = 1

            # Temporal probability
            temporal_prob = 1 - np.exp(-dt / self._landslides_return_time)

            # Combined probability
            combined_prob = temporal_prob * spatial_prob
            sliding = np.random.rand(combined_prob.size) < combined_prob

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

        # output variables
        suspended_sed = 0.0
        if self._verbose_landslides:
            print(f"nbSlides = {len(critical_landslide_nodes)}")

        store_cumul_volume = 0.0
        while critical_landslide_nodes.size > 0:
            crit_node = critical_landslide_nodes[0]  # start at first critical node
            crit_node_el = topo[crit_node]

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

            dist_to_initial_node = np.sqrt(
                np.add(
                    np.square(x_crit_node - self.grid.node_x[neighbors_up]),
                    np.square(y_crit_node - self.grid.node_y[neighbors_up]),
                )
            )
            slope_neighbors_to_crit_node = (
                topo[neighbors_up] - crit_node_el
            ) / dist_to_initial_node

            neighbors_up = neighbors_up[
                slope_neighbors_to_crit_node > self._angle_int_frict
            ]
            slope_neighbors_to_crit_node = slope_neighbors_to_crit_node[
                slope_neighbors_to_crit_node > self._angle_int_frict
            ]

            if slope_neighbors_to_crit_node.size > 0:
                slope_slide = max(slope_neighbors_to_crit_node)
                store_volume_bed = 0.0
                store_volume_sed = 0.0
                upstream_count = 0
                upstream_neighbors = neighbors_up
                if not self._landslides_on_boundary_nodes:
                    upstream_neighbors = upstream_neighbors[
                        ~self.grid.node_is_boundary(upstream_neighbors)
                    ]
                # Fix sliding angle of particular LS
                sliding_angle = (self._angle_int_frict + slope_slide) / 2.0
                nb_landslide_cells = 0

                # If landslides become unrealistically big, exit algorithm
                while upstream_neighbors.size > 0 and (
                    upstream_count <= self._max_pixelsize_landslide
                    and nb_landslide_cells < 1e5
                ):
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
                    new_el = crit_node_el + distance_to_crit_node * sliding_angle
                    nb_landslide_cells += 1
                    if new_el < topo[upstream_neighbors[0]]:
                        # Do actual slide
                        upstream_count += 1
                        sed_landslide_ero = np.clip(
                            min(
                                soil_d[upstream_neighbors[0]],
                                topo[upstream_neighbors[0]] - new_el,
                            ),
                            a_min=0.0,
                            a_max=None,
                        )
                        soil_d[upstream_neighbors[0]] -= sed_landslide_ero
                        topo[upstream_neighbors[0]] = new_el
                        bed_landslide_ero = np.clip(
                            bed[upstream_neighbors[0]]
                            - (new_el - soil_d[upstream_neighbors[0]]),
                            a_min=0.0,
                            a_max=None,
                        )
                        bed[upstream_neighbors[0]] -= bed_landslide_ero
                        topo[upstream_neighbors[0]] = new_el

                        vol_sed = (
                            sed_landslide_ero * (1 - self._phi) * (self.grid.dx**2)
                        )
                        vol_bed = bed_landslide_ero * (self.grid.dx**2)
                        store_volume_sed = store_volume_sed + vol_sed
                        store_volume_bed = store_volume_bed + vol_bed

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

                        landslide__ero[upstream_neighbors[0]] = (
                            sed_landslide_ero + bed_landslide_ero
                        )

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
                    self._landslides_volume_bed.append(store_volume_bed)

            if critical_landslide_nodes.size > 0:
                critical_landslide_nodes = np.delete(critical_landslide_nodes, 0, 0)

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
        """Advance BedrockLandslider component by one time step of size dt.

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
