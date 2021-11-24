# -*- coding: utf-8 -*-
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
    """
    Landlab component that calculates the location and magnitude of episodic
    bedrock landsliding following the Cullman criterion.
    See the publication:

    Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.
    (2020) HyLands 1.0: a hybrid landscape evolution model to simulate the
    impact of landslides and landslide-derived sediment on landscape evolution.
    Geosci Model Dev: 13(9):3863–86.
    `https://dx.doi.org/10.5194/esurf-6-1-2018 <https://dx.doi.org/10.5194/esurf-6-1-2018>`_


    Examples
    --------
    Make a raster model grid and create a plateau

    >>> import numpy as np
    >>> from numpy import testing
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import PriorityFloodFlowRouter, BedrockLandslider

    Make a raster model grid and create a plateau
    * 5x5 grid
    * Initial topography is set to plateau value of 10

    >>> nr = 5
    >>> nc = 5
    >>> dx = 1
    >>> mg = RasterModelGrid((nr, nc), xy_spacing=dx)
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> s = mg.add_zeros("soil__depth", at='node')
    >>> b = mg.add_zeros("bedrock__elevation", at='node')

    make plateau at 10m

    >>> b += 10

    Lower boundary cell to 0

    >>> b[2] = 0
    >>> z[:] = b + s

    Instantiate flow accumulation and BedrockLandslider

    >>> fd = PriorityFloodFlowRouter(
    ...     mg,
    ...     separate_hill_flow=True,
    ...     suppress_out=True,
    ... )
    >>> hy = BedrockLandslider(mg, landslides_return_time=1)

    run flow director and BedrockLandslider for one timestep
    >>> fd.run_one_step()
    >>> vol_SSY ,V_leaving = hy.run_one_step(dt=1)

    After one timestep, we can predict eactly where the landslide will occur.
    The return time is set to 1 year so that probability for sliding is 100%
    The angle of internal driction is 1m/m, the topographical gradient is 10 m/m
    At cardinal cells, the sliding plain will be at (1+10)/2 = 5.5 m/m.
    With a dx of 1, the cardial cell next to the critical sliding node must
    be 5.5 m and the diagonal one at 5.5 * sqrt(2) = 7.8 m

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
        "landslide__bed_erosion": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "erosion caused by bedrock landsliding ",
        },
        "landslide__deposition": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "deposition of bedrock derived sediment",
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

    _cite_as = """@Article{gmd-13-3863-2020,
                  AUTHOR = {Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.},
                  TITLE = {BedrockLandslider 1.0: a hybrid landscape evolution model to simulate the impact of landslides and landslide-derived sediment on landscape evolution.},
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
        store_landslides_size=True,
        store_landslides_volume=True,
        store_landslides_volume_sed=False,
        store_landslides_volume_bed=False,
    ):
        """Initialize the BedrockLandslider model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        angle_int_frict: float, optional
            Materials angle of internal friction in [m/m]
        cohesion_eff : float, optional
            Effective cohesion of material [m L^-1 T^-2].
        landslides_return_time  : float, optional
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
        store_landslides_size : bool  , optional
            Store the size of landslides through time. Values will be stored in
            a list as a property of a bedrock_landslider object.
        store_landslides_volume : bool , optional
            Store the volume of landslides through time. Values will be stored in
            a list as a property of a bedrock_landslider object.
        store_landslides_volume_sed  : bool , optional
            Store the volume of eroded bedrock through time. Values will be
            stored in a list as a property of a bedrock_landslider object.
        store_landslides_volume_bed : bool  , optional
            Store the volume of eroded sediment landslides through time.
            Values will be stored in a list as a property of a
            bedrock_landslider object.
        """
        super(BedrockLandslider, self).__init__(grid)

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

        # Data structures to store properties of simulated landslides.
        self._landslides_size = [] if store_landslides_size else None
        self._landslides_volume = [] if store_landslides_volume else None
        self._landslides_volume_sed = [] if store_landslides_volume_sed else None
        self._landslides_volume_bed = [] if store_landslides_volume_bed else None

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
    def landslides_size(self):
        return self._landslides_size

    @property
    def landslides_volume(self):
        return self._landslides_volume

    @property
    def landslides_volume_sed(self):
        return self._landslides_volume_sed

    @property
    def landslides_volume_bed(self):
        return self._landslides_volume_bed

    def _landslide_erosion(self, dt):
        """
        Calculate bedrock landsliding for a time period 'dt'.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.

        Returns
        -------
        suspended_Sed : float
            Volume of suspended sediment.

        """
        # Pointers
        topo = self.grid.at_node["topographic__elevation"]
        bed = self.grid.at_node["bedrock__elevation"]
        ss = self.grid.at_node["topographic__steepest_slope"]
        soil_d = self.grid.at_node["soil__depth"]
        ls_sed_in = self.grid.at_node["landslide_sediment_point_source"]
        lb_bed_ero = self.grid.at_node["landslide__bed_erosion"]

        # Reset LS Plains
        lb_bed_ero.fill(0.0)
        # Reset landslide sediment point source field
        ls_sed_in.fill(0.0)

        # Identify flooded nodes
        flood_status = self.grid.at_node["flood_status_code"]
        flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]

        """In the following section the location of critical nodes where
        landsldies are initatated is calcualted, unless these critical nodes
        are provided as critical_sliding_nodes"""
        if self._critical_sliding_nodes is None:
            # Calculate gradients
            H_el = topo - topo[self.grid.at_node["flow__receiver_node"]]

            H_el[flooded_nodes] = 0
            H_el[H_el > MAX_HEIGHT_SLOPE] = MAX_HEIGHT_SLOPE

            sc_rad = np.arctan(self._angle_int_frict)
            Hc = np.divide(
                (4 * self._cohesion_eff / (self._grav * self._rho_r))
                * (np.sin(np.arctan(ss)) * np.cos(sc_rad)),
                1 - np.cos(np.arctan(ss) - sc_rad),
                where=(1 - np.cos(np.arctan(ss) - sc_rad)) > 0,
                out=np.zeros_like(ss),
            )
            p_S = np.divide(H_el, Hc, where=Hc > 0, out=np.zeros_like(Hc))
            p_S[np.arctan(ss) <= sc_rad] = 0
            p_S[p_S > 1] = 1

            # Temporal probability
            p_t = 1 - np.exp(-dt / self._landslides_return_time)

            # Combined probability
            p = p_t * p_S
            slides = np.random.rand(p.size) < p

            # Now, find the critical node, which is the receiver of i_slide
            # Critical nodes must be unique (a given node can have more receivers...)
            i_slide = np.unique(
                self.grid.at_node["flow__receiver_node"][np.where(slides)]
            )
            # Remove boundary nodes
            if not self._landslides_on_boundary_nodes:
                i_slide = i_slide[~self.grid.node_is_boundary(i_slide)]
        else:
            i_slide = np.array(self._critical_sliding_nodes)

        # output variables
        suspended_Sed = 0.0
        if self._verbose_landslides:
            print(f"nbSlides = {len(i_slide)}")

        storeV_cum = 0.0
        while i_slide.size > 0:
            ind = 0
            cP = i_slide[ind]  # Critical Node
            cP_el = topo[cP]

            # get 8 neighbors and only keep those to active nodes which are upstream
            nb = np.concatenate(
                (
                    self.grid.active_adjacent_nodes_at_node[cP],
                    self.grid.diagonal_adjacent_nodes_at_node[cP],
                )
            )
            nb = nb[nb != -1]
            nb_up = nb[topo[nb] > cP_el]

            X_cP = self.grid.node_x[cP]
            Y_cP = self.grid.node_y[cP]

            distToIni_all = np.sqrt(
                np.add(
                    np.square(X_cP - self.grid.node_x[nb_up]),
                    np.square(Y_cP - self.grid.node_y[nb_up]),
                )
            )
            all_iP_el = topo[nb_up]
            s_slide_all = (all_iP_el - cP_el) / distToIni_all

            nb_up = nb_up[s_slide_all > self._angle_int_frict]
            s_slide_all = s_slide_all[s_slide_all > self._angle_int_frict]

            if s_slide_all.size > 0:
                s_slide = max(s_slide_all)
                storeV_bed = 0.0
                storeV_sed = 0.0
                upstream = 0
                uP = nb_up
                if not self._landslides_on_boundary_nodes:
                    uP = uP[~self.grid.node_is_boundary(uP)]
                # Fix sliding angle of particular LS
                ang_sl = (self._angle_int_frict + s_slide) / 2.0
                stall = 0

                while uP.size > 0 and (
                    upstream <= self._max_pixelsize_landslide and stall < 1e4
                ):
                    distToIni = np.sqrt(
                        np.add(
                            np.square(X_cP - self.grid.node_x[uP[0]]),
                            np.square(Y_cP - self.grid.node_y[uP[0]]),
                        )
                    )
                    newEl = cP_el + distToIni * ang_sl
                    stall += 1
                    if newEl < topo[uP[0]]:
                        # Do actual slide
                        upstream = upstream + 1
                        sed_LS_E = np.clip(
                            min(soil_d[uP[0]], topo[uP[0]] - newEl),
                            a_min=0.0,
                            a_max=None,
                        )
                        soil_d[uP[0]] -= sed_LS_E
                        topo[uP[0]] = newEl
                        bed_LS_E = np.clip(
                            bed[uP[0]] - (newEl - soil_d[uP[0]]), a_min=0.0, a_max=None
                        )
                        bed[uP[0]] -= bed_LS_E
                        topo[uP[0]] = newEl

                        vol_sed = sed_LS_E * (1 - self._phi) * (self.grid.dx ** 2)
                        vol_bed = bed_LS_E * (self.grid.dx ** 2)
                        storeV_sed = storeV_sed + vol_sed
                        storeV_bed = storeV_bed + vol_bed

                        nb = np.concatenate(
                            (
                                self.grid.active_adjacent_nodes_at_node[uP[0]],
                                self.grid.diagonal_adjacent_nodes_at_node[uP[0]],
                            )
                        )
                        nb = nb[nb != -1]
                        nb_up = nb[topo[nb] > cP_el]
                        uP = [*uP, *nb_up]

                        temp, idx = np.unique(uP, return_index=True)
                        uP = np.array(uP)
                        uP = uP[np.sort(idx)]
                        if not self._landslides_on_boundary_nodes:
                            uP = uP[~self.grid.node_is_boundary(uP)]
                        # if one of the LS pixels also appears in i_slide list,
                        # remove it there so that no new landslide is initialized
                        i_slide = i_slide[np.where((i_slide != uP[0]))]

                        lb_bed_ero[uP[0]] = sed_LS_E + bed_LS_E

                    uP = np.delete(uP, 0, 0)

                storeV = storeV_sed + storeV_bed
                storeV_cum += storeV
                if upstream > 0:
                    ls_sed_in[cP] += (storeV / dt) * (1.0 - self._fraction_fines_LS)
                    suspended_Sed += (storeV / dt) * self._fraction_fines_LS

                    if self._landslides_size is not None:
                        self._landslides_size.append(upstream)
                    if self._landslides_volume is not None:
                        self._landslides_volume.append(storeV)
                    if self._landslides_volume_sed is not None:
                        self._landslides_volume_sed.append(storeV_sed)
                    if self._landslides_volume_bed is not None:
                        self._landslides_volume_bed.append(storeV_bed)

            if i_slide.size > 0:
                i_slide = np.delete(i_slide, 0, 0)

        return suspended_Sed

    def _landslide_runout(self, dt):
        """
        calculate landslide runout using a non-local deposition algorithm based on:
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
            timestep.

        Returns
        -------
        dH_Hill : float
            Hillslope erosion over the simulated domain.
        V_leaving : float
            Total volume of sediment leaving the simulated domain.
        Qs_coreNodes : float
            Sediment flux over the simulated domain.

        """
        z = self.grid.at_node["topographic__elevation"]
        br = self.grid.at_node["bedrock__elevation"]
        H = self.grid.at_node["soil__depth"]
        sed_flux = self.grid.at_node["LS_sediment__flux"]
        stack_rev = np.flip(self.grid.at_node["flow__upstream_node_order"])
        ls_depo = self.grid.at_node["landslide__deposition"]
        ls_sed_in = self.grid.at_node["landslide_sediment_point_source"]
        node_status = self.grid.status_at_node

        # Only process core nodes
        stack_rev_sel = stack_rev[node_status[stack_rev] == NodeStatus.CORE]
        receivers = self.grid.at_node["hill_flow__receiver_node"]
        fract = self.grid.at_node["hill_flow__receiver_proportions"]

        # keep only steepest slope
        slope = np.max(self.grid.at_node["hill_topographic__steepest_slope"], axis=1)
        slope[slope < 0] = 0.0

        Qs_in = ls_sed_in * dt  # Qs_in, in m3 per timestep

        # L following carretier 2016
        L_Hill = np.matlib.divide(
            self.grid.dx,
            (
                1
                - np.matlib.minimum(
                    np.matlib.square(np.matlib.divide(slope, self._angle_int_frict)),
                    0.999,
                )
            ),
        )
        Qs_out = np.zeros(z.shape)
        dH_Hill = np.zeros(z.shape)
        H_i_temp = np.array(z)
        max_D = np.zeros(z.shape)

        _landslide_runout(
            self.grid.dx,
            self._phi,
            stack_rev_sel,
            receivers,
            fract,
            Qs_in,
            L_Hill,
            Qs_out,
            dH_Hill,
            H_i_temp,
            max_D,
        )
        sed_flux[:] = Qs_out

        Qs_coreNodes = np.sum(Qs_in[self.grid.status_at_node == 0])
        V_leaving = np.sum(Qs_in)  # Qs_leaving # in m3 per timestep

        # Change sediment layer
        H[:] += dH_Hill
        z[:] = br + H

        # Reset Qs
        ls_sed_in.fill(0.0)
        # Update deposition field
        ls_depo[:] = dH_Hill

        return dH_Hill, V_leaving, Qs_coreNodes

    def run_one_step(self, dt):
        """Advance BedrockLandslider component by one time step of size dt.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.

        Returns
        -------
        vol_SSY : float
            volume of sediment evacuated as syspended sediment.
        V_leaving : float
            Volume of sediment leaving the domain.
        """
        dt = float(dt)

        if self.current_time is None:
            self.current_time = dt
        else:
            self.current_time += dt

        # Landslides
        vol_SSY = self._landslide_erosion(dt)
        dH_Hill, V_leaving, Qs_coreNodes = self._landslide_runout(dt)

        return vol_SSY, V_leaving
