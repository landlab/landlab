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


class Hylands(Component):
    """ Deep-seated bedrock landsliding removing overlying soil and bedrock
    following the Cullman criterion

    Landlab component that calculates the location and magnitude of episodic
    bedrock landsliding.
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
    >>> from landlab.components import PriorityFloodFlowRouter, Hylands

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

    Instantiate flow accumulation and hylands

    >>> fd = PriorityFloodFlowRouter(
    ...     mg,
    ...     separate_hill_flow=True,
    ...     suppress_out=True)
    >>> hy = Hylands(mg, landslides_return_time=1)

    run flow director and Hylands for one timestep
    >>> fd.run_one_step()
    >>> vol_SSY ,V_leaving = hy.run_one_step(dt=1)

    After one timestep, we can predict eactly where the landslide will occur.
    The return time is set to 1 year so that probability for sliding is 100%
    The angle of internal driction is 1m/m, the topographical gradient is 10 m/m
    At cardinal cells, the sliding plain will be at (1+10)/2 = 5.5 m/m.
    With a dx of 1, the cardial cell next to the critical sliding node must
    be 5.5 m and the diagonal one at 5.5 * sqrt(2) = 7.8 m

    >>> err_msg ='Error in the calculation of the sliding plain'
    >>> testing.assert_almost_equal([5.5 * np.sqrt(2), 5.5, 5.5 * np.sqrt(2)],\
                                    z[6:9],decimal=5, err_msg=err_msg)

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.
    (2020) HyLands 1.0: a hybrid landscape evolution model to simulate the
    impact of landslides and landslide-derived sediment on landscape evolution.
    Geosci Model Dev: 13(9):3863–86.
    `https://dx.doi.org/10.5194/esurf-6-1-2018 <https://dx.doi.org/10.5194/esurf-6-1-2018>`_

    **Additional References**

    None Listed

    """

    _name = "Hylands"

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
        "bedrock__elevation": {
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
        "hill_flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "hill_flow__receiver_proportions": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "hill_topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "sediment__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "node",
            "doc": "Sediment flux (volume per unit time of sediment entering each node)",
        },
        "landslide__bed_erosion": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "landslide__bed_erosion",
        },
        "landslide__deposition": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "landslide__bed_erosion",
        },
    }

    _cite_as = """@Article{gmd-13-3863-2020,
                  AUTHOR = {Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.},
                  TITLE = {HyLands 1.0: a hybrid landscape evolution model to simulate the impact of landslides and landslide-derived sediment on landscape evolution.},
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
        landslides_size=None,
        landslides_volume=None,
        landslides_volume_sed=None,
        landslides_volume_bed=None,
    ):
        """Initialize the HyLands model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        angle_int_frict: float, optional
            Materials angle of internal friction in [m/m]
            Default = 1.0
        cohesion_eff : float, optional
            Effective cohesion of material [m L^-1 T^-2].
        landslides_return_time  : float, optional
            Return time for stochastic landslide events to occur
            Default = 1e5
        rho_r : float, optional
            Bulk density rock [m L^-3].
        fraction_fines_LS : float
            Fraction of permanently suspendable fines in bedrock
            Value must be between 0 and 1 [-].
        phi : float, optional
            Sediment porosity, value must be between 0 and 1 [-].
        max_pixelsize_landslide : int , optional
            Maximum size for landslides in number of pixels
            Default = 1e9
        verbose_landslides : bool , optional
            Print output as number of simulated landslides per timestep
            Default = False
        seed : float , optional
            Provide seed to set stochastic model.
            If not provided, seed is set to 2021.
            Provide None to keep current seed.
            Default = 2021
        landslides_on_boundary_nodes : bool, optional
            Allow landslides to initiate (critical node) and extend over
            boundary nodes.
            Default = True
        critical_sliding_nodes : list, optional
            Provide list with critical nodes where landslides have to initiate
            This cancels the stochastic part of the algorithm and allows the
            user to form landslides at the provided critical nodes.
            Default = None
        landslides_size : list  , optional
            Store the size of landslides through time
            Pass empty list ('[]') to store sizes, pass None to not store any data
            Default : None
        landslides_volume : list , optional
            Store the volume of landslides through time.
            Pass empty list ('[]') to store volumes, pass None to not store any data
            Default = None
        landslides_volume_sed  : list , optional
            Store the volume of eroded bedrock through time
            Pass empty list ('[]') to store volumes, pass None to not store any data
            Default = None
        landslides_volume_bed : list  , optional
            Store the volume of eroded sediment landslides through time
            Pass empty list ('[]') to store volumes, pass None to not store any data
            Default = None

        Returns
        -------
        vol_SSY : float
            volume of sediment evacuated as syspended sediment.
        V_leaving : float
            Volume of sediment leaving the domain.

        """
        super(Hylands, self).__init__(grid)

        # Check consistency of bedrock, soil and topogarphic elevation fields
        err_msg = (
            "The sum of bedrock elevation and topographic elevation should be equal"
        )
        np.testing.assert_almost_equal(
            grid.at_node["bedrock__elevation"] + grid.at_node["soil__depth"],
            grid.at_node["topographic__elevation"],
            decimal=5,
            err_msg=err_msg,
        )

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
        self._Qs_change_LS = np.zeros(self.grid.number_of_nodes)

        # Set seed
        if seed is not None:
            np.random.seed = seed
        self.landslides_on_boundary_nodes = landslides_on_boundary_nodes

        self.critical_sliding_nodes = critical_sliding_nodes

        # Data structure to store properties of simulated landslides.
        # If None, data will not be stored
        self.landslides_size = landslides_size
        self.landslides_volume = landslides_volume
        self.landslides_volume_sed = landslides_volume_sed
        self.landslides_volume_bed = landslides_volume_bed

        # Check input values
        if phi >= 1.0:
            raise ValueError("Porosity must be < 1.0")

        if fraction_fines_LS > 1.0:
            raise ValueError("Fraction of fines must be <= 1.0")

        if phi < 0.0:
            raise ValueError("Porosity must be > 0.0")

        if fraction_fines_LS < 0.0:
            raise ValueError("Fraction of fines must be > 0.0")

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

        dx = self.grid.dx

        # Reset LS Plains
        self.grid.at_node["landslide__bed_erosion"] = np.zeros(
            self.grid.number_of_nodes
        )
        self.grid.at_node["landslide__deposition"] = np.zeros(self.grid.number_of_nodes)

        # identify flooded nodes
        flood_status = self.grid.at_node["flood_status_code"]
        flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]

        """In the following section the location of critical nodes where
        landsldies are initatated is calcualted, unless these critical nodes
        are provided as critical_sliding_nodes"""
        if self.critical_sliding_nodes is None:
            # Calculate gradients
            H_el = topo - topo[self.grid.at_node["flow__receiver_node"]]

            H_el[flooded_nodes] = 0
            H_el[H_el > MAX_HEIGHT_SLOPE] = MAX_HEIGHT_SLOPE

            sc_rad = np.arctan(self._angle_int_frict)
            Hc = (4 * self._cohesion_eff / (self._grav * self._rho_r)) * (
                (np.sin(np.arctan(ss)) * np.cos(sc_rad))
                / (1 - np.cos(np.arctan(ss) - sc_rad))
            )
            Hc[Hc == 0] = np.nan
            # Spatial probability
            p_S = H_el / Hc
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
            if not self.landslides_on_boundary_nodes:
                i_slide = i_slide[(self.grid.node_is_boundary(i_slide) is False)]
        else:
            i_slide = np.array(self.critical_sliding_nodes)

        # output variables
        suspended_Sed = np.float64(0)
        if self._verbose_landslides:
            print("nbSlides = " + str(len(i_slide)))

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
                storeV_bed = np.float64(0.0)
                storeV_sed = np.float64(0.0)
                upstream = 0
                uP = nb_up
                if not self.landslides_on_boundary_nodes:
                    uP = uP[(self.grid.node_is_boundary(uP) is False)]
                # Fix sliding angle of particular LS
                ang_sl = np.float64((self._angle_int_frict + s_slide) / 2)
                stall = 0

                while uP.size > 0 & (
                    upstream <= self._max_pixelsize_landslide and stall < 1e4
                ):
                    # print(uP.size)
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
                        sed_LS_E = np.float64(
                            max(0.0, min(soil_d[uP[0]], topo[uP[0]] - newEl))
                        )
                        soil_d[uP[0]] -= sed_LS_E
                        topo[uP[0]] = newEl
                        bed_LS_E = max(0, bed[uP[0]] - (newEl - soil_d[uP[0]]))
                        bed[uP[0]] -= bed_LS_E
                        topo[uP[0]] = newEl

                        vol_sed = sed_LS_E * (1 - self._phi) * (dx * dx)
                        vol_bed = bed_LS_E * (dx * dx)
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

                        temp, idx = np.unique(uP, "first")
                        uP = np.array(uP)
                        uP = uP[np.sort(idx)]
                        if not self.landslides_on_boundary_nodes:
                            uP = uP[(self.grid.node_is_boundary(uP) is False)]
                        # if one of the LS pixels also appears in i_slide list,
                        # remove it there so that no new landslide is initialized
                        i_slide = i_slide[np.where((i_slide != uP[0]))]

                        self.grid.at_node["landslide__bed_erosion"][uP[0]] = (
                            sed_LS_E + bed_LS_E
                        )

                    uP = np.delete(uP, 0, 0)

                storeV = storeV_sed + storeV_bed
                storeV_cum += storeV
                if upstream > 0:
                    self._Qs_change_LS[cP] += np.float64(
                        (storeV / dt) * (1 - self._fraction_fines_LS)
                    )
                    suspended_Sed += np.float64((storeV / dt) * self._fraction_fines_LS)

                    if self.landslides_size is not None:
                        self.landslides_size.append(upstream)
                    if self.landslides_volume is not None:
                        self.landslides_volume.append(storeV)
                    if self.landslides_volume_sed is not None:
                        self.landslides_volume_sed.append(storeV_sed)
                    if self.landslides_volume_bed is not None:
                        self.landslides_volume_bed.append(storeV_bed)

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
        stack_rev = np.flip(self.grid.at_node["flow__upstream_node_order"])
        node_status = self.grid.status_at_node
        # Only process core nodes
        stack_rev_sel = stack_rev[(node_status[stack_rev] == NodeStatus.CORE)]
        receivers = self.grid.at_node["hill_flow__receiver_node"]
        fract = self.grid.at_node["hill_flow__receiver_proportions"]

        # keep only steepest slope
        slope = np.max(self.grid.at_node["hill_topographic__steepest_slope"], axis=1)
        slope[slope < 0] = 0

        Qs_in = self._Qs_change_LS * dt  # Qs_in, in m3 per timestep

        dx = self.grid.dx
        # L following carretier
        L_Hill = np.matlib.divide(
            dx,
            (
                1
                - np.matlib.minimum(
                    np.matlib.square(np.matlib.divide(slope, self._angle_int_frict)),
                    0.999,
                )
            ),
        )
        Qs_out = np.zeros(z.shape)
        # Node that this results in zero _landslide_runout at watershed divides
        dH_Hill = np.zeros(z.shape)
        H_i_temp = np.array(z)
        max_D = np.zeros(z.shape)
        max_dH = np.ones(z.shape) + np.inf

        phi = self._phi

        _landslide_runout(
            dx,
            phi,
            stack_rev_sel,
            receivers,
            fract,
            Qs_in,
            L_Hill,
            Qs_out,
            dH_Hill,
            H_i_temp,
            max_D,
            max_dH,
        )

        Qs_coreNodes = np.sum(Qs_in[self.grid.status_at_node == 0])
        V_leaving = np.sum(Qs_in)  # Qs_leaving # in m3 per timestep

        # Change sediment layer
        H[:] += dH_Hill[:]
        z[:] = br[:] + H[:]

        # Reset Qs
        self._Qs_change_LS = np.zeros(self.grid.number_of_nodes)
        self.grid.at_node["landslide__deposition"][:] = dH_Hill

        return dH_Hill, V_leaving, Qs_coreNodes

    def run_one_step(self, dt):
        """Advance hylands component by one time step of size dt.

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

        if self.current_time is None:
            self.current_time = dt
        else:
            self.current_time += dt

        dt = np.float64(dt)

        # Landslides
        vol_SSY = self._landslide_erosion(dt)
        dH_Hill, V_leaving, Qs_coreNodes = self._landslide_runout(dt)

        return vol_SSY, V_leaving
