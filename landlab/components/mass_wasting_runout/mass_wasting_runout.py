import warnings

import numpy as np
import pandas as pd

from landlab import Component
from landlab.components import FlowDirectorMFD
from landlab.components.mass_wasting_runout.mass_wasting_saver import MassWastingSaver


class MassWastingRunout(Component):
    """a cellular-automata mass wasting runout model that routes an initial mass
    wasting body (e.g., a landslide) through a watershed, determines erosion and
    aggradation depths, evolves the terrain and regolith and tracks attributes of
    the regolith. This model is intended for modeling the runout extent, topographic
    change and sediment transport caused by a mapped landslide(s) or landslides
    inferred from a landslide hazard map.


    Examples
    ----------
    Import necessary packages and components

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorMFD
    >>> from landlab.components.mass_wasting_runout import MassWastingRunout

    Define the topographic__elevation field of a 7 columns by 7 rows, 10-meter
    raster model grid

    >>> dem = np.array(
    ...     [
    ...         [10, 8, 4, 3, 4, 7.5, 10],
    ...         [10, 9, 3.5, 4, 5, 8, 10],
    ...         [10, 9, 6.5, 5, 6, 8, 10],
    ...         [10, 9.5, 7, 6, 7, 9, 10],
    ...         [10, 10, 9.5, 8, 9, 9.5, 10],
    ...         [10, 10, 10, 10, 10, 10, 10],
    ...         [10, 10, 10, 10, 10, 10, 10],
    ...     ]
    ... )

    >>> dem = np.hstack(dem).astype(float)
    >>> mg = RasterModelGrid((7, 7), 10)
    >>> _ = mg.add_field("topographic__elevation", dem, at="node")

    Define boundary conditions

    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    Add multiflow direction fields, soil thickness (here set to 1 meter)

    >>> fd = FlowDirectorMFD(mg, diagonals=True, partition_method="slope")
    >>> fd.run_one_step()
    >>> nn = mg.number_of_nodes
    >>> depth = np.ones(nn) * 1
    >>> _ = mg.add_field("node", "soil__thickness", depth)

    Define the initial landslide. Any mass_wasting_id value >1 is considered a
    landslide. The landslide extent is defined by assigining all nodes withing
    the landslide the same mass_wasting_id value.
    Here, the landslide is represented by a single node (node 38), which assigned
    a mass_wasting_id value of 1:

    >>> mg.at_node["mass__wasting_id"] = np.zeros(nn).astype(int)
    >>> mg.at_node["mass__wasting_id"][np.array([38])] = 1

    Add attributes of the regolith as fields of the raster model grid that will
    be tracked by the model. These could be any attribute in which the tracking
    method used by MassWastingRunout reasonably represents movement of the
    attribute. Here we track the particle diameter and organic content
    of the regolith. Note, a particle__diameter field is required if shear stress
    is determined as a function of grain size.

    >>> np.random.seed(seed=7)
    >>> mg.at_node["particle__diameter"] = np.random.uniform(0.05, 0.25, nn)
    >>> mg.at_node["organic__content"] = np.random.uniform(0.01, 0.10, nn)

    Next define parameter values for MassWastingRunout and instantiate the model:

    >>> Sc = [0.03]  # Sc, note: defined as a list (see below)
    >>> qsc = 0.01  # qsc
    >>> k = 0.02  # k
    >>> h_max = 1
    >>> tracked_attributes = ["particle__diameter", "organic__content"]
    >>> example_square_MWR = MassWastingRunout(
    ...     mg,
    ...     critical_slope=Sc,
    ...     threshold_flux=qsc,
    ...     erosion_coefficient=k,
    ...     tracked_attributes=tracked_attributes,
    ...     effective_qsi=True,
    ...     max_flow_depth_observed_in_field=h_max,
    ...     save=True,
    ... )

    Run MassWastingRunout

    >>> example_square_MWR.run_one_step()

    By subtracting the initial DEM from the final DEM, which has evolvod as
    a consequence of the runout, we can see areas of aggradation (positive values)
    and erosion (negative values). Nodes with non-zero topographic change
    represent the runout extent.

    >>> DEM_initial = mg.at_node["topographic__initial_elevation"]
    >>> DEM_final = mg.at_node["topographic__elevation"]
    >>> DEM_final - DEM_initial
    array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.96579201,
            0.52734339, -0.00778869,  0.        ,  0.        ,  0.        ,
            0.        , -0.00594927, -0.12261762, -0.0027898 ,  0.        ,
            0.        ,  0.        ,  0.        , -0.04562554, -0.10973222,
           -0.05776526,  0.        ,  0.        ,  0.        ,  0.        ,
           -0.01225359, -0.07973101, -0.04888238,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        , -1.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ])

    See how the landslide removes all of the regolith at node 38 (the negative -1)

    Look at the final spatial distribution of regolith particle diameter.

    >>> mg.at_node["particle__diameter"]
    array([0.06526166, 0.20598376, 0.13768185, 0.19469304, 0.2455979 ,
           0.15769917, 0.15022409, 0.06441023, 0.1036878 , 0.15619144,
           0.17680799, 0.21074781, 0.12618823, 0.06318727, 0.10762912,
           0.23191871, 0.09295928, 0.14042479, 0.23545374, 0.05497985,
           0.17010978, 0.2400259 , 0.09606058, 0.15969798, 0.23182567,
           0.07663389, 0.15468252, 0.20008197, 0.18380265, 0.14355057,
           0.09096982, 0.14815318, 0.12447694, 0.14548023, 0.12317808,
           0.2175836 , 0.2037295 , 0.11279894, 0.        , 0.10520981,
           0.14056859, 0.12059567, 0.18147989, 0.12407022, 0.1418186 ,
           0.19386482, 0.13259837, 0.23128465, 0.08609032])

    Also note that the attribute value is set to zero at any node in which the regolith
    depth is 0.


    References
    ----------
    Keck, J., Istanbulluoglu, E., Campforts, B., Tucker G., Horner-Devine A.,
    A landslide runout model for sediment transport, landscape evolution and hazard
    assessment applications, submitted to Earth Surface Dynamics (2023)

    """

    _name = "MassWastingRunout"

    _unit_agnostic = False

    _info = {
        "mass__wasting_id": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "interger or float id of each mass wasting area is assigned \
                to all nodes representing the mass wasting area.",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "soil__thickness": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "soil depth to restrictive layer",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__receiver_proportions": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "particle__diameter": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "representative particle diameter at each node, this might \
            vary with underlying geology, contributing area or field observations",
        },
    }

    def __init__(
        self,
        grid,
        critical_slope=0.05,
        threshold_flux=0.25,
        erosion_coefficient=0.005,
        tracked_attributes=None,
        deposition_rule="critical_slope",
        grain_shear=True,
        effective_qsi=False,
        settle_deposit=False,
        E_constraint=True,
        save=False,
        typical_flow_thickness_of_erosion_zone=2,
        typical_slope_of_erosion_zone=0.15,
        erosion_exponent=0.2,
        max_flow_depth_observed_in_field=None,
        vol_solids_concentration=0.6,
        density_solids=2650,
        density_fluid=1000,
        gravity=9.81,
        dist_to_full_qsc_constraint=0,
        itL=1000,
        run_id=0,
    ):
        """
        Parameters
        ----------
        grid: landlab raster model grid

        critical_slope: list of floats
            critical slope (angle of repose if no cohesion) of mass
            wasting material , L/L list of length 1 for a basin uniform Sc value
            list of length 2 for hydraulic geometry defined Sc, where the first
            and second values in the list are the coefficient and exponent of a
            user defined function for crictical slope that varies with contributing
            area [m2] to a node (e.g., [0.146,0.051] sets Sc<=0.1 at
            contributing area > ~1100 m2).

        threshold_flux: float
            minimum volumetric flux per unit contour width, [L3/L2/iterataion] or
            [L/iteration]. Flux below this threshold stops at the cell as a deposit

        erosion_coefficient: float
            coefficient used to convert total basal shear stress [kPa] of the runout
            material to a scour depth [m]

        tracked_attributes : list of str or None
            A list of the attribute names (strings) that will be tracked by the
            runout model. Attributes in tracked_attributes must also be a field
            on the model grid and names in list must match the grid field names.
            Default is None.

        deposition_rule : str
            Can be either "critical_slope", "L_metric" or "both".
            "critical_slope" is deposition rule used in Keck et al. 2023.
            "L_metric" is a variation of rule described by Campforts et al. 2020.
            "both" uses the minimum value of both rules.
            Default value is "critical_slope".

        grain_shear : bool
            Indicate whether to define shear stress at the base of the runout material
            as a function of grain size using Equation 13 (True) or the depth-slope
            approximation using Equation 12 (False). Default is True.

        effective_qsi : bool
            Indicate whether to limit erosion and aggradation rates to <= the
            erosion and aggradation rates coorisponding to the maximum observed flow
            depth. All results in Keck et al. 2023 use this constraint. Default is True.

        E_constraint : bool
             Indicate if erosion can not simultaneously occur with aggradation. If True,
             aggradation > 0, then erosion = 0. This is True in Keck et al., 2023.
             Default is True.

        settle_deposit : bool
            Indicate whether to allow deposits to settle before the next model iteration
            is implemented. Settlement is determined the critical slope as evaluated from
            the lowest adjacent node to the deposit. This is not used in Keck et al. 2023
            but tends to allow model to better reproduce smooth, evenly sloped deposits.
            Default is False.

        save : bool
            Save topographic elevation of watershed after each model iteration?
            This uses a lot of memory but is helpful for illustrating runout.
            The default is False.


        Other Parameters
        ----------
        These parameters have a lesser impact on model behavior or may not be
        applicable, depending on model run options.

        typical_flow_thickness_of_erosion_zone: float
            field estimated flow thickness in the erosion-dominatedd reaches of
            the runout path [m], used to estimate erosion_coefficient k using
            erosion_coef_k function. Default value: 3

        typical slope_of_erosion_zone: float
            field or remote sensing estimated slope in the scour dominated reach
            of the runout path [L/L], used to estimate erosion_coefficient k using
            erosion_coef_k function. Default value: 0.4

        erosion_exponent: float
            The exponent of equation 11, that scales erosion depth as a function of
            shear stress. Default value: 0.5

        max_flow_depth_observed: float
            Maximum observed flow depth, over the entire
            runout path [m], h_max in equation 24. Only used effective_qsi is True.
            Default value: 4

        vol_solids_concentration: float
            The ratio of the volume of the solids to the total volume of the flow
            mixture. Default value: 0.6

        density solids: float
            The density of the solids [kg/m3]. Default value: 2650

        density fluid: float
            The density of the fluid [kg/m3]. Default value: 1000

        gravity: float
            Acceleration due to gravity [m2/s]. Default value: 9.81

        dist_to_full_qsc_constraint : float
            distance in meters at which qsc is applied to runout. If the landslide
            initiates on relatively flat terrain, it may be difficult to determine
            a qsc value that allows the model start and deposit in a way that matches
            the observed. In Keck et al. 2023, dist_to_full_qsc_constraint = 0, but
            other landslides may need dist_to_full_qsc_constraint = 20 to 50 meters.

        itL : int
            maximum number of iterations the model runs before it
            is forced to stop. The default is 1000. Ideally, if properly parameterized,
            the model should stop on its own. All modeled runout in Keck et al. 2023
            stopped on its own.

        run_id : float, int or str
            label for landslide run, can be the time or some other identifier. This
            can be updated each time model is implemnted with "run_one_step"

        Returns
        -------
        None
        """
        if isinstance(critical_slope, (float, int)):
            critical_slope = (critical_slope,)

        super().__init__(grid)

        if len(critical_slope) > 1:
            self.variable_slpc = True
            self.a = critical_slope[0]
            self.b = critical_slope[1]
        else:
            self.variable_slpc = False
            self.slpc = critical_slope[0]
        self.qsc = threshold_flux
        self.k = erosion_coefficient
        self._tracked_attributes = tracked_attributes
        self.deposition_rule = deposition_rule
        self.grain_shear = grain_shear
        self.effective_qsi = effective_qsi
        self.settle_deposit = settle_deposit
        self.E_constraint = E_constraint
        self.save = save
        self.h = typical_flow_thickness_of_erosion_zone
        self.s = typical_slope_of_erosion_zone
        self.f = erosion_exponent
        self.qsi_max = max_flow_depth_observed_in_field
        if self.effective_qsi and self.qsi_max is None:
            raise ValueError(
                "Need to define the 'max_flow_depth_observed_in_field'"
                " or set effective_qsi to False"
            )
        self.vs = vol_solids_concentration
        self.ros = density_solids
        self.rof = density_fluid
        self.g = gravity
        self.dist_to_full_qsc_constraint = dist_to_full_qsc_constraint
        self.itL = itL
        self.run_id = run_id

        if tracked_attributes:
            self.track_attributes = True

            # check attributes are included in grid
            for key in self._tracked_attributes:
                if not self._grid.has_field(key, at="node"):
                    raise ValueError(f"{key} not included as field in grid")

            # if using grain size dependent erosion, check
            # particle_diameter is included as an attribute
            if (
                self.grain_shear
                and "particle__diameter" not in self._tracked_attributes
            ):
                raise ValueError(
                    "'particle__diameter' not included as field in grid and/or"
                    " key in tracked_attributes"
                )
        else:
            self.track_attributes = False

        # flow routing option
        # 'square_root_of_slope', see flow director
        self.routing_partition_method = "slope"

        # density of runout mixture
        self.ro_mw = self.vs * self.ros + (1 - self.vs) * self.rof
        # number of model iterations needed to reach dist_to_full_qsc_constraint
        self.d_it = int(self.dist_to_full_qsc_constraint / self._grid.dx)

        # define initial topographic + mass wasting thickness topography
        self._grid.at_node["energy__elevation"] = self._grid.at_node[
            "topographic__elevation"
        ].copy()
        self._grid.at_node["topographic__initial_elevation"] = self._grid.at_node[
            "topographic__elevation"
        ].copy()
        # prepare data containers for saving model images and behavior statistics
        if self.save:
            self.saver = MassWastingSaver(self)
            self.saver.prep_data_containers()

    def run_one_step(self):
        """run MWR"""

        # get all nodes that define the mass wasting events
        mask = self._grid.at_node["mass__wasting_id"] > 0

        # separate the mass wasting event nodes into individual events
        self.mw_ids = np.unique(self._grid.at_node["mass__wasting_id"][mask])
        innL = []  # innL is list of lists of nodes in each mass wasting event
        for mw_id in self.mw_ids:
            ls_mask = self._grid.at_node["mass__wasting_id"] == mw_id
            innL.append(np.hstack(self._grid.nodes)[ls_mask])

        # For each mass wasting event in list:
        for mw_i, inn in enumerate(innL):
            mw_id = self.mw_ids[mw_i]
            self._lsvol = (
                self._grid.at_node["soil__thickness"][inn].sum()
                * self._grid.dx
                * self._grid.dy
            )

            # prepare temporary data containers for each mass wasting event mw_i
            if self.save:
                self.saver.prep_mw_data_containers(mw_i, mw_id)

            # Algorithm 1, prepare initial mass wasting material (debritons) for release
            self._prep_initial_mass_wasting_material(inn, mw_i)

            # self.arndn_r[mw_id].append(self.arndn)
            if self.save:
                # save first set of data to reflect scar created by landslide
                self.saver.save_conditions_before_runout(mw_i, mw_id)

            # Algorith 2, now loop through each receiving nodes,
            # determine next set of recieving nodes,
            # repeat until no more receiving nodes (material deposits)
            self.c = 0  # model iteration counter
            while len(self.arn) > 0 and self.c < self.itL:
                # set qsc: the qsc constraint does not fully apply until runout
                # has traveled dist_to_full_qsc_constraint
                if self.d_it == 0:
                    self.qsc_v = self.qsc
                else:
                    self.qsc_v = self.qsc * (min(self.c / self.d_it, 1))

                # temporary data containers for each iteration of the while loop,
                # that store receiving node, flux and attributes to become the
                # input for the next iteration
                self.arndn_ns = np.array([])  # next iteration donor nodes
                self.arn_ns = np.array([])  # next iteration receiver nodes
                self.arqso_ns = np.array([])  # next iteration flux to receiver nodes
                self.arnL = []  # list of receiver nodes
                self.arqsoL = []  # list of flux out
                self.arndnL = []  # list of donor nodes
                if self.track_attributes:
                    self.aratt_ns = dict.fromkeys(
                        self._tracked_attributes, np.array([])
                    )  #
                    self.arattL = dict.fromkeys(self._tracked_attributes, [])

                # for each unique node in receiving node list self.arn
                self.arn_u = np.unique(self.arn).astype(int)  # unique arn list

                # determine the incoming flux to each node in self.arn_u
                self._determine_qsi()

                # update node elevation plus incoming flow thickness
                # this happens even if using topographic__elevation to route so that
                # the thickness of the debris flow is tracked for plotting
                self._update_E_dem()

                # determine erosion, aggradation, qso and attributes,
                # arranged in array nudat
                self._E_A_qso_determine_attributes()

                # from qso and flow direction at node, determine flux and attributes
                # sent to each receiver node
                self._determine_rn_proportions_attributes()

                # update grid field: topographic__elevation with the values in
                # nudat. Do this after directing flow, because assume deposition
                # does not impact flow direction
                self._update_dem()

                # update topographic slope field
                self._update_topographic_slope()

                # update tracked attribute grid fields
                if self._tracked_attributes:
                    for key in self._tracked_attributes:
                        self._update_attribute_at_node(key)

                # optional settlment of deposits and redistribution of attributes
                if self.settle_deposit:
                    self._settle()
                    self._update_topographic_slope()

                # once all nodes in this iteration have been processed, the lists of receiving
                # nodes (arn), donor nodes (arndn, which are the recieving nodes of this step),
                # outgoing node flux (arqso) and node attributes (artt) are updated
                # for the next iteration
                self.arndn = self.arndn_ns.astype(int)
                self.arn = self.arn_ns.astype(int)
                self.arqso = self.arqso_ns  #
                if self.track_attributes:
                    self.aratt = self.aratt_ns

                if self.save:
                    self.saver.save_conditions_after_one_iteration(mw_i, mw_id)

                # update iteration counter
                self.c += 1

    def _prep_initial_mass_wasting_material(self, inn, mw_i):
        """Algorithm 1 - from an initial source area (landslide), prepare the
        initial lists of receiving nodes and incoming fluxes and attributes
        and remove the source material from the DEM

        Parameters
        ----------
        inn: np.array
             node id's that make up the area of the initial mass wasting area
        mw_i: int
            index of the initial mass wasting area (e.g., if there are two landslides
                                                    the first landslide will be mw_i = 0,
                                                    the second will be mw_i = 0)
        """
        # data containers for initial recieving node, outgoing flux and attributes
        rni = np.array([])
        rqsoi = np.array([])
        if self._tracked_attributes:
            att = dict.fromkeys(self._tracked_attributes, np.array([]))

        # order source area nodes from lowest to highest elevation
        node_z = self._grid.at_node.dataset["topographic__elevation"][inn]
        zdf = pd.DataFrame({"nodes": inn, "z": node_z})
        zdf = zdf.sort_values("z")

        for ci, ni in enumerate(zdf["nodes"].values):
            # regolith (soil) thickness at node. soil thickness in source area
            # represents landslide thickness
            s_t = self._grid.at_node.dataset["soil__thickness"].values[ni]

            # remove soil (landslide) thickness at node
            self._grid.at_node.dataset["topographic__elevation"][ni] = (
                self._grid.at_node.dataset["topographic__elevation"][ni] - s_t
            )

            # update soil thickness at node (now = 0)
            self._grid.at_node["soil__thickness"][ni] = (
                self._grid.at_node["soil__thickness"][ni] - s_t
            )

            if (
                ci > 0
            ):  # use surface slope for first node to start movement of landslide
                # for all other nodes, update slope to reflect material removed from DEM
                self._update_topographic_slope()

            # get receiving nodes of node ni in mw index mw_i
            rn = self._grid.at_node.dataset["flow__receiver_node"].values[ni]
            rn = rn[np.where(rn != -1)]

            # receiving proportion of qso from cell n to each downslope cell
            rp = self._grid.at_node.dataset["flow__receiver_proportions"].values[ni]
            rp = rp[np.where(rp > 0)]  # only downslope cells considered

            # initial mass wasting thickness
            imw_t = s_t
            # get flux out of node ni
            qso = imw_t
            # divide into proportions going to each receiving node
            rqso = rp * qso

            if self._tracked_attributes:
                # get initial mass wasting attributes moving (out) of node ni
                self.att_ar_out = {}
                for key in self._tracked_attributes:
                    att_val = self._grid.at_node.dataset[key].values[ni]
                    # particle diameter to each recieving node
                    self.att_ar_out[key] = np.ones(len(rqso)) * att_val

                    # attribute value is zero at node after reglith leaves
                    self._grid.at_node[key][ni] = 0

                    att[key] = np.concatenate((att[key], self.att_ar_out[key]), axis=0)

            # append receiving node ids, fluxes and attributes to initial lists
            rni = np.concatenate((rni, rn), axis=0)
            rqsoi = np.concatenate((rqsoi, rqso), axis=0)

        self.arndn = np.ones([len(rni)]) * np.nan
        self.arn = rni
        self.arqso = rqsoi
        if self._tracked_attributes:
            self.aratt = att

    def _E_A_qso_determine_attributes(self):
        """mass conservation at a grid cell, implemented using the EAqA
        function below.
        """
        self.D_L = []  # list of deposition depths, if the _settle function is

        # implemented, this is used to determine settlement after deposition
        def EAqA(qsi_r):
            """function for iteratively determining Erosion and Aggradiation depths,
            the outgoing flux (qso) and Attribute values at each affected node and
            attribute values of the outgoing flux.

            thoughts for speeding this up: use cython or try to restructure
            EAqA function so that it can be applied with vector operations,
            then filter all changes to correct cells where computations
            should not have occurred

            Parameters
            ----------
            qsi_r: np array
                a row of the self.qsi_dat array

            """

            n = qsi_r[0]
            qsi = qsi_r[1]

            # maximum slope at node n
            slpn = self._grid.at_node["topographic__steepest_slope"][n].max()

            if self.effective_qsi:
                qsi_ = min(qsi, self.qsi_max)
            else:
                qsi_ = qsi

            # look up critical slope at node n
            if (
                self.variable_slpc
            ):  # if option 1, critical slope is not constant but depends on location
                self.slpc = self.a * self._grid.at_node["drainage_area"][n] ** self.b

            # incoming attributes (weighted average)
            if self._tracked_attributes:
                att_in = self._attributes_in(n, qsi)
            else:
                att_in = None

            rn_g = self._grid.at_node.dataset["flow__receiver_node"].values[n]
            rn_g = rn_g[np.where(rn_g != -1)]

            # if qsi less qsc (equation 1)
            if qsi <= (self.qsc_v):
                A = qsi
                qso = 0
                E = 0
                deta = A
                if self._tracked_attributes:
                    att_up = dict.fromkeys(self._tracked_attributes, 0)
                else:
                    att_up = None
                Tau = 0
                u = 0

            else:
                # aggradation
                A = min(qsi, self._aggradation(qsi, slpn, n))
                if A > 0 and self.E_constraint:
                    E = 0
                    if self._tracked_attributes:
                        att_up = dict.fromkeys(self._tracked_attributes, 0)
                    else:
                        att_up = None
                    Tau = 0
                    u = 0
                else:
                    # erosion
                    E, att_up, Tau, u = self._erosion(n, qsi_, slpn, att_in=att_in)

                ## flux out
                qso = qsi - A + E
                # small qso are considered zero
                qso = np.round(qso, decimals=8)

                # chage elevation
                deta = A - E

            # model behavior tracking
            if self.save:
                self.saver.save_flow_stats(E, A, qsi, slpn, Tau, u)

            # updated attribute values at node
            if self._tracked_attributes:
                n_att = self._attributes_node(n, att_in, E, A)
            else:
                n_att = None

            # list of deposition depths at cells in iteration
            self.D_L.append(A)

            # n_att, att_up, att_in are dictionaries of values of each attribute at
            # the node after erosion and deposition, in the eroded material before depostion
            # and incoming material in qsi respectively
            return deta, qso, qsi, E, A, n_att, att_up, att_in

        # apply EAqA function to all unique nodes in arn (arn_u)
        # create nudat, an np.array of data for updating fields at each node
        ll = np.array([EAqA(r) for r in self.qsi_dat], dtype=object)
        arn_ur = np.reshape(self.qsi_dat[:, 0], (-1, 1))
        self.nudat = np.concatenate((arn_ur, ll), axis=1)

    def _determine_rn_proportions_attributes(self):
        """determine how outgoing flux is partitioned to downslope cells and
        attributes of each parition"""

        def rn_proportions_attributes(nudat_r):
            n = nudat_r[0]
            qso = nudat_r[2]
            qsi = nudat_r[3]
            E = nudat_r[4]
            A = nudat_r[5]

            att_up = nudat_r[7]
            att_in = nudat_r[8]

            # get donor node ids
            dn = self.arndn[self.arn == n]

            # get receiver node idds
            rn = self._grid.at_node.dataset["flow__receiver_node"].values[n]
            rn_ = rn.copy()

            # remove node n and delivery nodes from the receiver node list
            rn_ = rn_[np.where(rn_ != -1)]
            rn_ = rn_[~np.isin(rn_, dn)]

            # these constraints needed to avoid infinite loop
            if qso > 0 and n not in self._grid.boundary_nodes:
                if len(rn_) < 1:  # if no receiver nodes, flux stays at node n
                    rp = np.array([1])
                    rn = np.array([n])
                    rqso = rp * qso
                else:  # flux is proportioned to remaining receiver nodes
                    rp = self._grid.at_node.dataset[
                        "flow__receiver_proportions"
                    ].values[n]
                    rp = rp[np.isin(rn, rn_)]
                    rp = rp / rp.sum()
                    rn = rn_
                    rqso = rp * qso

                # create donor node array, all values are donor node,
                # length equal to number of receiver nodes
                rndn = (np.ones(len(rn)) * n).astype(int)

                # outgoing attributes
                if self._tracked_attributes:
                    att_out = self._attribute_out(att_up, att_in, qsi, E, A)

                    for key in self._tracked_attributes:
                        ratt = np.ones(len(rqso)) * att_out[key]
                        self.aratt_ns[key] = np.concatenate(
                            (self.aratt_ns[key], ratt), axis=0
                        )  # next step receiving node incoming particle diameter list
                        self.arattL[key].append(ratt)

                # store receiving nodes and fluxes in temporary arrays
                self.arndn_ns = np.concatenate(
                    (self.arndn_ns, rndn), axis=0
                )  # next iteration donor nodes
                self.arn_ns = np.concatenate(
                    (self.arn_ns, rn), axis=0
                )  # next iteration receiving nodes
                self.arqso_ns = np.concatenate(
                    (self.arqso_ns, rqso), axis=0
                )  # next iteration qsi
                self.arnL.append(rn)
                self.arqsoL.append(rqso)
                self.arndnL.append(rndn)

        [rn_proportions_attributes(r) for r in self.nudat]

    def _determine_qsi(self):
        """determine flux of incoming material (qsi) to a node.
        returns self.qsi_dat: np array of receiving nodes [column 0],
        and qsi to those nodes [column 1]
        """

        def _qsi(n):
            """sum the incoming flux to node n"""
            qsi = np.sum(self.arqso[self.arn == n])
            return qsi

        ll = np.array([_qsi(n) for n in self.arn_u], dtype=object)
        ll = np.reshape(ll, (-1, 1))
        arn_ur = np.reshape(self.arn_u, (-1, 1))
        self.qsi_dat = np.concatenate((arn_ur, ll), axis=1)

    def _update_E_dem(self):
        """update energy__elevation"""
        n = self.qsi_dat[:, 0].astype(int)
        qsi = self.qsi_dat[:, 1]
        # energy elevation is equal to the topographic elevation plus qsi
        self._grid.at_node["energy__elevation"] = self._grid.at_node[
            "topographic__elevation"
        ].copy()
        self._grid.at_node["energy__elevation"][n] = (
            self._grid.at_node["energy__elevation"].copy()[n] + qsi
        )

    def _update_energy_slope(self):
        """updates the topographic__slope and flow directions grid fields using the
        energy__elevation field. This function is presently not used but may be useful
        for future implementations of MWR"""
        fd = FlowDirectorMFD(
            self._grid,
            surface="energy__elevation",
            diagonals=True,
            partition_method=self.routing_partition_method,
        )
        fd.run_one_step()

    def _update_dem(self):
        """updates the topographic elevation of the landscape dem and soil
        thickness fields"""
        n = self.nudat[:, 0].astype(int)
        deta = self.nudat[:, 1]
        self._grid.at_node["soil__thickness"][n] = (
            self._grid.at_node["soil__thickness"][n] + deta
        )
        self._grid.at_node["topographic__elevation"][n] = (
            self._grid.at_node["topographic__elevation"][n] + deta
        )

    def _update_topographic_slope(self):
        """updates the topographic__slope and flow directions fields using the
        topographic__elevation field"""
        fd = FlowDirectorMFD(
            self._grid,
            surface="topographic__elevation",
            diagonals=True,
            partition_method=self.routing_partition_method,
        )
        fd.run_one_step()

    def _update_attribute_at_node(self, key):
        """for each unique node in receiving node list, update the attribute
        using attribute value determined in the _E_A_qso_determine_attributes method

        Parameters
        ----------
        key: string
            one of the tracked attributes
        """
        n = self.nudat[:, 0].astype(int)
        new_node_pd = np.array([d[key] for d in self.nudat[:, 6]])
        if np.isnan(np.sum(new_node_pd)):
            raise ValueError(f"{key} is {new_node_pd}")
        self._grid.at_node[key][n] = new_node_pd

    def _settle(self):
        """for each unique node in receiving node list, after erosion, aggradation
        and change in node elevation have been determined, check that the height of the node
        is not greater than permitted by angle of repose/critical slope as evaluated from
        the lowest cell. Note, slope is not updated in this function. It is updated
        simultaneously at a later stage during the iteration.
        """
        # for each node in the list, use the slope field, computed from the previous
        # iteration, to compute settlment and settlment direction to adjacent cells
        for ii, n in enumerate(self.arn_u):
            if self.D_L[ii] > 0:  # only settle if node has had deposition...use dif?
                rn = self._grid.at_node.dataset["flow__receiver_node"].values[n]
                # slope to all receiving cells
                slpn = self._grid.at_node["topographic__steepest_slope"][n]

                # only consider downslope cells
                slpn = slpn[np.where(rn != -1)]
                rn = rn[np.where(rn != -1)]

                # critical slope
                if self.variable_slpc:
                    self.slpc = (
                        self.a * self._grid.at_node["drainage_area"][n] ** self.b
                    )

                # only consider all cells that slope > Sc
                rn = rn[slpn > self.slpc]
                slpn = slpn[slpn > self.slpc]

                # if slope to downlsope nodes > Sc, adjust elevation of node n
                if len(rn) >= 1:
                    # destribute material to downslope nodes based on weighted
                    # average slope (same as multiflow direciton proportions,
                    # but here only determined for downslope nodes in which
                    # S  > Sc )
                    sslp = sum(slpn)
                    pp = slpn / sslp

                    # determine the total flux sent to S > Sc downslope cells
                    # mean/min downslope cell elevation
                    zo = self._grid.at_node["topographic__elevation"][rn].min()

                    # node n elevation
                    zi = self._grid.at_node["topographic__elevation"][n]

                    # height of node n using Sc*dx above min downslope elevation
                    slp_h = self.slpc * self._grid.dx

                    # out going sediment depth, determined as half the depth
                    # above the critical slope
                    qso_s = (zi - (zo + slp_h)) / 2

                    if qso_s < 0:  # no negative (this shouldn't be needed because only
                        # nodes greater than slpc considered)
                        qso_s = 0

                    if qso_s > self.D_L[ii]:  # settlement out can not exceed deposit
                        qso_s = self.D_L[ii]

                    qso_s_i = qso_s * pp  # proportion sent to each receiving cell

                    # update the topographic elevation
                    self._grid.at_node["topographic__elevation"][n] = (
                        self._grid.at_node["topographic__elevation"][n] - qso_s
                    )
                    self._grid.at_node["topographic__elevation"][rn] = (
                        self._grid.at_node["topographic__elevation"][rn] + qso_s_i
                    )

                    # update the soil thickness
                    self._grid.at_node["soil__thickness"][n] = (
                        self._grid.at_node["soil__thickness"][n] - qso_s
                    )
                    self._grid.at_node["soil__thickness"][rn] = (
                        self._grid.at_node["soil__thickness"][rn] + qso_s_i
                    )

                    # update tracked attributes for sediment movement during settlement
                    if self._tracked_attributes:
                        # create the att_in dict, this is the same for each of the rn
                        att_in = {}
                        for key in self._tracked_attributes:
                            att_in[key] = self._grid.at_node[key][n]
                        for v, n_ in enumerate(rn):
                            A = qso_s_i[v]
                            n_att_d_ = self._attributes_node(n_, att_in, 0, A)
                            for key in self._tracked_attributes:
                                self._grid.at_node[key][n_] = n_att_d_[key]

    def _erosion(self, n, depth, slpn, att_in=None):
        """if self.grain_shear is True, determines the erosion depth using
        equation (13), otherwise uses equation (12).

        Parameters
        ----------
        n : int
            node id
        depth : float
            erosion depth
        slpn : float
            slope in [l/L]
        att_in: dict
            dictionary of the value of each attribute, this function
            only uses particle__diameter

        Returns
        -------
        E : float
            erosion depth [L]
        att_up : dict
            dictionary of the value of each attribute at node n
        Tau : float
            basal shear stress [Pa]
        u : float
            flow velocity [m/s]
        """
        theta = np.arctan(slpn)  # convert tan(theta) to theta
        # attributes of eroded material
        if self._tracked_attributes:
            att_up = {}
            for key in self._tracked_attributes:
                att_up[key] = self._grid.at_node[key][n]
        else:
            att_up = None

        if self.grain_shear:
            # shear stress approximated as a power function of inertial shear stress
            Dp = att_in["particle__diameter"]
            if depth < Dp:  # grain size dependent erosion breaks if depth<Dp
                Dp = depth * 0.99
            u = flow_velocity(Dp, depth, slpn, self.g)

            Tau = shear_stress_grains(self.vs, self.ros, Dp, depth, slpn, self.g)

            Ec = self._grid.dx * erosion_rate(self.k, Tau, self.f, self._grid.dx)

        else:
            # quasi-static approximation
            Tau = self.ro_mw * self.g * depth * (np.sin(theta))
            Ec = self.k * (Tau) ** self.f
            u = np.nan

        dmx = self._grid.at_node["soil__thickness"][n]

        E = min(dmx, Ec)

        return (E, att_up, Tau, u)

    def _aggradation(self, qsi, slpn, n):
        """determine aggradation depth as a function of a threshold-slope,
        L*qsi or a function of both. Where the L metric is computed following
        Campforts, et al., 2020 but is expressed as 1-(slpn/slpc)**2 rather
        than dx/(1-(slpn/slpc)**2)

        Parameters
        ----------
        qsi : float
            incoming flux [l3/iteration/l2]
        slpn : float
            slope at node in [L/L]
        n : int
            node id

        Returns
        -------
        A : float
            aggradation depth [L]
        """

        if self.deposition_rule == "L_metric":
            A = self._deposit_L_metric(qsi, slpn)
        elif self.deposition_rule == "critical_slope":
            A = self._deposit_friction_angle(qsi, n)
        elif self.deposition_rule == "both":
            A_L = self._deposit_L_metric(qsi, slpn)
            A_f = self._deposit_friction_angle(qsi, n)
            A = min(A_L, A_f)

        return A

    def _determine_zo(self, n, zi, qsi):
        """determine the minimum elevation of the adjacent nodes. If all adjacent
        nodes are higher than the elevation of the node + qsi, zo is set to zi

        Parameters
        ----------
        n : int
            node id
        zi : float
            topographic elevation at node n (eta_n)
        qsi : float
            incoming flux [l3/iteration/l2]

        Returns
        -------
        zo : float
            topographic elevation of the lowest elevation node [l],
            adjacent to node n
        """

        # get adjacent nodes
        adj_n = np.hstack(
            (
                self._grid.adjacent_nodes_at_node[n],
                self._grid.diagonal_adjacent_nodes_at_node[n],
            )
        )

        # exclude closed boundary nodes
        adj_n = adj_n[~np.isin(adj_n, self._grid.closed_boundary_nodes)]

        # elevation of flow surface at node... may not need this
        ei = qsi + zi

        # nodes below elevation of node n
        rn_e = adj_n[self._grid.at_node["topographic__elevation"][adj_n] < ei]

        if len(rn_e) > 0:
            zo = self._grid.at_node["topographic__elevation"][rn_e].min()

        else:  # an obstruction in the DEM
            zo = zi
        return zo

    def _deposit_L_metric(self, qsi, slpn):
        """
        determine the L metric similar to Campforts et al. (2020)

        Parameters
        ----------
        qsi : float
            in coming flux per unit contour width
        slpn : float
            slope of node, measured in downslope direction (downslope is postive)

        Returns
        -------
        A_L : float
            aggradation depth [L]
        """
        Lnum = np.max([(1 - (slpn / self.slpc) ** 2), 0])
        A_L = qsi * Lnum

        return A_L

    def _deposit_friction_angle(self, qsi, n):
        """determine deposition depth following equations 4 though 9

        Parameters
        ----------
        qsi : float
            incoming flux [l3/iteration/l2]
        n : int
            node id

        Returns
        -------
        A_f : float
            aggradation depth [L]
        """
        slp_h = self.slpc * self._grid.dx
        zi = self._grid.at_node["topographic__elevation"][n]
        zo = self._determine_zo(n, zi, qsi)
        rule = (zi - zo) <= (slp_h)

        def eq(qsi, zo, zi, slp_h):
            dx = self._grid.dx
            sc = self.slpc
            s = (zi - zo) / dx
            sd = sc - s
            D1 = sc * dx / 2
            a = 0.5 * dx * sd
            b = D1 - 0.5 * dx * sd
            c = -qsi
            N1 = -b + (((b**2) - 4 * a * c) ** 0.5) / (2 * a)
            N2 = -b - (((b**2) - 4 * a * c) ** 0.5) / (2 * a)
            ndn = np.round(max([N1, N2, 1]))  # aggradation on at least one node
            A = min((1 / ndn) * qsi + ((ndn - 1) / 2) * dx * sd, qsi)
            return A

        if rule:
            A_f = eq(qsi, zo, zi, slp_h)
        else:
            A_f = 0

        if A_f < 0:
            warnings.warn(
                f"negative aggradation!! n={n}, qsi={qsi}, A_f={A_f}, zo={zo}, zi={zi}",
                stacklevel=2,
            )
            # raise(ValueError)
        return A_f

    def _attributes_in(self, n, qsi):
        """determine the weighted average attribute value of the incoming
        flow

        Parameters
        ----------
        n : int
            node id
        qsi : float
            incoming flux [l3/iteration/l2]

        Returns
        -------
        att_in: dict
            dictionary of each attribute value flowing into the node
        """
        if qsi == 0:
            att_in = dict.fromkeys(self._tracked_attributes, 0)
        elif (np.isnan(qsi)) or (np.isinf(qsi)):
            msg = "in-flowing flux is nan or inf"
            raise ValueError(msg)
        else:
            att_in = {}
            for key in self._tracked_attributes:
                att_in[key] = np.sum(
                    (self.aratt[key][self.arn == n]) * (self.arqso[self.arn == n]) / qsi
                )
        return att_in

    def _attributes_node(self, n, att_in, E, A):
        """determine the weighted average attributes of the newly aggraded material
        + the inplace regolith

        Parameters
        ----------
        n : int
            node id
        att_in: dict
            dictionary of the value of each attribute flowing into the node
        E : float
            erosion depth [L]
        A : float
            aggradation depth [L]

        Returns
        -------
        n_att_d: dict
            dictionary of each attribute value at the node after erosion
            and aggradation
        """

        def weighted_avg_at_node(key):
            """determine weighted average attribute at the node. If all soil
            is eroded, attribute value is zero"""
            if A + self._grid.at_node["soil__thickness"][n] - E > 0:
                inatt = self._grid.at_node[key][n]
                n_att = (
                    inatt * (self._grid.at_node["soil__thickness"][n] - E)
                    + att_in[key] * A
                ) / (A + self._grid.at_node["soil__thickness"][n] - E)
            else:
                n_att = 0
            if (n_att < 0) or (np.isnan(n_att)) or (np.isinf(n_att)):
                msg = "node particle diameter is negative, nan or inf"
                raise ValueError(msg)
            return n_att

        n_att_d = {}
        for key in self._tracked_attributes:
            n_att_d[key] = weighted_avg_at_node(key)

        return n_att_d

    def _attribute_out(self, att_up, att_in, qsi, E, A):
        """determine the weighted average attributes of the outgoing
        flux

        Parameters
        ----------
        att_up: dict
            dictionary of each attribute value at the node before erosion
            or aggradation
        att_in: dict
            dictionary of each attribute value flowing into the node
        qsi : float
            incoming flux [l3/iteration/l2]
        E : float
            erosion depth [L]
        A : float
            aggradation depth [L]

        Returns
        -------
        att_out: dict
            dictionary of each attribute value flowing out of the node
        """
        att_out = {}
        for key in self._tracked_attributes:
            att_out[key] = np.sum(
                (att_up[key] * E + att_in[key] * (qsi - A)) / (qsi - A + E)
            )
            check_val = att_out[key]
            if (check_val <= 0) or (np.isnan(check_val)) or (np.isinf(check_val)):
                msg = f"out-flowing {key} is zero, negative, nan or inf"
                raise ValueError(msg)
        return att_out


def flow_velocity(Dp, h, s, g):
    us = (g * h * s) ** 0.5
    u = us * 5.75 * np.log10(h / Dp)
    return u


def shear_stress_grains(vs, ros, Dp, h, s, g):
    theta = np.arctan(s)
    phi = np.arctan(0.32)
    u = flow_velocity(Dp, h, s, g)
    dudz = u / h
    Tcn = np.cos(theta) * vs * ros * (Dp**2) * (dudz**2)
    tau = Tcn * np.tan(phi)
    return tau


def shear_stress_static(vs, ros, rof, h, s, g):
    theta = np.arctan(s)
    rodf = vs * ros + (1 - vs) * rof
    tau = rodf * g * h * (np.sin(theta))
    return tau


def erosion_coef_k(E_l, tau, f, dx):
    k = E_l * dx / (tau**f)
    return k


def erosion_rate(k, tau, f, dx):
    E_l = (k * tau**f) / dx
    return E_l
