class MassWastingSaver:
    """This class is instantiated and called by MassWastingRunout. It saves
    MWR model ouput. It is only called in MassWastingRunout if save = True"""

    def __init__(self, MassWastingRunout):
        self.MWR = MassWastingRunout

    def prep_data_containers(self):
        # lists and dictionaries for tracking model behavior
        self.EL = []  # entrainment depth / regolith depth
        self.AL = []  # aggradation (deposition) depth
        self.qsiL = []  # incoming flux (qsi)
        self.TauL = []  # basal shear stress
        self.slopeL = []  # slope
        self.velocityL = []  # velocity (if computed)
        self.arqso_r = {}  # flux out
        self.arn_r = {}  # receiver nodes
        self.arndn_r = {}  # donar nodes
        self.aratt_r = {}  # arriving attributes
        self.flowing_volume = (
            {}
        )  # the total volume [m3] of the mobilized runout material
        # dictionaries, that save the entire model grid field each model iteration
        # for all fields listed below
        # usefull for creating movies of how the flow and terrain evolve
        self.runout_evo_maps = {}  # runout material + topographic__elevation
        self.topo_evo_maps = {}  # topographic__elevation
        self.att_r = {}  # attribute value
        self.st_r = {}  # soil__thickness
        self.tss_r = {}  # topographic__steepest_slope
        self.frn_r = {}  # flow__receiver_node
        self.frp_r = {}  # 'flow__receiver_proportions'

    def prep_mw_data_containers(self, mw_i, mw_id):
        # containeers for each unique mass wasting ID
        self.runout_evo_maps[mw_i] = {}
        self.topo_evo_maps[mw_i] = {}
        self.flowing_volume[mw_id] = []
        self.st_r[mw_id] = []
        self.tss_r[mw_id] = []
        self.frn_r[mw_id] = []
        self.frp_r[mw_id] = []
        self.arqso_r[mw_id] = []
        self.arn_r[mw_id] = []
        self.arndn_r[mw_id] = []
        if self.MWR.track_attributes:
            self.att_r[mw_id] = dict.fromkeys(
                self.MWR._tracked_attributes, []
            )  # this becomes the data container for each attribute
            self.aratt_r[mw_id] = dict.fromkeys(self.MWR._tracked_attributes, [])

    def save_conditions_before_runout(self, mw_i, mw_id):
        # save first set of data to reflect scar/depression in DEM created by
        # mass wasting source area
        self.runout_evo_maps[mw_i][0] = self.MWR._grid.at_node[
            "energy__elevation"
        ].copy()
        self.topo_evo_maps[mw_i][0] = self.MWR._grid.at_node[
            "topographic__elevation"
        ].copy()
        self.flowing_volume[mw_id].append(0)
        if self.MWR.track_attributes:
            for key in self.MWR._tracked_attributes:
                self.att_r[mw_id][key].append(
                    self.MWR._grid.at_node[key].copy()
                )  # for each attribute, a copy of entire grid
                self.aratt_r[mw_id][key].append(self.MWR.aratt)
        self.st_r[mw_id].append(self.MWR._grid.at_node["soil__thickness"].copy())
        self.tss_r[mw_id].append(
            self.MWR._grid.at_node["topographic__steepest_slope"].copy()
        )
        self.frn_r[mw_id].append(self.MWR._grid.at_node["flow__receiver_node"].copy())
        self.frp_r[mw_id].append(
            self.MWR._grid.at_node["flow__receiver_proportions"].copy()
        )
        self.arqso_r[mw_id].append(self.MWR.arqso)
        self.arn_r[mw_id].append(self.MWR.arn)
        self.arndn_r[mw_id].append(self.MWR.arndn)

    def save_conditions_after_one_iteration(self, mw_i, mw_id):
        DEMf = self.MWR._grid.at_node["topographic__elevation"].copy()
        DEMdf_r = DEMf - self.MWR._grid.at_node["topographic__initial_elevation"]
        self.flowing_volume[mw_id].append(
            DEMdf_r.sum() * self.MWR._grid.dx * self.MWR._grid.dy
        )
        self.runout_evo_maps[mw_i][self.MWR.c + 1] = self.MWR._grid.at_node[
            "energy__elevation"
        ].copy()
        self.runout_evo_maps[mw_i][self.MWR.c + 1] = self.MWR._grid.at_node[
            "energy__elevation"
        ].copy()
        self.topo_evo_maps[mw_i][self.MWR.c + 1] = self.MWR._grid.at_node[
            "topographic__elevation"
        ].copy()
        if self.MWR.track_attributes:
            for key in self.MWR._tracked_attributes:
                self.att_r[mw_id][key].append(self.MWR._grid.at_node[key].copy())
                self.aratt_r[mw_id][key].append(self.MWR.arattL)
        self.st_r[mw_id].append(self.MWR._grid.at_node["soil__thickness"].copy())
        self.tss_r[mw_id].append(
            self.MWR._grid.at_node["topographic__steepest_slope"].copy()
        )
        self.frn_r[mw_id].append(self.MWR._grid.at_node["flow__receiver_node"].copy())
        self.frp_r[mw_id].append(
            self.MWR._grid.at_node["flow__receiver_proportions"].copy()
        )
        self.arqso_r[mw_id].append(self.MWR.arqsoL)
        self.arn_r[mw_id].append(self.MWR.arnL)
        self.arndn_r[mw_id].append(self.MWR.arndn)

    def save_flow_stats(self, E, A, qsi, slpn, Tau, u):
        self.EL.append(E)
        self.AL.append(A)
        self.qsiL.append(qsi)
        self.slopeL.append(slpn)  # slope
        self.TauL.append(Tau)
        self.velocityL.append(u)  # velocity (if any)
