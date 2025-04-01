"""
Created on Fri Jun 16 15:06:50 2023

@author: LaurentRoberge
"""

import numpy as np

from landlab import Component
from landlab import NodeStatus
from landlab.utils.return_array import return_array_at_node


class ConcentrationTrackerForSpace(Component):
    """This component tracks the concentration of any user-defined property of
    sediment using a mass balance approach in which concentration :math:`C_s`
    is calculated as:

    .. math::

                        ∂C_sH/∂t = C_s_w*D_s_w + C_s*E_s

    where :math:`H` is sediment depth, :math:`C_s_w` is concentration in
    sediment suspended in the water column, :math:`D_s_w` is volumetric
    depositional flux of sediment from the water column per unit bed area, and
    :math:`E_s` is volumetric erosional flux of sediment from the bed per unit
    bed area.

    .. note::

        This component requires the "sediment__outflux", "bedrock__erosion_flux"
        "sediment__erosion_flux", and "sediment__deposition_flux" fields
        calculated by the :class:`~.SpaceLargeScaleEroder` component. This
        component does not use the typical run_one_step(dt) method. Instead,
        a start_tracking() method is implemented immediately before every
        :class:`~.SpaceLargeScaleEroder` step and a stop_tracking(dt) method
        immediately after every :class:`~.SpaceLargeScaleEroder` step.
        See the docstring examples below.

        The required inputs "phi", "fraction_fines", and "settling_velocity"
        must have the same value as those used for the instance of
        :class:`~.SpaceLargeScaleEroder`.

    Examples
    --------

    A 1-D stream channel:

    >>> import numpy as np
    >>> from landlab import NodeStatus, RasterModelGrid
    >>> from landlab.components import PriorityFloodFlowRouter
    >>> from landlab.components import SpaceLargeScaleEroder
    >>> from landlab.components import ConcentrationTrackerForSpace

    >>> mg = RasterModelGrid((3, 5), xy_spacing=10.0)

    >>> mg.set_status_at_node_on_edges(
    ...     right=NodeStatus.CLOSED,
    ...     top=NodeStatus.CLOSED,
    ...     left=NodeStatus.CLOSED,
    ...     bottom=NodeStatus.CLOSED,
    ... )
    >>> mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE
    >>> mg.status_at_node.reshape(mg.shape)
    array([[4, 4, 4, 4, 4],
           [1, 0, 0, 0, 4],
           [4, 4, 4, 4, 4]], dtype=uint8)

    >>> mg.at_node["sediment_property__concentration"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 1.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.at_node["soil__depth"] = [
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ... ]
    >>> mg.at_node["bedrock__elevation"] = mg.node_x / 100
    >>> mg.at_node["topographic__elevation"] = (
    ...     mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ... )

    >>> fr = PriorityFloodFlowRouter(mg)
    >>> fr.run_one_step()
    >>> sp = SpaceLargeScaleEroder(mg, phi=0, F_f=0, v_s=1)
    >>> ct = ConcentrationTrackerForSpace(
    ...     mg,
    ...     phi=0,
    ...     fraction_fines=0,
    ...     settling_velocity=1,
    ... )

    >>> for i in range(40):
    ...     fr.run_one_step()
    ...     ct.start_tracking()
    ...     sp.run_one_step(10.0)
    ...     ct.stop_tracking(10.0)
    ...

    Erosion has lowered the topography and reduced channel bed sediment depth.
    >>> np.allclose(
    ...     mg.at_node["topographic__elevation"][mg.core_nodes],
    ...     np.array([1.00292211, 1.00902572, 1.0258774]),
    ... )
    True
    >>> np.allclose(
    ...     mg.at_node["soil__depth"][mg.core_nodes],
    ...     np.array([0.90294696, 0.80909071, 0.72601329]),
    ... )
    True

    Some high-concentration sediment has been transported from upstream to be
    deposited on the channel bed further downstream.
    >>> np.allclose(
    ...     mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...     np.array([0.0496547, 0.0997232, 0.9999151]),
    ... )
    True


    Now, a 2-D landscape with stream channels. All boundaries are closed except
    for Node 0, which is the outlet of the catchment.

    >>> mg = RasterModelGrid((6, 6), xy_spacing=10.0)

    >>> mg.set_status_at_node_on_edges(
    ...     right=NodeStatus.CLOSED,
    ...     top=NodeStatus.CLOSED,
    ...     left=NodeStatus.CLOSED,
    ...     bottom=NodeStatus.CLOSED,
    ... )
    >>> mg.status_at_node[0] = mg.BC_NODE_IS_FIXED_VALUE
    >>> mg.status_at_node.reshape(mg.shape)
    array([[4, 4, 4, 4, 4, 4],
           [4, 0, 0, 0, 0, 4],
           [4, 0, 0, 0, 0, 4],
           [4, 0, 0, 0, 0, 4],
           [4, 0, 0, 0, 0, 4],
           [1, 4, 4, 4, 4, 4]], dtype=uint8)


    >>> mg.at_node["sediment_property__concentration"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.at_node["soil__depth"] = [
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    ... ]

    # Add noise to the bedrock to create some topographic structure.
    >>> np.random.seed(5)
    >>> mg.add_zeros("bedrock__elevation", at="node")
    >>> mg.at_node["bedrock__elevation"] += np.random.rand(mg.number_of_nodes) / 100
    >>> mg.at_node["bedrock__elevation"][0] = 0

    >>> mg.at_node["topographic__elevation"] = (
    ...     mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ... )

    # Instantiate components.
    >>> fr = PriorityFloodFlowRouter(mg)
    >>> fr.run_one_step()
    >>> sp = SpaceLargeScaleEroder(mg, phi=0, F_f=0, v_s=1)
    >>> ct = ConcentrationTrackerForSpace(
    ...     mg,
    ...     phi=0,
    ...     fraction_fines=0,
    ...     settling_velocity=1,
    ... )

    # Run SPACE for 1,000 years to generate a fluvial network.
    >>> for i in range(1000):
    ...     mg.at_node["bedrock__elevation"][mg.core_nodes] += 0.001
    ...     mg.at_node["topographic__elevation"][:] = (
    ...         mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ...     )
    ...     fr.run_one_step()
    ...     sp.run_one_step(1.0)
    ...

    # Set high concentration at a headwater node to trace sediment downstream.
    >>> mg.at_node["sediment_property__concentration"][22] += 1

    >>> for i in range(100):
    ...     mg.at_node["bedrock__elevation"][mg.core_nodes] += 0.001
    ...     mg.at_node["topographic__elevation"][:] = (
    ...         mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ...     )
    ...     fr.run_one_step()
    ...     ct.start_tracking()
    ...     sp.run_one_step(1.0)
    ...     ct.stop_tracking(1.0)
    ...

    Some high-concentration sediment has been transported from the headwaters
    to be deposited on the channel bed further downstream. We can trace this
    sediment and see where the channel lies within the landscape.
    >>> np.allclose(
    ...     mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...     np.array(
    ...         [
    ...             [0.0288311, 0.0447778, 0.0, 0.0],
    ...             [0.0, 0.0, 0.0598574, 0.0],
    ...             [0.0, 0.0, 0.0, 0.9548471],
    ...             [0.0, 0.0, 0.0, 0.0],
    ...         ]
    ...     ).flatten(),
    ... )
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    CITATION

    """

    _name = "ConcentrationTrackerForSpace"

    _unit_agnostic = True

    _cite_as = """
    CITATION
    """

    _info = {
        "soil__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "sediment__outflux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m^3/yr",
            "mapping": "node",
            "doc": "Sediment flux (volume per unit time of sediment leaving each node)",
        },
        "bedrock__erosion_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": (
                "Bedrock erosion flux from bedrock to water column (depth eroded per"
                " unit time)"
            ),
        },
        "sediment__erosion_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": (
                "Sediment erosion flux from bed to water column (depth eroded per"
                " unit time)"
            ),
        },
        "sediment__deposition_flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": (
                "Sediment deposition flux from water column to bed (depth deposited"
                " per unit time)"
            ),
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "sediment_property__concentration": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3",
            "mapping": "node",
            "doc": "Mass concentration of property per unit volume of sediment",
        },
        "bedrock_property__concentration": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3",
            "mapping": "node",
            "doc": "Mass concentration of property per unit volume of bedrock",
        },
    }

    def __init__(
        self,
        grid,
        phi: float | None = None,
        fraction_fines: float | None = None,
        settling_velocity: float | None = None,
        concentration_initial=0,
        concentration_in_bedrock=0,
    ):
        """
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        phi: float
            Sediment porosity, [-]
        fraction_fines: float
            Fraction of permanently suspendable fines in bedrock, [-].
        settling_velocity: float
            Net effective settling velocity for chosen grain size metric, [L/T]
        concentration_initial: positive float, array, or field name (optional)
            Initial concentration in soil/sediment, [-/L^3]
        concentration_in_bedrock: positive float, array, or field name (optional)
            Concentration in bedrock, [-/L^3]
        """

        if phi is None:
            raise ValueError(
                "`phi` is a required input parameter. "
                "It must have the same value as the `phi` input "
                "parameter used for the `SpaceLargeScaleEroder` "
                "instance to which this component is coupled. "
                "See the docstring in each component for details."
            )
        if fraction_fines is None:
            raise ValueError(
                "`fraction_fines` is a required input parameter. "
                "It must have the same value as the same input "
                "parameter used for the `SpaceLargeScaleEroder` "
                "instance to which this component is coupled. "
                "The parameter is named `F_f` in SpaceLargeScaleEroder. "
                "See the docstring in each component for details."
            )
        if settling_velocity is None:
            raise ValueError(
                "`settling_velocity` is a required input parameter. "
                "It must have the same value as the same input "
                "parameter used for the `SpaceLargeScaleEroder` "
                "instance to which this component is coupled. "
                "It is named `v_s` in SpaceLargeScaleEroder. "
                "See the docstring in each component for details."
            )

        super().__init__(grid)
        # Store grid and parameters

        # use setters for C_init, C_br defined below
        self.C_init = concentration_initial
        self.C_br = concentration_in_bedrock

        # get reference to inputs
        self._phi = phi
        self._fraction_fines = fraction_fines
        self._settling_velocity = settling_velocity
        self._soil__depth = self._grid.at_node["soil__depth"]
        self._soil__depth_old = self._soil__depth.copy()
        self._Qs_out = self._grid.at_node["sediment__outflux"]

        # Define variables used for internal calculations
        self._cell_area = self._grid.cell_area_at_node
        self._C_sw = np.zeros(self._grid.number_of_nodes)
        self._QsCsw_in = np.zeros(self._grid.number_of_nodes)
        self._QsCsw_out = np.zeros(self._grid.number_of_nodes)
        self._BED_ero_depo_term = np.zeros(self._grid.number_of_nodes)

        # create outputs if necessary and get reference.
        self.initialize_output_fields()

        # Define concentration field (if all zeros, then add C_init)
        if np.allclose(self._grid.at_node["sediment_property__concentration"], 0.0):
            self._grid.at_node["sediment_property__concentration"] += self.C_init
        self._concentration = self._grid.at_node["sediment_property__concentration"]

        if np.allclose(self._grid.at_node["bedrock_property__concentration"], 0.0):
            self._grid.at_node["bedrock_property__concentration"] += self.C_br
        self.C_br = self._grid.at_node["bedrock_property__concentration"]

        if phi >= 1.0:
            raise ValueError("Porosity must be < 1.0")

        if fraction_fines > 1.0:
            raise ValueError("Fraction of fines must be <= 1.0")

        if phi < 0.0:
            raise ValueError("Porosity must be > 0.0")

        if fraction_fines < 0.0:
            raise ValueError("Fraction of fines must be > 0.0")

    @property
    def C_init(self):
        """Initial concentration in soil/sediment (kg/m^3)."""
        return self._C_init

    @property
    def C_br(self):
        """Concentration in bedrock (kg/m^3)."""
        return self._C_br

    @C_init.setter
    def C_init(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration in sediment cannot be negative")
        self._C_init = return_array_at_node(self._grid, new_val)

    @C_br.setter
    def C_br(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration in bedrock cannot be negative")
        self._C_br = return_array_at_node(self._grid, new_val)

    def _copy_old_soil_depth(self):
        """Store a copy of soil depth. This is used as the old soil depth when
        calculating changes in concentration.
        """

        self._soil__depth_old = self._soil__depth.copy()

    def _calc_concentration_watercolumn_and_bed(self, dt):
        """Calculate change in concentration within sediment transported in
        the water column and within sediment on the bed for a time period 'dt'.

        Parameters
        ----------

        dt: float (time)
            The imposed timestep.
        """
        # Define values generated by SPACE/SpaceLargeScaleEroder
        flow_receivers = self._grid.at_node["flow__receiver_node"]
        q = self._grid.at_node["surface_water__discharge"]
        Er = self._grid.at_node["bedrock__erosion_flux"]
        Es = self._grid.at_node["sediment__erosion_flux"]
        D_sw = self._grid.at_node["sediment__deposition_flux"]

        # Calculate portions of equation that have soil depth as denominator
        is_soil = self._soil__depth > 0.0

        old_depth_over_new = np.divide(
            self._soil__depth_old, self._soil__depth, where=is_soil
        )
        old_depth_over_new[~is_soil] = 0.0

        dt_over_depth = np.divide(dt, self._soil__depth, where=is_soil)
        dt_over_depth[~is_soil] = 0.0

        # Calculate mass balance terms that don't need downstream iteration
        WC_Es_term = Es * self._cell_area
        WC_Er_term = (1 - self._fraction_fines) * Er * self._cell_area
        WC_denominator_term = np.ones(np.shape(q))
        WC_denominator_term[q != 0] = (
            1 + self._settling_velocity * self._cell_area[q != 0] / q[q != 0]
        )
        BED_C_local_term = self._concentration * old_depth_over_new

        # Get stack of node ids from top to bottom of channel network
        node_status = self._grid.status_at_node
        stack_flip_ud = np.flipud(self._grid.at_node["flow__upstream_node_order"])
        # Select core nodes where qs >0
        stack_flip_ud_sel = stack_flip_ud[
            (node_status[stack_flip_ud] == NodeStatus.CORE) & (q[stack_flip_ud] > 0.0)
        ]

        # zero out array values that were updated in the old stack
        self._C_sw[:] = 0
        self._QsCsw_in[:] = 0
        self._BED_ero_depo_term[:] = 0

        # Iterate concentration calc (first BED, then WC) at each node
        for node_id in stack_flip_ud_sel:
            # Calculate QsCsw_out (i.e., QsCs in the water column)
            self._QsCsw_out[node_id] = (
                self._QsCsw_in[node_id]
                + self._concentration[node_id] * WC_Es_term[node_id]
                + self.C_br[node_id] * WC_Er_term[node_id]
            ) / WC_denominator_term[node_id]

            # Send QsCsw_out values to flow receiver nodes
            self._QsCsw_in[flow_receivers[node_id]] += self._QsCsw_out[node_id]

            # Divide QsCsw_out (from above) by Qs_out to get C_sw
            if self._Qs_out[node_id] > 0:
                self._C_sw[node_id] = self._QsCsw_out[node_id] / self._Qs_out[node_id]
            else:
                self._C_sw[node_id] = 0.0

            # Calculate BED erosion/deposition term (requires C_sw from above)
            self._BED_ero_depo_term[node_id] = (
                self._C_sw[node_id] * D_sw[node_id]
                - self._concentration[node_id] * Es[node_id]
            ) / (1 - self._phi)

            # Calculate BED concentration
            self._concentration[node_id] = (
                BED_C_local_term[node_id]
                + dt_over_depth[node_id] * self._BED_ero_depo_term[node_id]
            )

            self._concentration[~is_soil] = 0.0

    def start_tracking(self):
        """Stores values necessary for calculating changes in concentration.
        This method must be called prior to running the sediment flux component
        that changes physical properties of bedrock, soil, and/or topography.
        """

        self._copy_old_soil_depth()

    def stop_tracking(self, dt):
        """Calculate changes in concentration that have occurred in the timestep
        since tracking was started. This method must be called after running the
        sediment flux component that changes physical properties of bedrock,
        soil, and/or topography.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """

        self._calc_concentration_watercolumn_and_bed(dt)

    def run_one_step(self):
        """run_one_step is not implemented for this component."""
        raise NotImplementedError("run_one_step()")
