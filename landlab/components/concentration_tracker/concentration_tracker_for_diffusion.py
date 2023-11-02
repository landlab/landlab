"""
Created on Wed May 31 11:41:20 2023

@author: LaurentRoberge
"""

import numpy as np

from landlab import Component, LinkStatus
from landlab.grid.mappers import map_value_at_max_node_to_link
from landlab.utils.return_array import return_array_at_node


class ConcentrationTrackerForDiffusion(Component):

    """This component tracks the concentration of any user-defined property of
    sediment using a mass balance approach in which the concentration :math:`C`
    is calculated as:

    .. math::

        ∂CH / ∂t = [-(∂q_x C_x) / ∂x - (∂q_y C_y) / ∂y] + C_br * H_brw + PH + DH

    where :math:`H` is sediment depth, :math:`q_x` and :math:`q_y` are sediment
    fluxed in the x and y directions, :math:`C_br` is concentration in parent
    bedrock, :math:`H_brw` is the height of bedrock weathered into soil,
    :math:`P` is the local production rate, :math:`D` is the local decay rate.

    .. note::

        This component requires a soil flux field calculated by a hillslope diffusion
        component and must be run after every diffusion step. Currently, this component
        WILL ONLY WORK IF COUPLED with the :class:`~.DepthDependentDiffuser` or the
        :class:`~.DepthDependentTaylorDiffuser` (without the dynamic timestep option).

    Examples
    --------

    A 1-D hillslope:

    >>> import numpy as np
    >>> from landlab import NodeStatus, RasterModelGrid
    >>> from landlab.components import DepthDependentDiffuser
    >>> from landlab.components import ConcentrationTrackerForDiffusion

    >>> mg = RasterModelGrid((3, 5), xy_spacing=2.0)

    >>> mg.set_status_at_node_on_edges(
    ...     right=NodeStatus.CLOSED,
    ...     top=NodeStatus.CLOSED,
    ...     left=NodeStatus.CLOSED,
    ...     bottom=NodeStatus.CLOSED,
    ... )
    >>> mg.status_at_node[5] = NodeStatus.FIXED_VALUE
    >>> mg.status_at_node.reshape(mg.shape)
    array([[4, 4, 4, 4, 4],
           [1, 0, 0, 0, 4],
           [4, 4, 4, 4, 4]], dtype=uint8)

    >>> mg.at_node["sediment_property__concentration"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 1.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.at_node["soil__depth"] = mg.node_x.copy()
    >>> mg.at_node["bedrock__elevation"] = mg.node_x.copy()
    >>> mg.at_node["topographic__elevation"] = (
    ...     mg.at_node["soil__depth"] + mg.at_node["bedrock__elevation"]
    ... )
    >>> _ = mg.add_zeros("soil_production__rate", at="node")

    >>> ddd = DepthDependentDiffuser(mg)
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> ddd.run_one_step(1.0)
    >>> ct.run_one_step(1.0)
    >>> np.allclose(
    ...     mg.at_node["topographic__elevation"].reshape(mg.shape),
    ...     [
    ...         [0.0, 4.0, 8.0, 12.0, 16.0],
    ...         [0.0, 4.11701964, 8.01583689, 11.00247875, 16.0],
    ...         [0.0, 4.0, 8.0, 12.0, 16.0],
    ...     ],
    ... )
    True
    >>> np.allclose(
    ...     mg.at_node["sediment_property__concentration"].reshape(mg.shape),
    ...     [
    ...         [0.0, 0.0, 0.0, 0.0, 0.0],
    ...         [0.0, 0.0, 0.24839685, 1.0, 0.0],
    ...         [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     ],
    ... )
    True

    Now, a 2-D pyramid-shaped hillslope.

    >>> mg = RasterModelGrid((5, 5), xy_spacing=2.0)
    >>> c = mg.add_zeros("sediment_property__concentration", at="node")
    >>> h = mg.add_zeros("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> _ = mg.add_zeros("soil_production__rate", at="node")
    >>> c[12] += 1
    >>> h += 2
    >>> z_br += 8
    >>> z_br -= abs(4 - mg.node_x)
    >>> z_br -= abs(4 - mg.node_y)
    >>> z += z_br + h
    >>> ddd = DepthDependentDiffuser(mg)
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> ddd.run_one_step(1.0)
    >>> ct.run_one_step(1.0)
    >>> np.allclose(
    ...     mg.at_node["topographic__elevation"][mg.core_nodes],
    ...     np.array(
    ...         [
    ...             6.0,
    ...             7.13533528,
    ...             6.0,
    ...             7.13533528,
    ...             8.27067057,
    ...             7.13533528,
    ...             6.0,
    ...             7.13533528,
    ...             6.0,
    ...         ]
    ...     ),
    ... )
    True
    >>> np.allclose(
    ...     mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...     np.array(
    ...         [
    ...             0.0,
    ...             0.38079708,
    ...             0.0,
    ...             0.38079708,
    ...             1.0,
    ...             0.38079708,
    ...             0.0,
    ...             0.38079708,
    ...             0.0,
    ...         ]
    ...     ),
    ... )
    True

    And running one more step.

    >>> ddd.run_one_step(1.0)
    >>> ct.run_one_step(1.0)
    >>> np.allclose(
    ...     mg.at_node["topographic__elevation"][mg.core_nodes],
    ...     np.array(
    ...         [
    ...             5.52060315,
    ...             6.62473963,
    ...             5.52060315,
    ...             6.62473963,
    ...             8.00144598,
    ...             6.62473963,
    ...             5.52060315,
    ...             6.62473963,
    ...             5.52060315,
    ...         ]
    ...     ),
    ... )
    True
    >>> np.allclose(
    ...     mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...     np.array(
    ...         [
    ...             0.09648071,
    ...             0.44750673,
    ...             0.09648071,
    ...             0.44750673,
    ...             1.0,
    ...             0.44750673,
    ...             0.09648071,
    ...             0.44750673,
    ...             0.09648071,
    ...         ]
    ...     ),
    ... )
    True

    Finally, the same 2D hillslope now using the DepthDependentTaylorDiffuser.
    Note that the timestep must be smaller than 1 to maintain stability in the
    diffusion calculation. Typically, one could use the dynamic timestepping
    option. However, here it will provide incorrect soil flux values to the
    ConcentrationTrackerForDiffusion, which cannot do sub-timestep calculations.
    Use the if_unstable="warn" flag when instantiating the Taylor diffuser and
    pick a timestep that is stable.

    >>> from landlab.components import DepthDependentTaylorDiffuser
    >>> mg = RasterModelGrid((5, 5), xy_spacing=2.0)
    >>> c = mg.add_zeros("sediment_property__concentration", at="node")
    >>> h = mg.add_zeros("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> _ = mg.add_zeros("soil_production__rate", at="node")
    >>> c[12] += 1
    >>> h += 2
    >>> z_br += 8
    >>> z_br -= abs(4 - mg.node_x)
    >>> z_br -= abs(4 - mg.node_y)
    >>> z += z_br + h
    >>> ddtd = DepthDependentTaylorDiffuser(mg, if_unstable="warn")
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> ddtd.run_one_step(0.4)
    >>> ct.run_one_step(0.4)
    >>> np.allclose(
    ...     mg.at_node["topographic__elevation"][mg.core_nodes],
    ...     np.array(
    ...         [
    ...             6.0,
    ...             7.30826823,
    ...             6.0,
    ...             7.30826823,
    ...             8.61653645,
    ...             7.30826823,
    ...             6.0,
    ...             7.30826823,
    ...             6.0,
    ...         ]
    ...     ),
    ... )
    True
    >>> np.allclose(
    ...     mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...     np.array(
    ...         [
    ...             0.0,
    ...             0.26436925,
    ...             0.0,
    ...             0.26436925,
    ...             1.0,
    ...             0.26436925,
    ...             0.0,
    ...             0.26436925,
    ...             0.0,
    ...         ]
    ...     ),
    ... )
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    CITATION

    """

    _name = "ConcentrationTrackerForDiffusion"

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
        "soil__flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m^2/yr",
            "mapping": "link",
            "doc": "flux of soil in direction of link",
        },
        "soil_production__rate": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": "rate of soil production at nodes",
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
        "sediment_property_production__rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3/yr",
            "mapping": "node",
            "doc": "Production rate of property per unit volume of sediment per time",
        },
        "sediment_property_decay__rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-/m^3/yr",
            "mapping": "node",
            "doc": "Decay rate of property per unit volume of sediment per time",
        },
    }

    def __init__(
        self,
        grid,
        concentration_initial=0,
        concentration_in_bedrock=0,
        concentration_from_weathering=None,
        local_production_rate=0,
        local_decay_rate=0,
    ):
        """
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        concentration_initial: positive float, array, or field name (optional)
            Initial concentration in soil/sediment, -/m^3
        concentration_in_bedrock: positive float, array, or field name (optional)
            Concentration in bedrock, -/m^3
        concentration_from_weathering: positive float or array (optional)
            Concentration generated during the weathering process, -/m^3.
            Defaults to None, which causes all weathered bedrock to retain its 
            original parent material concentration (concentration_in_bedrock)
            as it weathers to soil. 
            Use this parameter to differentiate between the concentration in 
            weathered material compared to its parent bedrock.
        local_production_rate: float, array, or field name (optional)
            Rate of local production, -/m^3/yr
        local_decay_rate: float, array, or field name (optional)
            Rate of local decay, -/m^3/yr
        """

        super().__init__(grid)
        # Store grid and parameters

        # use setters for C_init, C_br, P, and D defined below
        self.C_init = concentration_initial
        self.C_br = concentration_in_bedrock
        self.C_w = concentration_from_weathering
        self.P = local_production_rate
        self.D = local_decay_rate

        # get reference to inputs
        self._soil__depth = self._grid.at_node["soil__depth"]
        self._soil__depth_old = self._soil__depth.copy()
        self._soil_prod_rate = self._grid.at_node["soil_production__rate"]
        self._flux = self._grid.at_link["soil__flux"]

        # create outputs if necessary and get reference.
        self.initialize_output_fields()

        # Define concentration field (if all zeros, then add C_init)
        if np.allclose(self._grid.at_node["sediment_property__concentration"], 0.0):
            self._grid.at_node["sediment_property__concentration"] += self.C_init
        self._concentration = self._grid.at_node["sediment_property__concentration"]

        if np.allclose(self._grid.at_node["sediment_property__concentration"], 0.0):
            self._grid.at_node["bedrock_property__concentration"] += self.C_br
        self.C_br = self._grid.at_node["bedrock_property__concentration"]

        if not self._grid.at_node["sediment_property_production__rate"].any():
            self._grid.at_node["sediment_property_production__rate"] += self.P
        self.P = self._grid.at_node["sediment_property_production__rate"]

        if not self._grid.at_node["sediment_property_decay__rate"].any():
            self._grid.at_node["sediment_property_decay__rate"] += self.D
        self.D = self._grid.at_node["sediment_property_decay__rate"]

        # Sediment property concentration field (at links, to calculate dQCdx)
        self._C_links = np.zeros(self._grid.number_of_links)

        # Sediment property mass field (at links, to calculate dQCdx)
        self._QC_links = np.zeros(self._grid.number_of_links)

    @property
    def C_init(self):
        """Initial concentration in soil/sediment (kg/m^3)."""
        return self._C_init

    @property
    def C_br(self):
        """Concentration in bedrock (kg/m^3)."""
        return self._C_br
    
    @property
    def C_w(self):
        """Concentration from the weathering process (kg/m^3)."""
        return self._C_w

    @property
    def P(self):
        """Rate of local production (kg/m^3/yr)."""
        return self._P

    @property
    def D(self):
        """Rate of local decay (kg/m^3/yr)."""
        return self._D

    @C_init.setter
    def C_init(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration cannot be negative")
        self._C_init = return_array_at_node(self._grid, new_val)

    @C_br.setter
    def C_br(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration in bedrock cannot be negative")
        self._C_br = return_array_at_node(self._grid, new_val)
        
    @C_w.setter
    def C_w(self, new_val):
        if new_val == None:
            self._C_w = self._C_br
        if np.any(new_val < 0.0):
            raise ValueError("Concentration cannot be negative")
        self._C_w = new_val

    @P.setter
    def P(self, new_val):
        self._P = return_array_at_node(self._grid, new_val)

    @D.setter
    def D(self, new_val):
        self._D = return_array_at_node(self._grid, new_val)

    def concentration(self, dt):
        """Calculate change in concentration for a time period 'dt'.

        Parameters
        ----------

        dt: float (time)
            The imposed timestep.
        """

        # Define concentration at previous timestep
        C_old = self._concentration.copy()

        # Map concentration from nodes to links (following soil flux direction)
        # Does this overwrite fixed-value/gradient links?
        self._C_links = map_value_at_max_node_to_link(
            self._grid, "topographic__elevation", "sediment_property__concentration"
        )
        # Replace values with zero for all INACTIVE links
        self._C_links[self._grid.status_at_link == LinkStatus.INACTIVE] = 0.0

        # Calculate QC at links (sediment flux times concentration)
        self._qc_at_link = self._grid.at_link["soil__flux"] * self._C_links

        # Calculate flux concentration divergence
        dQCdx = self._grid.calc_flux_div_at_node(self._qc_at_link)

        # Calculate other components of mass balance equation
        is_soil = self._soil__depth > 0.0
        
        old_depth_over_new = np.divide(self._soil__depth_old, self._soil__depth, where=is_soil)
        old_depth_over_new[~is_soil] = 0.0
        
        dt_over_depth = np.divide(dt, self._soil__depth, where=is_soil)
        dt_over_depth[~is_soil] = 0.0
        
        C_local = C_old * old_depth_over_new
        C_from_weathering = np.divide(
            self._C_w * self._soil_prod_rate * dt, self._soil__depth, where=is_soil
        )
        Production = (dt * self._P / 2.0) * (old_depth_over_new + 1.0)
        Decay = (dt * self._D / 2.0) * (old_depth_over_new + 1.0)

        # Calculate concentration
        self._concentration[:] = (
            C_local
            + C_from_weathering
            + dt_over_depth * (-dQCdx)
            + Production
            - Decay
        )
        
        self._concentration[~is_soil] = 0.0

        # Update old soil depth to new value
        self._soil__depth_old = self._soil__depth.copy()

    def run_one_step(self, dt):
        """

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """

        self.concentration(dt)
