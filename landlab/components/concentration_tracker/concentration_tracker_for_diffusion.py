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

        ∂CH/∂t = [-(∂q_x C_x)/∂x - (∂q_y C_y)/∂y] + C_br*H_brw + PH + DH

    where :math:`H` is sediment depth, :math:`q_x` and :math:`q_y` are sediment
    fluxed in the x and y directions, :math:`C_br` is concentration in parent
    bedrock, :math:`H_brw` is the height of bedrock weathered into soil,
    :math:`P` is the local production rate, :math:`D` is the local decay rate.

    NOTE: This component requires a soil flux field calculated by a hillslope
    diffusion component and must be run after every diffusion step. Currently,
    this component WILL ONLY WORK IF COUPLED with the DepthDependentDiffuser or
    the DepthDependentTaylorDiffuser (without the dynamic timestep option).

    Examples
    --------
    A 1-D hillslope:

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import DepthDependentDiffuser
    >>> from landlab.components import ConcentrationTrackerForDiffusion
    >>> mg = RasterModelGrid((3, 5),xy_spacing=2.)
    >>> mg.set_status_at_node_on_edges(right=4, top=4, left=4, bottom=4)
    >>> mg.status_at_node[5] = mg.BC_NODE_IS_FIXED_VALUE
    >>> c = mg.add_zeros('sediment_property__concentration', at='node')
    >>> h = mg.add_zeros("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> _ = mg.add_zeros('soil_production__rate', at='node')
    >>> c[8] += 1
    >>> h += mg.node_x
    >>> z_br += mg.node_x
    >>> z += z_br + h
    >>> ddd = DepthDependentDiffuser(mg)
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> ddd.run_one_step(1.)
    >>> ct.run_one_step(1.)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([4.11701964, 8.01583689, 11.00247875]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0., 0.24839685, 1.]))
    True

    Now, a 2-D pyramid-shaped hillslope.

    >>> mg = RasterModelGrid((5, 5),xy_spacing=2.)
    >>> c = mg.add_zeros('sediment_property__concentration', at='node')
    >>> h = mg.add_zeros("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> _ = mg.add_zeros('soil_production__rate', at='node')
    >>> c[12] += 1
    >>> h += 2
    >>> z_br += 8
    >>> z_br -= abs(4 - mg.node_x)
    >>> z_br -= abs(4 - mg.node_y)
    >>> z += z_br + h
    >>> ddd = DepthDependentDiffuser(mg)
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> ddd.run_one_step(1.)
    >>> ct.run_one_step(1.)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([6.        ,  7.13533528,  6.        ,
    ...                       7.13533528,  8.27067057,  7.13533528,
    ...                       6.        ,  7.13533528,  6.         ]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0.        ,  0.38079708,  0.        ,
    ...                       0.38079708,  1.        ,  0.38079708,
    ...                       0.        ,  0.38079708,  0.         ]))
    True

    And running one more step.

    >>> ddd.run_one_step(1.)
    >>> ct.run_one_step(1.)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([5.52060315,  6.62473963,  5.52060315,
    ...                       6.62473963,  8.00144598,  6.62473963,
    ...                       5.52060315,  6.62473963,  5.52060315 ]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0.09648071,  0.44750673,  0.09648071,
    ...                       0.44750673,  1.        ,  0.44750673,
    ...                       0.09648071,  0.44750673,  0.09648071 ]))
    True

    Finally, the same 2D hillslope now using the DepthDependentTaylorDiffuser.
    Note that the timestep must be smaller than 1 to maintain stability in the
    diffusion calculation. Typically, one could use the dynamic timestepping
    option. However, here it will provide incorrect soil flux values to the
    ConcentrationTrackerForDiffusion, which cannot do sub-timestep calculations.
    Use the if_unstable="warn" flag when instantiating the Taylor diffuser and
    pick a timestep that is stable.

    >>> from landlab.components import DepthDependentTaylorDiffuser
    >>> mg = RasterModelGrid((5, 5),xy_spacing=2.)
    >>> c = mg.add_zeros('sediment_property__concentration', at='node')
    >>> h = mg.add_zeros("soil__depth", at="node")
    >>> z_br = mg.add_zeros("bedrock__elevation", at="node")
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> _ = mg.add_zeros('soil_production__rate', at='node')
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
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes],
    ...             np.array([6.        ,  7.30826823,  6.        ,
    ...                       7.30826823,  8.61653645,  7.30826823,
    ...                       6.        ,  7.30826823,  6.         ]))
    True
    >>> np.allclose(mg.at_node["sediment_property__concentration"][mg.core_nodes],
    ...             np.array([0.        ,  0.26436925,  0.        ,
    ...                       0.26436925,  1.        ,  0.26436925,
    ...                       0.        ,  0.26436925,  0.         ]))
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
        "sediment_property__mass_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-^2/y",
            "mapping": "link",
            "doc": "Mass of property fluxing along links",
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
        if not self._grid.at_node["sediment_property__concentration"].any():
            self._grid.at_node["sediment_property__concentration"] += self.C_init
        self._concentration = self._grid.at_node["sediment_property__concentration"]

        if not self._grid.at_node["bedrock_property__concentration"].any():
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
        self._QC_links = self._grid.at_link["sediment_property__mass_flux"]

        # Check that concentration values are within physical limits
        if isinstance(concentration_initial, np.ndarray):
            if concentration_initial.any() < 0:
                raise ValueError("Concentration cannot be negative.")
        else:
            if concentration_initial < 0:
                raise ValueError("Concentration cannot be negative.")

        if isinstance(concentration_in_bedrock, np.ndarray):
            if concentration_in_bedrock.any() < 0:
                raise ValueError("Concentration in bedrock cannot be negative.")
        else:
            if concentration_in_bedrock < 0:
                raise ValueError("Concentration in bedrock cannot be negative.")

    @property
    def C_init(self):
        """Initial concentration in soil/sediment (kg/m^3)."""
        return self._C_init

    @property
    def C_br(self):
        """Concentration in bedrock (kg/m^3)."""
        return self._C_br

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
        self._C_init = return_array_at_node(self._grid, new_val)

    @C_br.setter
    def C_br(self, new_val):
        self._C_br = return_array_at_node(self._grid, new_val)

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
        self._grid.at_link["QC"] = (
            self._grid.at_link["soil__flux"][:] * self._C_links[:]
        )

        # Calculate flux concentration divergence
        dQCdx = self._grid.calc_flux_div_at_node(self._grid.at_link["QC"])

        # Calculate other components of mass balance equation
        with np.errstate(divide="ignore", invalid="ignore"):
            C_local = C_old * (self._soil__depth_old / self._soil__depth)
            C_from_weathering = (
                self._C_br * (self._soil_prod_rate * dt) / self._soil__depth
            )
            Production = (dt * self._P / 2) * (
                self._soil__depth_old / self._soil__depth + 1
            )
            Decay = (dt * self._D / 2) * (self._soil__depth_old / self._soil__depth + 1)

            # Calculate concentration
            self._concentration[:] = (
                C_local
                + C_from_weathering
                + (dt / self._soil__depth) * (-dQCdx)
                + Production
                - Decay
            )

        # Replace nan values (from dividing by zero soil depth)
        np.nan_to_num(C_local, copy=False)
        np.nan_to_num(C_from_weathering, copy=False)
        np.nan_to_num(Production, copy=False)
        np.nan_to_num(Decay, copy=False)
        np.nan_to_num(self._concentration, copy=False)

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
