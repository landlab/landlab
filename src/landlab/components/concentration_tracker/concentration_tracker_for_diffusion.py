"""
Created on Wed May 31 11:41:20 2023

@author: LaurentRoberge
"""

import numpy as np

from landlab import Component
from landlab import LinkStatus
from landlab.grid.mappers import map_value_at_max_node_to_link
from landlab.utils.return_array import return_array_at_node


class ConcentrationTrackerForDiffusion(Component):
    """Track the concentration of any user-defined property.

    This component tracks the concentration of any user-defined property of
    sediment using a mass balance approach in which the concentration :math:`C`
    is calculated as:

    .. math::

        ∂CH / ∂t = [-(∂q_x C_x) / ∂x - (∂q_y C_y) / ∂y] + C_br * H_brw

    where :math:`H` is sediment depth, :math:`q_x` and :math:`q_y` are sediment
    fluxed in the x and y directions, :math:`C_br` is concentration in parent
    bedrock, and :math:`H_brw` is the height of bedrock weathered into soil.

    .. note::

        This component requires a soil flux field calculated by a hillslope diffusion
        component. This component is implemented by applying a :meth:`start_tracking`
        method immediately before every diffusion step and a :meth:`stop_tracking`
        method immediately after every diffusion step. These methods are applied instead
        of the typical :meth:`run_one_step` method. See the docstring examples below.

        Currently, this component will only work if coupled with the
        :class:`~.DepthDependentDiffuser` or the :class:`~.DepthDependentTaylorDiffuser`
        (without the dynamic timestep option).

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
    >>> ct.start_tracking()
    >>> ddd.run_one_step(1.0)
    >>> ct.stop_tracking(1.0)

    >>> mg.at_node["topographic__elevation"].reshape(mg.shape)
    array([[ 0. ,  4.        ,  8.        , 12.        , 16. ],
           [ 0. ,  4.11701964,  8.01583689, 11.00247875, 16. ],
           [ 0. ,  4.        ,  8.        , 12.        , 16. ]])
    >>> mg.at_node["sediment_property__concentration"].reshape(mg.shape)
    array([[0. , 0. , 0.        , 0. , 0. ],
           [0. , 0. , 0.24839685, 1. , 0. ],
           [0. , 0. , 0.        , 0. , 0. ]])

    Now, a 2-D pyramid-shaped hillslope.

    >>> mg = RasterModelGrid((5, 5), xy_spacing=2.0)

    >>> mg.at_node["sediment_property__concentration"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 1.0, 0.0, 0.0],
    ... ]
    >>> h = mg.add_full("soil__depth", 2.0, at="node")
    >>> z_br = mg.add_field(
    ...     "bedrock__elevation",
    ...     8.0 - abs(4.0 - mg.node_x) - abs(4.0 - mg.node_y),
    ...     at="node",
    ... )
    >>> z = mg.add_field("topographic__elevation", z_br + h, at="node")
    >>> _ = mg.add_zeros("soil_production__rate", at="node")

    >>> ddd = DepthDependentDiffuser(mg)
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> ct.start_tracking()
    >>> ddd.run_one_step(1.0)
    >>> ct.stop_tracking(1.0)

    >>> mg.at_node["topographic__elevation"][mg.core_nodes].reshape((3, 3))
    array([[6. ,        7.13533528, 6. ],
           [7.13533528, 8.27067057, 7.13533528],
           [6. ,        7.13533528, 6. ]])
    >>> mg.at_node["sediment_property__concentration"][mg.core_nodes].reshape((3, 3))
    array([[0.        , 0.38079708, 0. ],
           [0.38079708, 1.        , 0.38079708],
           [0.        , 0.38079708, 0. ]])

    And running one more step.

    >>> ct.start_tracking()
    >>> ddd.run_one_step(1.0)
    >>> ct.stop_tracking(1.0)

    >>> mg.at_node["topographic__elevation"][mg.core_nodes].reshape((3, 3))
    array([[5.52060315, 6.62473963, 5.52060315],
           [6.62473963, 8.00144598, 6.62473963],
           [5.52060315, 6.62473963, 5.52060315]])
    >>> mg.at_node["sediment_property__concentration"][mg.core_nodes].reshape((3, 3))
    array([[0.09648071, 0.44750673, 0.09648071],
           [0.44750673, 1.        , 0.44750673],
           [0.09648071, 0.44750673, 0.09648071]])

    Finally, the same 2D hillslope now using the
    :class:`~.DepthDependentTaylorDiffuser`. Note that the timestep must be smaller
    than 1 to maintain stability in the diffusion calculation. Typically, one could
    use the dynamic timestepping option. However, here it will provide incorrect
    soil flux values to the :class:`~.ConcentrationTrackerForDiffusion`, which
    cannot do sub-timestep calculations. Use the ``if_unstable="warn"`` flag when
    instantiating the Taylor diffuser and pick a timestep that is stable.

    >>> from landlab.components import DepthDependentTaylorDiffuser
    >>> mg = RasterModelGrid((5, 5), xy_spacing=2.0)

    >>> mg.at_node["sediment_property__concentration"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.0, 1.0, 0.0, 0.0],
    ... ]
    >>> h = mg.add_full("soil__depth", 2.0, at="node")
    >>> z_br = mg.add_field(
    ...     "bedrock__elevation",
    ...     8.0 - abs(4.0 - mg.node_x) - abs(4.0 - mg.node_y),
    ...     at="node",
    ... )
    >>> z = mg.add_field("topographic__elevation", z_br + h, at="node")
    >>> _ = mg.add_zeros("soil_production__rate", at="node")

    >>> ddtd = DepthDependentTaylorDiffuser(mg, if_unstable="warn")
    >>> ct = ConcentrationTrackerForDiffusion(mg)
    >>> ct.start_tracking()
    >>> ddtd.run_one_step(0.4)
    >>> ct.stop_tracking(0.4)

    >>> mg.at_node["topographic__elevation"][mg.core_nodes].reshape((3, 3))
    array([[6.        , 7.30826823, 6.        ],
           [7.30826823, 8.61653645, 7.30826823],
           [6.        , 7.30826823, 6.        ]])
    >>> mg.at_node["sediment_property__concentration"][mg.core_nodes].reshape((3, 3))
    array([[0.        , 0.26436925, 0.        ],
           [0.26436925, 1.        , 0.26436925],
           [0.        , 0.26436925, 0.        ]])
    """

    _name = "ConcentrationTrackerForDiffusion"

    _unit_agnostic = True

    _cite_as = ""

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
    }

    def __init__(
        self,
        grid,
        concentration_initial=0,
        concentration_in_bedrock=0,
        concentration_from_weathering=None,
    ):
        """
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        concentration_initial: float, array, or str, optional
            Initial concentration in soil/sediment as either a scalar, array,
            or the name of an existing node field, -/m^3.
        concentration_in_bedrock: float, array, or str, optional
            Concentration in bedrock as either a scalar, array, or the name
            of an existing node field, -/m^3.
        concentration_from_weathering: float or array, optional
            Concentration generated during the weathering process, -/m^3.
            Defaults to ``None``, which causes all weathered bedrock to retain its
            original parent material concentration (`concentration_in_bedrock`)
            as it weathers to soil. Use this parameter to differentiate between
            the concentration in weathered material compared to its parent bedrock.
        """

        super().__init__(grid)
        # Store grid and parameters

        # use setters for conc_init, conc_br, and conc_w defined below
        self.conc_init = concentration_initial
        self.conc_br = concentration_in_bedrock

        # get reference to inputs
        self._soil__depth = self._grid.at_node["soil__depth"]
        self._soil__depth_old = self._soil__depth.copy()
        self._soil_prod_rate = self._grid.at_node["soil_production__rate"]
        self._flux = self._grid.at_link["soil__flux"]

        # create outputs if necessary and get reference.
        self.initialize_output_fields()

        # Define concentration field (if all zeros, then add conc_init)
        if np.allclose(self._grid.at_node["sediment_property__concentration"], 0.0):
            self._grid.at_node["sediment_property__concentration"] += self.conc_init
        self._concentration = self._grid.at_node["sediment_property__concentration"]

        if np.allclose(self._grid.at_node["bedrock_property__concentration"], 0.0):
            self._grid.at_node["bedrock_property__concentration"] += self.conc_br
        self.conc_br = self._grid.at_node["bedrock_property__concentration"]

        # use setter for conc_w defined below
        self.conc_w = concentration_from_weathering

        # Sediment property concentration field (at links, to calculate dqconc_dx)
        self._conc_at_links = np.zeros(self._grid.number_of_links)

        # Sediment property mass field (at links, to calculate dqconc_dx)
        self._qconc_at_links = np.zeros(self._grid.number_of_links)

    @property
    def conc_init(self):
        """Initial concentration in soil/sediment (kg/m^3)."""
        return self._conc_init

    @property
    def conc_br(self):
        """Concentration in bedrock (kg/m^3)."""
        return self._conc_br

    @property
    def conc_w(self):
        """Concentration from the weathering process (kg/m^3)."""
        return self._conc_w

    @conc_init.setter
    def conc_init(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration cannot be negative")
        self._conc_init = return_array_at_node(self._grid, new_val)

    @conc_br.setter
    def conc_br(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Concentration in bedrock cannot be negative")
        self._conc_br = return_array_at_node(self._grid, new_val)

    @conc_w.setter
    def conc_w(self, new_val):
        if new_val is None:
            new_val = self._conc_br
        if np.any(new_val < 0.0):
            raise ValueError("Concentration cannot be negative")
        self._conc_w = new_val

    def _copy_old_soil_depth(self):
        """Store a copy of soil depth. This is used as the old soil depth when
        calculating changes in concentration.
        """

        self._soil__depth_old = self._soil__depth.copy()

    def _calc_concentration(self, dt):
        """Calculate change in concentration for a time period 'dt'.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """

        # Define concentration at previous timestep
        conc_old = self._concentration.copy()

        # Map concentration from nodes to links (following soil flux direction)
        self._conc_at_links = map_value_at_max_node_to_link(
            self._grid, "topographic__elevation", "sediment_property__concentration"
        )
        # Replace values with zero for all INACTIVE links
        self._conc_at_links[self._grid.status_at_link == LinkStatus.INACTIVE] = 0.0

        # Calculate qconc at links (sediment flux times concentration)
        self._qconc_at_links = self._grid.at_link["soil__flux"] * self._conc_at_links

        # Calculate flux concentration divergence
        dqconc_dx = self._grid.calc_flux_div_at_node(self._qconc_at_links)

        # Calculate other components of mass balance equation
        is_soil = self._soil__depth > 0.0

        old_depth_over_new = np.divide(
            self._soil__depth_old, self._soil__depth, where=is_soil
        )
        old_depth_over_new[~is_soil] = 0.0

        dt_over_depth = np.divide(dt, self._soil__depth, where=is_soil)
        dt_over_depth[~is_soil] = 0.0

        conc_local = conc_old * old_depth_over_new
        conc_from_weathering = np.divide(
            self._conc_w * self._soil_prod_rate * dt, self._soil__depth, where=is_soil
        )

        # Calculate concentration
        self._concentration[:] = (
            conc_local + conc_from_weathering + dt_over_depth * (-dqconc_dx)
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

        self._calc_concentration(dt)

    def run_one_step(self):
        """run_one_step is not implemented for this component."""
        raise NotImplementedError("run_one_step")
