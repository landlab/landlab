"""Landlab component that simulates detachment-limited river erosion.

This component calculates changes in elevation in response to
vertical incision.
"""

import numpy as np
import scipy.constants

from landlab import Component


class DepthSlopeProductErosion(Component):
    """Calculate erosion rate as a function of the depth-slope product.

    Erosion rate is calculated as, ``erosion_rate = k_e * ((tau ** a - tau_crit ** a))``

    *k_e*
        Erodibility coefficient
    *tau*
        Bed shear stress: ``tau = rho * g * h * S``
    *rho*
        Density of fluid
    *g*
        Gravitational acceleration
    *h*
        Water depths
    *S*
        Slope
    *tau_crit*
        Critical shear stress
    *a*
        Positive exponent


    Note this equation was presented in Tucker, G.T., 2004, Drainage basin
    sensitivity to tectonic and climatic forcing: Implications of a stochastic
    model for the role of entrainment and erosion thresholds,
    Earth Surface Processes and Landforms.

    More generalized than other erosion components, as it doesn't require the
    upstream node order, links to flow receiver and flow receiver fields. Instead,
    takes in the water depth and slope fields on NODES calculated by the
    OverlandFlow class and erodes the landscape in response to the hydrograph
    generted by that method.

    As of right now, this component relies on the OverlandFlow component
    for stability. There are no stability criteria implemented in this class.
    To ensure model stability, use StreamPowerEroder or FastscapeEroder
    components instead.

    .. codeauthor:: Jordan Adams

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import DepthSlopeProductErosion

    Create a grid on which to calculate detachment ltd sediment transport.

    >>> grid = RasterModelGrid((5, 5))

    The grid will need some data to provide the detachment limited sediment
    transport component. To check the names of the fields that provide input to
    the detachment ltd transport component, use the *input_var_names* class
    property.

    Create fields of data for each of these input variables.

    First create topography. This is a flat surface of elevation 10 m.

    >>> grid.at_node["topographic__elevation"] = np.ones(grid.number_of_nodes)
    >>> grid.at_node["topographic__elevation"] *= 10.0
    >>> grid.at_node["topographic__elevation"] = [
    ...     [10.0, 10.0, 10.0, 10.0, 10.0],
    ...     [10.0, 10.0, 10.0, 10.0, 10.0],
    ...     [10.0, 10.0, 10.0, 10.0, 10.0],
    ...     [10.0, 10.0, 10.0, 10.0, 10.0],
    ...     [10.0, 10.0, 10.0, 10.0, 10.0],
    ... ]

    Now we'll add an arbitrary water depth field on top of that topography.

    >>> grid.at_node["surface_water__depth"] = [
    ...     [5.0, 5.0, 5.0, 5.0, 5.0],
    ...     [4.0, 4.0, 4.0, 4.0, 4.0],
    ...     [3.0, 3.0, 3.0, 3.0, 3.0],
    ...     [2.0, 2.0, 2.0, 2.0, 2.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ... ]

    Using the set topography, now we will calculate slopes on all nodes.

    First calculating slopes on links

    >>> grid.at_link["water_surface__slope"] = grid.calc_grad_at_link(
    ...     "surface_water__depth"
    ... )

    Now putting slopes on nodes

    >>> grid.at_node["water_surface__slope"] = (
    ...     grid.at_link["water_surface__slope"][grid.links_at_node]
    ...     * grid.active_link_dirs_at_node
    ... ).max(axis=1)
    >>> grid.at_node["water_surface__slope"][grid.core_nodes]
    array([1., 1., 1., 1., 1., 1., 1., 1., 1.])


    Instantiate the `DepthSlopeProductErosion` component to work on this grid, and
    run it. In this simple case, we need to pass it a time step ('dt') and also
    an erodibility factor ('k_e').

    >>> dt = 1.0
    >>> dspe = DepthSlopeProductErosion(
    ...     grid, k_e=0.00005, g=9.81, slope="water_surface__slope"
    ... )
    >>> dspe.run_one_step(
    ...     dt=dt,
    ... )

    Now we test to see how the topography changed as a function of the erosion
    rate. First, we'll look at the erosion rate:

    >>> dspe.dz.reshape(grid.shape)
    array([[ 0.    , -2.4525, -2.4525, -2.4525,  0.    ],
           [ 0.    , -1.962 , -1.962 , -1.962 ,  0.    ],
           [ 0.    , -1.4715, -1.4715, -1.4715,  0.    ],
           [ 0.    , -0.981 , -0.981 , -0.981 ,  0.    ],
           [ 0.    ,  0.    ,  0.    ,  0.    ,  0.    ]])

    Now, our updated topography...

    >>> grid.at_node["topographic__elevation"].reshape(grid.shape)
    array([[10.    ,  7.5475,  7.5475,  7.5475, 10.    ],
           [10.    ,  8.038 ,  8.038 ,  8.038 , 10.    ],
           [10.    ,  8.5285,  8.5285,  8.5285, 10.    ],
           [10.    ,  9.019 ,  9.019 ,  9.019 , 10.    ],
           [10.    , 10.    , 10.    , 10.    , 10.    ]])
    """

    _name = "DepthSlopeProductErosion"

    _unit_agnostic = True

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__slope": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "gradient of the ground surface",
        },
    }

    def __init__(
        self,
        grid,
        k_e=0.001,
        fluid_density=1000.0,
        g=scipy.constants.g,
        a_exp=1.0,
        tau_crit=0.0,
        uplift_rate=0.0,
        slope="topographic__slope",
    ):
        """Calculate detachment limited erosion rate on nodes using the shear
        stress equation, solved using the depth slope product.

        Landlab component that generalizes the detachment limited erosion
        equation, primarily to be coupled to the the Landlab OverlandFlow
        component.

        This component adjusts topographic elevation and is contained in the
        landlab.components.detachment_ltd_erosion folder.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab grid.
        k_e : float
            Erodibility parameter, (m^(1+a_exp)*s^(2*a_exp-1)/kg^a_exp)
        fluid_density : float, optional
            Density of fluid, default set to water density of 1000 kg / m^3
        g : float, optional
            Acceleration due to gravity (m/s^2).
        a_exp : float, optional
            exponent on shear stress, positive, unitless
        tau_crit : float, optional
            threshold for sediment movement, (kg/m/s^2)
        uplift_rate : float, optional
            uplift rate applied to the topographic surface, m/s
        slope : str
            Field name of an at-node field that contains the slope.
        """
        super().__init__(grid)

        assert slope in grid.at_node

        self._slope = slope
        self._a = a_exp
        self._g = g
        self._rho = fluid_density
        self._E = self._grid.zeros(at="node")
        self._uplift_rate = uplift_rate
        self._tau_crit = tau_crit
        self._k_e = k_e

        self._dz = self._grid.zeros(at="node")

    def run_one_step(self, dt):
        """Erode into grid topography.

        For one time step, this erodes into the grid topography using
        the water discharge and topographic slope.

        The grid field 'topographic__elevation' is altered each time step.

        Parameters
        ----------
        dt : float
            Time step.
        """
        S = self._grid.at_node[self._slope]
        h = self._grid.at_node["surface_water__depth"]

        self._tau = self._rho * self._g * h * S

        (greater_than_tc,) = np.where(self._tau >= self._tau_crit)
        (less_than_tc,) = np.where(self._tau < self._tau_crit)

        self._E[less_than_tc] = 0.0

        self._E[greater_than_tc] = self._k_e * (
            (self._tau[greater_than_tc] ** self._a) - (self._tau_crit**self._a)
        )

        self._E[self._E < 0.0] = 0.0

        self._dz = (self._uplift_rate - self._E) * dt

        self._grid["node"]["topographic__elevation"] += self._dz

    @property
    def dz(self):
        """Magnitude of change of the topographic__elevation due to erosion
        [L]."""
        return self._dz
