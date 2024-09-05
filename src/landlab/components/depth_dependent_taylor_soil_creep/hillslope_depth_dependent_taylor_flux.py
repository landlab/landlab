"""DepthDependentTaylorNonLinearDiffuser Component.

@author: R Glade
@author: K Barnhart
@author: G Tucker
"""

import numpy as np

from landlab import Component
from landlab import LinkStatus
from landlab.core.messages import deprecation_message


class DepthDependentTaylorDiffuser(Component):
    r"""
    This component implements a depth-dependent Taylor series diffusion rule,
    combining concepts of Ganti et al. (2012) and Johnstone and Hilley (2014).

    Hillslope sediment flux uses a Taylor Series expansion of the Andrews-
    Bucknam formulation of nonlinear hillslope flux derived following following
    Ganti et al., 2012 with a depth dependent component inspired Johnstone and
    Hilley (2014). The flux :math:`q_s` is given as:

    .. math::

        q_s = - K H_* \nabla \eta (
                1 + (S/S_c)^2 + (S/S_c)^4 + .. + (S/S_c)^2(n-1)
            ) (1 - exp( - H / H_*)

    where :math:`K` is a transport velocity coefficient, :math:`\eta` is land
    surface elevation, :math:`S` is the slope gradient (defined as
    positive downward), :math:`S_c` is the critical slope, :math:`n` is the
    number of terms, :math:`H` is the soil depth on links, and :math:`H_*` is
    the soil transport decay depth.

    The default behavior uses two terms to produce a slope dependence as
    described by Equation 6 of Ganti et al. (2012).

    This component will ignore soil thickness located at non-core nodes.

    Examples
    --------
    First lets make a simple example with flat topography.

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeatherer
    >>> from landlab.components import DepthDependentTaylorDiffuser
    >>> mg = RasterModelGrid((5, 5))
    >>> soilTh = mg.add_zeros("node", "soil__depth")
    >>> z = mg.add_zeros("node", "topographic__elevation")
    >>> BRz = mg.add_zeros("node", "bedrock__elevation")
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> np.allclose(mg.at_node["soil_production__rate"][mg.core_nodes], 1.0)
    True
    >>> DDdiff.run_one_step(2.0)
    >>> np.allclose(mg.at_node["topographic__elevation"][mg.core_nodes], 0.0)
    True
    >>> np.allclose(mg.at_node["bedrock__elevation"][mg.core_nodes], -2.0)
    True
    >>> np.allclose(mg.at_node["soil__depth"][mg.core_nodes], 2.0)
    True

    Now a more complicated example with a slope.

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros("node", "soil__depth")
    >>> z = mg.add_zeros("node", "topographic__elevation")
    >>> BRz = mg.add_zeros("node", "bedrock__elevation")
    >>> z += mg.node_x.copy()
    >>> BRz += mg.node_x / 2.0
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> np.allclose(
    ...     mg.at_node["soil_production__rate"][mg.core_nodes],
    ...     np.array([0.60653066, 0.36787944, 0.22313016]),
    ... )
    True
    >>> DDdiff.run_one_step(0.1)
    >>> np.allclose(
    ...     mg.at_node["topographic__elevation"][mg.core_nodes],
    ...     np.array([1.04773024, 2.02894986, 3.01755898]),
    ... )
    True
    >>> np.allclose(
    ...     mg.at_node["bedrock__elevation"][mg.core_nodes],
    ...     np.array([0.43934693, 0.96321206, 1.47768698]),
    ... )
    True
    >>> np.allclose(mg.at_node["soil__depth"], z - BRz)
    True

    The DepthDependentTaylorDiffuser makes and moves soil at a rate proportional
    to slope, this means that there is a characteristic time scale for soil
    transport and an associated stability criteria for the timestep. The
    maximum characteristic time scale, :math:`De_{max}`, is given as a function of the
    hillslope diffustivity, :math:`D`, the maximum slope, :math:`S_{max}`,
    and the critical slope :math:`S_c`.

    .. math::

        De_{max} = D
            \left(
            1 +
            \left( \frac{S_{max}{S_c}\right )^2 +
            \left( \frac{S_{max}{S_c}\right )^4 +
            \dots +
            \left( \frac{S_{max}{S_c}\right )^{( 2 * ( n - 1 ))}
            \right)

    The maximum stable time step is given by

    .. math::

        dtmax = courant_factor * dx * dx / Demax

    Where the courant factor is a user defined scale (default is 0.2), and
    dx is the length of the shortest link in the grid.

    The DepthDependentTaylorDiffuser has a boolean flag that permits a user
    to be warned if timesteps are too large for the slopes in the model grid
    (if_unstable = 'warn') and a boolean flag that turns on dynamic timestepping
    (dynamic_dt = False).

    >>> DDdiff = DepthDependentTaylorDiffuser(mg, if_unstable="warn")
    >>> DDdiff.run_one_step(2.0)
    Topographic slopes are high enough such that the Courant condition is
    exceeded AND you have not selected dynamic timestepping with
    dynamic_dt=True. This may lead to infinite and/or nan values for slope,
    elevation, and soil depth. Consider using a smaller time step or dynamic
    timestepping. The Courant condition recommends a timestep of
    0.09534076073069653 or smaller.

    Alternatively you can specify if_unstable='raise', and a Runtime Error will
    be raised if this condition is not met.

    Next, lets do an example with dynamic timestepping.

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros("node", "soil__depth")
    >>> z = mg.add_zeros("node", "topographic__elevation")
    >>> BRz = mg.add_zeros("node", "bedrock__elevation")

    We'll use a steep slope and very little soil.

    >>> z += mg.node_x.copy() ** 2
    >>> BRz = z.copy() - 1.0
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)

    Lets try to move the soil with a large timestep. Without dynamic time
    steps, this gives a warning that we've exceeded the dynamic timestep size
    and should use a smaller timestep. We could either use the smaller timestep,
    or specify that we want to use a dynamic timestep.

    >>> DDdiff = DepthDependentTaylorDiffuser(mg, if_unstable="warn", dynamic_dt=False)
    >>> expweath.calc_soil_prod_rate()
    >>> DDdiff.run_one_step(10)
    Topographic slopes are high enough such that the Courant condition is
    exceeded AND you have not selected dynamic timestepping with
    dynamic_dt=True. This may lead to infinite and/or nan values for slope,
    elevation, and soil depth. Consider using a smaller time step or dynamic
    timestepping. The Courant condition recommends a timestep of
    0.004 or smaller.

    Now, we'll re-build the grid and do the same example with dynamic timesteps.

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros("node", "soil__depth")
    >>> z = mg.add_zeros("node", "topographic__elevation")
    >>> BRz = mg.add_zeros("node", "bedrock__elevation")
    >>> z += mg.node_x.copy() ** 2
    >>> BRz = z.copy() - 1.0
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg, if_unstable="warn", dynamic_dt=True)
    >>> expweath.calc_soil_prod_rate()
    >>> DDdiff.run_one_step(10)
    >>> np.any(np.isnan(z))
    False

    Now, we'll test that changing the transport decay depth behaves as expected.

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros("node", "soil__depth")
    >>> z = mg.add_zeros("node", "topographic__elevation")
    >>> BRz = mg.add_zeros("node", "bedrock__elevation")
    >>> z += mg.node_x.copy() ** 0.5
    >>> BRz = z.copy() - 1.0
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg, soil_transport_decay_depth=0.1)
    >>> DDdiff.run_one_step(1)
    >>> soil_decay_depth_point1 = mg.at_node["topographic__elevation"][mg.core_nodes]
    >>> z[:] = 0
    >>> z += mg.node_x.copy() ** 0.5
    >>> BRz = z.copy() - 1.0
    >>> soilTh[:] = z - BRz
    >>> DDdiff = DepthDependentTaylorDiffuser(mg, soil_transport_decay_depth=1.0)
    >>> DDdiff.run_one_step(1)
    >>> soil_decay_depth_1 = mg.at_node["topographic__elevation"][mg.core_nodes]
    >>> np.greater(soil_decay_depth_1[1], soil_decay_depth_point1[1])
    False

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Barnhart, K., Glade, R., Shobe, C., Tucker, G. (2019). Terrainbento 1.0: a
    Python package for multi-model analysis in long-term drainage basin
    evolution. Geoscientific Model Development  12(4), 1267--1297.
    https://dx.doi.org/10.5194/gmd-12-1267-2019

    **Additional References**

    Ganti, V., Passalacqua, P., Foufoula-Georgiou, E. (2012). A sub-grid scale
    closure for nonlinear hillslope sediment transport models Journal of
    Geophysical Research: Earth Surface  117(F2).
    https://dx.doi.org/10.1029/2011jf002181

    Johnstone, S., Hilley, G. (2015). Lithologic control on the form of
    soil-mantled hillslopes Geology  43(1), 83-86.
    https://doi.org/10.1130/G36052.1

    """

    _name = "DepthDependentTaylorDiffuser"

    _unit_agnostic = True

    _cite_as = """
    @article{barnhart2019terrain,
      author = {Barnhart, Katherine R and Glade, Rachel C and Shobe, Charles M
                and Tucker, Gregory E},
      title = {{Terrainbento 1.0: a Python package for multi-model analysis in
                long-term drainage basin evolution}},
      doi = {10.5194/gmd-12-1267-2019},
      pages = {1267---1297},
      number = {4},
      volume = {12},
      journal = {Geoscientific Model Development},
      year = {2019},
    }
    """

    _info = {
        "bedrock__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of the bedrock surface",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "soil__flux": {
            "dtype": float,
            "intent": "out",
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
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__slope": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "gradient of the ground surface",
        },
    }

    def __init__(
        self,
        grid,
        linear_diffusivity=None,
        slope_crit=1.0,
        soil_transport_decay_depth=1.0,
        nterms=2,
        dynamic_dt=False,
        if_unstable="pass",
        courant_factor=0.2,
        soil_transport_velocity=1.0,
    ):
        """Initialize the DepthDependentTaylorDiffuser.

        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        linear_diffusivity: float, optional, DEPRECATED
            Hillslope diffusivity / decay depth, m/yr
            Default = 1.0
        slope_crit: float, optional
            Critical gradient parameter, m/m
            Default = 1.0
        soil_transport_decay_depth: float, optional
            characteristic transport soil depth, m
            Default = 1.0
        nterms: int, optional. default = 2
            number of terms in the Taylor expansion.
            Two terms (default) gives the behavior
            described in Ganti et al. (2012).
        dynamic_dt : bool, optional, default  = False
            Whether internal timestepping is used.
        if_unstable : str, optional, default = "pass"
            What to do if unstable (options are "pass",
            "raise", "warn")
        courant_factor : float, optional, default = 0.2
            Courant factor for timestep calculation.
        soil_transport_velocity : float, optional, default = 1.0
            Velocity parameter for soil transport, m/yr. Diffusivity is the
            product of this parameter and soil_transport_decay_depth.
        """
        super().__init__(grid)

        # Handle now-deprecated diffusivity argument
        if linear_diffusivity is None:
            self._K = soil_transport_velocity
        else:
            message = """Use of linear_diffusivity is deprecated, because the
                         name is misleading: it is actually a velocity;
                         diffusivity is obtained by multiplying by soil
                         transport decay depth. Use soil_transport_velocity
                         instead."""
            print(deprecation_message(message))
            self._K = linear_diffusivity
        self._soil_transport_decay_depth = soil_transport_decay_depth
        self._slope_crit = slope_crit
        self._nterms = nterms

        self._dynamic_dt = dynamic_dt
        self._if_unstable = if_unstable
        self._courant_factor = courant_factor
        self._shortest_link = np.amin(grid.length_of_link)  # for Courant

        # get reference to inputs
        self._elev = self._grid.at_node["topographic__elevation"]
        self._soil_prod_rate = self._grid.at_node["soil_production__rate"]
        self._depth = self._grid.at_node["soil__depth"]

        # create outputs if necessary and get reference.
        self.initialize_output_fields()
        self._slope = self._grid.at_link["topographic__slope"]
        self._flux = self._grid.at_link["soil__flux"]
        self._bedrock = self._grid.at_node["bedrock__elevation"]

    def soilflux(self, dt):
        """Calculate soil flux for a time period 'dt'.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        # establish time left as all of dt
        time_left = dt

        # begin while loop for time left
        while time_left > 0.0:
            # calculate soil__depth
            self._grid.at_node["soil__depth"][:] = (
                self._grid.at_node["topographic__elevation"]
                - self._grid.at_node["bedrock__elevation"]
            )

            # Calculate soil depth at links.
            self._H_link = self._grid.map_value_at_max_node_to_link(
                "topographic__elevation", "soil__depth"
            )

            # Calculate gradients
            self._slope = self._grid.calc_grad_at_link(self._elev)
            self._slope[self._grid.status_at_link == LinkStatus.INACTIVE] = 0.0

            # Test for time stepping courant condition
            # Test for time stepping courant condition
            courant_slope_term = 0.0
            courant_s_over_scrit = self._slope.max() / self._slope_crit
            for i in range(0, 2 * self._nterms, 2):
                courant_slope_term += courant_s_over_scrit**i
                if np.any(np.isinf(courant_slope_term)):
                    message = (
                        "Soil flux term is infinite in Courant condition "
                        "calculation. This is likely due to "
                        "using too many terms in the Taylor expansion."
                    )
                    raise RuntimeError(message)
            # Calculate De Max
            De_max = self._K * (courant_slope_term)
            # Calculate longest stable timestep
            self._dt_max = self._courant_factor * (self._shortest_link**2) / De_max

            # Test for the Courant condition and print warning if user intended
            # for it to be printed.
            if (
                (self._dt_max < dt)
                and (not self._dynamic_dt)
                and (self._if_unstable != "pass")
            ):
                message = (
                    "Topographic slopes are high enough such that the "
                    "Courant condition is exceeded AND you have not "
                    "selected dynamic timestepping with dynamic_dt=True. "
                    "This may lead to infinite and/or nan values for "
                    "slope, elevation, and soil depth. Consider using a "
                    "smaller time step or dynamic timestepping. The "
                    "Courant condition recommends a timestep of "
                    "" + str(self._dt_max) + " or smaller."
                )
                if self._if_unstable == "raise":
                    raise RuntimeError(message)
                if self._if_unstable == "warn":
                    print(message)

            # if dynamic dt is selected, use it, otherwise, use the entire time
            if self._dynamic_dt:
                self._sub_dt = np.min([dt, self._dt_max])
                time_left -= self._sub_dt
            else:
                self._sub_dt = dt
                time_left = 0

            # update sed flux, topography, soil, and bedrock based on the
            # current self._sub_dt
            self._update_flux_topography_soil_and_bedrock()

    def _update_flux_topography_soil_and_bedrock(self):
        """Calculate soil flux and update topography."""
        # Calculate flux
        slope_term = 0.0
        s_over_scrit = self._slope / self._slope_crit
        for i in range(0, 2 * self._nterms, 2):
            slope_term += s_over_scrit**i
            if np.any(np.isinf(slope_term)):
                message = (
                    "Soil flux term is infinite. This is likely due to "
                    "using too many terms in the Taylor expansion."
                )
                raise RuntimeError(message)

        self._flux[:] = -(
            (self._K * self._slope * self._soil_transport_decay_depth)
            * (slope_term)
            * (1.0 - np.exp(-self._H_link / self._soil_transport_decay_depth))
        )

        # Calculate flux divergence
        dqdx = self._grid.calc_flux_div_at_node(self._flux)

        # Calculate change in soil depth
        dhdt = self._soil_prod_rate - dqdx

        # Calculate soil depth at nodes
        self._depth[self._grid.core_nodes] += dhdt[self._grid.core_nodes] * self._sub_dt

        # prevent negative soil thickness
        self._depth[self._depth < 0.0] = 0.0

        # Calculate bedrock elevation
        self._bedrock[self._grid.core_nodes] -= (
            self._soil_prod_rate[self._grid.core_nodes] * self._sub_dt
        )

        # Update topography
        self._elev[self._grid.core_nodes] = (
            self._depth[self._grid.core_nodes] + self._bedrock[self._grid.core_nodes]
        )

    def run_one_step(self, dt):
        """

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        self.soilflux(dt)
