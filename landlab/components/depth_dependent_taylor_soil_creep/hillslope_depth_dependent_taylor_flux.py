# -*- coding: utf-8 -*-
"""
DepthDependentTaylorNonLinearDiffuser Component

@author: R Glade
@author: K Barnhart
@author: G Tucker
"""

import numpy as np

from landlab import INACTIVE_LINK, Component


class DepthDependentTaylorDiffuser(Component):

    """
    This component implements a depth-dependent Taylor series diffusion rule,
    combining concepts of Ganti et al. (2012) and Johnstone and Hilley (2014).

    Hillslope sediment flux uses a Taylor Series expansion of the Andrews-
    Bucknam formulation of nonlinear hillslope flux derived following following
    Ganti et al., 2012 with a depth dependent component inspired Johnstone and
    Hilley (2014). The flux is given as:

        qs = DS ( 1 + (S/Sc)**2 + (S/Sc)**4 + .. + (S/Sc)**2(n-1) ) * (1.0 - exp( H / Hstar )

    where D is is the diffusivity, S is the slope, Sc is the critical slope, n
    is the number of terms, H is the soil depth on links, and Hstar is the soil
    transport decay depth.

    The default behavior uses two terms to produce a slope dependence as
    described by Equation 6 of Ganti et al., (2012).

    This component will ignore soil thickness located at non-core nodes.

    Parameters
    ----------
    grid: ModelGrid
        Landlab ModelGrid object
    linear_diffusivity: float, optional.
        Hillslope diffusivity, m**2/yr
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

    Examples
    --------
    First lets make a simple example with flat topography.

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeatherer
    >>> from landlab.components import DepthDependentTaylorDiffuser
    >>> mg = RasterModelGrid((5, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> np.allclose(mg.at_node['soil_production__rate'][mg.core_nodes], 1.)
    True
    >>> DDdiff.soilflux(2.)
    >>> np.allclose(mg.at_node['topographic__elevation'][mg.core_nodes], 0.)
    True
    >>> np.allclose(mg.at_node['bedrock__elevation'][mg.core_nodes], -2.)
    True
    >>> np.allclose(mg.at_node['soil__depth'][mg.core_nodes], 2.)
    True

    Now a more complicated example with a slope.

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> z += mg.node_x.copy()
    >>> BRz += mg.node_x / 2.
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> np.allclose(
    ...     mg.at_node['soil_production__rate'][mg.core_nodes],
    ...     np.array([ 0.60653066, 0.36787944, 0.22313016]))
    True
    >>> DDdiff.soilflux(0.1)
    >>> np.allclose(
    ...     mg.at_node['topographic__elevation'][mg.core_nodes],
    ...     np.array([ 1.04773024, 2.02894986, 3.01755898]))
    True
    >>> np.allclose(mg.at_node['bedrock__elevation'][mg.core_nodes],
    ...     np.array([ 0.43934693, 0.96321206, 1.47768698]))
    True
    >>> np.allclose(mg.at_node['soil__depth'], z - BRz)
    True

    The DepthDependentTaylorDiffuser makes and moves soil at a rate proportional
    to slope, this means that there is a characteristic time scale for soil \
    transport and an associated stability criteria for the timestep. The
    maximum characteristic time scale, Demax, is given as a function of the
    hillslope diffustivity, D, the maximum slope, Smax, and the critical slope
    Sc.

        Demax = D ( 1 + ( Smax / Sc )**2 ( Smax / Sc )**4 + .. + ( Smax / Sc )**( 2 * ( n - 1 )) )

    The maximum stable time step is given by

        dtmax = courant_factor * dx * dx / Demax

    Where the courant factor is a user defined scale (default is 0.2)

    The DepthDependentTaylorDiffuser has a boolean flag that permits a user
    to be warned if timesteps are too large for the slopes in the model grid
    (if_unstable = 'warn') and a boolean flag that turns on dynamic timestepping
    (dynamic_dt = False).

    >>> DDdiff.soilflux(2., if_unstable='warn')
    Topographic slopes are high enough such that the Courant condition is
    exceeded AND you have not selected dynamic timestepping with
    dynamic_dt=True. This may lead to infinite and/or nan values for slope,
    elevation, and soil depth. Consider using a smaller time step or dynamic
    timestepping. The Courant condition recommends a timestep of
    0.0953407607307 or smaller.

    Alternatively you can specify if_unstable='raise', and a Runtime Error will
    be raised if this condition is not met.

    Next, lets do an example with dynamic timestepping.

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')

    We'll use a steep slope and very little soil.

    >>> z += mg.node_x.copy()**2
    >>> BRz = z.copy() - 1.0
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()

    Lets try to move the soil with a large timestep. Without dynamic time
    steps, this gives a warning that we've exceeded the dynamic timestep size
    and should use a smaller timestep. We could either use the smaller timestep,
    or specify that we want to use a dynamic timestep.

    >>> DDdiff.soilflux(10, if_unstable='warn', dynamic_dt=False)
    Topographic slopes are high enough such that the Courant condition is
    exceeded AND you have not selected dynamic timestepping with
    dynamic_dt=True. This may lead to infinite and/or nan values for slope,
    elevation, and soil depth. Consider using a smaller time step or dynamic
    timestepping. The Courant condition recommends a timestep of
    0.004 or smaller.

    Now, we'll re-build the grid and do the same example with dynamic timesteps.

    >>> mg = RasterModelGrid((3, 5))
    >>> soilTh = mg.add_zeros('node', 'soil__depth')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> BRz = mg.add_zeros('node', 'bedrock__elevation')
    >>> z += mg.node_x.copy()**2
    >>> BRz = z.copy() - 1.0
    >>> soilTh[:] = z - BRz
    >>> expweath = ExponentialWeatherer(mg)
    >>> DDdiff = DepthDependentTaylorDiffuser(mg)
    >>> expweath.calc_soil_prod_rate()
    >>> DDdiff.soilflux(10, if_unstable='warn', dynamic_dt=True)
    >>> np.any(np.isnan(z))
    False
    """

    _name = "DepthDependentTaylorDiffuser"

    _input_var_names = (
        "topographic__elevation",
        "soil__depth",
        "soil_production__rate",
    )

    _output_var_names = (
        "soil__flux",
        "topographic__slope",
        "topographic__elevation",
        "bedrock__elevation",
        "soil__depth",
    )

    _var_units = {
        "topographic__elevation": "m",
        "topographic__slope": "m/m",
        "soil__depth": "m",
        "soil__flux": "m^2/yr",
        "soil_production__rate": "m/yr",
        "bedrock__elevation": "m",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "topographic__slope": "link",
        "soil__depth": "node",
        "soil__flux": "link",
        "soil_production__rate": "node",
        "bedrock__elevation": "node",
    }

    _var_doc = {
        "topographic__elevation": "elevation of the ground surface",
        "topographic__slope": "gradient of the ground surface",
        "soil__depth": "depth of soil/weather bedrock",
        "soil__flux": "flux of soil in direction of link",
        "soil_production__rate": "rate of soil production at nodes",
        "bedrock__elevation": "elevation of the bedrock surface",
    }

    def __init__(
        self,
        grid,
        linear_diffusivity=1.0,
        slope_crit=1.0,
        soil_transport_decay_depth=1.0,
        nterms=2,
    ):
        """Initialize the DepthDependentTaylorDiffuser.

        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        linear_diffusivity: float, optional.
            Hillslope diffusivity, m**2/yr
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
        """
        # Store grid and parameters
        self._grid = grid
        self.K = linear_diffusivity
        self.soil_transport_decay_depth = soil_transport_decay_depth
        self.slope_crit = slope_crit
        self.nterms = nterms

        # create fields
        # elevation
        if "topographic__elevation" in self.grid.at_node:
            self.elev = self.grid.at_node["topographic__elevation"]
        else:
            self.elev = self.grid.add_zeros("node", "topographic__elevation")

        # slope
        if "topographic__slope" in self.grid.at_link:
            self.slope = self.grid.at_link["topographic__slope"]
        else:
            self.slope = self.grid.add_zeros("link", "topographic__slope")

        # soil depth
        if "soil__depth" in self.grid.at_node:
            self.depth = self.grid.at_node["soil__depth"]
        else:
            self.depth = self.grid.add_zeros("node", "soil__depth")

        # soil flux
        if "soil__flux" in self.grid.at_link:
            self.flux = self.grid.at_link["soil__flux"]
        else:
            self.flux = self.grid.add_zeros("link", "soil__flux")

        # weathering rate
        if "soil_production__rate" in self.grid.at_node:
            self.soil_prod_rate = self.grid.at_node["soil_production__rate"]
        else:
            self.soil_prod_rate = self.grid.add_zeros("node", "soil_production__rate")

        # bedrock elevation
        if "bedrock__elevation" in self.grid.at_node:
            self.bedrock = self.grid.at_node["bedrock__elevation"]
        else:
            self.bedrock = self.grid.add_zeros("node", "bedrock__elevation")

    def soilflux(self, dt, dynamic_dt=False, if_unstable="pass", courant_factor=0.2):
        """Calculate soil flux for a time period 'dt'.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        dynamic_dt : boolean (optional, default is False)
            Keyword argument to turn on or off dynamic time-stepping
        if_unstable : string (optional, default is "pass")
            Keyword argument to determine how potential instability due to
            slopes that are too high is handled. Options are "pass", "warn",
            and "raise".
        courant_factor : float (optional, default = 0.2)
            Factor to identify stable time-step duration when using dynamic
            timestepping.
        """
        # establish time left as all of dt
        time_left = dt

        # begin while loop for time left
        while time_left > 0.0:

            # calculate soil__depth
            self.grid.at_node["soil__depth"][:] = (
                self.grid.at_node["topographic__elevation"]
                - self.grid.at_node["bedrock__elevation"]
            )

            # Calculate soil depth at links.
            self.H_link = self.grid.map_value_at_max_node_to_link(
                "topographic__elevation", "soil__depth"
            )

            # Calculate gradients
            self.slope = self.grid.calc_grad_at_link(self.elev)
            self.slope[self.grid.status_at_link == INACTIVE_LINK] = 0.0

            # Test for time stepping courant condition
            # Test for time stepping courant condition
            courant_slope_term = 0.0
            courant_s_over_scrit = self.slope.max() / self.slope_crit
            for i in range(0, 2 * self.nterms, 2):
                courant_slope_term += courant_s_over_scrit ** i
                if np.any(np.isinf(courant_slope_term)):
                    message = (
                        "Soil flux term is infinite in Courant condition "
                        "calculation. This is likely due to "
                        "using too many terms in the Taylor expansion."
                    )
                    raise RuntimeError(message)
            # Calculate De Max
            De_max = self.K * (courant_slope_term)
            # Calculate longest stable timestep
            self.dt_max = courant_factor * (self.grid.dx ** 2) / De_max

            # Test for the Courant condition and print warning if user intended
            # for it to be printed.
            if (self.dt_max < dt) and (not dynamic_dt) and (if_unstable != "pass"):
                message = (
                    "Topographic slopes are high enough such that the "
                    "Courant condition is exceeded AND you have not "
                    "selected dynamic timestepping with dynamic_dt=True. "
                    "This may lead to infinite and/or nan values for "
                    "slope, elevation, and soil depth. Consider using a "
                    "smaller time step or dynamic timestepping. The "
                    "Courant condition recommends a timestep of "
                    "" + str(self.dt_max) + " or smaller."
                )
                if if_unstable == "raise":
                    raise RuntimeError(message)
                if if_unstable == "warn":
                    print(message)

            # if dynamic dt is selected, use it, otherwise, use the entire time
            if dynamic_dt:
                self.sub_dt = np.min([dt, self.dt_max])
                time_left -= self.sub_dt
            else:
                self.sub_dt = dt
                time_left = 0

            # update sed flux, topography, soil, and bedrock based on the
            # current self.sub_dt
            self._update_flux_topography_soil_and_bedrock()

    def _update_flux_topography_soil_and_bedrock(self):
        """Calculate soil flux and update topography. """
        # Calculate flux
        slope_term = 0.0
        s_over_scrit = self.slope / self.slope_crit
        for i in range(0, 2 * self.nterms, 2):
            slope_term += s_over_scrit ** i
            if np.any(np.isinf(slope_term)):
                message = (
                    "Soil flux term is infinite. This is likely due to "
                    "using too many terms in the Taylor expansion."
                )
                raise RuntimeError(message)

        self.flux[:] = -(
            (self.K * self.slope)
            * (slope_term)
            * (1.0 - np.exp(-self.H_link / self.soil_transport_decay_depth))
        )

        # Calculate flux divergence
        dqdx = self.grid.calc_flux_div_at_node(self.flux)

        # Calculate change in soil depth
        dhdt = self.soil_prod_rate - dqdx

        # Calculate soil depth at nodes
        self.depth[self.grid.core_nodes] += dhdt[self.grid.core_nodes] * self.sub_dt

        # prevent negative soil thickness
        self.depth[self.depth < 0.0] = 0.0

        # Calculate bedrock elevation
        self.bedrock[self.grid.core_nodes] -= (
            self.soil_prod_rate[self.grid.core_nodes] * self.sub_dt
        )

        # Update topography
        self.elev[self.grid.core_nodes] = (
            self.depth[self.grid.core_nodes] + self.bedrock[self.grid.core_nodes]
        )

    def run_one_step(self, dt, **kwds):
        """

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        self.soilflux(dt, **kwds)
