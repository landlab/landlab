# -*- coding: utf-8 -*-
"""TaylorNonLinearDiffuser Component.

@author: R Glade
@author: K Barnhart
@author: G Tucker
"""

# Cubic hillslope flux component

import numpy as np

from landlab import Component, LinkStatus


class TaylorNonLinearDiffuser(Component):
    """Hillslope evolution using a Taylor Series expansion of the Andrews-
    Bucknam formulation of nonlinear hillslope flux derived following following
    Ganti et al., 2012. The flux is given as:

        qs = KS ( 1 + (S/Sc)**2 + (S / Sc)**4 + .. + (S / Sc)**2(n - 1) )

    where K is is the diffusivity, S is the slope, Sc is the critical slope, and
    n is the number of terms.

    The default behavior uses two terms to produce a flux law as described by
    Equation 6 of Ganti et al., (2012).

    Parameters
    ----------
    grid: ModelGrid
        Landlab ModelGrid object
    linear_diffusivity: float, optional
        Hillslope diffusivity, m**2/yr
        Default = 1.0
    slope_crit: float, optional
        Critical slope
        Default = 1.0
    nterms: int, optional
        number of terms in the Taylor expansion.
        Two terms (Default) gives the behavior
        described in Ganti et al. (2012).

    Examples
    --------
    >>> import numpy as np
    >>> import decimal
    >>> from landlab import RasterModelGrid
    >>> from landlab.plot.imshow import imshow_grid
    >>> mg = RasterModelGrid((3, 3))
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> initial_slope=1.0
    >>> leftmost_elev=1000.
    >>> z[:] = leftmost_elev
    >>> z[:] += (
    ...    (initial_slope * np.amax(mg.x_of_node)) -
    ...    (initial_slope * mg.x_of_node))
    >>> mg.set_closed_boundaries_at_grid_edges(False, True, False, True)
    >>> cubicflux = TaylorNonLinearDiffuser(mg, slope_crit=0.1)
    >>> cubicflux.run_one_step(1.)
    >>> np.allclose(
    ...     mg.at_node['topographic__elevation'],
    ...     np.array([ 1002.,  1001.,  1000.,  1002.,  1001.,  1000.,  1002.,  1001.,
    ...     1000.]))
    True

    The TaylorNonLinearDiffuser makes and moves soil at a rate proportional
    to slope, this means that there is a characteristic time scale for soil
    transport and an associated stability criteria for the timestep. The
    maximum characteristic time scale, Demax, is given as a function of the
    hillslope diffustivity, D, the maximum slope, Smax, and the critical slope
    Sc.

        Demax = D ( 1 + ( Smax / Sc )**2 ( Smax / Sc )**4 + .. + ( Smax / Sc )**( 2 * ( n - 1 )) )

    The maximum stable time step is given by

        dtmax = courant_factor * dx * dx / Demax

    Where the courant factor is a user defined scale (default is 0.2)

    The TaylorNonLinearDiffuser has a boolean flag that permits a user to be
    warned if timesteps are too large for the slopes in the model grid
    (if_unstable = 'warn') and a boolean flag that turns on dynamic
    timesteppping (dynamic_dt = False).

    >>> cubicflux = TaylorNonLinearDiffuser(mg, slope_crit=0.1, if_unstable='warn')
    >>> cubicflux.run_one_step(2.)
    Topographic slopes are high enough such that the Courant condition is
    exceeded AND you have not selected dynamic timestepping with
    dynamic_dt=True. This may lead to infinite and/or nan values for slope,
    elevation, and soil depth. Consider using a smaller time step or dynamic
    timestepping. The Courant condition recommends a timestep of  0.2 or
    smaller.

    Alternatively you can specify if_unstable='raise', and a Runtime Error will
    be raised if this condition is not met.

    Next, lets do an example with dynamic timestepping.

    >>> mg = RasterModelGrid((5, 5))
    >>> z = mg.add_zeros("topographic__elevation", at="node")

    We'll use a steep slope.

    >>> z += mg.node_x.copy()**2
    >>> cubicflux = TaylorNonLinearDiffuser(
    ...     mg,
    ...     if_unstable='warn',
    ...     dynamic_dt=False)

    Lets try to move the soil with a large timestep. Without dynamic time
    steps, this gives a warning that we've exceeded the dynamic timestep size
    and should use a smaller timestep. We could either use the smaller timestep,
    or specify that we want to use a dynamic timestep.

    >>> cubicflux.run_one_step(10.)
    Topographic slopes are high enough such that the Courant condition is
    exceeded AND you have not selected dynamic timestepping with
    dynamic_dt=True. This may lead to infinite and/or nan values for slope,
    elevation, and soil depth. Consider using a smaller time step or dynamic
    timestepping. The Courant condition recommends a timestep of
    0.004 or smaller.

    Now, we'll re-build the grid and do the same example with dynamic
    timesteps.

    >>> mg = RasterModelGrid((5, 5))
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> z += mg.node_x.copy()**2
    >>> cubicflux = TaylorNonLinearDiffuser(
    ...     mg,
    ...     if_unstable='warn',
    ...     dynamic_dt=True)
    >>> cubicflux.run_one_step(10.)
    >>> np.any(np.isnan(z))
    False
    """

    _name = "TaylorNonLinearDiffuser"

    _info = {
        "soil__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m^2/yr",
            "mapping": "link",
            "doc": "flux of soil in direction of link",
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
        linear_diffusivity=1.0,
        slope_crit=1.0,
        nterms=2,
        dynamic_dt=False,
        if_unstable="pass",
        courant_factor=0.2,
    ):
        """Initialize the TaylorNonLinearDiffuser.

        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        linear_diffusivity: float, optional
            Hillslope diffusivity, m**2/yr
            Default = 1.0
        slope_crit: float, optional
            Critical slope
            Default = 1.0
        nterms: int, optional. default = 2
            number of terms in the Taylor expansion.
            Two terms (Default) gives the behavior
            described in Ganti et al. (2012).
        dynamic_dt : boolean (optional, default is False)
            Keyword argument to turn on or off dynamic time-stepping.
        if_unstable : string (optional, default is "pass")
            Keyword argument to determine how potential instability due to
            slopes that are too high is handled. Options are "pass", "warn",
            and "raise".
        courant_factor : float (optional, default = 0.2)
            Factor to identify stable time-step duration when using dynamic
            timestepping.
        """
        super(TaylorNonLinearDiffuser, self).__init__(grid)

        # Store grid and parameters

        self._K = linear_diffusivity
        self._slope_crit = slope_crit
        self._nterms = nterms
        self._dynamic_dt = dynamic_dt
        self._courant_factor = courant_factor
        self._if_unstable = if_unstable

        # Create fields:

        # elevation
        self._elev = self._grid.at_node["topographic__elevation"]

        # slope gradient
        if "topographic__slope" in self._grid.at_link:
            self._slope = self._grid.at_link["topographic__slope"]
        else:
            self._slope = self._grid.add_zeros("topographic__slope", at="link")

        # soil flux
        if "soil__flux" in self._grid.at_link:
            self._flux = self._grid.at_link["soil__flux"]
        else:
            self._flux = self._grid.add_zeros("soil__flux", at="link")

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

            # Calculate gradients
            self._slope[:] = self._grid.calc_grad_at_link(self._elev)
            self._slope[self._grid.status_at_link == LinkStatus.INACTIVE] = 0.0

            # Test for time stepping courant condition
            courant_slope_term = 0.0
            courant_s_over_scrit = self._slope.max() / self._slope_crit
            for i in range(0, 2 * self._nterms, 2):
                courant_slope_term += courant_s_over_scrit ** i
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
            self._dt_max = self._courant_factor * (self._grid.dx ** 2) / De_max

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

            # Calculate flux
            slope_term = 0.0
            s_over_scrit = self._slope / self._slope_crit
            for i in range(0, 2 * self._nterms, 2):
                slope_term += s_over_scrit ** i
                if np.any(np.isinf(slope_term)):
                    message = (
                        "Soil flux term is infinite. This is likely due to "
                        "using too many terms in the Taylor expansion."
                    )
                    raise RuntimeError(message)

            self._flux[:] = -((self._K * self._slope) * (slope_term))

            # Calculate flux divergence
            dqdx = self._grid.calc_flux_div_at_node(self._flux)

            # Update topography
            self._elev[self._grid.core_nodes] -= (
                dqdx[self._grid.core_nodes] * self._sub_dt
            )

    def run_one_step(self, dt):
        """Advance cubic soil flux component by one time step of size dt.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        self.soilflux(dt)
