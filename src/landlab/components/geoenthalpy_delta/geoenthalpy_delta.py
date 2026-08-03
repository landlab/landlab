"""Landlab component for 2D enthalpy-based sediment diffusion transport.

This version, v1.0, only deals with a fixed basement. A future version
will account for an input (erodible) topography.
"""

import numpy as np

from landlab import Component
from landlab import RasterModelGrid


def _enthalpy_diffusion_step(
    basement,
    H,
    feeder_mask,
    qdens,
    Z,
    Ctop_x,
    Ctop_y,
    Cfore_x,
    Cfore_y,
    Dtop_x,
    Dtop_y,
    Dfore_x,
    Dfore_y,
    dx,
    dy,
    dt,
):
    """Advance the sediment thickness ``H`` by one explicit time step.

    This is a free function (rather than a method) so it can be called
    directly, outside of a :class:`GeoEnthalpyDelta` instance, by users who
    want to drive the same finite-volume update from their own loop (e.g.
    to modify ``H`` between steps for a different project).
    """
    nx, ny = basement.shape
    eta = basement + H
    Hnew = np.zeros((nx, ny))
    qx = np.zeros((nx + 1, ny))
    qy = np.zeros((nx, ny + 1))
    fac_y = np.ones((nx, ny))

    for j in range(ny):
        if feeder_mask[j]:
            qx[0, j] = qdens

    # Candidate signed lateral fluxes, positive in +y.
    for i in range(nx):
        for j in range(ny - 1):
            slope = (eta[i, j] - eta[i, j + 1]) / dy
            donor_top = eta[i, j] >= Z if slope >= 0.0 else eta[i, j + 1] >= Z
            D = Dtop_y if donor_top else Dfore_y
            C = Ctop_y if donor_top else Cfore_y
            excess = abs(slope) - C
            if excess > 0.0:
                qy[i, j + 1] = D * excess if slope >= 0.0 else -D * excess

    # Lateral donor limiter.
    tiny = 1.0e-30
    for i in range(nx):
        for j in range(ny):
            out_rate = 0.0
            if qy[i, j] < 0.0:
                out_rate += -qy[i, j] / dy
            if qy[i, j + 1] > 0.0:
                out_rate += qy[i, j + 1] / dy
            if out_rate > 0.0:
                f = (H[i, j] / dt) / (out_rate + tiny)
                if f < 1.0:
                    fac_y[i, j] = max(f, 0.0)
    for i in range(nx):
        for j in range(ny - 1):
            face = qy[i, j + 1]
            qy[i, j + 1] = fac_y[i, j] * face if face >= 0.0 else fac_y[i, j + 1] * face

    # Sequential downstream limiter, including incoming and lateral flux.
    for i in range(nx - 1):
        for j in range(ny):
            slope = (eta[i, j] - eta[i + 1, j]) / dx
            donor_top = eta[i, j] >= Z
            D = Dtop_x if donor_top else Dfore_x
            C = Ctop_x if donor_top else Cfore_x
            candidate = D * (slope - C) if slope > C else 0.0
            lateral_net = (qy[i, j] - qy[i, j + 1]) / dy
            available = qx[i, j] + dx * (H[i, j] / dt + lateral_net)
            qx[i + 1, j] = min(candidate, max(available, 0.0))

    # Conservative update. Basement is fixed and H is clipped only at zero.
    for i in range(nx):
        for j in range(ny):
            rate = (qx[i, j] - qx[i + 1, j]) / dx + (qy[i, j] - qy[i, j + 1]) / dy
            Hnew[i, j] = max(H[i, j] + dt * rate, 0.0)

    return Hnew


class GeoEnthalpyDelta(Component):
    r"""
    Simulate 2D sediment diffusion, transport, and deposition using an
    enthalpy formulation of a topset/foreset delta model.

    This component is a structured Landlab wrapper around the core physics
    of a manuscript on 2D GeoEnthalpy-Delta modeling.

    Sediment is supplied through a feeder with assigned width and is transported as a nonlinear,
    slope-threshold diffusive flux. At every node the transport diffusivity and slope threshold
    depend on whether that node is a "topset" (subaerial or shallow, ``eta >= Z``) or
    "foreset" (below sea level, ``eta < Z``) node, where ``eta`` is the land
    surface elevation and ``Z`` is the (possibly time varying) sea level.

    The land surface elevation is the sum of a fixed, non-erodible
    ``basement__elevation`` and a mobile ``sediment__thickness``:

    .. math::

       \eta = \eta_b + H

    Lateral (grid-y) fluxes between neighboring nodes are calculated first,
    limited so that no node can lose more sediment than it holds. Downstream
    (grid-x) fluxes are then calculated node-by-node from west to east,
    accounting for the upstream flux and the (already limited) net lateral
    flux, which guarantees a non-negative thickness update everywhere without
    an iterative solve.

    Because the transport scheme distinguishes an upstream (grid-x) and a
    cross-stream (grid-y) direction and assumes a uniform, structured grid,
    this component requires a :class:`~.RasterModelGrid`.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import GeoEnthalpyDelta
    >>> nrows, ncols = 50, 50
    >>> dx = dy = 0.2
    >>> grid = RasterModelGrid((nrows, ncols), xy_spacing=dx)
    >>> x = grid.x_of_node.reshape(grid.shape)
    >>> x[0, 1], x[0, -1]
    (0.2, 9.8)
    >>> basement = grid.add_zeros("basement__elevation", at="node")
    >>> basement[:] = (-x).reshape(-1)  # planar, non-erodible basement E(x) = -x
    >>> basement.min(), basement.max()
    (-9.8, -0.0)
    >>> sea_level = -5.0
    >>> esd = GeoEnthalpyDelta(
    ...     grid,
    ...     sediment_flux=0.5,
    ...     feeder_width=1.0,
    ...     sea_level_start=sea_level,
    ...     sea_level_rise_rate=0.1,
    ...     topset_threshold_x=0.1,
    ...     topset_threshold_y=0.1,
    ...     foreset_threshold_x=2.0,
    ...     foreset_threshold_y=2.0,
    ...     topset_diffusivity_x=1.0,
    ...     topset_diffusivity_y=1.0,
    ...     foreset_diffusivity_x=1.0,
    ...     foreset_diffusivity_y=1.0,
    ... )
    >>> dt = esd.calc_stable_time_step()
    >>> nsteps = 50
    >>> for _ in range(nsteps):
    ...     esd.run_one_step(dt)
    ...
    >>> model_volume = np.sum(grid.at_node["sediment__thickness"]) * grid.dx * grid.dy
    >>> expected_volume = esd.sediment_flux * nsteps * dt
    >>> np.isclose(model_volume, expected_volume, rtol=1e-6)
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Lorenzo-Trueba, J., Anderson, W., Bui, V., and Voller, V. R.:
    GeoEnthalpy-Delta v1.0: an enthalpy-based model for coupled subaerial and
    subaqueous delta evolution with diagnostic moving boundaries,
    manuscript in preparation for Geoscientific Model Development.

    **Additional References**

    https://github.com/GeoJorge/GeoEnthalpy-Delta/

    """

    _name = "GeoEnthalpyDelta"

    _unit_agnostic = True  # unitless

    _info = {
        "basement__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Fixed, non-erodible basement elevation",
        },
        "sediment__thickness": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Thickness of the mobile sediment deposit",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Land surface elevation (basement__elevation + sediment__thickness)",
        },
    }

    def __init__(
        self,
        grid,
        sediment_flux=0.5,
        feeder_width=1.0,
        feeder_y_center=None,
        sea_level_start=0.0,
        sea_level_rise_rate=0.0,
        topset_threshold_x=0.0,
        topset_threshold_y=0.0,
        foreset_threshold_x=2.0,
        foreset_threshold_y=2.0,
        topset_diffusivity_x=1.0,
        topset_diffusivity_y=1.0,
        foreset_diffusivity_x=1.0,
        foreset_diffusivity_y=1.0,
        cfl=0.4,
    ):
        """Initialize the GeoEnthalpyDelta component.

        Parameters
        ----------
        grid : RasterModelGrid
        sediment_flux : float, optional
            Total volumetric sediment input rate delivered through the
            feeder (volume per time).
        feeder_width : float, optional
            Width, in grid-y units, of the sediment feeder positioned along
            the grid's west (minimum-x) edge.
        feeder_y_center : float, optional
            Center, in grid-y coordinates, of the feeder. Defaults to the
            midpoint of the grid's y-extent.
        sea_level_start : float, optional
            Sea level elevation at the time the component is instantiated.
        sea_level_rise_rate : float, optional
            Rate of change of sea level with time.
        topset_threshold_x, topset_threshold_y : float, optional
            Critical slope thresholds above which topset (subaerial) transport
            occurs, in the grid-x and grid-y directions, respectively.
        foreset_threshold_x, foreset_threshold_y : float, optional
            Critical slope thresholds above which foreset (subaqueous) transport
            occurs, in the grid-x and grid-y directions, respectively.
        topset_diffusivity_x, topset_diffusivity_y : float, optional
            Diffusivities for topset transport.
        foreset_diffusivity_x, foreset_diffusivity_y : float, optional
            Diffusivities for foreset transport.
        cfl : float, optional
            Courant-Friedrichs-Lewy stability factor (0, 1] used by
            :meth:`calc_stable_time_step`.
        """
        if not isinstance(grid, RasterModelGrid):
            raise TypeError(
                "GeoEnthalpyDelta requires a RasterModelGrid because "
                "its explicit finite-volume scheme assumes a uniform, "
                "structured (x, y) grid."
            )

        super().__init__(grid)

        self.initialize_output_fields()
        self.initialize_optional_output_fields()

        self._basement = grid.at_node["basement__elevation"]
        self._thickness = grid.at_node["sediment__thickness"]
        self._topo = grid.at_node["topographic__elevation"]
        self._topo[:] = self._basement + self._thickness

        self._nrows, self._ncols = grid.shape
        self._dx = float(grid.dx)
        self._dy = float(grid.dy)

        y_of_row = grid.y_of_node.reshape(grid.shape)[:, 0]
        if feeder_y_center is None:
            feeder_y_center = 0.5 * (y_of_row.min() + y_of_row.max())
        self._feeder_mask = (
            np.abs(y_of_row - feeder_y_center) <= 0.5 * feeder_width + 1.0e-12
        )
        if not np.any(self._feeder_mask):
            raise ValueError(
                "feeder_width is smaller than the grid-y spacing: no nodes "
                "fall within the feeder."
            )

        self.sediment_flux = sediment_flux
        self.sea_level_start = sea_level_start
        self.sea_level_rise_rate = sea_level_rise_rate
        self.topset_threshold_x = topset_threshold_x
        self.topset_threshold_y = topset_threshold_y
        self.foreset_threshold_x = foreset_threshold_x
        self.foreset_threshold_y = foreset_threshold_y
        self.topset_diffusivity_x = topset_diffusivity_x
        self.topset_diffusivity_y = topset_diffusivity_y
        self.foreset_diffusivity_x = foreset_diffusivity_x
        self.foreset_diffusivity_y = foreset_diffusivity_y
        self.cfl = cfl

        self._time_elapsed = 0.0

    @property
    def time_elapsed(self):
        """Cumulative model time advanced by :meth:`run_one_step`."""
        return self._time_elapsed

    @property
    def sea_level(self):
        """Sea level elevation at the current model time."""
        return self.sea_level_start + self.sea_level_rise_rate * self._time_elapsed

    def calc_stable_time_step(self):
        """Calculate a Courant-Friedrichs-Lewy limited stable time step.
        https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

        The candidate diffusivity is floored at 1.0, the paper's
        beta-normalized reference topset diffusivity (Sect. 4.3, where
        ``D_top^x = 1`` by construction). This keeps the estimate safe
        even if all four diffusivities are configured below that
        reference scale, at the cost of a more conservative (smaller)
        ``dt`` in that case.

        If the stable time step returned by `calc_stable_time_step()` is smaller than the time step you want
        to advance by in your own project, subcycle: run several inner steps per outer step.
        Compute the number of inner steps as `n_inner = ceil(outer_dt / calculated_dt)`,
        then call `run_one_step` that many times, each with `dt = outer_dt / n_inner`,
        which is guaranteed to be `<= calculated_dt`.

        Returns
        -------
        float
            A time step that satisfies the explicit stability criterion for
            the current diffusivities and grid spacing.
        """
        d_max = max(
            self.topset_diffusivity_x,
            self.topset_diffusivity_y,
            self.foreset_diffusivity_x,
            self.foreset_diffusivity_y,
            1.0,
        )
        return self.cfl / (2.0 * d_max * (1.0 / self._dx**2 + 1.0 / self._dy**2))

    def run_one_step(self, dt, feeder_mask=None):
        """Advance the sediment diffusion model by one time step.

        Parameters
        ----------
        dt : float
            Time step duration. For stability, ``dt`` should not exceed the
            value returned by :meth:`calc_stable_time_step`.
        feeder_mask : array_like of bool, optional
            Per-row boolean mask (length equal to the number of grid rows)
            overriding which rows act as the sediment feeder for this and
            subsequent steps. Defaults to the mask most recently set (from
            ``feeder_width``/``feeder_y_center`` at construction, or from
            an earlier call to this method).
        """
        if feeder_mask is not None:
            self._feeder_mask = feeder_mask
        H_xy = self._thickness.reshape(self._nrows, self._ncols).T.copy()
        basement_xy = self._basement.reshape(self._nrows, self._ncols).T

        qdens = self.sediment_flux / (np.count_nonzero(self._feeder_mask) * self._dy)

        Hnew_xy = _enthalpy_diffusion_step(
            basement=basement_xy,
            H=H_xy,
            feeder_mask=self._feeder_mask,
            qdens=qdens,
            Z=self.sea_level,
            Ctop_x=self.topset_threshold_x,
            Ctop_y=self.topset_threshold_y,
            Cfore_x=self.foreset_threshold_x,
            Cfore_y=self.foreset_threshold_y,
            Dtop_x=self.topset_diffusivity_x,
            Dtop_y=self.topset_diffusivity_y,
            Dfore_x=self.foreset_diffusivity_x,
            Dfore_y=self.foreset_diffusivity_y,
            dx=self._dx,
            dy=self._dy,
            dt=dt,
        )

        self._thickness[:] = Hnew_xy.T.reshape(-1)
        self._topo[:] = self._basement + self._thickness
        self._time_elapsed += dt
