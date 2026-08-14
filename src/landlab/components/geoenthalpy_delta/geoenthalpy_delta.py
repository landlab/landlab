"""Landlab component for 2D enthalpy-based sediment diffusion transport.

This version, v1.0, only deals with a fixed basement. A future version
will account for an input (erodible) topography.
"""

import math

import numpy as np

from landlab import Component
from landlab import RasterModelGrid


class GeoEnthalpyDelta(Component):
    """
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
    ``basement__elevation`` and a mobile sediment thickness ``H``:
    ``eta = eta_b + H``. Rather than tracking ``H`` as its own field, it is
    exposed as the derived :attr:`sediment_thickness` property,
    ``topographic__elevation - basement__elevation``.

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
    ...     topset_threshold=(0.1, 0.1),
    ...     foreset_threshold=(2.0, 2.0),
    ...     topset_diffusivity=(1.0, 1.0),
    ...     foreset_diffusivity=(1.0, 1.0),
    ... )
    >>> nsteps = 50
    >>> for _ in range(nsteps):
    ...     esd.run_one_step()  # dt chosen automatically for CFL stability
    ...
    >>> model_volume = np.sum(esd.sediment_thickness) * grid.dx * grid.dy
    >>> expected_volume = esd.sediment_flux * esd.time_elapsed
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
        "topographic__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Land surface elevation (basement__elevation + sediment thickness)",
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
        topset_threshold=(0.0, 0.0),
        foreset_threshold=(2.0, 2.0),
        topset_diffusivity=(1.0, 1.0),
        foreset_diffusivity=(1.0, 1.0),
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
        topset_threshold : (float, float), optional
            Critical slope thresholds, in the grid-x and grid-y directions
            respectively, above which topset (subaerial) transport occurs.
        foreset_threshold : (float, float), optional
            Critical slope thresholds, in the grid-x and grid-y directions
            respectively, above which foreset (subaqueous) transport occurs.
        topset_diffusivity : (float, float), optional
            Diffusivities for topset transport, in the grid-x and grid-y
            directions respectively.
        foreset_diffusivity : (float, float), optional
            Diffusivities for foreset transport, in the grid-x and grid-y
            directions respectively.
        cfl : float, optional
            Courant-Friedrichs-Lewy stability factor (0, 1] used to pick a
            stable time step automatically in :meth:`run_one_step`.
        """
        if not isinstance(grid, RasterModelGrid):
            raise TypeError(
                "GeoEnthalpyDelta requires a RasterModelGrid because "
                "its explicit finite-volume scheme assumes a uniform, "
                "structured (x, y) grid."
            )

        super().__init__(grid)

        self.initialize_output_fields()

        grid.at_node["topographic__elevation"][:] = grid.at_node["basement__elevation"]

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

        self._sediment_flux = sediment_flux
        self._sea_level_start = sea_level_start
        self._sea_level_rise_rate = sea_level_rise_rate
        self._topset_threshold = topset_threshold
        self._foreset_threshold = foreset_threshold
        self._topset_diffusivity = topset_diffusivity
        self._foreset_diffusivity = foreset_diffusivity
        self._cfl = cfl

        self._time_elapsed = 0.0

    @property
    def sediment_flux(self):
        """Total volumetric sediment input rate delivered through the feeder.

        May be changed between calls to :meth:`run_one_step` to vary the
        sediment supply over the course of a run.
        """
        return self._sediment_flux

    @sediment_flux.setter
    def sediment_flux(self, sediment_flux):
        self._sediment_flux = sediment_flux

    @property
    def time_elapsed(self):
        """Cumulative model time advanced by :meth:`run_one_step`."""
        return self._time_elapsed

    @property
    def sea_level(self):
        """Sea level elevation at the current model time."""
        return self._sea_level_start + self._sea_level_rise_rate * self._time_elapsed

    @property
    def sediment_thickness(self):
        """Thickness of the mobile sediment deposit at each node.

        Calculated as ``topographic__elevation - basement__elevation``
        rather than tracked as its own field, so it is always consistent
        with the two elevation fields even if another component modifies
        them between calls to :meth:`run_one_step`.
        """
        return (
            self.grid.at_node["topographic__elevation"]
            - self.grid.at_node["basement__elevation"]
        )

    def _calc_stable_time_step(self):
        """Calculate a Courant-Friedrichs-Lewy limited stable time step.
        https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

        The candidate diffusivity is floored at 1.0, the paper's
        beta-normalized reference topset diffusivity (Sect. 4.3, where
        ``D_top^x = 1`` by construction). This keeps the estimate safe
        even if all four diffusivities are configured below that
        reference scale, at the cost of a more conservative (smaller)
        ``dt`` in that case.

        Returns
        -------
        float
            A time step that satisfies the explicit stability criterion for
            the current diffusivities and grid spacing.
        """
        d_max = max(*self._topset_diffusivity, *self._foreset_diffusivity, 1.0)
        return self._cfl / (
            2.0 * d_max * (1.0 / self.grid.dx**2 + 1.0 / self.grid.dy**2)
        )

    def _advance_substep(self, basement, thickness, dt):
        """Advance sediment thickness ``thickness`` by one explicit time step.

        Parameters
        ----------
        basement : ndarray of shape (n_grid_rows, n_grid_cols)
            Fixed basement elevation, in (x, y) grid order.
        thickness : ndarray of shape (n_grid_rows, n_grid_cols)
            Sediment thickness at the start of the step, in (x, y) grid
            order.
        dt : float
            Duration of this step.

        Returns
        -------
        ndarray
            The updated sediment thickness, in (x, y) grid order.
        """
        H = thickness
        feeder_mask = self._feeder_mask
        qdens = self.sediment_flux / (np.count_nonzero(feeder_mask) * self.grid.dy)
        Z = self.sea_level
        Ctop_x, Ctop_y = self._topset_threshold
        Cfore_x, Cfore_y = self._foreset_threshold
        Dtop_x, Dtop_y = self._topset_diffusivity
        Dfore_x, Dfore_y = self._foreset_diffusivity

        nx, ny = basement.shape
        eta = basement + H
        qx = np.zeros((nx + 1, ny))
        qy = np.zeros((nx, ny + 1))

        qx[0, feeder_mask] = qdens

        # Candidate signed lateral fluxes, positive in +y.
        slope = (eta[:, :-1] - eta[:, 1:]) / self.grid.dy
        donor_top = np.where(slope >= 0.0, eta[:, :-1] >= Z, eta[:, 1:] >= Z)
        D = np.where(donor_top, Dtop_y, Dfore_y)
        C = np.where(donor_top, Ctop_y, Cfore_y)
        excess = np.abs(slope) - C
        qy[:, 1:-1] = np.where(
            excess > 0.0, np.where(slope >= 0.0, D * excess, -D * excess), 0.0
        )

        # Lateral donor limiter.
        tiny = 1.0e-30
        out_rate = np.where(
            qy[:, :-1] < 0.0, -qy[:, :-1] / self.grid.dy, 0.0
        ) + np.where(qy[:, 1:] > 0.0, qy[:, 1:] / self.grid.dy, 0.0)
        f = (H / dt) / (out_rate + tiny)
        fac_y = np.where(out_rate > 0.0, np.clip(f, 0.0, 1.0), 1.0)
        face = qy[:, 1:-1]
        qy[:, 1:-1] = np.where(face >= 0.0, fac_y[:, :-1] * face, fac_y[:, 1:] * face)

        # Sequential downstream limiter, including incoming and lateral flux.
        # qx[i + 1] depends on qx[i], so the loop over grid-x rows is
        # inherently sequential, but each row is vectorized over grid-y.
        for i in range(nx - 1):
            slope = (eta[i, :] - eta[i + 1, :]) / self.grid.dx
            donor_top = eta[i, :] >= Z
            D = np.where(donor_top, Dtop_x, Dfore_x)
            C = np.where(donor_top, Ctop_x, Cfore_x)
            candidate = np.where(slope > C, D * (slope - C), 0.0)
            lateral_net = (qy[i, :-1] - qy[i, 1:]) / self.grid.dy
            available = qx[i, :] + self.grid.dx * (H[i, :] / dt + lateral_net)
            qx[i + 1, :] = np.minimum(candidate, np.maximum(available, 0.0))

        # Conservative update. Basement is fixed and H is clipped only at zero.
        rate = (qx[:-1, :] - qx[1:, :]) / self.grid.dx + (
            qy[:, :-1] - qy[:, 1:]
        ) / self.grid.dy
        return np.maximum(H + dt * rate, 0.0)

    def run_one_step(self, dt=None):
        """Advance the sediment diffusion model by a time step ``dt``.

        Internally, ``dt`` is divided into one or more substeps that
        satisfy the CFL stability criterion (see
        :meth:`_calc_stable_time_step`), so any ``dt`` produces a
        numerically stable result without the caller having to manage
        substepping.

        Parameters
        ----------
        dt : float, optional
            Time step duration. If not given, a single CFL-stable substep
            is taken.
        """
        stable_dt = self._calc_stable_time_step()
        if dt is None:
            dt = stable_dt
        n_steps = math.ceil(dt / stable_dt)
        substep_dt = dt / n_steps

        # _advance_substep takes 2D arrays, so reshape the flat node fields
        # into explicit H_xy, basement_xy parameters.
        H_xy = (
            (
                self.grid.at_node["topographic__elevation"]
                - self.grid.at_node["basement__elevation"]
            )
            .reshape(self.grid.shape)
            .T.copy()
        )
        basement_xy = self.grid.at_node["basement__elevation"].reshape(
            self.grid.shape
        ).T

        # Internal loop for running with user input dt
        for _ in range(n_steps):
            H_xy = self._advance_substep(basement_xy, H_xy, substep_dt)
            self._time_elapsed += substep_dt

        self.grid.at_node["topographic__elevation"][:] = (
            self.grid.at_node["basement__elevation"] + H_xy.T.reshape(-1)
        )
