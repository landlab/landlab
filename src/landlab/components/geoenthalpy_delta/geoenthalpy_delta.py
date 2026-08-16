"""Landlab component for 2D enthalpy-based sediment diffusion transport.

This version, v1.0, only deals with a fixed basement. A future version
will account for an input (erodible) topography.
"""

import math

import numpy as np
from requireit import require_between
from requireit import require_nonnegative
from requireit import require_positive

from landlab import Component
from landlab import RasterModelGrid


def _normalize_as_xy(value):
    return tuple(float(v) for v in np.broadcast_to(value, (2,)))


class GeoEnthalpyDelta(Component):
    """
    Simulate 2D sediment diffusion, transport, and deposition using an
    enthalpy formulation of a topset/foreset delta model.

    This component is a structured Landlab wrapper around the core physics
    of a manuscript on 2D GeoEnthalpy-Delta modeling.

    Sediment is supplied at the grid's west (minimum-x) edge according to the
    ``sediment__influx`` field and is transported as a nonlinear,
    slope-threshold diffusive flux. At every node the transport diffusivity and slope threshold
    depend on whether that node is a "topset" (subaerial or shallow, ``eta >= Z``) or
    "foreset" (below sea level, ``eta < Z``) node, where ``eta`` is the land
    surface elevation and ``Z`` is the (possibly time varying) sea level.

    The land surface elevation ``eta`` (``topographic__elevation``) is the
    sum of a fixed, non-erodible basement elevation ``eta_b`` and a mobile
    sediment thickness ``H``: ``eta = eta_b + H``. Neither the basement nor
    the sediment thickness is a field; both are tracked internally and
    held fixed or updated, respectively, only by this component, so
    another component sharing the grid can't corrupt them. The basement is
    computed once at construction as
    ``topographic__elevation - sediment_thickness`` (see
    :attr:`sediment_thickness` to read the current thickness) and
    held fixed for the lifetime of the component.

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
    >>> topo = grid.add_zeros("topographic__elevation", at="node")
    >>> topo[:] = (-x).reshape(-1)  # planar surface, e.g. from a DEM
    >>> topo.min(), topo.max()
    (-9.8, -0.0)
    >>> sea_level = grid.add_field("sea_level__elevation", -5.0, at="grid")
    >>> influx = grid.add_zeros("sediment__influx", at="node")
    >>> influx.reshape(grid.shape)[20:25, 0] = 0.1  # feeder on the west edge
    >>> component = GeoEnthalpyDelta(
    ...     grid,
    ...     topset_threshold=(0.1, 0.1),
    ...     foreset_threshold=(2.0, 2.0),
    ...     topset_diffusivity=(1.0, 1.0),
    ...     foreset_diffusivity=(1.0, 1.0),
    ... )
    >>> nsteps = 50
    >>> for _ in range(nsteps):
    ...     component.run_one_step()  # dt chosen automatically for CFL stability
    ...     component.sea_level += 0.1  # sea level rises by 0.1/step, set externally
    ...
    >>> model_volume = np.sum(component.sediment_thickness) * grid.dx * grid.dy
    >>> expected_volume = np.sum(influx) * component.time_elapsed
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
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "sea_level__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "grid",
            "doc": "Sea level elevation",
        },
        "sediment__influx": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Sediment flux (volume per unit time of sediment entering each node)"
            ),
        },
    }

    def __init__(
        self,
        grid,
        sediment_thickness=0.0,
        topset_threshold=0.0,
        foreset_threshold=2.0,
        topset_diffusivity=1.0,
        foreset_diffusivity=1.0,
        cfl=0.4,
    ):
        """Initialize the GeoEnthalpyDelta component.

        ``topographic__elevation``, ``sea_level__elevation``, and
        ``sediment__influx`` must already exist on the grid before this
        component is constructed. The component has no
        ``basement__elevation`` or ``sediment__thickness`` field: instead,
        it tracks a fixed, non-erodible basement and the mobile sediment
        thickness internally (see :attr:`sediment_thickness`), so
        neither can be corrupted by another component writing to a shared
        grid field. The basement is computed once at construction as
        ``topographic__elevation - sediment_thickness`` and held
        internally, never updated afterward. This means the initial
        ``topographic__elevation`` you supply (e.g. from a DEM) should
        already include any sediment thickness you pass in via
        *sediment_thickness*.

        Sea level and sediment supply are both external forcing:

        - Sea level: read and, if it varies with time, update it yourself
          via the :attr:`sea_level` property (or ``grid.at_grid``
          directly) between calls to :meth:`run_one_step`.
        - Sediment supply: set ``sediment__influx`` values (volume per
          time) at nodes along the grid's west (minimum-x) edge before
          construction. Those values define both the feeder's location
          and its spatial distribution; their sum is the total sediment
          supply rate. Update the field yourself between calls to
          :meth:`run_one_step` for a time-varying or nonuniform supply.

        Parameters
        ----------
        grid : RasterModelGrid
        sediment_thickness : float or array_like, optional
            Initial sediment thickness at each node. A scalar applies
            everywhere. Must be non-negative.
        topset_threshold : float or (float, float), optional
            Critical slope thresholds, in the grid-x and grid-y directions
            respectively, above which topset (subaerial) transport occurs.
            A scalar applies to both directions. Must be non-negative.
        foreset_threshold : float or (float, float), optional
            Critical slope thresholds, in the grid-x and grid-y directions
            respectively, above which foreset (subaqueous) transport occurs.
            A scalar applies to both directions. Must be non-negative.
        topset_diffusivity : float or (float, float), optional
            Diffusivities for topset transport, in the grid-x and grid-y
            directions respectively. A scalar applies to both directions.
            Must be positive.
        foreset_diffusivity : float or (float, float), optional
            Diffusivities for foreset transport, in the grid-x and grid-y
            directions respectively. A scalar applies to both directions.
            Must be positive.
        cfl : float, optional
            Courant-Friedrichs-Lewy stability factor used to pick a stable
            time step automatically in :meth:`run_one_step`. Must be in
            the interval (0, 1].
        """
        if not isinstance(grid, RasterModelGrid):
            raise TypeError(
                "GeoEnthalpyDelta requires a RasterModelGrid because "
                "its explicit finite-volume scheme assumes a uniform, "
                "structured (x, y) grid."
            )

        super().__init__(grid)

        sediment_thickness = require_nonnegative(
            sediment_thickness, name="sediment_thickness"
        )
        self._thickness = np.broadcast_to(
            sediment_thickness, (grid.number_of_nodes,)
        ).astype(float)
        self._basement = grid.at_node["topographic__elevation"] - self._thickness

        require_nonnegative(topset_threshold, name="topset_threshold")
        require_nonnegative(foreset_threshold, name="foreset_threshold")
        require_positive(topset_diffusivity, name="topset_diffusivity")
        require_positive(foreset_diffusivity, name="foreset_diffusivity")

        self._topset_threshold = _normalize_as_xy(topset_threshold)
        self._foreset_threshold = _normalize_as_xy(foreset_threshold)
        self._topset_diffusivity = _normalize_as_xy(topset_diffusivity)
        self._foreset_diffusivity = _normalize_as_xy(foreset_diffusivity)

        self._cfl = require_between(
            cfl, 0.0, 1.0, inclusive_min=False, inclusive_max=True, name="cfl"
        )

        self._time_elapsed = 0.0

    @property
    def time_elapsed(self):
        """Cumulative model time advanced by :meth:`run_one_step`."""
        return self._time_elapsed

    @property
    def sea_level(self):
        """Sea level elevation, read from the ``sea_level__elevation`` grid field.

        This is external forcing: update ``grid.at_grid["sea_level__elevation"]``
        (or set this property) between calls to :meth:`run_one_step` if you
        want sea level to vary with time.
        """
        return self.grid.at_grid["sea_level__elevation"]

    @sea_level.setter
    def sea_level(self, sea_level):
        self.grid.at_grid["sea_level__elevation"] = float(sea_level)

    @property
    def sediment_thickness(self):
        """Thickness of the mobile sediment deposit at each node.

        Tracked internally rather than as a field, so it can't be
        corrupted by another component writing to a shared grid field.
        """
        return self._thickness

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

    def _calc_lateral_flux(self, eta):
        """Calculate candidate signed lateral (grid-y) flux, positive in +y."""
        Ctop_y = self._topset_threshold[1]
        Cfore_y = self._foreset_threshold[1]
        Dtop_y = self._topset_diffusivity[1]
        Dfore_y = self._foreset_diffusivity[1]
        Z = self.sea_level

        nx, ny = eta.shape
        qy = np.zeros((nx, ny + 1))

        slope = (eta[:, :-1] - eta[:, 1:]) / self.grid.dy
        donor_top = np.where(slope >= 0.0, eta[:, :-1] >= Z, eta[:, 1:] >= Z)
        D = np.where(donor_top, Dtop_y, Dfore_y)
        C = np.where(donor_top, Ctop_y, Cfore_y)
        excess = np.abs(slope) - C
        qy[:, 1:-1] = np.where(
            excess > 0.0, np.where(slope >= 0.0, D * excess, -D * excess), 0.0
        )
        return qy

    def _limit_lateral_flux(self, qy, thickness, dt):
        """Limit lateral flux so no node loses more sediment than it holds."""
        tiny = 1.0e-30
        out_rate = np.where(
            qy[:, :-1] < 0.0, -qy[:, :-1] / self.grid.dy, 0.0
        ) + np.where(qy[:, 1:] > 0.0, qy[:, 1:] / self.grid.dy, 0.0)
        f = (thickness / dt) / (out_rate + tiny)
        fac_y = np.where(out_rate > 0.0, np.clip(f, 0.0, 1.0), 1.0)
        face = qy[:, 1:-1]
        limited = qy.copy()
        limited[:, 1:-1] = np.where(
            face >= 0.0, fac_y[:, :-1] * face, fac_y[:, 1:] * face
        )
        return limited

    def _calc_downstream_flux(self, eta, thickness, qy, influx, dt):
        """Calculate downstream (grid-x) flux node-by-node from west to east.

        qx[i + 1] depends on qx[i], so the loop over grid-x rows is
        inherently sequential, but each row is vectorized over grid-y.
        """
        Ctop_x = self._topset_threshold[0]
        Cfore_x = self._foreset_threshold[0]
        Dtop_x = self._topset_diffusivity[0]
        Dfore_x = self._foreset_diffusivity[0]
        Z = self.sea_level

        nx, ny = eta.shape
        qx = np.zeros((nx + 1, ny))
        qx[0, :] = influx[0, :] / self.grid.dy

        for i in range(nx - 1):
            slope = (eta[i, :] - eta[i + 1, :]) / self.grid.dx
            donor_top = eta[i, :] >= Z
            D = np.where(donor_top, Dtop_x, Dfore_x)
            C = np.where(donor_top, Ctop_x, Cfore_x)
            candidate = np.where(slope > C, D * (slope - C), 0.0)
            lateral_net = (qy[i, :-1] - qy[i, 1:]) / self.grid.dy
            available = qx[i, :] + self.grid.dx * (thickness[i, :] / dt + lateral_net)
            qx[i + 1, :] = np.minimum(candidate, np.maximum(available, 0.0))

        return qx

    def _apply_conservative_update(self, thickness, qx, qy, dt):
        """Update sediment thickness from the net flux divergence at each node."""
        rate = (qx[:-1, :] - qx[1:, :]) / self.grid.dx + (
            qy[:, :-1] - qy[:, 1:]
        ) / self.grid.dy
        return np.maximum(thickness + dt * rate, 0.0)

    def _advance_substep(self, basement, thickness, influx, dt):
        """Advance sediment thickness ``thickness`` by one explicit time step.

        Parameters
        ----------
        basement : ndarray of shape (n_grid_rows, n_grid_cols)
            Fixed basement elevation, in (x, y) grid order.
        thickness : ndarray of shape (n_grid_rows, n_grid_cols)
            Sediment thickness at the start of the step, in (x, y) grid
            order.
        influx : ndarray of shape (n_grid_rows, n_grid_cols)
            Volumetric sediment influx at each node, in (x, y) grid order.
            Only the west-boundary row (``influx[0, :]``) is used.
        dt : float
            Duration of this step.

        Returns
        -------
        ndarray
            The updated sediment thickness, in (x, y) grid order.
        """
        eta = basement + thickness
        qy = self._calc_lateral_flux(eta)
        qy = self._limit_lateral_flux(qy, thickness, dt)
        qx = self._calc_downstream_flux(eta, thickness, qy, influx, dt)
        return self._apply_conservative_update(thickness, qx, qy, dt)

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

        # _advance_substep takes 2D arrays, so reshape the flat internal
        # state into explicit H_xy, basement_xy parameters.
        H_xy = self._thickness.reshape(self.grid.shape).T.copy()
        basement_xy = self._basement.reshape(self.grid.shape).T
        influx_xy = self.grid.at_node["sediment__influx"].reshape(self.grid.shape).T

        # Internal loop for running with user input dt
        for _ in range(n_steps):
            H_xy = self._advance_substep(basement_xy, H_xy, influx_xy, substep_dt)
            self._time_elapsed += substep_dt

        self._thickness = H_xy.T.reshape(-1)
        self.grid.at_node["topographic__elevation"][:] = (
            self._basement + self._thickness
        )
