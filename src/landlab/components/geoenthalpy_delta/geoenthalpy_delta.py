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
    ``sediment__thickness`` ``H``: ``eta = eta_b + H``. The basement is not
    a field; it is computed once at construction as
    ``topographic__elevation - sediment__thickness`` and held fixed
    internally for the lifetime of the component.

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
    >>> esd = GeoEnthalpyDelta(
    ...     grid,
    ...     topset_threshold=(0.1, 0.1),
    ...     foreset_threshold=(2.0, 2.0),
    ...     topset_diffusivity=(1.0, 1.0),
    ...     foreset_diffusivity=(1.0, 1.0),
    ... )
    >>> nsteps = 50
    >>> for _ in range(nsteps):
    ...     esd.run_one_step()  # dt chosen automatically for CFL stability
    ...     esd.sea_level += 0.1  # sea level rises by 0.1 per step, set externally
    ...
    >>> model_volume = np.sum(grid.at_node["sediment__thickness"]) * grid.dx * grid.dy
    >>> expected_volume = np.sum(influx) * esd.time_elapsed
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
        "sediment__thickness": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Thickness of the mobile sediment deposit",
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
        topset_threshold=0.0,
        foreset_threshold=2.0,
        topset_diffusivity=1.0,
        foreset_diffusivity=1.0,
        cfl=0.4,
    ):
        """Initialize the GeoEnthalpyDelta component.

        ``topographic__elevation``, ``sea_level__elevation``, and
        ``sediment__influx`` must already exist on the grid before this
        component is constructed. ``sediment__thickness`` is optional; if
        not supplied, it defaults to zero everywhere. The component has no
        ``basement__elevation`` field: a fixed, non-erodible basement is
        computed once at construction as
        ``topographic__elevation - sediment__thickness`` and held
        internally, never updated afterward. This means the initial
        ``topographic__elevation`` you supply (e.g. from a DEM) is only
        ever the basement if ``sediment__thickness`` is zero; if you seed
        a nonzero initial thickness, ``topographic__elevation`` should
        already include it.

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

        self.initialize_output_fields()

        self._basement = (
            grid.at_node["topographic__elevation"]
            - grid.at_node["sediment__thickness"]
        )

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
        H = thickness
        qdens = influx[0, :] / self.grid.dy
        Z = self.sea_level
        Ctop_x, Ctop_y = self._topset_threshold
        Cfore_x, Cfore_y = self._foreset_threshold
        Dtop_x, Dtop_y = self._topset_diffusivity
        Dfore_x, Dfore_y = self._foreset_diffusivity

        nx, ny = basement.shape
        eta = basement + H
        qx = np.zeros((nx + 1, ny))
        qy = np.zeros((nx, ny + 1))

        qx[0, :] = qdens

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
            self.grid.at_node["sediment__thickness"].reshape(self.grid.shape).T.copy()
        )
        basement_xy = self._basement.reshape(self.grid.shape).T
        influx_xy = (
            self.grid.at_node["sediment__influx"].reshape(self.grid.shape).T
        )

        # Internal loop for running with user input dt
        for _ in range(n_steps):
            H_xy = self._advance_substep(basement_xy, H_xy, influx_xy, substep_dt)
            self._time_elapsed += substep_dt

        thickness = H_xy.T.reshape(-1)
        self.grid.at_node["sediment__thickness"][:] = thickness
        self.grid.at_node["topographic__elevation"][:] = self._basement + thickness
