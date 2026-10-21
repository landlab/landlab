import numpy as np
from requireit import require_nonnegative
from requireit import require_positive

from landlab import Component


class SoilDepthEvolver(Component):
    """
    Evolve soil depth through soil production and hillslope transport.

    SoilDepthEvolver couples depth-dependent soil production
    with an externally supplied Landlab hillslope-transport component.
    The component tracks changes in soil thickness through time while
    allowing the user to choose the transport formulation independently.

    Soil is produced according to an exponential soil-production function::

      production_rate = rock_to_soil_density_ratio
                      * maximum_production_rate
                      * exp(-soil_depth / decay_depth)

    where ``P`` is the soil-production rate, ``rho_r`` is rock density,
    ``rho_s`` is soil density, ``P_0`` is the maximum soil-production
    rate, ``H`` is the current soil depth, and ``H_0`` is the
    soil-production decay depth.

    During each call to :meth:`run_one_step`, the component first
    calculates soil production from the soil depth at the beginning of
    the timestep. It then advances the supplied hillslope-transport
    component, which modifies ``topographic__elevation`` in place.
    Surface-elevation change is calculated from the difference between
    the pre- and post-transport elevations.

    Soil depth is updated from the combined effects of surface-elevation
    change and soil production. Positive elevation change increases soil
    thickness, whereas negative elevation change removes soil. If modeled
    surface lowering exceeds the soil thickness available at the beginning
    of the timestep, the pre-existing soil is treated as exhausted and
    only soil produced during that timestep remains.

    The hillslope-transport component is supplied externally rather than
    constructed internally. This allows different compatible Landlab
    transport components to be used without modifying
    SoilDepthEvolver. The supplied component must operate on
    the same grid, update ``topographic__elevation``, and provide a
    ``run_one_step(dt)`` method.

    Notes
    --------

    This component was developed in support of the RESET (Recurring Soil
    Evacuation in Topographic Hollows) modeling framework, which couples
    soil production, hillslope transport, and slope-stability modeling to
    investigate soil infilling and shallow-landslide recurrence in
    topographic hollows of the Oregon Coast Range. The associated manuscript
    is currently under review at *Journal of Geophysical Research: Earth
    Surface*.


    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import TaylorNonLinearDiffuser
    >>> grid = RasterModelGrid((5, 5), xy_spacing=1.0)
    >>> elevation = grid.add_zeros("topographic__elevation", at="node")
    >>> elevation[:] = grid.node_x * 0.1
    >>> soil_depth = grid.add_zeros("soil__depth", at="node")
    >>> soil_depth[:] = 0.5
    >>> diffuser = TaylorNonLinearDiffuser(
    ...     grid,
    ...     linear_diffusivity=0.0042,
    ...     slope_crit=1.25,
    ...     nterms=2,
    ...     dynamic_dt=True,
    ... )
    >>> component = SoilDepthEvolver(
    ...     grid,
    ...     diffuser=diffuser,
    ... )
    >>> result = component.run_one_step(1.0)
    >>> component.current_time
    1.0
    >>> result["soil_depth"].shape
    (25,)
    >>> bool((result["soil_depth"] >= 0.0).all())
    True
    """

    _name = "SoilDepthEvolver"

    _unit_agnostic = True

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Surface elevation.",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Current soil thickness.",
        },
    }

    def __init__(
        self,
        grid,
        *,
        diffuser=None,
        soil_production_rate: float = 0.0003,
        soil_production_decay_depth: float = 0.5,
        rock_to_soil_density_ratio: float = 1.25,
    ) -> None:
        """SoilDepthEvolver.

        Notes
        -----
        ``topographic__elevation`` and ``soil__depth`` must already exist as
        node fields on the grid before this component is constructed.

        The component does not define initial soil-depth distributions,
        geomorphic regions, masks, or model boundary conditions. These should
        be configured by the caller before initialization.

        Timestep diagnostics are returned by :meth:`run_one_step` rather than
        stored as additional Landlab fields. This includes soil depth,
        elevation change, soil-depth change, soil produced, and soil-production
        rate.

        Parameters
        ----------
        grid : ModelGrid
            Landlab grid containing ``topographic__elevation`` and
            ``soil__depth`` at nodes.
        diffuser : Landlab Component or None, optional
            Hillslope-transport component operating on the same grid. If provided,
            the component must provide ``run_one_step(dt)`` and update
            ``topographic__elevation``. If None, no sediment transport is applied.
            Default is None.
        soil_production_rate : float, optional
            Maximum soil-production rate for zero soil thickness, in meters
            per year. Must be nonnegative. Default is 0.0003 m/yr.
        soil_production_decay_depth : float, optional
            Characteristic soil depth controlling the exponential decline in
            soil-production rate, in meters. Must be positive. Default is
            0.5 m.
        rock_to_soil_density_ratio : float, optional
            Ratio of parent-rock density to produced-soil bulk density.
            Must be positive. Default is 1.25.
        """

        super().__init__(grid)

        self._maximum_production_rate = float(
            require_nonnegative(
                soil_production_rate,
                name="soil_production_rate",
            )
        )

        self._production_decay_depth = float(
            require_positive(
                soil_production_decay_depth,
                name="soil_production_decay_depth",
            )
        )

        self._rock_to_soil_density_ratio = float(
            require_positive(
                rock_to_soil_density_ratio, name="rock_to_soil_density_ratio"
            )
        )

        # Diffuser is created in the driver and passed in.
        self._diffuser = diffuser
        self._current_time = 0.0

    def run_one_step(self, dt: float) -> dict:
        """Advance soil production and transport by one timestep.

        Parameters
        ----------
        dt : float
            Model timestep in years. Must be positive and finite.

        Returns
        -------
        dict
            Timestep diagnostics containing ``soil_depth``,
            ``elevation_change``, ``soil_depth_change``,
            ``soil_produced``, and ``production_rate``.
        """

        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be a positive, finite number.")

        elevation = self.grid.at_node["topographic__elevation"]
        soil_depth = self.grid.at_node["soil__depth"]

        if np.any(soil_depth < 0.0):
            raise ValueError("soil__depth cannot contain negative values.")

        elevation_before = elevation.copy()
        soil_depth_before = soil_depth.copy()

        production_rate = (
            self._rock_to_soil_density_ratio
            * self._maximum_production_rate
            * np.exp(-soil_depth / self._production_decay_depth)
        )

        soil_produced = production_rate * dt

        if self._diffuser is not None:
            self._diffuser.run_one_step(dt)

        elevation_change = elevation - elevation_before

        self._update_soil_depth(
            soil_depth=soil_depth,
            soil_depth_before=soil_depth_before,
            elevation_change=elevation_change,
            soil_produced=soil_produced,
        )

        soil_depth_change = soil_depth - soil_depth_before

        self._current_time += dt

        return {
            "soil_depth": soil_depth.copy(),
            "elevation_change": elevation_change,
            "soil_depth_change": soil_depth_change,
            "soil_produced": soil_produced,
            "production_rate": production_rate,
        }

    def _update_soil_depth(
        self,
        *,
        soil_depth_before: np.ndarray,
        elevation_change: np.ndarray,
        soil_produced: np.ndarray,
    ) -> None:
        """
        Update soil depth after transport and soil production.

        Parameters
        ----------
        soil_depth_before : ndarray
            Soil thickness at the beginning of the timestep.
        elevation_change : ndarray
            Surface-elevation change caused by hillslope transport.
        soil_produced : ndarray
            Soil thickness produced during the timestep.
        """

        updated_depth = soil_depth_before + elevation_change + soil_produced

        exhausted = (elevation_change < 0.0) & (-elevation_change >= soil_depth_before)

        updated_depth[exhausted] = soil_produced[exhausted]

        np.maximum(
            updated_depth,
            0.0,
            out=updated_depth,
        )

        self._soil_depth[:] = updated_depth

    @property
    def current_time(self) -> float:
        """Elapsed model time in years."""
        return self._current_time
