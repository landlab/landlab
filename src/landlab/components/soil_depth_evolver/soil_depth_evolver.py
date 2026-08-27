# -*- coding: utf-8 -*-

import numpy as np

from requireit import require_nonnegative
from requireit import require_positive

from landlab import Component


class SoilDepthEvolver(Component):
    """
    Evolve soil depth through soil production and hillslope transport.

    SoilProductionAndTransporter couples depth-dependent soil production
    with an externally supplied Landlab hillslope-transport component.
    The component tracks changes in soil thickness through time while
    allowing the user to choose the transport formulation independently.
    
    Soil is produced according to an exponential soil-production function,
    
    .. math::
    
        P = \\frac{\\rho_r}{\\rho_s} P_0
            \\exp\\left(-\\frac{H}{H_0}\\right),
    
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
    SoilProductionAndTransporter. The supplied component must operate on
    the same grid, update ``topographic__elevation``, and provide a
    ``run_one_step(dt)`` method.

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
        diffuser,
        soil_production_rate: float = 0.0003,
        soil_production_decay_depth: float = 0.5,
        rock_density: float = 2000.0,
        soil_density: float = 1600.0,
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
        diffuser : Landlab Component
            Hillslope-transport component operating on the same grid. The
            component must provide ``run_one_step(dt)`` and update
            ``topographic__elevation``.
        soil_production_rate : float, optional
            Maximum soil-production rate for zero soil thickness, in meters
            per year. Must be nonnegative. Default is 0.0003 m/yr.
        soil_production_decay_depth : float, optional
            Characteristic soil depth controlling the exponential decline in
            soil-production rate, in meters. Must be positive. Default is
            0.5 m.
        rock_density : float, optional
            Density of parent rock, in kilograms per cubic meter. Must be
            positive. Default is 2000 kg/m3.
        soil_density : float, optional
            Bulk density of produced soil, in kilograms per cubic meter.
            Must be positive. Default is 1600 kg/m3.
        """

        super().__init__(grid)

        soil_production_rate = require_nonnegative(
            soil_production_rate,
            name="soil_production_rate",
        )
    
        soil_production_decay_depth = require_positive(
            soil_production_decay_depth,
            name="soil_production_decay_depth",
        )
    
        rock_density = require_positive(
            rock_density,
            name="rock_density",
        )
    
        soil_density = require_positive(
            soil_density,
            name="soil_density",
        )

        self._maximum_production_rate = float(soil_production_rate)
        self._production_decay_depth = float(soil_production_decay_depth)

        self._rock_density = float(rock_density)
        self._soil_density = float(soil_density)

        self._density_ratio = (self._rock_density / self._soil_density)

        # Diffuser is created in the driver and passed in.
        self._diffuser = diffuser

        # Landlab fields
        self._elevation = self.grid.at_node[
            "topographic__elevation"
        ]

        self._soil_depth = self.grid.at_node[
            "soil__depth"
        ]

        if np.any(self._soil_depth < 0.0):
            raise ValueError(
                "soil__depth cannot contain negative values."
            )

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
            raise ValueError(
                "dt must be a positive, finite number."
            )

        elevation_before = self._elevation.copy()
        soil_depth_before = self._soil_depth.copy()

        production_rate = (
            self._density_ratio
            * self._maximum_production_rate
            * np.exp(
                -self._soil_depth
                / self._production_decay_depth
            )
        )

        soil_produced = production_rate * dt


        self._diffuser.run_one_step(dt)

        elevation_change = (
            self._elevation - elevation_before
        )


        self._update_soil_depth(
            soil_depth_before=soil_depth_before,
            elevation_change=elevation_change,
            soil_produced=soil_produced,
        )


        soil_depth_change = (
            self._soil_depth - soil_depth_before
        )

        self._current_time += dt

        return {
            "soil_depth": self._soil_depth.copy(),
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

        updated_depth = (
            soil_depth_before
            + elevation_change
            + soil_produced
        )

        exhausted = (
            (elevation_change < 0.0)
            & (-elevation_change >= soil_depth_before)
        )

        updated_depth[exhausted] = (
            soil_produced[exhausted]
        )

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