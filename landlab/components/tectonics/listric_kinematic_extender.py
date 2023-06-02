#!/usr/bin/env python3
"""
Apply tectonic extension and subsidence kinematically.

Landlab component that simulates development of an asymmetric rift on a listric
fault plane.

See notebook tutorial for complete examples.

@author: gtucker
"""

import numpy as np

from landlab import Component, RasterModelGrid, HexModelGrid
from landlab.components import AdvectionSolverTVD


def dist_to_line(Px, Py, x0, y0, alpha):
    """Calculate and return the distance of point(x) (Px, Py) to the
    line described by x = x0 + t cos alpha, y = y0 + t sin alpha.

    Parameters
    ----------
    Px : float
        x-coordinate of point(s)
    Py : float
        y-coordinate of point(s)
    x0 : float
        x intercept of line
    y0 : float
        y intercept of line
    alpha : float, degrees
        angle of line, counter-clockwise from positive x-axis

    Examples
    --------
    >>> np.round(dist_to_line(1, 1, 0, 0, 90), 6)
    1.0
    >>> np.round(dist_to_line(0, 1, 1, 0, 90), 6)
    -1.0
    >>> np.round(dist_to_line(1, 1, 0, 0, 0), 6)
    -1.0
    >>> np.round(dist_to_line(2.0**0.5, 0, 0, 0, 45), 6)
    1.0
    >>> np.round(dist_to_line(0, 2.0**0.5, 0, 0, 45), 6)
    -1.0
    """
    alpha_r = np.radians(alpha)
    return np.sin(alpha_r) * (Px - x0) - np.cos(alpha_r) * (Py - y0)


class ListricKinematicExtender(Component):
    """Apply tectonic extension and subsidence kinematically to a raster or
    hex grid.

    Extension is east-west, with a north-south fault in the case of a raster
    grid, and a N30E fault in the case of a hex grid (obique extension).

    The component requires a Raster or Hex grid.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ListricKinematicExtender
    >>> grid = RasterModelGrid((3, 7), xy_spacing=2500.0)
    >>> topo = grid.add_zeros('topographic__elevation', at='node')
    >>> lke = ListricKinematicExtender(grid)

    ...TO BE CONTINUED
    """

    _name = "ListricKinematicExtender"

    _time_units = "y"

    _unit_agnostic = True

    _info = {
        "advection__velocity": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m/y",
            "mapping": "link",
            "doc": "Link-parallel advection velocity magnitude",
        },
        "fault_plane__elevation": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Elevation of fault plane",
        },
        "hangingwall__thickness": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Thickness of material in hangingwall block",
        },
        "hangingwall__velocity": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "link",
            "doc": "Horizontal velocity of material in hangingwall block",
        },
        "topographic__elevation": {
            "dtype": "float",
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        extension_rate_x=0.001,
        extension_rate_y=0.0,
        fault_dip=60.0,
        fault_x0=0.0,
        fault_y0=0.0,
        fault_strike=45.0,
        detachment_depth=1.0e4,
        track_crustal_thickness=False,
        fields_to_advect=None,
        advection_direction_is_steady=False,
    ):
        """Deform vertically and horizontally to represent tectonic extension.

        Parameters
        ----------
        grid: RasterModelGrid
            A landlab grid.
        extension_rate_x: float, optional
            Rate of x-directed horizontal motion of hangingwall relative to footwall
            (m / y), default 0.001 m/y.
        extension_rate_y: float, optional
            Rate of y-directed horizontal motion of hangingwall relative to footwall
            (m / y), default 0.
        fault_x0: float, optional
            x intercept of zero-surface fault trace, m (default 0).
        fault_y0: float, optional
            y intercept of zero-surface fault trace, m (default 0).
        fault_strike: float, optional
            Strike of zero-surface fault trace, degrees (default 45).
        detachment_depth: float, optional
            Depth to horizontal detachment (m), default 10 km.
        track_crustal_thickness: bool, optional
            Option to keep track of changes in crustal thickness (default False)
        fields_to_advect: list of str, optional
            List of names of fields, in addition to 'hangingwall__thickness'
            and (if track_crustal_thickness==True) 'upper_crust_thickness',
            that should be advected horizontally.
        advection_direction_is_steady : bool (default False)
            Indicates whether the directions of advection are expected to remain
            steady throughout a run. If True, some computation time is saved
            by calculating upwind links only once.
        """
        if not (isinstance(grid, RasterModelGrid) or isinstance(grid, HexModelGrid)):
            raise (TypeError, "grid must be a RasterModelGrid or HexModelGrid")

        fields_to_advect = [] if fields_to_advect is None else fields_to_advect

        super().__init__(grid)
        self.initialize_output_fields()

        self._elev = grid.at_node["topographic__elevation"]
        self._fault_plane_elev = grid.at_node["fault_plane__elevation"]
        self._hw_thick = grid.at_node["hangingwall__thickness"]

        self._fields_to_advect = fields_to_advect.copy()
        self._fields_to_advect.append("hangingwall__thickness")

        self._setup_fault_plane_elevation_and_hangingwall_thickness(
            grid, fault_x0, fault_y0, fault_strike, fault_dip, detachment_depth
        )
        self._setup_advection_component(
            grid, extension_rate_x, extension_rate_y, advection_direction_is_steady
        )
        self._track_thickness = track_crustal_thickness
        if self._track_thickness:
            try:
                self._thickness = grid.at_node["upper_crust_thickness"]
            except KeyError as exc:
                raise KeyError(
                    "When handle_thickness is True you must provide an"
                    "'upper_crust_thickness' node field."
                ) from exc
            self._fields_to_advect.append("upper_crust_thickness")

    def _setup_fault_plane_elevation_and_hangingwall_thickness(
        self, grid, fault_x0, fault_y0, fault_strike, fault_dip, detachment_depth
    ):
        """Initialize fields fault_plane__elevation and hangingwall__thickness.

        Calculate and store the fault plane elevation at grid nodes using an
        exponential function of (signed) distance to fault, with topographic
        elevation as the minimum. Calculate the thickness of the hangingwall
        block at grid nodes by subtracting fault plane elevation from
        topographic elevation.

        Parameters
        ----------
        fault_x0 : float
            x-intercept of zero-surface fault trace, m
        fault_y0 : float
            y-intercept of zero-surface fault trace, m
        fault_strike : float
            strike angle of fault trace, degrees ccw from +x
        fault_dip : float
            dip angle of fault at the zero elevation point, degrees
        detachment_depth : float
            depth to the point where the detachment is horizontal, m

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import ListricKinematicExtender
        >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
        >>> _ = grid.add_zeros("topographic__elevation", at="node")
        >>> extender = ListricKinematicExtender(grid, fault_strike=90.0)
        >>> round(grid.at_node["fault_plane__elevation"][4])
        -1590
        >>> round(grid.at_node["hangingwall__thickness"][4])
        1590
        """
        fault_grad = np.tan(np.deg2rad(fault_dip))
        dist_to_fault = dist_to_line(
            grid.x_of_node, grid.y_of_node, fault_x0, fault_y0, fault_strike
        )
        self._fault_plane_elev[:] = np.minimum(
            -detachment_depth
            * (1.0 - np.exp(-dist_to_fault * fault_grad / detachment_depth)),
            self._elev,
        )
        self._hw_thick[:] = self._elev - self._fault_plane_elev

    def _setup_advection_component(
        self, grid, extension_rate_x, extension_rate_y, advection_direction_is_steady
    ):
        """Instantiate and initialize AdvectionSolverTVD.

        If the link field advection__velocity already exists and contains
        non-zero values, these values are used for the advection field.
        Otherwise, the field is created if needed, and initialized by mapping
        the vector components extension_rate_x and extension_rate_y onto the
        grid links.
        """
        if "advection__velocity" not in grid.at_link.keys():
            grid.add_zeros("advection__velocity", at="link")
        self._advec_velocity = grid.at_link["advection__velocity"]
        if np.amax(np.abs(self._advec_velocity)) == 0.0:  # if no nonzero values
            grid.map_vectors_to_links(
                extension_rate_x, extension_rate_y, out=self._advec_velocity
            )
        self.advector = AdvectionSolverTVD(
            grid,
            field_to_advect="hangingwall__thickness",
            advection_direction_is_steady=advection_direction_is_steady,
        )

    def run_one_step(self, dt):
        """Apply extensional motion to grid for one time step."""
        self.advector.run_one_step(dt)
        self._elev[self.grid.core_nodes] = (
            self._fault_plane_elev[self.grid.core_nodes]
            + self._hw_thick[self.grid.core_nodes]
        )
