#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 09:59:02 2020

@author: gtucker
"""

import numpy as np
from landlab import Component, RasterModelGrid, HexModelGrid
#from landlab.utils.structured_grid import neighbor_node_array


class ListricKinematicExtender(Component):
    """Apply tectonic extension and subsidence kinematically to a raster or
    hex grid.

    Extension in  raster grid can be east-west (default) or north-south.
    In a hex grid, extension is oblique, with a fault at N30E (default) and
    east-west extension, or N60E (vertically oriented grid) and north-south
    extension.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ListricKinematicExtender
    >>> grid = RasterModelGrid((3, 7), xy_spacing=2500.0)
    >>> topo = grid.add_zeros('topographic__elevation', at='node')
    >>> lke = ListricKinematicExtender(grid, fault_location=2500.0)
    >>> lke.update_subsidence_rate()
    >>> np.round(grid.at_node["subsidence_rate"][8:13], 6)
    array([ 0.      ,  0.001123,  0.000729,  0.000472,  0.000306])
    """

    _name = "ListricKinematicExtender"

    _time_units = "y"

    _info = {
        "topographic__elevation": {
            "dtype": "float",
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "subsidence_rate": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Rate of tectonic subsidence in hangingwall area",
        },
        "upper_crust_thickness": {
            "dtype": "float",
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Thickness of upper crust (arbitrary datum)",
        },
        "cumulative_subsidence_depth": {
            "dtype": "float",
            "intent": "out",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Cumulative depth of tectonic subsidence",
        },
    }

    def __init__(
        self,
        grid,
        extension_rate=0.001,
        fault_dip=60.0,
        fault_location=None,
        detachment_depth=1.0e4,
        track_crustal_thickness=False,
        fields_to_shift=[],
        extension_direction="east-west",
    ):
        """Deform vertically and horizontally to represent tectonic extension.

        Parameters
        ----------
        grid: RasterModelGrid
            A landlab grid.
        extension_rate: float, optional
            Rate of horizontal motion of hangingwall relative to footwall
            (m / y), default 1 mm/y.
        fault_dip: float, optional
            Dip of the fault, degrees (default 45).
        fault_location: float, optional
            Distance of fault trace from x=0 (m) (default = grid width / 2).
        elastic_decay_length: float, optional
            Decay length scale for vertical motion (m), default 100 km.
        extension_direction: string, optional
            'east-west' (default), 'north-south' (raster only)
        """
        is_raster = isinstance(grid, RasterModelGrid)
        is_hex = isinstance(grid, HexModelGrid)
        if not (is_raster or is_hex):
            raise TypeError("Grid must be RasterModelGrid or HexModelGrid")
        elif isinstance(grid, HexModelGrid) and (
            grid.node_layout == "hex" or grid.number_of_node_columns % 2 == 0
        ):
            raise TypeError(
                "Hex grid must have rectangular layout & odd number of columns"
            )

        super().__init__(grid)
        self.initialize_output_fields()

        self._extension_rate = extension_rate
        self._fault_grad = np.tan(np.deg2rad(fault_dip))
        self._detachment_depth = detachment_depth
        self._decay_length = detachment_depth / self._fault_grad

        self._elev = grid.at_node["topographic__elevation"]
        self._subs_rate = grid.at_node["subsidence_rate"]

        # Handle fields to shift right as hangingwall moves
        self._fields_to_shift = fields_to_shift.copy()
        self._fields_to_shift.append("topographic__elevation")

        # Handle option for tracking crustal thickness field
        self._track_thickness = track_crustal_thickness
        if self._track_thickness:
            try:
                self._thickness = grid.at_node["upper_crust_thickness"]
            except KeyError:
                raise KeyError(
                    "When handle_thickness is True you must provide an 'upper_crust_thickness' node field."
                )
            self._cum_subs = grid.add_zeros("cumulative_subsidence_depth", at="node")
            self._fields_to_shift.append("upper_crust_thickness")

        # handle grid type and extension direction
        if isinstance(grid, HexModelGrid):
            if fault_location is None:
                self._fault_loc = 0.0
            else:
                self._fault_loc = fault_location
            self._ds = grid.spacing
            if grid.orientation[0] == "h":
                phi = np.deg2rad(-30.0)
            else:
                raise NotImplementedError("vertical orientation hex grids not currently handled")
            self._fault_normal_coord = grid.x_of_node + grid.y_of_node * np.tan(phi)
        else:
            if extension_direction == "north-south":
                raise NotImplementedError("north-south extension not currently handled")
            else:
                self._fault_normal_coord = grid.x_of_node
                self._ds = grid.dy
            if fault_location is None:
                self._fault_loc = 0.5 * (
                    np.amax(self._fault_normal_coord)
                    - np.amin(self._fault_normal_coord)
                )
            else:
                self._fault_loc = fault_location

        self._elev = grid.at_node["topographic__elevation"]

        # shorthand to make the next block of code easier to read
        s = self._fault_normal_coord
        fault = self._fault_loc

        # set up data structure for horizontal shift of elevation values
        self._horiz_displacement = 0.0  # horiz displ since last grid shift
        self._current_time = (
            0.0  # cumulative time (needed for adding new nodes during shift)
        )
        self._hangwall_edge = fault
        self._footwall = np.where(s <= fault)[0]
        self._hangwall = np.where(s > fault)[0]
        self._update_hangingwall_nodes()

    def _update_hangingwall_nodes(self):
        """Update the hw_downwind and hw_upwind arrays."""
        self._hangwall = np.where(self._fault_normal_coord > self._hangwall_edge)[0]
        self._hw_downwind = np.where(
            self._fault_normal_coord > self._hangwall_edge + self._ds
        )[0]
        self._hw_upwind = self._hw_downwind - 1
        self._fault_nodes = np.logical_and(
            self._fault_normal_coord > self._hangwall_edge,
            self._fault_normal_coord <= self._hangwall_edge + self._ds,
        )

    def update_subsidence_rate(self):
        """Update subsidence rate array."""
        dist_to_fault = np.abs(self._fault_normal_coord - self._fault_loc)
        self._subs_rate[self._hangwall] = (
            self._extension_rate
            * self._fault_grad
            * np.exp(-dist_to_fault[self._hangwall] / self._decay_length)
        )

    def fault_plane_elev(self, dist_from_fault_trace):
        """Return elevation of fault plane at given distance from fault trace.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import ListricKinematicExtender
        >>> grid = RasterModelGrid((3, 7), xy_spacing=5000.0)
        >>> topo = grid.add_zeros('topographic__elevation', at='node')
        >>> lke = ListricKinematicExtender(grid, fault_location=10000.0)
        >>> np.round(lke.fault_plane_elev(np.array([0.0, 5000.0, 10000.0])))
        array([   -0., -5794., -8231.])
        """
        return -(
            self._detachment_depth
            * (
                1.0
                - np.exp(
                    -dist_from_fault_trace * self._fault_grad / self._detachment_depth
                )
            )
        )

    def debug_print(self):
        if self.grid.number_of_node_columns == 51:
            print(self._elev[61:72])

    def run_one_step(self, dt):
        """Apply extensional motion to grid for one time step."""
        self.update_subsidence_rate()
        self._elev[self._hangwall] -= self._subs_rate[self._hangwall] * dt
        self._horiz_displacement += self._extension_rate * dt
        if self._track_thickness:
            self._cum_subs[self._hangwall] += self._subs_rate[self._hangwall] * dt

        self._current_time += dt

        # Shift footwall nodes
        if self._horiz_displacement >= self._ds:

            # Shift node fields in hangingwall, including elevation, by one cell
            for fieldname in self._fields_to_shift:
                self.grid.at_node[fieldname][self._hw_downwind] = self.grid.at_node[
                    fieldname
                ][self._hw_upwind]

            # Pull elevation along the current fault nodes down to either the
            # fault plane, or their original values, whichever is lower.
            self._elev[self._fault_nodes] = np.minimum(
                self._elev[self._fault_nodes],
                self.fault_plane_elev(
                    (self._hangwall_edge - self._fault_loc) + self._ds
                ),
            )

            # Subtract cumulative subsidence from thickness, and zero out
            # cumulative subsidence (it is subsidence since latest horiz shift)
            if self._track_thickness:
                self._thickness -= self._cum_subs
                self._cum_subs[:] = 0.0

            # Update horizontal displacement (since latest shift) and
            # hangingwall locations
            self._horiz_displacement -= self._ds
            self._hangwall_edge += self._ds
            self._update_hangingwall_nodes()
