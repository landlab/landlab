#!/usr/bin/env python3
"""Solve advection numerically using Total Variation Diminishing method."""

import numpy as np

from landlab import Component
from landlab.utils.return_array import return_array_at_node


def find_upwind_link_at_link(grid, u):
    """
    For all links, return ID of upwind link, defined based on the sign of u, or
    -1 if upwind direction has a boundary node, or if link is inactive.

    For instance (see examples below), consider a 3x4 raster grid with link
    numbering:

    .-14-.-15-.-16-.
    |    |    |    |
    10  11   12   13
    |    |    |    |
    .--7-.--8-.--9-.
    |    |    |    |
    3    4    5    6
    |    |    |    |
    .--0-.--1-.--2-.

    There are at most 7 active links (4, 5, 7, 8, 9, 11, 12). If u is positive
    everywhere, then the upwind links are:

    .----.----.----.
    |    |    |    |
    |    4    5    |
    |    |    |    |
    .----.--7-.--8-.
    |    |    |    |
    |    |    |    |
    |    |    |    |
    .----.----.----.

    If u is negative everywhere, then the upwind links are:

    .----.----.----.
    |    |    |    |
    |    |    |    |
    |    |    |    |
    .--8-.--9-.----.
    |    |    |    |
    |   11   12    |
    |    |    |    |
    .----.----.----.


    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 4))
    >>> uwl = find_upwind_link_at_link(grid, 1.0)
    >>> uwl[grid.active_links]
    array([-1, -1, -1,  7,  8,  4,  5])
    >>> uwl = find_upwind_link_at_link(grid, -1.0)
    >>> uwl[grid.active_links]
    array([11, 12,  8,  9, -1, -1, -1])
    >>> u = np.zeros(grid.number_of_links)
    >>> u[4:6] = -1
    >>> u[7] = -1
    >>> u[8:10] = 1
    >>> u[11:13] = 1
    >>> uwl = find_upwind_link_at_link(grid, u)
    >>> uwl[grid.active_links]
    array([11, 12,  8,  7,  8,  4,  5])
    """
    pll = grid.parallel_links_at_link
    uwl = np.zeros(grid.number_of_links, dtype=int)
    if isinstance(u, np.ndarray):
        uwl[:] = pll[:, 0]
        uneg = np.where(u < 0)[0]
        uwl[uneg] = pll[uneg, 1]
    elif u > 0:
        uwl[:] = pll[:, 0]
    else:
        uwl[:] = pll[:, 1]
    return uwl


def upwind_to_local_grad_ratio(grid, v, uwll, out=None):
    """
    Calculate and return ratio of upwind to local gradient in v.

    Gradients are defined on links. Upwind is pre-determined via
    parameter uwll (upwind link at link), which can be obtained
    using the find_upwind_link_at_link function.

    In Total Variation Diminishing (TVD) numerical schemes, this
    ratio is input to a flux limiter to calculate the weighting factor
    for higher-order vs. lower-order terms.
    """
    if out is None:
        out = np.ones(grid.number_of_links)
    else:
        out[:] = 1.0

    h = grid.node_at_link_head
    t = grid.node_at_link_tail
    local_diff = v[h] - v[t]
    upwind_exists_and_nonzero_local = np.logical_and(uwll != -1, local_diff != 0.0)
    out[upwind_exists_and_nonzero_local] = (
        local_diff[uwll[upwind_exists_and_nonzero_local]]
        / local_diff[upwind_exists_and_nonzero_local]
    )
    return out


class AdvectionSolverTVD(Component):
    r"""Component that implements numerical solution for advection using a
    Total Variation Diminishing method.

    The component is restricted to regular grids (e.g., Raster or Hex).

    Parameters
    ----------
    grid : RasterModelGrid or HexModelGrid
        A Landlab grid object.

    Examples
    --------

    References
    ----------
    """

    _name = "AdvectionSolverTVD"

    _unit_agnostic = True

    _info = {
        "advection__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m2/y",
            "mapping": "link",
            "doc": "Link-parallel advection flux",
        },
        "advection__velocity": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/y",
            "mapping": "link",
            "doc": "Link-parallel advection velocity magnitude",
        },
    }

    def __init__(
        self,
        grid,
        field_to_advect,
        advection_direction_is_steady=False,
    ):
        """Initialize AdvectionSolverTVD."""

        # Call base class methods to check existence of input fields,
        # create output fields, etc.
        super().__init__(grid)
        self.initialize_output_fields()

        self._scalar = return_array_at_node(field_to_advect)
        self._vel = self.grid.at_link["advection__velocity"]
        self._flux_at_link = self.grid.at_link["advection__flux"]

        self._advection_direction_is_steady = advection_direction_is_steady
        if advection_direction_is_steady:  # if so, only need to do this once
            self._upwind_link_at_link = find_upwind_link_at_link(self.grid, self._vel)

    def calc_rate_of_change_at_nodes(self):
        if not self._advection_direction_is_steady:
            self._upwind_link_at_link = find_upwind_link_at_link(self.grid, self._vel)
        s_link_low = self.grid.map_node_to_link_linear_upwind(self._scalar, self._vel)
        s_link_high = self.grid.map_node_to_link_lax_wendroff(self._scalar, self._vel)
        r = upwind_to_local_grad_ratio(
            self.grid, self._scalar, self._upwind_link_at_link
        )
        psi = flux_lim_vanleer(r)
        s_at_link = psi * s_link_high + (1.0 - psi) * s_link_low
        self._flux_at_link[self.grid.active_links] = (
            self._vel[self.grid.active_links] * s_at_link[self.grid.active_links]
        )
        return -self.grid.calc_flux_div_at_node(self.flux_at_link)

    def update(self, dt):
        roc = self.calc_rate_of_change_at_nodes()
        self._scalar[self.grid.core_nodes] += roc[self.grid.core_nodes] * dt

    def run_one_step(self, dt):
        self.update(dt)
