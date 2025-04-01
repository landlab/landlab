#!/usr/bin/env python3
"""Solve advection numerically using Total Variation Diminishing method."""

import numpy as np

from landlab import Component
from landlab import LinkStatus
from landlab.components.advection.flux_limiters import flux_lim_vanleer
from landlab.field.errors import FieldError
from landlab.utils.return_array import return_array_at_link
from landlab.utils.return_array import return_array_at_node


def find_upwind_link_at_link(grid, u):
    """Return the upwind link at every link.

    For all links, return ID of upwind link, defined based on the sign of `u`.
    If `u` is zero, the upwind link is found as though `u` were positive.

    For instance (see examples below), consider a 3x4 raster grid with link
    numbering::

        .-14-.-15-.-16-.
        |    |    |    |
        10  11   12   13
        |    |    |    |
        .--7-.--8-.--9-.
        |    |    |    |
        3    4    5    6
        |    |    |    |
        .--0-.--1-.--2-.

    There are at most 7 active links (4, 5, 7, 8, 9, 11, 12).
    If `u` is positive everywhere, then the upwind links are::

        .----.-14-.-15-.
        |    |    |    |
        3    4    5    6
        |    |    |    |
        .----.--7-.--8-.
        |    |    |    |
        |    |    |    |
        |    |    |    |
        .----.--0-.--1-.

    If `u` is negative everywhere, then the upwind links are::

        .-15-.-16-.----.
        |    |    |    |
        |    |    |    |
        |    |    |    |
        .--8-.--9-.----.
        |    |    |    |
        10  11   12   13
        |    |    |    |
        .--1-.--2-.----.

    Parameters
    ----------
    grid : RasterModelGrid or HexModelGrid
        A landlab grid.
    u : float or (n_links,) ndarray
        Array of *at-link* values used to determine which node is
        upwind.

    Returns
    -------
    (n_links,) ndarray of int
        The upwind links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 4))

    >>> uwl = find_upwind_link_at_link(grid, 1.0)
    >>> uwl[grid.vertical_links].reshape((2, 4))
    array([[-1, -1, -1, -1],
           [ 3,  4,  5,  6]])
    >>> uwl[grid.horizontal_links].reshape((3, 3))
    array([[-1,  0,  1],
           [-1,  7,  8],
           [-1, 14, 15]])

    >>> uwl = find_upwind_link_at_link(grid, -1.0)
    >>> uwl[grid.vertical_links].reshape((2, 4))
    array([[10, 11, 12, 13],
           [-1, -1, -1, -1]])
    >>> uwl[grid.horizontal_links].reshape((3, 3))
    array([[ 1,  2, -1],
           [ 8,  9, -1],
           [15, 16, -1]])

    >>> u = np.zeros(grid.number_of_links)
    >>> u[4:6] = -1
    >>> u[7] = -1
    >>> u[8:10] = 1
    >>> u[11:13] = 1
    >>> u[grid.vertical_links].reshape((2, 4))
    array([[ 0., -1., -1.,  0.],
           [ 0.,  1.,  1.,  0.]])
    >>> u[grid.horizontal_links].reshape((3, 3))
    array([[ 0.,  0.,  0.],
           [-1.,  1.,  1.],
           [ 0.,  0.,  0.]])
    >>> uwl = find_upwind_link_at_link(grid, u)
    >>> uwl[grid.vertical_links].reshape((2, 4))
    array([[-1, 11, 12, -1],
           [ 3,  4, 5,   6]])
    >>> uwl[grid.horizontal_links].reshape((3, 3))
    array([[-1,  0,  1],
           [ 8,  7,  8],
           [-1, 14, 15]])
    """
    pll = grid.parallel_links_at_link

    cols = np.choose(np.broadcast_to(u, len(pll)) >= 0, [1, 0])
    uwl = pll[np.arange(len(pll)), cols]

    return uwl


def upwind_to_local_grad_ratio(grid, v, uwll, out=None):
    """Calculate and return ratio of upwind to local gradient in v.

    Gradients are defined on links. Upwind is pre-determined via
    parameter uwll (upwind link at link), which can be obtained
    using the find_upwind_link_at_link function.

    In Total Variation Diminishing (TVD) numerical schemes, this
    ratio is input to a flux limiter to calculate the weighting factor
    for higher-order vs. lower-order terms.

    Parameters
    ----------
    grid : RasterModelGrid or HexModelGrid
        A landlab grid.
    v : (n_links,) ndarray
        Array of *at-link* values of which to calculate the gradient.
    uwll : (n_links,) ndarray
        Array of upwind links for every link (as returned, for example, by
        :func:`~.find_upwind_link_at_link`).
    out : (n_links,) ndarray, optional
        If provided, place output into this array. Otherwise, create a new array.

    Returns
    -------
    (n_links,) ndarray of int
        The ratio of the gradients. For links that have a gradient of zero, the ratio
        is set to one. For links that do not have an upwind link, the ratio is also
        set to one.
    """
    if out is None:
        out = np.ones(grid.number_of_links)
    else:
        out[:] = 1.0

    local_diff = v[grid.node_at_link_head] - v[grid.node_at_link_tail]

    np.divide(
        local_diff[uwll], local_diff, where=(uwll != -1) & (local_diff != 0.0), out=out
    )

    return out


class AdvectionSolverTVD(Component):
    """Numerical solution for advection using a Total Variation Diminishing method.

    The component is restricted to regular grids (e.g., Raster or Hex).
    If multiple fields are advected, the advection__flux field will apply to
    the last one listed.

    Parameters
    ----------
    grid : RasterModelGrid or HexModelGrid
        A Landlab grid object.
    fields_to_advect : field name or list or (n_nodes,) array (default None)
        A node field of scalar values that will be advected, or list of fields.
        If not given, the component creates a generic field, initialized to zeros,
        called advected__quantity. If list >1 element given, advection will be
        applied to each field in it.
    advection_direction_is_steady : bool (default False)
        Indicates whether the directions of advection are expected to remain
        steady throughout a run. If True, some computation time is saved
        by calculating upwind links only once.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import AdvectionSolverTVD
    >>> grid = RasterModelGrid((3, 7))
    >>> s = grid.add_zeros("advected__quantity", at="node")
    >>> s[9:12] = np.array([1.0, 2.0, 1.0])
    >>> u = grid.add_zeros("advection__velocity", at="link")
    >>> u[grid.horizontal_links] = 1.0
    >>> advec = AdvectionSolverTVD(grid, fields_to_advect="advected__quantity")
    >>> for _ in range(5):
    ...     advec.update(0.2)
    ...
    >>> np.argmax(s[7:14])
    4
    """

    _name = "AdvectionSolverTVD"

    _unit_agnostic = True

    _info = {
        "advected__quantity": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Scalar quantity advected",
        },
        "advection__flux": {
            "dtype": float,
            "intent": "out",
            "optional": True,
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
        fields_to_advect=None,
        advection_direction_is_steady=False,
    ):
        """Initialize AdvectionSolverTVD."""

        # Call base class methods to check existence of input fields,
        # create output fields, etc.
        super().__init__(grid)
        self.initialize_output_fields()

        self._scalars = []  # list of fields to advect
        self._fluxes = []  # list of flux fields
        if fields_to_advect is None:
            try:
                self._scalars.append(self.grid.at_node["advected__quantity"])
            except KeyError:
                self._scalars.append(
                    self.grid.add_zeros("advected__quantity", at="node")
                )
            try:
                self._fluxes.append(self.grid.at_link["advection__flux"])
            except KeyError:
                self._fluxes.append(self.grid.add_zeros("advection__flux", at="link"))
        elif isinstance(fields_to_advect, list):
            flux_counter = 0
            for field in fields_to_advect:
                self._scalars.append(return_array_at_node(self.grid, field))
                if isinstance(field, str):
                    flux_name = "flux_of_" + field
                else:
                    flux_name = "advection__flux_" + str(flux_counter)
                    flux_counter += 1
                try:
                    flux = return_array_at_link(self.grid, flux_name)
                except FieldError:
                    flux = grid.add_zeros(flux_name, at="link")
                self._fluxes.append(flux)
        else:
            self._scalars.append(return_array_at_node(self.grid, fields_to_advect))
            if isinstance(fields_to_advect, str):
                flux_name = "flux_of_" + fields_to_advect
            else:
                flux_name = "advection__flux"
            try:
                flux = return_array_at_link(self.grid, flux_name)
            except FieldError:
                flux = grid.add_zeros(flux_name, at="link")
            self._fluxes.append(flux)

        self._vel = self.grid.at_link["advection__velocity"]

        self._advection_direction_is_steady = advection_direction_is_steady
        if advection_direction_is_steady:  # if so, only need to do this once
            self._upwind_link_at_link = find_upwind_link_at_link(self.grid, self._vel)
            self._upwind_link_at_link[
                self.grid.status_at_link == LinkStatus.INACTIVE
            ] = -1

    def calc_rate_of_change_at_nodes(self, scalar, flux, dt, update_upwind_links=False):
        """Calculate time rate of change in the advected quantity at nodes.

        Parameters
        ----------
        scalar : (n_nodes, ) array
            Scalar at-node field of values to be advected.
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        update_upwind_links : bool (optional; default False)
            If True, upwind links will be updated (set to True if the
            direction of advection is changing through time; if there are
            multiple advected quantities, it only needs to be True for the
            first one updated, and the update will be used for the others)
        """
        if update_upwind_links:
            self._upwind_link_at_link = find_upwind_link_at_link(self.grid, self._vel)
            self._upwind_link_at_link[
                self.grid.status_at_link == LinkStatus.INACTIVE
            ] = -1
        s_link_low = self.grid.map_node_to_link_linear_upwind(scalar, self._vel)
        s_link_high = self.grid.map_node_to_link_lax_wendroff(
            scalar, dt * self._vel / self.grid.length_of_link
        )
        r = upwind_to_local_grad_ratio(self.grid, scalar, self._upwind_link_at_link)
        psi = flux_lim_vanleer(r)
        s_at_link = psi * s_link_high + (1.0 - psi) * s_link_low
        flux[self.grid.active_links] = (
            self._vel[self.grid.active_links] * s_at_link[self.grid.active_links]
        )
        return -self.grid.calc_flux_div_at_node(flux)

    def update(self, dt):
        """Update the solution by one time step dt.

        Same as :meth:`~.run_one_step`.

        Parameters
        ----------
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        """
        update_upwinds = not self._advection_direction_is_steady
        for i in range(len(self._scalars)):  # update each of the advected scalars
            scalar = self._scalars[i]
            flux = self._fluxes[i]
            roc = self.calc_rate_of_change_at_nodes(
                scalar, flux, dt, update_upwind_links=update_upwinds
            )
            scalar[self.grid.core_nodes] += roc[self.grid.core_nodes] * dt
            update_upwinds = False  # should always be False after 1st scalar done

    def run_one_step(self, dt):
        """Update the solution by one time step dt.

        Same as :meth:`~.update`.

        Parameters
        ----------
        dt : float
            Time-step duration. Needed to calculate the Courant number.
        """
        self.update(dt)
