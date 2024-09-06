"""This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router. Credit: Voller, Hobley, Paola.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

# ##in the diagonal case, even closed edges can produce "drag". Is this right?
# Could suppress by mirroring the diagonals

import numpy as np

from landlab import Component
from landlab import LinkStatus
from landlab import RasterModelGrid


class PotentialityFlowRouter(Component):
    """Multidirectional flow routing using a novel method.

    This class implements Voller, Hobley, and Paola's experimental matrix
    solutions for flow routing. The method works by solving for a potential
    field at all nodes on the grid, which enforces both mass conservation
    and flow downhill along topographic gradients. It is order n and highly
    efficient, but does not return any information about flow connectivity.

    Options are permitted to allow "abstract" routing (flow enforced downslope,
    but no particular assumptions are made about the governing equations), or
    routing according to the Chezy or Manning equations. This routine assumes
    that water is distributed evenly over the surface of the cell in deriving
    the depth, and does not assume channelization. You will need to back-
    calculate channel depths for yourself using known widths at each node
    if that is what you want.

    It is VITAL you initialize this component AFTER setting boundary
    conditions.

    If Manning or Chezy specified, the surface_water__depth is the depth of
    flow in the cell, calculated assuming flow occurs over the whole surface.

    Note that this component offers the property `discharges_at_links`. This
    returns the discharges at all links. If method=='D8', this list will
    include diagonal links after the orthogonal links, which is why this
    information is not returned as a field.

    Discharges at nodes are recorded as the outgoing total discharge (i.e.,
    including any contribution from 'water__unit_flux_in').

    The primary method of this class is :func:`run_one_step`.

    Notes
    -----
    This is a "research grade" component, and is subject to dramatic change
    with little warning. No guarantees are made regarding its accuracy or
    utility. It is not recommended for user use yet!

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> import numpy as np
    >>> mg = HexModelGrid(
    ...     (4, 6), spacing=2.0, node_layout="rect", orientation="vertical"
    ... )
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> Q_in = mg.add_ones("water__unit_flux_in", at="node")
    >>> z += mg.node_y.copy()
    >>> potfr = PotentialityFlowRouter(mg)
    >>> potfr.run_one_step()
    >>> mg.at_node["surface_water__discharge"][mg.core_nodes]
    array([11.70706863, 11.5709712 , 10.41329927,  9.24959728,
            6.65448576,  6.39262702,  5.71410162,  5.04743495])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed

    """

    _name = "PotentialityFlowRouter"

    _unit_agnostic = False

    _info = {
        "flow__potential": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": (
                "Value of the hypothetical field 'K', used to force water "
                "flux to flow downhill"
            ),
        },
        "surface_water__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water__unit_flux_in": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": (
                "External volume water per area per time input to each node "
                "(e.g., rainfall rate)"
            ),
        },
    }

    _min_slope_thresh = 1.0e-24
    # if your flow isn't connecting up, this probably needs to be reduced

    def __init__(
        self, grid, method="D8", flow_equation="default", Chezys_C=30.0, Mannings_n=0.03
    ):
        """
        Parameters
        ----------
        grid : ModelGrid
            A grid.
        method : {'D8', 'D4'}, optional
            Routing method ('D8' is the default). This keyword has no effect for a
            Voronoi-based grid.
        flow_equation : {'default', 'Manning', 'Chezy'}, optional
            If Manning or Chezy, flow is routed according to the Manning or Chezy
            equation; discharge is allocated to multiple downslope nodes
            proportional to the square root of discharge; and a water__depth field
            is returned. If default, flow is allocated to multiple nodes linearly
            with slope; and the water__depth field is not calculated.
        Chezys_C : float (optional)
            Required if flow_equation == 'Chezy'.
        Mannings_n : float (optional)
            Required if flow_equation == 'Manning'.
        """
        super().__init__(grid)

        if isinstance(grid, RasterModelGrid):
            assert grid.number_of_node_rows >= 3
            assert grid.number_of_node_columns >= 3
            self._raster = True
        else:
            self._raster = False

        self._equation = flow_equation
        assert self._equation in ("default", "Chezy", "Manning")
        if self._equation == "Chezy":
            self._chezy_C = Chezys_C
        elif self._equation == "Manning":
            self._manning_n = Mannings_n
        assert method in ("D8", "D4")
        if method == "D8":
            self._route_on_diagonals = True
        else:
            self._route_on_diagonals = False

        self.initialize_output_fields()

        if self._raster:
            self._equiv_circ_diam = 2.0 * np.sqrt(grid.dx * grid.dy / np.pi)
        else:
            for_cell_areas = 2.0 * np.sqrt(grid.area_of_cell / np.pi)
            mean_A = for_cell_areas.mean()
            self._equiv_circ_diam = for_cell_areas[grid.cell_at_node]
            self._equiv_circ_diam[grid.cell_at_node == -1] = mean_A
        # ^this is the equivalent seen CSWidth of a cell for a flow in a
        # generic 360 direction
        if self._route_on_diagonals and self._raster:
            self._discharges_at_link = np.empty(grid.number_of_d8)
        else:
            self._discharges_at_link = self._grid.empty("link")

    def run_one_step(self):
        """Route surface-water flow over a landscape.

        Both convergent and divergent flow can occur.
        """
        grid = self._grid
        self._K = grid.at_node["flow__potential"]
        self._Qw = grid.at_node["surface_water__discharge"]
        z = grid.at_node["topographic__elevation"]
        qwater_in = grid.at_node["water__unit_flux_in"].copy()
        qwater_in[grid.node_at_cell] *= grid.area_of_cell
        prev_K = self._K.copy()
        mismatch = 10000.0
        # do the ortho nodes first, in isolation
        g = grid.calc_grad_at_link(z)
        if self._equation != "default":
            g = np.sign(g) * np.sqrt(np.fabs(g))
            # ^...because both Manning and Chezy actually follow sqrt
            # slope, not slope
        # weight by face width - NO, because diags
        # g *= grid.length_of_face[grid.face_at_link]
        link_grad_at_node_w_dir = g[grid.links_at_node] * grid.active_link_dirs_at_node
        # active_link_dirs strips "wrong" face widths

        # now outgoing link grad sum
        outgoing_sum = (
            np.sum((link_grad_at_node_w_dir).clip(0.0), axis=1) + self._min_slope_thresh
        )
        pos_incoming_link_grads = (-link_grad_at_node_w_dir).clip(0.0)

        if not self._route_on_diagonals or not self._raster:
            while mismatch > 1.0e-6:
                K_link_ends = self._K[grid.adjacent_nodes_at_node]
                incoming_K_sum = (pos_incoming_link_grads * K_link_ends).sum(
                    axis=1
                ) + self._min_slope_thresh
                self._K[:] = (incoming_K_sum + qwater_in) / outgoing_sum
                mismatch = np.sum(np.square(self._K - prev_K))
                prev_K = self._K.copy()

            upwind_K = grid.map_value_at_max_node_to_link(z, self._K)
            self._discharges_at_link[:] = upwind_K * g
            self._discharges_at_link[grid.status_at_link == LinkStatus.INACTIVE] = 0.0
        else:
            # grad on diags:
            gwd = np.empty(grid.number_of_d8, dtype=float)
            gd = gwd[grid.number_of_links :]
            gd[:] = z[grid.nodes_at_diagonal[:, 1]] - z[grid.nodes_at_diagonal[:, 0]]
            gd /= grid.length_of_d8[grid.number_of_links :]
            if self._equation != "default":
                gd[:] = np.sign(gd) * np.sqrt(np.fabs(gd))
            diag_grad_at_node_w_dir = (
                gwd[grid.d8s_at_node[:, 4:]] * self._grid.active_diagonal_dirs_at_node
            )

            outgoing_sum += np.sum(diag_grad_at_node_w_dir.clip(0.0), axis=1)
            pos_incoming_diag_grads = (-diag_grad_at_node_w_dir).clip(0.0)
            while mismatch > 1.0e-6:
                K_link_ends = self._K[grid.adjacent_nodes_at_node]
                K_diag_ends = self._K[grid.diagonal_adjacent_nodes_at_node]
                incoming_K_sum = (
                    (pos_incoming_link_grads * K_link_ends).sum(axis=1)
                    + (pos_incoming_diag_grads * K_diag_ends).sum(axis=1)
                    + self._min_slope_thresh
                )
                self._K[:] = (incoming_K_sum + qwater_in) / outgoing_sum
                mismatch = np.sum(np.square(self._K - prev_K))
                prev_K = self._K.copy()

            # ^this is necessary to suppress stupid apparent link Qs at flow
            # edges, if present.
            upwind_K = grid.map_value_at_max_node_to_link(z, self._K)
            upwind_diag_K = np.where(
                z[grid.nodes_at_diagonal[:, 1]] > z[grid.nodes_at_diagonal[:, 0]],
                self._K[grid.nodes_at_diagonal[:, 1]],
                self._K[grid.nodes_at_diagonal[:, 0]],
            )
            self._discharges_at_link[: grid.number_of_links] = upwind_K * g
            self._discharges_at_link[grid.number_of_links :] = upwind_diag_K * gd
            self._discharges_at_link[grid.status_at_d8 == LinkStatus.FIXED] = 0.0

        np.multiply(self._K, outgoing_sum, out=self._Qw)
        # there is no sensible way to save discharges at links, if we route
        # on diagonals.
        # for now, let's make a property

        # now process uval and vval to give the depths, if Chezy or Manning:
        if self._equation == "Chezy":
            # Chezy: Q = C*Area*sqrt(depth*slope)
            grid.at_node["surface_water__depth"][:] = (
                grid.at_node["flow__potential"] / self._chezy_C / self._equiv_circ_diam
            ) ** (2.0 / 3.0)
        elif self._equation == "Manning":
            # Manning: Q = w/n*depth**(5/3)
            grid.at_node["surface_water__depth"][:] = (
                grid.at_node["flow__potential"]
                * self._manning_n
                / self._equiv_circ_diam
            ) ** 0.6
        else:
            pass

    @property
    def discharges_at_links(self):
        """Return the discharges at links.

        Note that if diagonal routing, this will return number_of_d8.
        Otherwise, it will be number_of_links.
        """
        return self._discharges_at_link
