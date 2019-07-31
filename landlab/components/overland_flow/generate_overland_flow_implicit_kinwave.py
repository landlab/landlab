# -*- coding: utf-8 -*-
"""
Landlab component for overland flow using a local implicit solution to the
kinematic-wave approximation.

Created on Fri May 27 14:26:13 2016

@author: gtucker
"""


import numpy as np
from scipy.optimize import newton

from landlab import Component
from landlab.components import FlowAccumulator


def water_fn(x, a, b, c, d, e):
    r"""Evaluates the solution to the water-depth equation.

    Called by scipy.newton() to find solution for :math:`x`
    using Newton's method.

    Parameters
    ----------
    x : float
        Water depth at new time step.
    a : float
        "alpha" parameter (see below)
    b : float
        Weighting factor on new versus old time step. :math:`b=1` means purely
        implicit solution with all weight on :math:`H` at new time
        step. :math:`b=0` (not recommended) would mean purely explicit.
    c : float
        Water depth at old time step (time step :math:`t` instead
        of :math:`t+1`)
    d : float
        Depth-discharge exponent; normally either 5/3 (Manning) or 3/2 (Chezy)
    e : float
        Water inflow volume per unit cell area in one time step.

    This equation represents the implicit solution for water depth
    :math:`H` at the next time step. In the code below, it is
    formulated in a generic way.  Written using more familiar
    terminology, the equation is:

    .. math::

        H - H_0 + \alpha ( w H + (w-1) H_0)^d - \Delta t (R + Q_{in} / A)

    .. math::

        \alpha = \frac{\Delta t \sum S^{1/2}}{C_f A}

    where :math:`H` is water depth at the given node at the new
    time step, :math:`H_0` is water depth at the prior time step,
    :math:`w` is a weighting factor, :math:`d` is the depth-discharge
    exponent (2/3 or 1/2), :math:`\Delta t` is time-step duration,
    :math:`R` is local runoff rate, :math:`Q_{in}` is inflow
    discharge, :math:`A` is cell area, :math:`C_f` is a
    dimensional roughness coefficient, and :math:`\sum S^{1/2}`
    represents the sum of square-root-of-downhill-gradient over
    all outgoing (downhill) links.
    """
    return x - c + a * (b * x + (b - 1.0) * c) ** d - e


class KinwaveImplicitOverlandFlow(Component):
    r"""Calculate shallow water flow over topography.

    Landlab component that implements a two-dimensional kinematic wave model.
    This is a form of the 2D shallow-water equations in which energy slope is
    assumed to equal bed slope. The solution method is locally implicit, and
    works as follows. At each time step, we iterate from upstream to downstream
    over the topography. Because we are working downstream, we can assume that
    we know the total water inflow to a given cell. We solve the following mass
    conservation equation at each cell:

    .. math::

        (H^{t+1} - H^t)/\Delta t = Q_{in}/A - Q_{out}/A + R

    where :math:`H` is water depth, :math:`t` indicates time step
    number, :math:`\Delta t` is time step duration, :math:`Q_{in}` is
    total inflow discharge, :math:`Q_{out}` is total outflow
    discharge, :math:`A` is cell area, and :math:`R` is local
    runoff rate (precipitation minus infiltration; could be
    negative if runon infiltration is occurring).

    The specific outflow discharge leaving a cell along one of its faces is:

    .. math::

        q = (1/C_r) H^\alpha S^{1/2}

    where :math:`C_r` is a roughness coefficient (such as
    Manning's n), :math:`\alpha` is an exponent equal to :math:`5/3`
    for the Manning equation and :math:`3/2` for the Chezy family,
    and :math:`S` is the downhill-positive gradient of the link
    that crosses this particular face. Outflow discharge is zero
    for links that are flat or "uphill" from the given node.
    Total discharge out of a cell is then the sum of (specific
    discharge x face width) over all outflow faces

    .. math::

        Q_{out} = \sum_{i=1}^N (1/C_r) H^\alpha S_i^{1/2} W_i

    where :math:`N` is the number of outflow faces (i.e., faces
    where the ground slopes downhill away from the cell's node),
    and :math:`W_i` is the width of face :math:`i`.

    We use the depth at the cell's node, so this simplifies to:

    .. math::

        Q_{out} = (1/C_r) H'^\alpha \sum_{i=1}^N S_i^{1/2} W_i

    We define :math:`H` in the above as a weighted sum of
    the "old" (time step :math:`t`) and "new" (time step :math:`t+1`)
    depth values:

    .. math::

        H' = w H^{t+1} + (1-w) H^t

    If :math:`w=1`, the method is fully implicit. If :math:`w=0`,
    it is a simple forward explicit method.

    When we combine these equations, we have an equation that includes the
    unknown :math:`H^{t+1}` and a bunch of terms that are known.
    If :math:`w\ne 0`, it is a nonlinear equation in :math:`H^{t+1}`,
    and must be solved iteratively. We do this using a root-finding
    method in the scipy.optimize library.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid((4, 5), xy_spacing=10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> kw = KinwaveImplicitOverlandFlow(rg)
    >>> round(kw.runoff_rate * 1.0e7, 2)
    2.78
    >>> kw.vel_coef  # default value
    100.0
    >>> rg.at_node['surface_water__depth'][6:9]
    array([ 0.,  0.,  0.])
    """

    _name = "KinwaveImplicitOverlandFlow"

    _input_var_names = ("topographic__elevation",)

    _output_var_names = (
        "topographic__gradient",
        "surface_water__depth",
        # 'water__velocity',
        # 'water__specific_discharge',
        "surface_water_inflow__discharge",
    )

    _var_units = {
        "topographic__elevation": "m",
        "topographic__slope": "m/m",
        "surface_water__depth": "m",
        "water__velocity": "m/s",
        "water__specific_discharge": "m2/s",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "topographic__gradient": "link",
        "surface_water__depth": "node",
        # 'water__velocity': 'link',
        # 'water__specific_discharge': 'link',
        "surface_water_inflow__discharge": "node",
    }

    _var_doc = {
        "topographic__elevation": "elevation of the ground surface relative to some datum",
        "topographic__gradient": "gradient of the ground surface",
        "surface_water__depth": "depth of water",
        #        'water__velocity':
        #            'flow velocity component in the direction of the link',
        #        'water__specific_discharge':
        #            'flow discharge component in the direction of the link',
        "surface_water_inflow__discharge": "water volume inflow rate to the cell around each node",
    }

    def __init__(
        self,
        grid,
        runoff_rate=1.0,
        roughness=0.01,
        changing_topo=False,
        depth_exp=1.5,
        weight=1.0,
        **kwds
    ):
        """Initialize the KinwaveOverlandFlowModel.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        runoff_rate : float, optional (defaults to 1 mm/hr)
            Precipitation rate, mm/hr
        roughnes : float, defaults to 0.01
            Manning roughness coefficient, s/m^1/3
        changing_topo : boolean, optional (defaults to False)
            Flag indicating whether topography changes between time steps
        depth_exp : float (defaults to 1.5)
            Exponent on water depth in velocity equation (3/2 for Darcy/Chezy,
            5/3 for Manning)
        weight : float (defaults to 1.0)
            Weighting on depth at new time step versus old time step (1 = all
            implicit; 0 = explicit)
        """

        # Store grid and parameters and do unit conversion
        self._grid = grid
        self.runoff_rate = runoff_rate / 3600000.0  # convert to m/s
        self.vel_coef = 1.0 / roughness  # do division now to save time
        self.changing_topo = changing_topo
        self.depth_exp = depth_exp
        self.weight = weight

        # Get elevation field
        self.elev = grid.at_node["topographic__elevation"]

        # Create fields...
        #  Water depth
        if "surface_water__depth" in grid.at_node:
            self.depth = grid.at_node["surface_water__depth"]
        else:
            self.depth = grid.add_zeros("node", "surface_water__depth")
        #   Slope
        if "topographic__gradient" in grid.at_link:
            self.slope = grid.at_link["topographic__gradient"]
        else:
            self.slope = grid.add_zeros("link", "topographic__gradient")
        #  Velocity
        #        if 'water__velocity' in grid.at_link:
        #            self.vel = grid.at_link['water__velocity']
        #        else:
        #            self.vel = grid.add_zeros('link', 'water__velocity')
        #  Discharge
        #        if 'surface_water__specific_discharge' in grid.at_link:
        #            self.disch = grid.at_link['surface_water__specific_discharge']
        #        else:
        #            self.disch = grid.add_zeros('link',
        #                                        'surface_water__specific_discharge')
        #  Inflow discharge at nodes
        if "surface_water_inflow__discharge" in grid.at_node:
            self.disch_in = grid.at_node["surface_water_inflow__discharge"]
        else:
            self.disch_in = grid.add_zeros("node", "surface_water_inflow__discharge")

        # This array holds, for each node, the sum of sqrt(slope) x face width
        # for each link/face.
        self.grad_width_sum = grid.zeros("node")

        # This array holds the prefactor in the algebraic equation that we
        # will find a solution for.
        self.alpha = grid.zeros("node")

        # Instantiate flow router
        self.flow_accum = FlowAccumulator(
            grid,
            "topographic__elevation",
            flow_director="MFD",
            partition_method="square_root_of_slope",
        )

        # Flag to let us know whether this is our first iteration
        self.first_iteration = True

    def run_one_step(self, dt, current_time=0.0, runoff_rate=None, **kwds):
        """Calculate water flow for a time period `dt`.
        """
        # Handle runoff rate
        if runoff_rate is None:
            runoff_rate = self.runoff_rate

        # If it's our first iteration, or if the topography may be changing,
        # do flow routing and calculate square root of slopes at links
        if self.changing_topo or self.first_iteration:

            # Calculate the ground-surface slope
            self.slope[self.grid.active_links] = self._grid.calc_grad_at_link(
                self.elev
            )[self._grid.active_links]

            # Take square root of slope magnitude for use in velocity eqn
            self.sqrt_slope = np.sqrt(np.abs(self.slope))

            # Re-route flow, which gives us the downstream-to-upstream
            # ordering
            self.flow_accum.run_one_step()
            self.nodes_ordered = self.grid.at_node["flow__upstream_node_order"]
            self.flow_lnks = self.grid.at_node["flow__link_to_receiver_node"]

            # (Re)calculate, for each node, sum of sqrt(gradient) x width
            self.grad_width_sum[:] = 0.0
            for i in range(self.flow_lnks.shape[1]):
                self.grad_width_sum[:] += (
                    self.sqrt_slope[self.flow_lnks[:, i]]
                    * self._grid.width_of_face[
                        self.grid.face_at_link[self.flow_lnks[:, i]]
                    ]
                )

            # Calculate values of alpha, which is defined as
            #
            #   $\alpha = \frac{\Sigma W S^{1/2} \Delta t}{A C_r}$
            cores = self.grid.core_nodes
            self.alpha[cores] = (
                self.vel_coef
                * self.grad_width_sum[cores]
                * dt
                / (self.grid.area_of_cell[self.grid.cell_at_node[cores]])
            )

        # Zero out inflow discharge
        self.disch_in[:] = 0.0

        # Upstream-to-downstream loop
        for i in range(len(self.nodes_ordered) - 1, -1, -1):
            n = self.nodes_ordered[i]
            if self.grid.status_at_node[n] == 0:

                # Solve for new water depth
                aa = self.alpha[n]
                cc = self.depth[n]
                ee = (dt * runoff_rate) + (
                    dt
                    * self.disch_in[n]
                    / self.grid.area_of_cell[self.grid.cell_at_node[n]]
                )
                self.depth[n] = newton(
                    water_fn,
                    self.depth[n],
                    args=(aa, self.weight, cc, self.depth_exp, ee),
                )

                # Calc outflow
                Heff = self.weight * self.depth[n] + (1.0 - self.weight) * cc
                outflow = (
                    self.vel_coef * (Heff ** self.depth_exp) * self.grad_width_sum[n]
                )  # this is manning/chezy/darcy

                # Send flow downstream. Here we take total inflow discharge
                # and partition it among the node's neighbors. For this, we use
                # the flow director's "proportions" array, which contains, for
                # each node, the proportion of flow that heads out toward each
                # of its N neighbors. The proportion is zero if the neighbor is
                # uphill; otherwise, it is S^1/2 / sum(S^1/2). If for example
                # we have a raster grid, there will be four neighbors and four
                # proportions, some of which may be zero and some between 0 and
                # 1.
                self.disch_in[self.grid.adjacent_nodes_at_node[n]] += (
                    outflow * self.flow_accum.flow_director.proportions[n]
                )

                # TODO: the above is enough to implement the solution for flow
                # depth, but it does not provide any information about flow
                # velocity or discharge on links. This could be added as an
                # optional method, perhaps done just before output.


if __name__ == "__main__":
    import doctest

    doctest.testmod()
