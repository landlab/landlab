#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Tue Feb 25, 2022.

@author: BenjaminCampforts

PLAN, GT & CLB May 2025:
- make this into two components, one for depth-dep and one for non-depth-dep, to mirror others
- deprecate the original TLD, which has unit/dimension errors and could benefit from up-to-downstream integration
- get input from BC re naming (CarrierDiffuser? TransportLength? NonLocal? EroDep?)
"""


import numpy as np

from landlab import Component

from landlab.grid.nodestatus import NodeStatus

from .cfuncs import non_local_depo


class NonlocalDiffuser(Component):
    r"""Non local diffusion.

    #TODO
    Correct for rho soil vs bedrock

        Hillslope diffusion component in the style of Carretier et al. (2016,
    ESurf), and Davy and Lague (2009)

    .. math::

        \frac{dz}{dt} = -E + D (+ U)

        D = \frac{q_s}{L}

        E = k S (1 - \exp(-H / H_*))

        L = \frac{dx}{(1 - (S / S_c)^2}

    Works on regular raster-type grid (RasterModelGrid, dx=dy).
    To be coupled with FlowDirectorSteepest for the calculation of steepest
    slope at each timestep.

    Component written by Benjamin Campforts, 2022 based on earlier
    non-depth-limited version by Margaux Mouchene. Modified 2025
    by Caroline Le Bouteiller, Greg Tucker, and Yuval Shmilovitz
    for *only* non-depth-dependent functionality 
    (see also DepthDependentNonLocalDiffuser).

    Parameters
    ----------

    grid : ModelGrid
        Landlab ModelGrid object
    erodibility: float
        Erodibility coefficient [L/T]
    slope_crit: float (default=1.)
        Critical slope [L/L]
    depositOnBoundaries: boolean (default=False)

    H_star=1.0 : float (default=1.)
        Decay scale for depth dependent transport rate (m)

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorSteepest
    >>> from landlab.components import NonLocalDiffuser

    Define grid and initial topography:

        - 3x5 grid
        - east and west boundaries are open, north and south are closed
        - Initial topography is plane at base level on the boundaries and
          1m of elevation elsewhere (core)

    >>> mg = RasterModelGrid((5, 5))
    >>> mg.set_closed_boundaries_at_grid_edges(False, True, False, True)
    >>> z = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 1.0, 1.0, 1.0, 0.0],
    ...     [0.0, 1.0, 1.0, 1.0, 0.0],
    ...     [0.0, 1.0, 1.0, 1.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> _ = mg.add_field("topographic__elevation", z, at="node")

    Instantiate Flow director (steepest slope type) and NL hillslope diffuser

    >>> fdir = FlowDirectorSteepest(mg)
    >>> mg.at_node.keys()
    >>> nl_diff = NonLocalDiffuser(
    ...     mg, erodibility=0.001, slope_crit=0.6
    ... )

    Run the components for ten short timepsteps

    >>> for t in range(10):
    ...     fdir.run_one_step()
    ...     nl_diff.run_one_step(1.0)
    ...

    Check final topography

    >>> np.allclose(
    ...     mg.at_node["topographic__elevation"].reshape(mg.shape),
    ...     [
    ...         [0.0, 0.0, 0.0, 0.0, 0.0],
    ...         [0.0, 0.96175283, 0.99982519, 0.96175283, 0.0],
    ...         [0.0, 0.96175283, 0.99982519, 0.96175283, 0.0],
    ...         [0.0, 0.96175283, 0.99982519, 0.96175283, 0.0],
    ...         [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     ],
    ... )
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Carretier, S., Martinod, P., Reich, M., Godderis, Y. (2016). Modelling
    sediment clasts transport during landscape evolution. Earth Surface Dynamics
    4(1), 237-251. https://dx.doi.org/10.5194/esurf-4-237-2016

    Davy, P., Lague, D. (2009). Fluvial erosion/transport equation of landscape
    evolution models revisited. Journal of Geophysical Research  114(F3),
    F03007. https://dx.doi.org/10.1029/2008jf001146

    """

    _name = "NonLocalDiffuser"

    _unit_agnostic = True

    _info = {
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/m",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "bedrock__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Bedrock elevation",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Sediment thickness",
        },
        "sediment_flux_out": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/y",
            "mapping": "node",
            "doc": "Sediment flux at boundary nodes in m3/y",
        },
    }

    def __init__(
        self,
        grid,
        erodibility=0.001,
        slope_crit=1.0,
        depositOnBoundaries=False,
        transportLengthCoefficient=None,
    ):
        """Initialize Diffuser.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        erodibility: float
            Erodibility coefficient [L/T]
        slope_crit: float (default=1.)
            Critical slope [L/L]
        depositOnBoundaries: boolean (default=False)

        transportLengthCoefficient [default = dx]

        """
        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that TransportLengthHillslopeDiffuser is compatible "
                "with route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        # Store grid and parameters

        self._k = erodibility
        self._slope_crit = slope_crit
        if transportLengthCoefficient is None:
            self._transportLengthCoefficient = grid.dx
        else:
            if transportLengthCoefficient < grid.dx:
                raise ValueError(
                    "The value for transportLengthCoefficient must be larger than the grid resolution"
                )
            else:
                self._transportLengthCoefficient = transportLengthCoefficient
        # Create fields:
        # Elevation
        self._el = self._grid.at_node["topographic__elevation"]

        self._steepest = self._grid.at_node["topographic__steepest_slope"]
        self._r = self._grid.at_node["flow__receiver_node"]
        self._lk_rcvr = self.grid.at_node["flow__link_to_receiver_node"]
        self._stack = self.grid.at_node["flow__upstream_node_order"]
        self._link_lengths = self.grid.length_of_d8

        self.initialize_output_fields()
        self._depositOnBoundaries = depositOnBoundaries

        self._flux = self._grid.at_node["sediment_flux_out"]

    def nldiffusion(self, dt):
        """Calculate hillslope diffusion for a time period 'dt'.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        dt: float (time)
            The imposed timestep.
        """

        dx = self._grid.dx
        lakes = self._steepest < 0

        # Calcualte erosion -- in comparison to v1, not curring off at Sc
        ero = self._k * self._steepest

        ero[lakes] = 0

        L = np.where(
            self._steepest < self._slope_crit,
            self._transportLengthCoefficient
            / (1 - (self._steepest / self._slope_crit) ** 2),
            1e9,
        )

        qs_out = np.zeros_like(self._el)
        depo = np.zeros_like(self._el)

        stack_flip_ud = np.flipud(self._stack)
        if self._depositOnBoundaries:
            stack_flip_ud_sel = stack_flip_ud
        else:
            node_status = self.grid.status_at_node
            stack_flip_ud_sel = stack_flip_ud[
                (node_status[stack_flip_ud] == NodeStatus.CORE)
            ]

        # C-code
        non_local_depo(dx, stack_flip_ud_sel, self._r, qs_out, L, ero, depo)
        # Non C-Code
        # for node in stack_flip_ud_sel:
        #     # L has to be larger than dx
        #     depo[node] = qs_out[node]/L[node]
        #     qs_out[self._r[node]] += qs_out[node]+ (ero[node]- depo[node])*dx

        # Update flux
        self._flux[:] = qs_out * dx  # in m3
        # Update elevation
        self._el += (-ero + depo) * dt

    def run_one_step(self, dt):
        """Advance one timestep.

        Advance transport length-model hillslope diffusion component
        by one time step of size dt and tests for timestep stability.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        elev_dif_before = self._el - self._el[self._r]
        flow__sink_flag = elev_dif_before < 0
        self.nldiffusion(dt)

        # Test code stability for timestep dt
        # Raise unstability error if local slope is reversed by erosion
        # and deposition during a timestep dt
        elev_dif = self._el - self._el[self._r]
        s = elev_dif[np.where(flow__sink_flag == 0)]
        if np.any(s < -1) is True:
            raise ValueError(
                "The component is unstable" " for such a large timestep " "on this grid"
            )
        else:
            pass
