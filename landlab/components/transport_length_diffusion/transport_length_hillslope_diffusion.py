#!/usr/bin/env python3
"""Created on Tue Apr 11 10:13:38 2017.

@author: margauxmouchene
"""


import numpy as np

from landlab import Component


class TransportLengthHillslopeDiffuser(Component):
    r"""Transport length hillslope diffusion.

    Hillslope diffusion component in the style of Carretier et al. (2016,
    ESurf), and Davy and Lague (2009)

    .. math::

        \frac{dz}{dt} = -E + D (+ U)

        D = \frac{q_s}{L}

        E = k S

        L = \frac{dx}{(1 - (S / S_c)^2}

    Works on regular raster-type grid (RasterModelGrid, dx=dy).
    To be coupled with FlowDirectorSteepest for the calculation of steepest
    slope at each timestep.

    Component written by Margaux Mouchene, 2017

    Parameters
    ----------
    grid : ModelGrid
        Landlab ModelGrid object
    erodibility: float
        Erodibility coefficient [L/T]
    slope_crit: float (default=1.)
        Critical slope [L/L]

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorSteepest
    >>> from landlab.components import TransportLengthHillslopeDiffuser

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

    Instantiate Flow director (steepest slope type) and TL hillslope diffuser

    >>> fdir = FlowDirectorSteepest(mg)
    >>> tl_diff = TransportLengthHillslopeDiffuser(
    ...     mg, erodibility=0.001, slope_crit=0.6
    ... )

    Run the components for ten short timepsteps

    >>> for t in range(10):
    ...     fdir.run_one_step()
    ...     tl_diff.run_one_step(1.0)
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

    _name = "TransportLengthHillslopeDiffuser"

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
        "sediment__deposition_coeff": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Fraction of incoming sediment that is deposited on the node",
        },
        "sediment__deposition_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": "Deposition rate on node",
        },
        "sediment__erosion_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": "Erosion rate on node",
        },
        "sediment__flux_in": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": "Incoming sediment rate on node (=qs/dx)",
        },
        "sediment__flux_out": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": (
                "Outgoing sediment rate on node = sediment eroded on "
                "node + sediment transported across node from upstream"
            ),
        },
        "sediment__transfer_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": (
                "Rate of transferred sediment across a node (incoming "
                "sediment - deposited sediment on node)"
            ),
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/m",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
    }

    def __init__(self, grid, erodibility=0.001, slope_crit=1.0):
        """Initialize Diffuser.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        erodibility: float
            Erodibility coefficient [L/T]
        slope_crit: float (default=1.)
            Critical slope [L/L]
        """
        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that TransportLengthHillslopeDiffuser is compatible "
                "with route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        # Store grid and parameters

        self._k = erodibility
        self._slope_crit = slope_crit

        # Create fields:
        # Elevation
        self._elev = self._grid.at_node["topographic__elevation"]

        # Downstream steepest slope at node:
        self._steepest = self._grid.at_node["topographic__steepest_slope"]
        # On each node, node ID of downstream receiver node
        # (on node (i), ID of node that receives flow from node (i)):
        self._receiver = self._grid.at_node["flow__receiver_node"]

        self.initialize_output_fields()
        # Deposition
        self._depo = self._grid.at_node["sediment__deposition_rate"]

        # Transferred sediments (crossing over node)
        self._trans = self._grid.at_node["sediment__transfer_rate"]

        # Transport coefficient
        self._d_coeff = self._grid.at_node["sediment__deposition_coeff"]

        # Flux in
        self._flux_in = self._grid.at_node["sediment__flux_in"]

        # Flux out
        self._flux_out = self._grid.at_node["sediment__flux_out"]

        # Erosion
        self._erosion = self._grid.at_node["sediment__erosion_rate"]

    def tldiffusion(self, dt):
        """Calculate hillslope diffusion for a time period 'dt'.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        dt: float (time)
            The imposed timestep.
        """

        # Reset erosion, depo, trans and flux_in to 0
        self._erosion[:] = 0.0
        self._depo[:] = 0.0
        self._trans[:] = 0.0
        self._flux_in[:] = 0.0

        dx = self._grid.dx
        cores = self._grid.core_nodes

        # Calculate influx rate on node i  = outflux of nodes
        # whose receiver is i
        for i in self._grid.core_nodes:
            self._flux_in[self._receiver[i]] += self._flux_out[i]

            # Calculate transport coefficient
            # When S ~ Scrit, d_coeff is set to "infinity", for stability and
            # so that there is no deposition
            if self._steepest[i] >= self._slope_crit:
                self._d_coeff[i] = 1000000000.0
            else:
                self._d_coeff[i] = 1 / (
                    1 - (np.power(((self._steepest[i]) / self._slope_crit), 2))
                )

        # Calculate deposition rate on node
        self._depo[cores] = self._flux_in[cores] / self._d_coeff[cores]

        # Calculate erosion rate on node (positive value)
        # If S > Scrit, erosion is simply set for the slope to return to Scrit
        # Otherwise, erosion is slope times erodibility coefficent
        for i in self._grid.core_nodes:
            if self._steepest[i] > self._slope_crit:
                self._erosion[i] = (
                    dx * (self._steepest[i] - self._slope_crit) / (100 * dt)
                )
            else:
                self._erosion[i] = self._k * self._steepest[i]

            # Update elevation
            self._elev[i] += (-self._erosion[i] + self._depo[i]) * dt

        # Calculate transfer rate over node
        self._trans[cores] = self._flux_in[cores] - self._depo[cores]

        # Calculate outflux rate
        self._flux_out[:] = self._erosion + self._trans

    def run_one_step(self, dt):
        """Advance one timestep.

        Advance transport length-model hillslope diffusion component
        by one time step of size dt and tests for timestep stability.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        self.tldiffusion(dt)

        # Test code stability for timestep dt
        # Raise unstability error if local slope is reversed by erosion
        # and deposition during a timestep dt
        elev_dif = self._elev - self._elev[self._receiver]
        s = elev_dif[np.where(self._grid.at_node["flow__sink_flag"] == 0)]
        if np.any(s < -1) is True:
            raise ValueError(
                "The component is unstable" " for such a large timestep " "on this grid"
            )
        else:
            pass
