r"""
stream_power_smooth_threshold.py: Defines the StreamPowerSmoothThresholdEroder,
which is a version of the FastscapeEroder (derived from it).

StreamPowerSmoothThresholdEroder uses a mathematically smooth threshold
formulation, rather than one with a singularity. The erosion rate is defined as
follows:

$\omega = K A^m S$

$E = \omega - \omega_c \left[ 1 - \exp ( -\omega / \omega_c ) \right]$

Created on Sat Nov 26 08:36:49 2016

@author: gtucker
"""

import numpy as np

from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import smooth_stream_power_eroder_solver
from .fastscape_stream_power import FastscapeEroder


class StreamPowerSmoothThresholdEroder(FastscapeEroder):
    """Stream erosion component with smooth threshold function.

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    K_sp : float, array, or field name
        K in the stream power equation (units vary with other parameters).
    m_sp : float, optional
        m in the stream power equation (power on drainage area).
    n_sp : float, optional, ~ 0.5<n_sp<4.
        n in the stream power equation (power on slope). NOTE: NOT PRESENTLY
        HONORED BY StreamPowerSmoothThresholdEroder (TODO)
    threshold_sp : float (TODO: array, or field name)
        The threshold stream power.
    discharge_field : float, field name, or array, optional
        Discharge [L^2/T]. The default is to use the grid field
        'drainage_area'. To use custom spatially/temporally varying
        rainfall, use 'water__unit_flux_in' to specify water input to the
        FlowAccumulator and use "surface_water__discharge" for this
        keyword argument.
    erode_flooded_nodes : bool (optional)
        Whether erosion occurs in flooded nodes identified by a
        depression/lake mapper (e.g., DepressionFinderAndRouter). When set
        to false, the field *flood_status_code* must be present on the grid
        (this is created by the DepressionFinderAndRouter). Default True.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid((3, 4))
    >>> rg.set_closed_boundaries_at_grid_edges(False, True, True, True)
    >>> z = rg.add_zeros("node", "topographic__elevation")
    >>> z[5] = 2.0
    >>> z[6] = 1.0
    >>> from landlab.components import FlowAccumulator
    >>> fr = FlowAccumulator(rg, flow_director="D4")
    >>> fr.run_one_step()
    >>> from landlab.components import StreamPowerSmoothThresholdEroder
    >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0)
    >>> sp.thresholds
    1.0
    >>> sp.run_one_step(dt=1.0)
    >>> import numpy as np
    >>> np.round(z[5:7], 3)
    array([1.646, 0.667])
    >>> z[5] = 2.0
    >>> z[6] = 1.0
    >>> import numpy as np
    >>> q = np.zeros(rg.number_of_nodes) + 0.25
    >>> q[6] = 100.0
    >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0, discharge_field=q)
    >>> sp.run_one_step(dt=1.0)
    >>> np.round(z[5:7], 3)
    array([1.754, 0.164])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Barnhart, K., Glade, R., Shobe, C., Tucker, G. (2019). Terrainbento 1.0: a
    Python package for multi-model analysis in long-term drainage basin
    evolution. Geoscientific Model Development  12(4), 1267--1297.
    https://dx.doi.org/10.5194/gmd-12-1267-2019

    **Additional References**

    Braun, J., Willett, S. (2013). A very efficient O(n), implicit and parallel
    method to solve the stream power equation governing fluvial incision and
    landscape evolution. Geomorphology  180-181(C), 170-179.
    https://dx.doi.org/10.1016/j.geomorph.2012.10.008

    """

    _name = "StreamPowerSmoothThresholdEroder"

    _unit_agnostic = True

    _cite_as = """
    @article{barnhart2019terrain,
      author = {Barnhart, Katherine R and Glade, Rachel C and Shobe, Charles M
                and Tucker, Gregory E},
      title = {{Terrainbento 1.0: a Python package for multi-model analysis in
                long-term drainage basin evolution}},
      doi = {10.5194/gmd-12-1267-2019},
      pages = {1267---1297},
      number = {4},
      volume = {12},
      journal = {Geoscientific Model Development},
      year = {2019},
    }
    """

    _info = {
        "drainage_area": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "topographic__elevation": {
            "dtype": float,
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
        K_sp=None,
        m_sp=0.5,
        n_sp=1.0,
        threshold_sp=1.0,
        discharge_field="drainage_area",
        erode_flooded_nodes=True,
    ):
        """Initialize StreamPowerSmoothThresholdEroder."""
        if "flow__receiver_node" in grid.at_node and grid.at_node[
            "flow__receiver_node"
        ].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that StreamPowerSmoothThresholdEroder is compatible "
                "with route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        if not erode_flooded_nodes and "flood_status_code" not in grid.at_node:
            raise ValueError(
                "In order to not erode flooded nodes another component "
                "must create the field *flood_status_code*. You want to "
                "run a lake mapper/depression finder."
            )

        if n_sp != 1.0:
            raise ValueError(
                "StreamPowerSmoothThresholdEroder currently only " "supports n_sp = 1"
            )

        # Call base-class init
        super().__init__(
            grid,
            K_sp=K_sp,
            m_sp=m_sp,
            n_sp=n_sp,
            threshold_sp=threshold_sp,
            discharge_field=discharge_field,
        )

        # Arrays with parameters for use in implicit solver
        self._gamma = grid.empty(at="node")
        self._delta = grid.empty(at="node")

    @property
    def alpha(self):
        """Erosion term divided by link length.

        Alpha is given as::

            alpha = K A^m dt / L

        where K is the erodibility, A is the drainage area, m is the
        drainage area exponent, dt is the timestep, and L is the link length.
        """
        return self._alpha

    @property
    def gamma(self):
        """Erosion threshold times timestep."""
        return self._gamma

    @property
    def thresholds(self):
        """Erosion thresholds."""
        return self._thresholds

    @property
    def delta(self):
        """Erosion term divided by link length and erosion threshold.

        delta is given as::

            delta = K A^m dt / (L * omega_c)

        where K is the erodibility, A is the drainage area, m is the
        drainage area exponent, dt is the timestep, L is the link length, and
        omega_c is the erosion threshold.
        """
        return self._delta

    def run_one_step(self, dt, runoff_rate=None):
        """Run one forward iteration of duration dt.

        Parameters
        ----------
        dt : float
            Time step size
        runoff_rate : (not used yet)
            (to be added later)

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((3, 3))
        >>> rg.set_closed_boundaries_at_grid_edges(False, True, True, True)
        >>> z = rg.add_zeros("node", "topographic__elevation")
        >>> z[4] = 1.0
        >>> from landlab.components import FlowAccumulator
        >>> fr = FlowAccumulator(rg, flow_director="D4")
        >>> fr.run_one_step()
        >>> from landlab.components import StreamPowerSmoothThresholdEroder
        >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0)
        >>> sp.run_one_step(dt=1.0)
        >>> sp.alpha
        array([0., 0., 0., 0., 1., 0., 0., 0., 0.])
        >>> sp.gamma
        array([0., 0., 0., 0., 1., 0., 0., 0., 0.])
        >>> sp.delta
        array([0., 0., 0., 0., 1., 0., 0., 0., 0.])
        """
        if not self._erode_flooded_nodes:
            flood_status = self._grid.at_node["flood_status_code"]
            flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]
        else:
            flooded_nodes = []

        # Set up needed arrays
        #
        # Get shorthand for elevation field ("z"), and for up-to-downstream
        # ordered array of node IDs ("upstream_order_IDs")
        node_id = np.arange(self._grid.number_of_nodes)
        upstream_order_IDs = self._grid["node"]["flow__upstream_node_order"]
        z = self._grid["node"]["topographic__elevation"]
        flow_receivers = self._grid["node"]["flow__receiver_node"]

        # Get an array of flow-link length for each node that has a defined
        # receiver (i.e., that drains to another node).
        defined_flow_receivers = np.not_equal(
            self._grid["node"]["flow__link_to_receiver_node"], self._grid.BAD_INDEX
        )
        defined_flow_receivers[flow_receivers == node_id] = False

        if flooded_nodes is not None:
            defined_flow_receivers[flooded_nodes] = False
        flow_link_lengths = self._grid.length_of_d8[
            self._grid.at_node["flow__link_to_receiver_node"][defined_flow_receivers]
        ]

        # Set up alpha, beta, delta arrays
        #
        #   First, compute drainage area or discharge raised to the power m.
        np.power(self._A, self._m, out=self._A_to_the_m)

        #   Alpha
        self._alpha[~defined_flow_receivers] = 0.0
        self._alpha[defined_flow_receivers] = (
            self._K[defined_flow_receivers]
            * dt
            * self._A_to_the_m[defined_flow_receivers]
            / flow_link_lengths
        )

        #   Gamma
        self._gamma[~defined_flow_receivers] = 0.0

        #   Delta
        self._delta[~defined_flow_receivers] = 0.0

        if isinstance(self._thresholds, np.ndarray):
            self._gamma[defined_flow_receivers] = (
                dt * self._thresholds[defined_flow_receivers]
            )

            self._delta[defined_flow_receivers] = (
                self._K[defined_flow_receivers]
                * self._A_to_the_m[defined_flow_receivers]
            ) / (self._thresholds[defined_flow_receivers] * flow_link_lengths)

            self._delta[defined_flow_receivers][
                self._thresholds[defined_flow_receivers] == 0.0
            ] = 0.0
        else:
            self._gamma[defined_flow_receivers] = dt * self._thresholds
            if self._thresholds == 0:
                self._delta[defined_flow_receivers] = 0.0
            else:
                self._delta[defined_flow_receivers] = (
                    self._K[defined_flow_receivers]
                    * self._A_to_the_m[defined_flow_receivers]
                ) / (self._thresholds * flow_link_lengths)

        # Iterate over nodes from downstream to upstream, using scipy's
        # 'newton' function to find new elevation at each node in turn.
        smooth_stream_power_eroder_solver(
            upstream_order_IDs, flow_receivers, z, self._alpha, self._gamma, self._delta
        )
