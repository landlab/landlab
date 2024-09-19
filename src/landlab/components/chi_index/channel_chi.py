"""Created March 2016.

@author: dejh
"""

import numpy as np

from landlab import Component
from landlab import RasterModelGrid

try:
    from itertools import izip
except ImportError:
    izip = zip


class ChiFinder(Component):
    """Calculate Chi Indices.

    This component calculates chi indices, sensu Perron & Royden, 2013,
    for a Landlab landscape.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, FastscapeEroder
    >>> from landlab.components import ChiFinder

    >>> mg = RasterModelGrid((3, 4))
    >>> for nodes in (
    ...     mg.nodes_at_right_edge,
    ...     mg.nodes_at_bottom_edge,
    ...     mg.nodes_at_top_edge,
    ... ):
    ...     mg.status_at_node[nodes] = mg.BC_NODE_IS_CLOSED
    >>> _ = mg.add_field("topographic__elevation", mg.node_x, at="node")
    >>> fr = FlowAccumulator(mg, flow_director="D8")
    >>> cf = ChiFinder(mg, min_drainage_area=1.0, reference_concavity=1.0)
    >>> fr.run_one_step()
    >>> cf.calculate_chi()
    >>> mg.at_node["channel__chi_index"].reshape(mg.shape)[1, :]
    array([0.5, 1. , 2. , 0. ])

    >>> mg2 = RasterModelGrid((5, 5), xy_spacing=100.0)
    >>> for nodes in (
    ...     mg2.nodes_at_right_edge,
    ...     mg2.nodes_at_bottom_edge,
    ...     mg2.nodes_at_top_edge,
    ... ):
    ...     mg2.status_at_node[nodes] = mg2.BC_NODE_IS_CLOSED
    >>> _ = mg2.add_zeros("node", "topographic__elevation")
    >>> mg2.at_node["topographic__elevation"][mg2.core_nodes] = (
    ...     mg2.node_x[mg2.core_nodes] / 1000.0
    ... )
    >>> np.random.seed(0)
    >>> mg2.at_node["topographic__elevation"][mg2.core_nodes] += np.random.rand(
    ...     mg2.number_of_core_nodes
    ... )
    >>> fr2 = FlowAccumulator(mg2, flow_director="D8")
    >>> sp2 = FastscapeEroder(mg2, K_sp=0.01)
    >>> cf2 = ChiFinder(mg2, min_drainage_area=0.0, reference_concavity=0.5)
    >>> for i in range(10):
    ...     mg2.at_node["topographic__elevation"][mg2.core_nodes] += 10.0
    ...     fr2.run_one_step()
    ...     sp2.run_one_step(1000.0)
    ...
    >>> fr2.run_one_step()
    >>> cf2.calculate_chi()
    >>> mg2.at_node["channel__chi_index"].reshape(mg2.shape)
    array([[0.        , 0.        , 0.        , 0.        , 0.        ],
           [0.77219416, 1.54438833, 2.63643578, 2.61419437, 0.        ],
           [1.09204746, 2.18409492, 1.52214691, 2.61419437, 0.        ],
           [0.44582651, 0.89165302, 1.66384718, 2.75589464, 0.        ],
           [0.        , 0.        , 0.        , 0.        , 0.        ]])

    >>> cf3 = ChiFinder(
    ...     mg2,
    ...     min_drainage_area=20000.0,
    ...     use_true_dx=True,
    ...     reference_concavity=0.5,
    ...     reference_area=mg2.at_node["drainage_area"].max(),
    ...     clobber=True,
    ... )
    >>> cf3.calculate_chi()
    >>> cf3.chi_indices.reshape(mg2.shape)
    array([[  0. ,   0.        ,   0.        ,   0. ,  0. ],
           [  0. , 173.20508076,   0.        ,   0. ,  0. ],
           [  0. ,   0.        , 270.71067812,   0. ,  0. ],
           [  0. , 100.        , 236.60254038,   0. ,  0. ],
           [  0. ,   0.        ,   0.        ,   0. ,  0. ]])
    >>> cf3.hillslope_mask.reshape(mg2.shape)
    array([[ True,  True,  True,  True,  True],
           [False, False,  True,  True,  True],
           [ True,  True, False,  True,  True],
           [False, False, False,  True,  True],
           [ True,  True,  True,  True,  True]])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Perron, J., Royden, L. (2012). An integral approach to bedrock river
    profile analysis Earth Surface Processes and Landforms  38(6), 570-576.
    https://dx.doi.org/10.1002/esp.3302

    """

    _name = "ChiFinder"

    _unit_agnostic = True

    _info = {
        "channel__chi_index": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "variable",
            "mapping": "node",
            "doc": "the local steepness index",
        },
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
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
    }

    def __init__(
        self,
        grid,
        reference_concavity=0.5,
        min_drainage_area=1.0e6,
        reference_area=1.0,
        use_true_dx=False,
        clobber=False,
    ):
        """
        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        reference_concavity : float
            The reference concavity to use in the calculation.
        min_drainage_area : float (m**2)
            The drainage area down to which to calculate chi.
        reference_area : float or None (m**2)
            If None, will default to the mean core cell area on the grid.
            Else, provide a value to use. Essentially becomes a prefactor on the
            value of chi.
        use_true_dx : bool (default False)
            If True, integration to give chi is performed using each value of node
            spacing along the channel (which can lead to a quantization effect,
            and is not preferred by Taylor & Royden). If False, the mean value of
            node spacing along the all channels is assumed everywhere.
        clobber : bool (default False)
            Raise an exception if adding an already existing field.

        """
        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that ChiFinder is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        if isinstance(self._grid, RasterModelGrid):
            self._link_lengths = self._grid.length_of_d8
        else:
            self._link_lengths = self._grid.length_of_link  # not tested

        self._reftheta = reference_concavity
        self._min_drainage = min_drainage_area

        self._set_up_reference_area(reference_area)

        self._use_true_dx = use_true_dx
        self._chi = self._grid.add_zeros("node", "channel__chi_index", clobber=clobber)
        self._mask = self._grid.ones("node", dtype=bool)
        self._elev = self._grid.at_node["topographic__elevation"]

    def _set_up_reference_area(self, reference_area):
        """Set up and validate reference_area."""
        if reference_area <= 0.0:
            raise ValueError(
                "ChiFinder: reference_area must be positive."
            )  # not tested
        self._A0 = reference_area

    def calculate_chi(self):
        """Calculate local chi indices.

        This is the main method. Call it to calculate local chi indices
        at all points with drainage areas greater than `min_drainage_area`.

        Chi of any node without a defined value is reported as 0. These nodes
        are also identified in the mask retrieved with :func:`hillslope_mask`.
        """
        self._mask.fill(True)
        self._chi.fill(0.0)

        reftheta = self._reftheta
        min_drainage = self._min_drainage
        reference_area = self._A0
        self._set_up_reference_area(reference_area)

        use_true_dx = self._use_true_dx

        upstr_order = self._grid.at_node["flow__upstream_node_order"]
        # get an array of only nodes with A above threshold:
        valid_upstr_order = upstr_order[
            self._grid.at_node["drainage_area"][upstr_order] >= min_drainage
        ]
        valid_upstr_areas = self._grid.at_node["drainage_area"][valid_upstr_order]
        if not use_true_dx:
            chi_integrand = (self._A0 / valid_upstr_areas) ** reftheta
            mean_dx = self.mean_channel_node_spacing(valid_upstr_order)
            self.integrate_chi_avg_dx(
                valid_upstr_order, chi_integrand, self._chi, mean_dx
            )
        else:
            chi_integrand = self._grid.zeros("node")
            chi_integrand[valid_upstr_order] = (
                self._A0 / valid_upstr_areas
            ) ** reftheta
            self.integrate_chi_each_dx(valid_upstr_order, chi_integrand, self._chi)
        # stamp over the closed nodes, as it's possible they can receive infs
        # if min_drainage_area < grid.cell_area_at_node
        self._chi[self._grid.status_at_node == self._grid.BC_NODE_IS_CLOSED] = 0.0
        self._mask[valid_upstr_order] = False

    def integrate_chi_avg_dx(
        self, valid_upstr_order, chi_integrand, chi_array, mean_dx
    ):
        """Calculates chi at each channel node by summing chi_integrand.

        This method assumes a uniform, mean spacing between nodes. Method is
        deliberately split out for potential cythonization at a later stage.

        Parameters
        ----------
        valid_upstr_order : array of ints
            nodes in the channel network in upstream order.
        chi_integrand : array of floats
            The value (A0/A)**concavity, in upstream order.
        chi_array : array of floats
            Array in which to store chi.
        mean_dx : float
            The mean node spacing in the network.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import ChiFinder
        >>> mg = RasterModelGrid((5, 4))
        >>> for nodes in (
        ...     mg.nodes_at_right_edge,
        ...     mg.nodes_at_bottom_edge,
        ...     mg.nodes_at_top_edge,
        ... ):
        ...     mg.status_at_node[nodes] = mg.BC_NODE_IS_CLOSED
        >>> z = mg.node_x.copy()
        >>> z[[5, 13]] = z[6]  # guard nodes
        >>> _ = mg.add_field("topographic__elevation", z, at="node")
        >>> fr = FlowAccumulator(mg, flow_director="D8")
        >>> cf = ChiFinder(mg)
        >>> fr.run_one_step()
        >>> ch_nodes = np.array([4, 8, 12, 5, 9, 13, 6, 10, 14])
        >>> ch_integrand = 3.0 * np.ones(9, dtype=float)  # to make calc clearer
        >>> chi_array = np.zeros(mg.number_of_nodes, dtype=float)
        >>> cf.integrate_chi_avg_dx(ch_nodes, ch_integrand, chi_array, 0.5)
        >>> chi_array.reshape(mg.shape)
        array([[0. , 0. , 0. , 0. ],
               [1.5, 3. , 4.5, 0. ],
               [1.5, 3. , 4.5, 0. ],
               [1.5, 3. , 4.5, 0. ],
               [0. , 0. , 0. , 0. ]])
        """
        receivers = self._grid.at_node["flow__receiver_node"]
        # because chi_array is all zeros, BC cases where node is receiver
        # resolve themselves
        for node, integrand in izip(valid_upstr_order, chi_integrand):
            dstr_node = receivers[node]
            chi_array[node] = chi_array[dstr_node] + integrand
        chi_array *= mean_dx

    def integrate_chi_each_dx(
        self, valid_upstr_order, chi_integrand_at_nodes, chi_array
    ):
        """Calculates chi at each channel node by summing chi_integrand*dx.

        This method accounts explicitly for spacing between each node. Method
        is deliberately split out for potential cythonization at a later
        stage. Uses a trapezium integration method.

        Parameters
        ----------
        valid_upstr_order : array of ints
            nodes in the channel network in upstream order.
        chi_integrand_at_nodes : array of floats
            The value (A0/A)**concavity, in *node* order.
        chi_array : array of floats
            Array in which to store chi.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import ChiFinder
        >>> mg = RasterModelGrid((5, 4), xy_spacing=3.0)
        >>> for nodes in (
        ...     mg.nodes_at_right_edge,
        ...     mg.nodes_at_bottom_edge,
        ...     mg.nodes_at_top_edge,
        ... ):
        ...     mg.status_at_node[nodes] = mg.BC_NODE_IS_CLOSED
        >>> z = mg.node_x.copy()
        >>> z[[5, 13]] = z[6]  # guard nodes
        >>> _ = mg.add_field("topographic__elevation", z, at="node")
        >>> fr = FlowAccumulator(mg, flow_director="D8")
        >>> cf = ChiFinder(mg)
        >>> fr.run_one_step()
        >>> ch_nodes = np.array([4, 8, 12, 5, 9, 13, 6, 10, 14])
        >>> ch_integrand = 2.0 * np.ones(
        ...     mg.number_of_nodes, dtype=float
        ... )  # to make calc clearer
        >>> chi_array = np.zeros(mg.number_of_nodes, dtype=float)
        >>> cf.integrate_chi_each_dx(ch_nodes, ch_integrand, chi_array)
        >>> chi_array.reshape(mg.shape)
        array([[ 0.        ,  0.        ,  0.        ,  0.        ],
               [ 0.        ,  6.        , 14.48528137,  0.        ],
               [ 0.        ,  6.        , 12.        ,  0.        ],
               [ 0.        ,  6.        , 14.48528137,  0.        ],
               [ 0.        ,  0.        ,  0.        ,  0.        ]])
        >>> from landlab.components import FastscapeEroder
        >>> mg2 = RasterModelGrid((5, 5), xy_spacing=100.0)
        >>> for nodes in (
        ...     mg2.nodes_at_right_edge,
        ...     mg2.nodes_at_bottom_edge,
        ...     mg2.nodes_at_top_edge,
        ... ):
        ...     mg2.status_at_node[nodes] = mg2.BC_NODE_IS_CLOSED
        >>> _ = mg2.add_zeros("node", "topographic__elevation")
        >>> mg2.at_node["topographic__elevation"][mg2.core_nodes] = (
        ...     mg2.node_x[mg2.core_nodes] / 1000.0
        ... )
        >>> np.random.seed(0)
        >>> mg2.at_node["topographic__elevation"][mg2.core_nodes] += np.random.rand(
        ...     mg2.number_of_core_nodes
        ... )
        >>> fr2 = FlowAccumulator(mg2, flow_director="D8")
        >>> sp2 = FastscapeEroder(mg2, K_sp=0.01)
        >>> cf2 = ChiFinder(
        ...     mg2, min_drainage_area=1.0, reference_concavity=0.5, use_true_dx=True
        ... )
        >>> for i in range(10):
        ...     mg2.at_node["topographic__elevation"][mg2.core_nodes] += 10.0
        ...     fr2.run_one_step()
        ...     sp2.run_one_step(1000.0)
        ...
        >>> fr2.run_one_step()
        >>> output_array = np.zeros(25, dtype=float)
        >>> cf2.integrate_chi_each_dx(
        ...     mg2.at_node["flow__upstream_node_order"],
        ...     np.ones(25, dtype=float),
        ...     output_array,
        ... )
        >>> output_array.reshape(mg2.shape)
        array([[  0. ,   0. ,   0.        ,   0.        ,   0. ],
               [  0. , 100. , 200.        , 382.84271247,   0. ],
               [  0. , 100. , 241.42135624, 341.42135624,   0. ],
               [  0. , 100. , 200.        , 300.        ,   0. ],
               [  0. ,   0. ,   0.        ,   0.        ,   0. ]])
        """
        receivers = self._grid.at_node["flow__receiver_node"]
        links = self._grid.at_node["flow__link_to_receiver_node"]

        # because chi_array is all zeros, BC cases where node is receiver
        # resolve themselves
        half_integrand = 0.5 * chi_integrand_at_nodes
        for node in valid_upstr_order:
            dstr_node = receivers[node]
            dstr_link = links[node]
            if dstr_link != self._grid.BAD_INDEX:
                dstr_length = self._link_lengths[dstr_link]
                half_head_val = half_integrand[node]
                half_tail_val = half_integrand[dstr_node]
                mean_val = half_head_val + half_tail_val
                chi_to_add = mean_val * dstr_length
                chi_array[node] = chi_array[dstr_node] + chi_to_add

    def mean_channel_node_spacing(self, ch_nodes):
        """Calculates the mean spacing between all adjacent channel nodes.

        Parameters
        ----------
        ch_nodes : array of ints
            The nodes within the defined channel network.

        Returns
        -------
        mean_spacing : float (m)
            The mean spacing between all nodes in the network.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import ChiFinder
        >>> mg = RasterModelGrid((5, 4), xy_spacing=2.0)
        >>> for nodes in (
        ...     mg.nodes_at_right_edge,
        ...     mg.nodes_at_bottom_edge,
        ...     mg.nodes_at_top_edge,
        ... ):
        ...     mg.status_at_node[nodes] = mg.BC_NODE_IS_CLOSED
        >>> z = mg.node_x.copy()
        >>> z[[5, 13]] = z[6]  # guard nodes
        >>> _ = mg.add_field("topographic__elevation", z, at="node")
        >>> fr = FlowAccumulator(mg, flow_director="D8")
        >>> cf = ChiFinder(mg)
        >>> fr.run_one_step()
        >>> ch_nodes = np.array([4, 8, 12, 5, 9, 13, 6, 10, 14])
        >>> cf.mean_channel_node_spacing(ch_nodes)
        2.2761423749153966
        """
        ch_links = self._grid.at_node["flow__link_to_receiver_node"][ch_nodes]
        ch_links_valid = ch_links[ch_links != self._grid.BAD_INDEX]

        valid_link_lengths = self._link_lengths[ch_links_valid]
        return valid_link_lengths.mean()

    @property
    def chi_indices(self):
        """Return the array of channel steepness indices.

        Nodes not in the channel receive zeros.
        """
        return self._chi

    @property
    def hillslope_mask(self):
        """Return a boolean array, False where steepness indices exist."""
        return self._mask

    def best_fit_chi_elevation_gradient_and_intercept(self, ch_nodes=None):
        """Returns least squares best fit for a straight line through a chi
        plot.

        Parameters
        ----------
        ch_nodes : array of ints or None
            Nodes at which to consider chi and elevation values. If None,
            will use all nodes in grid with area greater than the component
            min_drainage_area.

        Returns
        -------
        coeffs : array(gradient, intercept)
            A len-2 array containing the m then z0, where z = z0 + m * chi.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import ChiFinder
        >>> mg = RasterModelGrid((3, 4))
        >>> for nodes in (
        ...     mg.nodes_at_right_edge,
        ...     mg.nodes_at_bottom_edge,
        ...     mg.nodes_at_top_edge,
        ... ):
        ...     mg.status_at_node[nodes] = mg.BC_NODE_IS_CLOSED
        >>> z = mg.add_field("topographic__elevation", mg.node_x.copy(), at="node")
        >>> z[4:8] = np.array([0.5, 1.0, 2.0, 0.0])
        >>> fr = FlowAccumulator(mg, flow_director="D8")
        >>> cf = ChiFinder(mg, min_drainage_area=1.0, reference_concavity=1.0)
        >>> fr.run_one_step()
        >>> cf.calculate_chi()
        >>> mg.at_node["channel__chi_index"].reshape(mg.shape)[1, :]
        array([0.5, 1. , 2. , 0. ])
        >>> coeffs = cf.best_fit_chi_elevation_gradient_and_intercept()
        >>> np.allclose(np.array([1.0, 0.0]), coeffs)
        True
        """
        if ch_nodes is None:
            good_vals = np.logical_not(self._mask)
        else:
            good_vals = np.array(ch_nodes)  # not tested
        chi_vals = self._chi[good_vals]
        elev_vals = self._grid.at_node["topographic__elevation"][good_vals]
        coeffs = np.polyfit(chi_vals, elev_vals, 1)
        return coeffs

    def nodes_downstream_of_channel_head(self, channel_head):
        """Find and return an array with nodes downstream of channel_head.

        Parameters
        ----------
        channel_head : int
            Node ID of channel head from which to get downstream nodes.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator, ChiFinder
        >>> mg = RasterModelGrid((3, 4))
        >>> for nodes in (
        ...     mg.nodes_at_right_edge,
        ...     mg.nodes_at_bottom_edge,
        ...     mg.nodes_at_top_edge,
        ... ):
        ...     mg.status_at_node[nodes] = mg.BC_NODE_IS_CLOSED
        >>> z = mg.add_field("topographic__elevation", mg.node_x.copy(), at="node")
        >>> z[4:8] = np.array([0.5, 1.0, 2.0, 0.0])
        >>> fr = FlowAccumulator(mg, flow_director="D8")
        >>> fr.run_one_step()
        >>> mg.at_node["flow__receiver_node"]
        array([ 0,  1,  2,  3,  4,  4,  5,  7,  8,  9, 10, 11])
        >>> cf = ChiFinder(mg, min_drainage_area=0.0, reference_concavity=1.0)
        >>> cf.calculate_chi()
        >>> cf.nodes_downstream_of_channel_head(6)
        [6, 5, 4]
        """
        ch_nodes = []
        current_node = channel_head
        while True:
            ch_A = self._grid.at_node["drainage_area"][current_node]
            if ch_A > self._min_drainage:
                ch_nodes.append(current_node)
            next_node = self._grid.at_node["flow__receiver_node"][current_node]
            if next_node == current_node:
                break
            else:
                current_node = next_node
        return ch_nodes

    def create_chi_plot(
        self,
        channel_heads=None,
        label_axes=True,
        symbol="kx",
        plot_line=False,
        line_symbol="r-",
    ):
        """Plots a "chi plot" (chi vs elevation for points in channel network).

        If channel_heads is provided, only the channel nodes downstream of
        the provided points (and with area > min_drainage_area) will be
        plotted.

        Parameters
        ----------

        channel_heads : int, list or array of ints, or None
            Node IDs of channel heads to from which plot downstream.
        label_axes : bool
            If True, labels the axes as "Chi" and "Elevation (m)".
        symbol : str
            A matplotlib-style string for the style to use for the points.
        plot_line : bool
            If True, will plot a linear best fit line through the data cloud.
        line_symbol : str
            A matplotlib-style string for the style to use for the line, if
            plot_line.
        """
        from matplotlib.pyplot import clf
        from matplotlib.pyplot import figure
        from matplotlib.pyplot import plot
        from matplotlib.pyplot import xlabel
        from matplotlib.pyplot import ylabel

        figure("Chi plot")
        clf()
        if channel_heads is not None:
            if plot_line:
                good_nodes = set()
            if isinstance(channel_heads, int):
                channel_heads = [channel_heads]
            for head in channel_heads:
                ch_nodes = self.nodes_downstream_of_channel_head(head)
                plot(
                    self._chi[ch_nodes],
                    self._grid.at_node["topographic__elevation"][ch_nodes],
                    symbol,
                )
                if plot_line:
                    good_nodes.update(ch_nodes)
        else:
            ch_nodes = np.logical_not(self._mask)
            plot(
                self._chi[ch_nodes],
                self._grid.at_node["topographic__elevation"][ch_nodes],
                symbol,
            )
            good_nodes = ch_nodes
        if plot_line:
            coeffs = self.best_fit_chi_elevation_gradient_and_intercept(good_nodes)
            p = np.poly1d(coeffs)
            chirange = np.linspace(
                self._chi[good_nodes].min(), self._chi[good_nodes].max(), 100
            )
            plot(chirange, p(chirange), line_symbol)
        if label_axes:
            ylabel("Elevation (m)")
            xlabel("Chi")

    @property
    def masked_chi_indices(self):
        """Returns a masked array version of the 'channel__chi_index' field.
        This enables easier plotting of the values with.

        :func:`landlab.imshow_grid_at_node` or similar.

        Examples
        --------
        Make a topographic map with an overlay of chi values:

        >>> from landlab import imshow_grid_at_node
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator, FastscapeEroder
        >>> from landlab.components import ChiFinder
        >>> mg = RasterModelGrid((5, 5), xy_spacing=100.0)
        >>> for nodes in (
        ...     mg.nodes_at_right_edge,
        ...     mg.nodes_at_bottom_edge,
        ...     mg.nodes_at_top_edge,
        ... ):
        ...     mg.status_at_node[nodes] = mg.BC_NODE_IS_CLOSED
        >>> _ = mg.add_zeros("node", "topographic__elevation")
        >>> mg.at_node["topographic__elevation"][mg.core_nodes] = (
        ...     mg.node_x[mg.core_nodes] / 1000.0
        ... )
        >>> np.random.seed(0)
        >>> mg.at_node["topographic__elevation"][mg.core_nodes] += np.random.rand(
        ...     mg.number_of_core_nodes
        ... )
        >>> fr = FlowAccumulator(mg, flow_director="D8")
        >>> sp = FastscapeEroder(mg, K_sp=0.01)
        >>> cf = ChiFinder(mg, min_drainage_area=20000.0)
        >>> for i in range(10):
        ...     mg.at_node["topographic__elevation"][mg.core_nodes] += 10.0
        ...     fr.run_one_step()
        ...     sp.run_one_step(1000.0)
        ...
        >>> fr.run_one_step()
        >>> cf.calculate_chi()

        >>> imshow_grid_at_node(mg, "topographic__elevation", allow_colorbar=False)
        >>> imshow_grid_at_node(
        ...     mg, cf.masked_chi_indices, color_for_closed=None, cmap="winter"
        ... )
        """
        return np.ma.array(self._chi, mask=self._mask)
