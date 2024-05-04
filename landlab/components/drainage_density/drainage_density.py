"""Landlab component to calculate drainage density."""

from warnings import warn

import numpy as np

from landlab import Component


class DrainageDensity(Component):
    r"""
    Calculate drainage density over a DEM.

    Landlab component that implements the distance to channel algorithm of
    Tucker et al., 2001.

    This component requires EITHER a ``channel__mask array`` with 1's
    where channels exist and 0's elsewhere, OR a set of coefficients
    and exponents for a slope-area relationship and a
    channelization threshold to compare against that relationship.

    If an array is provided it MUST be of type ``np.uint8``. See the example
    below for how to make such an array.

    The ``channel__mask`` array will be assigned to an at-node field with the
    name ``channel__mask``. If the channel__mask was originaly created from a
    passed array, a user can update this array to change the mask.

    If the ``channel__mask`` is created using an area coefficent,
    slope coefficient, area exponent, slope exponent, and channelization
    threshold, the location of the mask will be re-update when
    calculate_drainage_density is called.

    If an area coefficient, :math:`C_A`, a slope coefficent, :math:`C_S`, an
    area exponent, :math:`m_r`, a slope exponent, :math:`n_r`, and
    channelization threshold :math:`T_C` are provided, nodes that meet the
    criteria

    .. math::

       C_A A^{m_r} C_s S^{n_r} > T_c

    where :math:`A` is the drainage density and :math:`S` is the local slope,
    will be marked as channel nodes.

    The ``calculate_drainage_density`` function returns drainage density for the
    model domain. This function calculates the distance from every node to the
    nearest channel node :math:`L` along the flow line of steepest descent
    (assuming D8 routing if the grid is a RasterModelGrid).

    This component stores this distance a field, called:
    ``surface_to_channel__minimum_distance``. The drainage density is then
    calculated (after Tucker et al., 2001):

    .. math::

       D_d = \frac{1}{2\overline{L}}

    where :math:`\overline{L}` is the mean L for the model domain.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, FastscapeEroder
    >>> mg = RasterModelGrid((10, 10))
    >>> _ = mg.add_zeros("node", "topographic__elevation")
    >>> np.random.seed(50)
    >>> noise = np.random.rand(100)
    >>> mg.at_node["topographic__elevation"] += noise
    >>> mg.at_node["topographic__elevation"].reshape(mg.shape)
    array([[0.49460165, 0.2280831 , 0.25547392, 0.39632991, 0.3773151 ,
            0.99657423, 0.4081972 , 0.77189399, 0.76053669, 0.31000935],
           [0.3465412 , 0.35176482, 0.14546686, 0.97266468, 0.90917844,
            0.5599571 , 0.31359075, 0.88820004, 0.67457307, 0.39108745],
           [0.50718412, 0.5241035 , 0.92800093, 0.57137307, 0.66833757,
            0.05225869, 0.3270573 , 0.05640164, 0.17982769, 0.92593317],
           [0.93801522, 0.71409271, 0.73268761, 0.46174768, 0.93132927,
            0.40642024, 0.68320577, 0.64991587, 0.59876518, 0.22203939],
           [0.68235717, 0.8780563 , 0.79671726, 0.43200225, 0.91787822,
            0.78183368, 0.72575028, 0.12485469, 0.91630845, 0.38771099],
           [0.29492955, 0.61673141, 0.46784623, 0.25533891, 0.83899589,
            0.1786192 , 0.22711417, 0.65987645, 0.47911625, 0.07344734],
           [0.13896007, 0.11230718, 0.47778497, 0.54029623, 0.95807105,
            0.58379231, 0.52666409, 0.92226269, 0.91925702, 0.25200886],
           [0.68263261, 0.96427612, 0.22696165, 0.7160172 , 0.79776011,
            0.9367512 , 0.8537225 , 0.42154581, 0.00543987, 0.03486533],
           [0.01390537, 0.58890993, 0.3829931 , 0.11481895, 0.86445401,
            0.82165703, 0.73749168, 0.84034417, 0.4015291 , 0.74862   ],
           [0.55962945, 0.61323757, 0.29810165, 0.60237917, 0.42567684,
            0.53854438, 0.48672986, 0.49989164, 0.91745948, 0.26287702]])
    >>> fr = FlowAccumulator(mg, flow_director="D8")
    >>> fsc = FastscapeEroder(mg, K_sp=0.01, m_sp=0.5, n_sp=1)
    >>> for x in range(100):
    ...     fr.run_one_step()
    ...     fsc.run_one_step(dt=10.0)
    ...     mg.at_node["topographic__elevation"][mg.core_nodes] += 0.01
    ...
    >>> channels = np.array(mg.at_node["drainage_area"] > 5, dtype=np.uint8)
    >>> dd = DrainageDensity(mg, channel__mask=channels)
    >>> mean_drainage_density = dd.calculate_drainage_density()
    >>> np.isclose(mean_drainage_density, 0.3831100571)
    True

    Alternatively you can pass a set of coefficients to identify the channel
    mask. Next shows the same example as above, but with these coefficients
    provided.

    >>> mg = RasterModelGrid((10, 10))
    >>> _ = mg.add_zeros("node", "topographic__elevation")
    >>> np.random.seed(50)
    >>> noise = np.random.rand(100)
    >>> mg.at_node["topographic__elevation"] += noise
    >>> fr = FlowAccumulator(mg, flow_director="D8")
    >>> fsc = FastscapeEroder(mg, K_sp=0.01, m_sp=0.5, n_sp=1)
    >>> for x in range(100):
    ...     fr.run_one_step()
    ...     fsc.run_one_step(dt=10.0)
    ...     mg.at_node["topographic__elevation"][mg.core_nodes] += 0.01
    ...
    >>> channels = np.array(mg.at_node["drainage_area"] > 5, dtype=np.uint8)
    >>> dd = DrainageDensity(
    ...     mg,
    ...     area_coefficient=1.0,
    ...     slope_coefficient=1.0,
    ...     area_exponent=1.0,
    ...     slope_exponent=0.0,
    ...     channelization_threshold=5,
    ... )
    >>> mean_drainage_density = dd.calculate_drainage_density()
    >>> np.isclose(mean_drainage_density, 0.3831100571)
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Tucker, G., Catani, F., Rinaldo, A., Bras, R. (2001). Statistical analysis
    of drainage density from digital terrain data. Geomorphology 36(3-4),
    187-202. https://dx.doi.org/10.1016/s0169-555x(00)00056-8

    """

    _name = "DrainageDensity"

    _unit_agnostic = True

    _info = {
        "area_coefficient": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Area coefficient to define channels.",
        },
        "area_exponent": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Area exponent to define channels.",
        },
        "channel__mask": {
            "dtype": np.uint8,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Logical map of at which grid nodes channels are present",
        },
        "channelization_threshold": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Channelization threshold for use with area and slope "
                "coefficients and exponents."
            ),
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
        "slope_coefficient": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Slope coefficient to define channels.",
        },
        "slope_exponent": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Slope exponent to define channels.",
        },
        "surface_to_channel__minimum_distance": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Distance from each node to the nearest channel",
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
        channel__mask=None,
        area_coefficient=None,
        slope_coefficient=None,
        area_exponent=None,
        slope_exponent=None,
        channelization_threshold=None,
    ):
        """Initialize the DrainageDensity component.

        Parameters
        ----------
        grid : ModelGrid
        channel__mask : Array that holds 1's where
            channels exist and 0's elsewhere
        area_coefficient : coefficient to multiply drainage area by,
            for calculating channelization threshold
        slope_coefficient : coefficient to multiply slope by,
            for calculating channelization threshold
        area_exponent : exponent to raise drainage area to,
            for calculating channelization threshold
        slope_exponent : exponent to raise slope to,
            for calculating channelization threshold
        channelization_threshold : threshold value above
            which channels exist
        """
        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that DrainageDensity is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        if channel__mask is not None:
            if area_coefficient is not None:
                warn(
                    "Channel mask and area "
                    "coefficient supplied. Defaulting "
                    "to channel mask, ignoring area "
                    "coefficient."
                )
            if slope_coefficient is not None:
                warn(
                    "Channel mask and slope "
                    "coefficient supplied. Defaulting "
                    "to channel mask, ignoring slope "
                    "coefficient."
                )
            if area_exponent is not None:
                warn(
                    "Channel mask and area "
                    "exponent supplied. Defaulting "
                    "to channel mask, ignoring area "
                    "exponent."
                )
            if slope_exponent is not None:
                warn(
                    "Channel mask and slope "
                    "exponent supplied. Defaulting "
                    "to channel mask, ignoring slope "
                    "exponent."
                )
            if channelization_threshold is not None:
                warn(
                    "Channel mask and channelization "
                    "threshold supplied. Defaulting "
                    "to channel mask, ignoring "
                    "threshold."
                )
            if grid.number_of_nodes != len(channel__mask):
                raise ValueError(
                    "Length of channel mask is not equal to " "number of grid nodes"
                )

            if "channel__mask" in grid.at_node:
                warn("Existing channel__mask grid field was overwritten.")

            if channel__mask.dtype.type is not np.uint8:
                raise ValueError("mask must by np.uint8")

            self._mask_as_array = True
            self._update_channel_mask = self._update_channel_mask_array
            grid.at_node["channel__mask"] = channel__mask

        if channel__mask is None:
            if area_coefficient is None:
                raise ValueError(
                    "No channel mask and no area "
                    "coefficient supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if slope_coefficient is None:
                raise ValueError(
                    "No channel mask and no slope "
                    "coefficient supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if area_exponent is None:
                raise ValueError(
                    "No channel mask and no area "
                    "exponent supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if slope_exponent is None:
                raise ValueError(
                    "No channel mask and no slope "
                    "exponent supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )
            if channelization_threshold is None:
                raise ValueError(
                    "No channel mask and no channelization "
                    "threshold supplied. Either "
                    "a channel mask or all 5 threshold "
                    "parameters are needed."
                )

            self._mask_as_array = False
            self._update_channel_mask = self._update_channel_mask_values
            self._area_coefficient = area_coefficient
            self._slope_coefficient = slope_coefficient
            self._area_exponent = area_exponent
            self._slope_exponent = slope_exponent
            self._channelization_threshold = channelization_threshold

            self._update_channel_mask()

        # for this component to work with Cython acceleration,
        # the channel_network must be uint8, not bool...
        self._channel_network = grid.at_node["channel__mask"]

        # Flow receivers
        self._flow_receivers = grid.at_node["flow__receiver_node"]

        # Links to receiver nodes
        self._stack_links = grid.at_node["flow__link_to_receiver_node"]

        # Upstream node order
        self._upstream_order = grid.at_node["flow__upstream_node_order"]

        # Distance to channel
        if "surface_to_channel__minimum_distance" not in grid.at_node:
            grid.add_zeros(
                "surface_to_channel__minimum_distance", at="node", dtype=float
            )
        self._distance_to_channel = grid.at_node["surface_to_channel__minimum_distance"]

        # Use the appropriate array for link or d8 lengths
        try:
            self._length_of_link = self._grid.length_of_d8
        except AttributeError:
            self._length_of_link = self._grid.length_of_link

    def _update_channel_mask_array(self):
        raise NotImplementedError(
            "If you provided a channel mask to "
            "DrainageDensity, update it by updating the "
            "model grid field channel__mask"
        )

    def _update_channel_mask_values(self):
        channel__mask = (
            self._area_coefficient
            * np.power(self._grid.at_node["drainage_area"], self._area_exponent)
            * self._slope_coefficient
            * np.power(
                self._grid.at_node["topographic__steepest_slope"], self._slope_exponent
            )
        ) > self._channelization_threshold
        self._grid.at_node["channel__mask"] = channel__mask.astype(np.uint8)

    def calculate_drainage_density(self):
        """Calculate drainage density.

        If the channel mask is defined based on slope and area coefficients,
        it will be update based on the current drainage area and slope fields.

        Returns
        -------
        landscape_drainage_density : float (1/m)
            Drainage density over the model domain.
        """
        from .cfuncs import _calc_dists_to_channel

        if self._mask_as_array is False:
            self._update_channel_mask()

        _calc_dists_to_channel(
            self._channel_network,
            self._flow_receivers,
            self._upstream_order,
            self._length_of_link,
            self._stack_links,
            self._distance_to_channel,
            self._grid.number_of_nodes,
        )
        landscape_drainage_density = 1.0 / (
            2.0
            * np.mean(
                self._grid.at_node["surface_to_channel__minimum_distance"][
                    self._grid.core_nodes
                ]
            )
        )
        # this is THE drainage density
        return landscape_drainage_density
