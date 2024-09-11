"""Landlab component that simulates detachment-limited river erosion.

This component calculates changes in elevation in response to
vertical incision.
"""

import numpy as np

from landlab import Component


class DetachmentLtdErosion(Component):
    """Simulate detachment limited sediment transport.

    Landlab component that simulates detachment limited sediment transport is more
    general than the stream power component. Doesn't require the upstream node
    order, links to flow receiver and flow receiver fields. Instead, takes in
    the discharge values on NODES calculated by the OverlandFlow class and
    erodes the landscape in response to the output discharge.

    As of right now, this component relies on the OverlandFlow component
    for stability. There are no stability criteria implemented in this class.
    To ensure model stability, use StreamPowerEroder or FastscapeEroder
    components instead.

    .. codeauthor:: Jordan Adams

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import DetachmentLtdErosion

    Create a grid on which to calculate detachment ltd sediment transport.

    >>> grid = RasterModelGrid((4, 5))

    The grid will need some data to provide the detachment limited sediment
    transport component. To check the names of the fields that provide input to
    the detachment ltd transport component, use the *input_var_names* class
    property.

    Create fields of data for each of these input variables.

    >>> grid.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [1.0, 1.0, 1.0, 1.0, 1.0],
    ...     [2.0, 2.0, 2.0, 2.0, 2.0],
    ...     [3.0, 3.0, 3.0, 3.0, 3.0],
    ... ]

    Using the set topography, now we will calculate slopes on all nodes.

    >>> grid.at_node["topographic__slope"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.70710678, 1.0, 1.0, 1.0, 0.70710678],
    ...     [0.70710678, 1.0, 1.0, 1.0, 0.70710678],
    ...     [0.70710678, 1.0, 1.0, 1.0, 0.70710678],
    ... ]


    Now we will arbitrarily add water discharge to each node for simplicity.

    >>> grid.at_node["surface_water__discharge"] = [
    ...     [30.0, 30.0, 30.0, 30.0, 30.0],
    ...     [20.0, 20.0, 20.0, 20.0, 20.0],
    ...     [10.0, 10.0, 10.0, 10.0, 10.0],
    ...     [5.0, 5.0, 5.0, 5.0, 5.0],
    ... ]

    Instantiate the `DetachmentLtdErosion` component to work on this grid, and
    run it. In this simple case, we need to pass it a time step ('dt')

    >>> dt = 10.0
    >>> dle = DetachmentLtdErosion(grid)
    >>> dle.run_one_step(dt=dt)

    After calculating the erosion rate, the elevation field is updated in the
    grid. Use the *output_var_names* property to see the names of the fields that
    have been changed.

    >>> dle.output_var_names
    ('topographic__elevation',)

    The `topographic__elevation` field is defined at nodes.

    >>> dle.var_loc("topographic__elevation")
    'node'


    Now we test to see how the topography changed as a function of the erosion
    rate.

    >>> grid.at_node["topographic__elevation"].reshape(grid.shape)
    array([[0.        , 0.        , 0.        , 0.        , 0.        ],
           [0.99936754, 0.99910557, 0.99910557, 0.99910557, 0.99936754],
           [1.99955279, 1.99936754, 1.99936754, 1.99936754, 1.99955279],
           [2.99968377, 2.99955279, 2.99955279, 2.99955279, 2.99968377]])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Howard, A. (1994). A detachment-limited model of drainage basin evolution. Water
    Resources Research  30(7), 2261-2285. https://dx.doi.org/10.1029/94wr00757
    """

    _name = "DetachmentLtdErosion"

    _unit_agnostic = True

    _info = {
        "surface_water__discharge": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__slope": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "gradient of the ground surface",
        },
    }

    def __init__(
        self,
        grid,
        K_sp=0.00002,
        m_sp=0.5,
        n_sp=1.0,
        uplift_rate=0.0,
        entrainment_threshold=0.0,
        slope="topographic__slope",
    ):
        """Calculate detachment limited erosion rate on nodes.

        Landlab component that generalizes the detachment limited erosion
        equation, primarily to be coupled to the the Landlab OverlandFlow
        component.

        This component adjusts topographic elevation.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab grid.
        K_sp : float, optional
            K in the stream power equation (units vary with other parameters -
            if used with the de Almeida equation it is paramount to make sure
            the time component is set to *seconds*, not *years*!)
        m_sp : float, optional
            Stream power exponent, power on discharge
        n_sp : float, optional
            Stream power exponent, power on slope
        uplift_rate : float, optional
            changes in topographic elevation due to tectonic uplift
        entrainment_threshold : float, optional
            threshold for sediment movement
        slope : str
            Field name of an at-node field that contains the slope.
        """
        super().__init__(grid)

        assert slope in grid.at_node

        self._K = K_sp
        self._m = m_sp
        self._n = n_sp

        self._I = self._grid.zeros(at="node")  # noqa: E741
        self._uplift_rate = uplift_rate
        self._entrainment_threshold = entrainment_threshold

        self._dzdt = self._grid.zeros(at="node")

    def run_one_step(self, dt):
        """Erode into grid topography.

        For one time step, this erodes into the grid topography using
        the water discharge and topographic slope.

        The grid field 'topographic__elevation' is altered each time step.

        Parameters
        ----------
        dt : float
            Time step.
        """

        S = self._grid.at_node["topographic__slope"]
        Q = self._grid.at_node["surface_water__discharge"]

        Q_to_m = np.power(Q, self._m)

        S_to_n = np.power(S, self._n)

        self._I = (
            self._K * Q_to_m * S_to_n
        ) - self._entrainment_threshold  # noqa: E741

        self._I[self._I < 0.0] = 0.0

        self._dz = (self._uplift_rate - self._I) * dt

        self._grid["node"]["topographic__elevation"] += self._dz
