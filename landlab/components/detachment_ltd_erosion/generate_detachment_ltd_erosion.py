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

>>> grid.at_node['topographic__elevation'] = np.array([
...     0., 0., 0., 0., 0.,
...     1., 1., 1., 1., 1.,
...     2., 2., 2., 2., 2.,
...     3., 3., 3., 3., 3.])

Using the set topography, now we will calculate slopes on all nodes.


>>> grid.at_node['topographic__slope'] = np.array([
...     -0.        , -0.        , -0.        , -0.        , -0,
...      0.70710678,  1.        ,  1.        ,  1.        ,  0.70710678,
...      0.70710678,  1.        ,  1.        ,  1.        ,  0.70710678,
...      0.70710678,  1.        ,  1.        ,  1.        ,  0.70710678])


Now we will arbitrarily add water discharge to each node for simplicity.
>>> grid.at_node['surface_water__discharge'] = np.array([
...     30., 30., 30., 30., 30.,
...     20., 20., 20., 20., 20.,
...     10., 10., 10., 10., 10.,
...      5., 5., 5., 5., 5.])

Instantiate the `DetachmentLtdErosion` component to work on this grid, and
run it. In this simple case, we need to pass it a time step ('dt')

>>> dt = 10.0
>>> dle = DetachmentLtdErosion(grid)
>>> dle.erode(dt=dt)

After calculating the erosion rate, the elevation field is updated in the
grid. Use the *output_var_names* property to see the names of the fields that
have been changed.

>>> dle.output_var_names
('topographic__elevation',)

The `topographic__elevation` field is defined at nodes.

>>> dle.var_loc('topographic__elevation')
'node'


Now we test to see how the topography changed as a function of the erosion
rate.

>>> grid.at_node['topographic__elevation'] # doctest: +NORMALIZE_WHITESPACE
array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
        0.99936754,  0.99910557,  0.99910557,  0.99910557,  0.99936754,
        1.99955279,  1.99936754,  1.99936754,  1.99936754,  1.99955279,
        2.99968377,  2.99955279,  2.99955279,  2.99955279,  2.99968377])

"""

import numpy as np

from landlab import Component
from landlab.field.scalar_data_fields import FieldError


class DetachmentLtdErosion(Component):

    """Landlab component that simulates detachment-limited river erosion.

    This component calculates changes in elevation in response to vertical
    incision.
    """

    _name = "DetachmentLtdErosion"

    _input_var_names = (
        "topographic__elevation",
        "topographic__slope",
        "surface_water__discharge",
    )

    _output_var_names = ("topographic__elevation",)

    _var_units = {
        "topographic__elevation": "m",
        "topographic__slope": "-",
        "surface_water__discharge": "m^3/s",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "topographic__slope": "node",
        "surface_water__discharge": "node",
    }

    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "topographic__slope": "Slope of ",
        "surface_water__discharge": "node",
    }

    def __init__(
        self,
        grid,
        K_sp=0.00002,
        m_sp=0.5,
        n_sp=1.0,
        uplift_rate=0.0,
        entrainment_threshold=0.0,
        **kwds
    ):
        """Calculate detachment limited erosion rate on nodes.

        Landlab component that generalizes the detachment limited erosion
        equation, primarily to be coupled to the the Landlab OverlandFlow
        component.

        This component adjusts topographic elevation and is contained in the
        landlab.components.detachment_ltd_sed_trp folder.

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
        """
        super(DetachmentLtdErosion, self).__init__(grid, **kwds)

        self.K = K_sp
        self.m = m_sp
        self.n = n_sp

        self.I = self._grid.zeros(at="node")  # noqa: E741
        self.uplift_rate = uplift_rate
        self.entrainment_threshold = entrainment_threshold

        self.dzdt = self._grid.zeros(at="node")

    def erode(
        self,
        dt,
        elevs="topographic__elevation",
        discharge_cms="surface_water__discharge",
        slope="topographic__slope",
    ):
        """Erode into grid topography.

        For one time step, this erodes into the grid topography using
        the water discharge and topographic slope.

        The grid field 'topographic__elevation' is altered each time step.

        Parameters
        ----------
        dt : float
            Time step.
        discharge_cms : str, optional
            Name of the field that represents discharge on the nodes, if
            from the de Almeida solution have units of cubic meters per second.
        slope : str, optional
            Name of the field that represent topographic slope on each node.
        """
        try:
            S = self._grid.at_node[slope]
        except FieldError:
            raise ValueError("missing field for slope")

        if type(discharge_cms) is str:
            Q = self._grid.at_node[discharge_cms]
        else:
            Q = discharge_cms

        Q_to_m = np.power(Q, self.m)

        S_to_n = np.power(S, self.n)

        self.I = (self.K * Q_to_m * S_to_n) - self.entrainment_threshold  # noqa: E741

        self.I[self.I < 0.0] = 0.0

        self.dz = (self.uplift_rate - self.I) * dt

        self._grid["node"][elevs] += self.dz
