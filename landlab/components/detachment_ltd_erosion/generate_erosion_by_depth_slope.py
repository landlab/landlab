"""Calculates erosion rate as a function of the depth-slope product

Erosion rate = k_e * ((Tau**a - Tau_crit**a))

k_e = erodibility coefficient
Tau = bed shear stress
    = density of fluid (rho) * gravitational acceleration (g) * water depths (h) * slopes (S)
Tau_crit = critical shear stress
a = positive exponent

Note this equation was presented in Tucker, G.T., 2004, Drainage basin
sensitivityto tectonic and climatic forcing: Implications of a stochastic
model for the role of entrainment and erosion thresholds,
Earth Surface Processes and Landforms.

More generalized than other erosion components, as it doesn't require the
upstream node order, links to flow receiver and flow receiver fields. Instead,
takes in the water depth and slope fields on NODES calculated by the
OverlandFlow class and erodes the landscape in response to the hydrograph
generted by that method.

As of right now, this component relies on the OverlandFlow component
for stability. There are no stability criteria implemented in this class.
To ensure model stability, use StreamPowerEroder or FastscapeEroder
components instead.

.. codeauthor:: Jordan Adams

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components import DepthSlopeProductErosion

Create a grid on which to calculate detachment ltd sediment transport.

>>> grid = RasterModelGrid((5, 5))

The grid will need some data to provide the detachment limited sediment
transport component. To check the names of the fields that provide input to
the detachment ltd transport component, use the *input_var_names* class
property.

Create fields of data for each of these input variables.

First create toopgraphy. This is a flat surface of elevation 10 m.
>>> grid.at_node['topographic__elevation'] = np.ones(grid.number_of_nodes)
>>> grid.at_node['topographic__elevation'] *= 10.
>>> grid.at_node['topographic__elevation'] = np.array([
...      10., 10., 10., 10., 10.,
...      10., 10., 10., 10., 10.,
...      10., 10., 10., 10., 10.,
...      10., 10., 10., 10., 10.,
...      10., 10., 10., 10., 10.])

Now we'll add an arbitrary water depth field on top of that topography.
>>> grid.at_node['surface_water__depth'] = np.array([
...      5., 5., 5., 5., 5.,
...      4., 4., 4., 4., 4.,
...      3., 3., 3., 3., 3.,
...      2., 2., 2., 2., 2.,
...      1., 1., 1., 1., 1.])

Using the set topography, now we will calculate slopes on all nodes.

First calculating slopes on links
>>> grid.at_link['water_surface__slope'] = grid.calc_grad_at_link('surface_water__depth')

Now putting slopes on nodes

>>> grid['node']['water_surface__slope'] = (grid['link']['water_surface__slope'][grid.links_at_node] * grid.active_link_dirs_at_node).max(axis=1) # doctest: +NORMALIZE_WHITESPACE
>>> grid.at_node['water_surface__slope']
array([ 0.,  1.,  1.,  1.,  0., -0.,  1.,  1.,  1.,  0., -0.,  1.,  1.,
        1.,  0., -0.,  1.,  1.,  1.,  0.,  0.,  0.,  0.,  0.,  0.])


Instantiate the `DepthSlopeProductErosion` component to work on this grid, and
run it. In this simple case, we need to pass it a time step ('dt') and also
an erodibility factor ('k_e').

>>> dt = 1.
>>> dspe = DepthSlopeProductErosion(grid, k_e=0.00005)
>>> dspe.erode(dt=dt, slope='water_surface__slope')

Now we test to see how the topography changed as a function of the erosion
rate. First, we'll look at the erosion rate:

>>> dspe.dz   # doctest: +NORMALIZE_WHITESPACE
array([ 0.    , -2.4525, -2.4525, -2.4525,  0.    ,  0.    , -1.962 ,
       -1.962 , -1.962 ,  0.    ,  0.    , -1.4715, -1.4715, -1.4715,
        0.    ,  0.    , -0.981 , -0.981 , -0.981 ,  0.    ,  0.    ,
        0.    ,  0.    ,  0.    ,  0.    ])

Now, our updated topography...
>>> grid.at_node['topographic__elevation'] # doctest: +NORMALIZE_WHITESPACE
array([ 10.    ,   7.5475,   7.5475,   7.5475,  10.    ,  10.    ,
         8.038 ,   8.038 ,   8.038 ,  10.    ,  10.    ,   8.5285,
         8.5285,   8.5285,  10.    ,  10.    ,   9.019 ,   9.019 ,
         9.019 ,  10.    ,  10.    ,  10.    ,  10.    ,  10.    ,  10.    ])
"""

import numpy as np

from landlab import Component
from landlab.field.scalar_data_fields import FieldError


class DepthSlopeProductErosion(Component):

    """Landlab component that simulates detachment-limited river erosion.

    This component calculates changes in elevation in response to vertical
    incision.
    """

    _name = "DetachmentLtdErosion"

    _input_var_names = (
        "topographic__elevation",
        "topographic__slope",
        "surface_water__depth",
    )

    _output_var_names = ("topographic__elevation",)

    _var_units = {
        "topographic__elevation": "m",
        "topographic__slope": "-",
        "surface_water__depth": "m",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "topographic__slope": "node",
        "surface_water__depth": "node",
    }

    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "topographic__slope": "Slope of the land surface",
        "surface_water__depth": "Depth of water on the surface",
    }

    def __init__(
        self,
        grid,
        k_e,
        fluid_density=1000.0,
        g=9.81,
        a_exp=1.0,
        tau_crit=0.0,
        uplift_rate=0.0,
        **kwds
    ):
        """Calculate detachment limited erosion rate on nodes using the shear
        stress equation, solved using the depth slope product.

        Landlab component that generalizes the detachment limited erosion
        equation, primarily to be coupled to the the Landlab OverlandFlow
        component.

        This component adjusts topographic elevation and is contained in the
        landlab.components.detachment_ltd_erosion folder.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab grid.
        k_e : float
            Erodibility parameter, (m^(1+a_exp)*s^(2*a_exp-1)/kg^a_exp)
        fluid_density : float, optional
            Density of fluid, default set to water density of 1000 kg / m^3
        g : float, optional
            Acceleration due to gravity (m/s^2).
        a_exp : float, optional
            exponent on shear stress, positive, unitless
        tau_crit : float, optional
            threshold for sediment movement, (kg/m/s^2)
        uplift_rate : float, optional
            uplift rate applied to the topographic surface, m/s
        """
        super(DepthSlopeProductErosion, self).__init__(grid, **kwds)

        self.a = a_exp
        self.g = g
        self.rho = fluid_density
        self.E = self._grid.zeros(at="node")
        self.uplift_rate = uplift_rate
        self.tau_crit = tau_crit
        self.k_e = k_e

        self.dz = self._grid.zeros(at="node")

    def erode(
        self,
        dt,
        elevs="topographic__elevation",
        depth="surface_water__depth",
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
        elevs : str, optional
            Name of the field that represents topographic elevation on nodes.
        depth : str, optional
            Name of the field that represents water depths on nodes.
        slope : str, optional
            Name of the field that represent topographic slope on each node.
        """
        try:
            S = self._grid.at_node[slope]
        except FieldError:
            raise ValueError("Slope field is missing!")

        try:
            h = self._grid.at_node[depth]
        except FieldError:
            raise ValueError("Depth field is missing!")

        self.tau = self.rho * self.g * h * S

        greater_than_tc, = np.where(self.tau >= self.tau_crit)
        less_than_tc, = np.where(self.tau < self.tau_crit)

        self.E[less_than_tc] = 0.0

        self.E[greater_than_tc] = self.k_e * (
            (self.tau[greater_than_tc] ** self.a) - (self.tau_crit ** self.a)
        )

        self.E[self.E < 0.0] = 0.0

        self.dz = (self.uplift_rate - self.E) * dt

        self._grid["node"][elevs] += self.dz

    def run_one_step(
        self,
        dt,
        elevs="topographic__elevation",
        depth="surface_water__depth",
        slope="topographic__slope",
    ):

        self.erode(dt=dt, elevs=elevs, depth=depth, slope=slope)
