"""Landlab component that simulates overland flow.

This component simulates overland flow using the 2-D numerical model of
shallow-water flow over topography using the de Almeida et al., 2012
algorithm for storage-cell inundation modeling.

.. codeauthor:: Jordan Adams

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.overland_flow import OverlandFlow

Create a grid on which to calculate overland flow.

>>> grid = RasterModelGrid((4, 5))

The grid will need some data to provide the overland flow component. To
check the names of the fields that provide input to the overland flow
component use the *input_var_names* class property.

>>> OverlandFlow.input_var_names
('surface_water__depth', 'topographic__elevation')

Create fields of data for each of these input variables.

>>> grid.at_node["topographic__elevation"] = [
...     [0.0, 0.0, 0.0, 0.0, 0.0],
...     [1.0, 1.0, 1.0, 1.0, 1.0],
...     [2.0, 2.0, 2.0, 2.0, 2.0],
...     [3.0, 3.0, 3.0, 3.0, 3.0],
... ]
>>> grid.at_node["surface_water__depth"] = [
...     [0.0, 0.0, 0.0, 0.0, 0.0],
...     [0.0, 0.0, 0.0, 0.0, 0.0],
...     [0.0, 0.0, 0.0, 0.0, 0.0],
...     [0.1, 0.1, 0.1, 0.1, 0.1],
... ]

Instantiate the `OverlandFlow` component to work on this grid, and run it.

>>> of = OverlandFlow(grid, steep_slopes=True)
>>> of.run_one_step()

After calculating the overland flow, new fields have been added to the
grid. Use the *output_var_names* property to see the names of the fields that
have been changed.

>>> of.output_var_names
('surface_water__depth', 'surface_water__discharge', 'water_surface__gradient')

The `surface_water__depth` field is defined at nodes.

>>> of.var_loc("surface_water__depth")
'node'
>>> grid.at_node["surface_water__depth"]
array([1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05,
       1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05,
       1.0000e-05, 2.0010e-02, 2.0010e-02, 2.0010e-02, 1.0000e-05,
       1.0001e-01, 1.0001e-01, 1.0001e-01, 1.0001e-01, 1.0001e-01])

The `surface_water__discharge` field is defined at links. Because our initial
topography was a dipping plane, there is no water discharge in the horizontal
direction, only toward the bottom of the grid.

>>> of.var_loc("surface_water__discharge")
'link'
>>> q = grid.at_link["surface_water__discharge"]
>>> np.all(q[grid.horizontal_links] == 0.0)
True
>>> np.all(q[grid.vertical_links] <= 0.0)
True

The *water_surface__gradient* is also defined at links.

>>> of.var_loc("water_surface__gradient")
'link'
>>> grid.at_link["water_surface__gradient"]
array([0. , 0. , 0. , 0. ,
       0. , 1. , 1. , 1. , 0. ,
       0. , 0. , 0. , 0. ,
       0. , 1. , 1. , 1. , 0. ,
       0. , 0. , 0. , 0. ,
       0. , 1.1, 1.1, 1.1, 0. ,
       0. , 0. , 0. , 0. ])
"""

import numpy as np
import scipy.constants
from line_profiler import profile

from landlab import Component
from landlab.components.overland_flow._neighbors_at_link import (
    adjust_supercritial_discharge,
)
from landlab.components.overland_flow._neighbors_at_link import (
    adjust_unstable_discharge,
)
from landlab.components.overland_flow._neighbors_at_link import (
    calc_bates_flow_height_at_some_links,
)
from landlab.components.overland_flow._neighbors_at_link import (
    calc_discharge_at_some_links,
)
from landlab.components.overland_flow._neighbors_at_link import update_water_depths
from landlab.components.overland_flow._neighbors_at_link import (
    weighted_mean_of_parallel_links,
)
from landlab.grid.linkstatus import is_inactive_link


class OverlandFlow(Component):
    """Simulate overland flow using de Almeida approximations.

    Landlab component that simulates overland flow using the de Almeida
    et al., 2012 approximations of the 1D shallow water equations to be used
    for 2D flood inundation modeling.

    This component calculates discharge, depth and shear stress after some
    precipitation event across any raster grid. Default input file is named
    "overland_flow_input.txt' and is contained in the
    landlab.components.overland_flow folder.

    The primary method of this class is :func:`run_one_step`.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Adams, J., Gasparini, N., Hobley, D., Tucker, G., Hutton, E., Nudurupati,
    S., Istanbulluoglu, E. (2017). The Landlab v1. 0 OverlandFlow component:
    a Python tool for computing shallow-water flow across watersheds.
    Geoscientific Model Development  10(4), 1645.
    https://dx.doi.org/10.5194/gmd-10-1645-2017

    **Additional References**

    de Almeida, G., Bates, P., Freer, J., Souvignet, M. (2012). Improving the
    stability of a simple formulation of the shallow water equations for 2-D
    flood modeling. Water Resources Research 48(5)
    https://dx.doi.org/10.1029/2011wr011570

    """

    _name = "OverlandFlow"

    _unit_agnostic = False

    _cite_as = """@article{adams2017landlab,
        title={The Landlab v1. 0 OverlandFlow component: a Python
            tool for computing shallow-water flow across watersheds},
        author={Adams, Jordan M and Gasparini, Nicole M and
            Hobley, Daniel EJ and Tucker, Gregory E and
            Hutton, Eric WH and Nudurupati, Sai S and
            Istanbulluoglu, Erkan},
        journal={Geoscientific Model Development},
        volume={10},
        number={4},
        pages={1645},
        year={2017},
        publisher={Copernicus GmbH}
        }
    """

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "link",
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
        "water_surface__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "link",
            "doc": "Downstream gradient of the water surface.",
        },
    }

    def __init__(
        self,
        grid,
        h_init=0.00001,
        alpha=0.7,
        mannings_n=0.03,
        g=scipy.constants.g,
        theta=0.8,
        rainfall_intensity=0.0,
        steep_slopes=False,
    ):
        """Create an overland flow component.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab grid.
        h_init : float, optional
            Thicknes of initial thin layer of water to prevent divide by zero
            errors (m).
        alpha : float, optional
            Time step coeffcient, described in Bates et al., 2010 and
            de Almeida et al., 2012.
        mannings_n : float, optional
            Manning's roughness coefficient.
        g : float, optional
            Acceleration due to gravity (m/s^2).
        theta : float, optional
            Weighting factor from de Almeida et al., 2012.
        rainfall_intensity : float or array of float, optional
            Rainfall intensity. Default is zero.
        steep_slopes : bool, optional
            Modify the algorithm to handle steeper slopes at the expense of
            speed. If model runs become unstable, consider setting to True.
        """
        super().__init__(grid)

        self._h_init = h_init
        self._alpha = alpha

        if isinstance(mannings_n, str):
            self._mannings_n = self._grid.at_link[mannings_n]
        else:
            self._mannings_n = np.broadcast_to(mannings_n, self._grid.number_of_links)

        self._g = g
        self._theta = theta
        self.rainfall_intensity = rainfall_intensity
        self._steep_slopes = steep_slopes

        # Now setting up fields at the links...
        # For water discharge
        if "surface_water__discharge" not in grid.at_link:
            grid.add_empty(
                "surface_water__discharge",
                at="link",
                units=self._info["surface_water__discharge"]["units"],
            )
        grid.at_link["surface_water__discharge"].fill(0.0)

        # For water depths calculated at links
        if "surface_water__depth" not in grid.at_link:
            grid.add_empty(
                "surface_water__depth",
                at="link",
                units=self._info["surface_water__depth"]["units"],
            )
        grid.at_link["surface_water__depth"].fill(self._h_init)
        grid.at_node["surface_water__depth"] += self._h_init

        # For water surface slopes at links
        if "water_surface__gradient" not in grid.at_link:
            grid.add_empty("water_surface__gradient", at="link")
        grid.at_link["water_surface__gradient"].fill(0.0)

        self._dt = None

        # When we instantiate the class we recognize that neighbors have not
        # been found. After the user either calls self.set_up_neighbor_array
        # or self.overland_flow this will be set to True. This is done so
        # that every iteration of self.overland_flow does NOT need to
        # reinitalize the neighbors and saves computation time.
        self._active_link_flag = False

    @property
    def h(self):
        """The depth of water at each node."""
        return self._grid.at_node["surface_water__depth"]

    @property
    def dt(self):
        """dt: Component timestep."""
        return self._dt

    @dt.setter
    def dt(self, dt):
        if dt <= 0:
            raise ValueError("timestep dt must be positive")

        self._dt = dt

    @property
    def rainfall_intensity(self):
        """rainfall_intensity: the rainfall rate [m/s]

        Must be positive.
        """
        return self._rainfall_intensity

    @rainfall_intensity.setter
    def rainfall_intensity(self, new_val):
        if np.any(new_val < 0.0):
            raise ValueError("Rainfall intensity must be positive")
        self._rainfall_intensity = np.broadcast_to(new_val, self.grid.number_of_nodes)

    def calc_time_step(self):
        """Calculate time step.

        Adaptive time stepper from Bates et al., 2010 and de Almeida et
        al., 2012
        """
        max_water_depth = np.amax(self._grid.at_node["surface_water__depth"])
        if max_water_depth <= 0.0:
            raise RuntimeError("no water")

        self._dt = self._alpha * self._grid.dx / np.sqrt(self._g * max_water_depth)

        return self._dt

    def set_up_active_link_arrays(self):
        self._link_is_inactive = is_inactive_link(
            self.grid.status_at_node[self.grid.nodes_at_link]
        )

        self._active_links = np.where(~self._link_is_inactive)[0]

        self._active_link_flag = True

    @profile
    def overland_flow(self, dt=None):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.
        """
        if dt is None:
            dt = self.calc_time_step()

        h_at_node = self._grid.at_node["surface_water__depth"]
        z_at_node = self._grid.at_node["topographic__elevation"]
        q_at_link = self._grid.at_link["surface_water__discharge"]
        h_at_link = self._grid.at_link["surface_water__depth"]
        water_surface_slope = self._grid.at_link["water_surface__gradient"]
        q_mean_at_link = self.grid.empty(at="link")
        q_at_node = self.grid.empty(at="node")

        core_nodes = self._grid.core_nodes

        if not self._active_link_flag:
            self.set_up_active_link_arrays()

        time_remaining = dt
        while time_remaining > 0.0:
            dt_local = min(self.calc_time_step(), time_remaining)
            time_remaining -= dt_local

            # Per Bates et al., 2010, this solution needs to find difference
            # between the highest water surface in the two cells and the
            # highest bed elevation
            calc_bates_flow_height_at_some_links(
                z_at_node,
                h_at_node,
                self.grid.nodes_at_link,
                self._active_links,
                h_at_link,
            )

            # Now we calculate the slope of the water surface elevation
            self._grid.calc_grad_at_link(h_at_node + z_at_node, out=water_surface_slope)

            q_at_link[self._link_is_inactive] = 0.0

            weighted_mean_of_parallel_links(
                self.grid.shape,
                self._theta,
                q_at_link,
                q_mean_at_link,
            )

            calc_discharge_at_some_links(
                q_at_link,
                q_mean_at_link,
                h_at_link,
                water_surface_slope,
                self._mannings_n,
                self._active_links,
                self._g,
                dt_local,
            )
            q_at_link[h_at_link <= 0] = 0.0

            if self._steep_slopes:
                # To prevent water from draining too fast for our time steps...
                # Our Froude number.
                froude = 1.0

                # Our two limiting factors, the froude number and courant number.
                # Looking a calculated q to be compared to our Fr number.
                # Where does our calculated q exceed the Froude number? If q
                # does exceed the Froude number, we are getting supercritical
                # flow and discharge needs to be reduced to maintain stability.
                adjust_supercritial_discharge(
                    q_at_link,
                    h_at_link,
                    self._active_links,
                    self._g,
                    froude,
                )

                # Looking at our calculated q and comparing it to Courant no.,
                # Where does our calculated q exceed the Courant number and
                # water depth divided amongst 4 links? If the calculated q
                # exceeds the Courant number and is greater than the water
                # depth divided by 4 links, we reduce discharge to maintain
                # stability.
                adjust_unstable_discharge(
                    q_at_link,
                    h_at_link,
                    self._active_links,
                    self._grid.dx,
                    dt_local,
                )

            # Once stability has been restored, we calculate the change in
            # water depths on all core nodes by finding the difference between
            # inputs (rainfall) and the inputs/outputs (flux divergence of
            # discharge)
            self._grid.calc_flux_div_at_node(q_at_link, out=q_at_node)

            update_water_depths(
                q_at_node,
                h_at_node,
                self._rainfall_intensity,
                core_nodes,
                dt_local,
            )

            # To prevent divide by zero errors, a minimum threshold water depth
            # must be maintained. To reduce mass imbalances, this is set to
            # find locations where water depth is smaller than h_init (default
            # is 0.001) and the new value is self._h_init * 10^-3. This was set
            # as it showed the smallest amount of mass creation in the grid
            # during testing.
            # if self._steep_slopes:
            #     h_at_node[h_at_node < self._h_init] = self._h_init * 1e-3

            # And reset our field values with the newest water depth and
            # discharge.
            # self._grid.at_node["surface_water__depth"][:] = h_at_node
            # self._grid.at_link["surface_water__discharge"][:] = q_at_link

    def run_one_step(self, dt=None):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.
        """
        self.overland_flow(dt=dt)

    def discharge_mapper(self, input_discharge, convert_to_volume=False):
        """Maps discharge value from links onto nodes.

        This method takes the discharge values on links and determines the
        links that are flowing INTO a given node. The fluxes moving INTO a
        given node are summed.

        This method ignores all flow moving OUT of a given node.

        This takes values from the OverlandFlow component (by default) in
        units of [L^2/T]. If the convert_to_cms flag is raised as True, this
        method converts discharge to units [L^3/T] - as of Aug 2016, only
        operates for square RasterModelGrid instances.

        The output array is of length grid.number_of_nodes and can be used
        with the Landlab imshow_grid plotter.

        Returns a numpy array (discharge_vals)
        """

        discharge_vals = np.zeros(self._grid.number_of_links)
        discharge_vals[:] = input_discharge[:]

        if convert_to_volume:
            discharge_vals *= self._grid.dx

        discharge_vals = (
            discharge_vals[self._grid.links_at_node] * self._grid.link_dirs_at_node
        )

        discharge_vals = discharge_vals.flatten()

        discharge_vals[np.where(discharge_vals < 0)] = 0.0

        discharge_vals = discharge_vals.reshape(self._grid.number_of_nodes, 4)

        discharge_vals = discharge_vals.sum(axis=1)

        return discharge_vals
