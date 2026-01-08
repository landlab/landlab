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
from numpy.typing import NDArray

from landlab import Component
from landlab.components.overland_flow._calc import adjust_supercritical_discharge
from landlab.components.overland_flow._calc import adjust_unstable_discharge
from landlab.components.overland_flow._calc import calc_bates_flow_height
from landlab.components.overland_flow._calc import calc_discharge_at_links
from landlab.components.overland_flow._calc import calc_grad_at_link
from landlab.components.overland_flow._calc import calc_weighted_mean_of_parallel_links
from landlab.components.overland_flow._calc import zero_out_dry_links
from landlab.components.overland_flow._links import active_link_ids
from landlab.components.overland_flow._links import horizontal_active_link_ids
from landlab.components.overland_flow._links import horizontal_east_link_neighbor
from landlab.components.overland_flow._links import horizontal_west_link_neighbor
from landlab.components.overland_flow._links import is_horizontal_link
from landlab.components.overland_flow._links import is_vertical_link
from landlab.components.overland_flow._links import nth_horizontal_link
from landlab.components.overland_flow._links import nth_vertical_link
from landlab.components.overland_flow._links import vertical_active_link_ids
from landlab.components.overland_flow._links import vertical_north_link_neighbor
from landlab.components.overland_flow._links import vertical_south_link_neighbor
from landlab.core._validate import require_between
from landlab.core._validate import require_nonnegative
from landlab.core._validate import require_positive
from landlab.core.errors import Error
from landlab.grid.linkstatus import LinkStatus

_SEVEN_OVER_THREE = 7.0 / 3.0


class NoWaterError(Error):
    pass


def _active_links_at_node(grid, *args):
    """_active_links_at_node([node_ids]) Active links of a node.

    .. note::

        This function returns links that are in *clockwise* order,
        rather than the standard *counterclockwise* ordering that
        landlab uses everywhere else.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    node_ids : int or list of ints
        ID(s) of node(s) for which to find connected active links

    Returns
    -------
    (4, N) ndarray
        The ids of active links attached to grid nodes with
        *node_ids*. If *node_ids* is not given, return links for all of the
        nodes in the grid. Link ids are listed in clockwise order starting
        with the south link. Diagonal links are never returned.

    Examples
    --------
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 4))
    >>> grid.links_at_node[5]
    array([ 8, 11,  7,  4])
    >>> _active_links_at_node(grid, (5, 6))
    array([[ 4,  5],
           [ 7,  8],
           [11, 12],
           [ 8,  9]])
    >>> _active_links_at_node(grid)
    array([[-1, -1, -1, -1, -1,  4,  5, -1, -1, 11, 12, -1],
           [-1, -1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1],
           [-1,  4,  5, -1, -1, 11, 12, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1, -1]])

    :meta landlab: deprecated, info-link, info-node
    """
    active_links_at_node = grid.links_at_node.copy()
    active_links_at_node[grid.active_link_dirs_at_node == 0] = -1
    active_links_at_node = active_links_at_node[:, (3, 2, 1, 0)]

    if len(args) == 0:
        return active_links_at_node.T
    elif len(args) == 1:
        node_ids = np.broadcast_arrays(args[0])[0]
        return active_links_at_node[node_ids, :].T
    else:
        raise ValueError("only zero or one arguments accepted")


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
            "optional": True,
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
            "optional": True,
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
        mannings_n : array_like or str, optional
            Manning's roughness coefficient. If a `str`, use the corresponding
            *at-link* field of the provided grid.
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
        for var in ("surface_water__discharge", "water_surface__gradient"):
            at = self._info[var]["mapping"]
            units = self._info[var]["units"]
            if not grid.has_field(var, at=at):
                grid.add_empty(var, at=at, units=units)

        super().__init__(grid)

        self._h_init = require_nonnegative(h_init)
        self._alpha = require_positive(alpha)

        if isinstance(mannings_n, str):
            self._mannings_n = mannings_n
        else:
            self._mannings_n = np.require(
                np.broadcast_to(mannings_n, grid.number_of_links),
                requirements=("C_CONTIGUOUS",),
            )

        self._g = require_positive(g)
        self._theta = require_between(
            theta, 0.0, 1.0, inclusive_min=True, inclusive_max=True
        )
        self.rainfall_intensity = require_nonnegative(rainfall_intensity)
        self._steep_slopes = steep_slopes

        grid.at_link["surface_water__discharge"].fill(0.0)

        # For water depths calculated at links
        grid.at_node["surface_water__depth"] += self._h_init
        grid.at_link["water_surface__gradient"].fill(0.0)
        self._h_links = grid.empty(at="link")
        self._h_links.fill(self._h_init)

        # Start time of simulation is at 1.0 s
        self._elapsed_time = 1.0

        # When we instantiate the class we recognize that neighbors have not
        # been found. After the user either calls self.set_up_neighbor_array
        # or self.overland_flow this will be set to True. This is done so
        # that every iteration of self.overland_flow does NOT need to
        # reinitalize the neighbors and saves computation time.
        self._neighbor_flag = False

    @property
    def h(self):
        """The depth of water at each node."""
        return self._grid.at_node["surface_water__depth"]

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

        self._rainfall_intensity = new_val

    def calc_time_step(self):
        """Calculate time step.

        Adaptive time stepper from Bates et al., 2010 and de Almeida et
        al., 2012.

        Returns
        -------
        time_step : float
            A stable time step for the current landscape.

        Raises
        ------
        NoWaterError
            If there is no water on active nodes (i.e., the maximum water depth is
            zero), this exception is raised to indicate that a time step cannot be
            determined.
        """
        active_nodes = self._update_active_nodes()

        h: NDArray = self.grid.at_node["surface_water__depth"][active_nodes]
        if h.size == 0:
            raise NoWaterError("no active links, unable to determine time step")

        max_water_depth: float = np.max(h)
        if max_water_depth <= 0.0:
            raise NoWaterError("no water on landscape, unable to determine time step")

        return (
            self._alpha
            * min(self._grid.dx, self._grid.dy)
            / np.sqrt(self._g * max_water_depth)
        )

    def _update_active_nodes(self, *, clear_cache: bool = False) -> NDArray[np.bool_]:
        """Compute and optionally re-compute the active-node mask.

        Parameters
        ----------
        clear_cache : bool, optional
            If ``True`` clear any cached value and recompute the active
            node mask.

        Returns
        -------
        active_nodes : ndarray of bool, size (n_nodes,)
            The active nodes.
        """
        if clear_cache:
            self._cached_is_active_node = None

        if getattr(self, "_cached_is_active_node", None) is None:
            self._cached_is_active_node: NDArray[np.bool_] = np.any(
                self.grid.link_status_at_node == LinkStatus.ACTIVE, axis=1
            )

        return self._cached_is_active_node

    def set_up_neighbor_arrays(self):
        """Create and initialize link neighbor arrays.

        Set up arrays of neighboring horizontal and vertical links that
        are needed for the de Almeida solution.
        """
        # First we identify all active links

        active_ids = active_link_ids(self._grid.shape, self._grid.status_at_node)

        active_links_at_open_bdy = _active_links_at_node(
            self.grid, self.grid.open_boundary_nodes
        ).transpose()

        active_links_at_open_bdy = active_links_at_open_bdy[
            np.where(active_links_at_open_bdy > -1)
        ]

        # Find all horizontal active link ids
        _horizontal_active_link_ids = horizontal_active_link_ids(
            self._grid.shape, active_ids
        )

        # Find the *active* verical link ids
        _vertical_active_link_ids = vertical_active_link_ids(
            self._grid.shape, active_ids
        )

        vert_bdy_ids = active_links_at_open_bdy[
            is_vertical_link(self._grid.shape, active_links_at_open_bdy)
        ]

        vert_bdy_ids = nth_vertical_link(self._grid.shape, vert_bdy_ids)

        horiz_bdy_ids = active_links_at_open_bdy[
            is_horizontal_link(self._grid.shape, active_links_at_open_bdy)
        ]

        horiz_bdy_ids = nth_horizontal_link(self._grid.shape, horiz_bdy_ids)

        # Using the active vertical link ids we can find the north
        # and south vertical neighbors
        self._north_neighbors = vertical_north_link_neighbor(
            self._grid.shape, _vertical_active_link_ids
        )
        self._south_neighbors = vertical_south_link_neighbor(
            self._grid.shape, _vertical_active_link_ids
        )

        # Using the horizontal active link ids, we can find the west and
        # east neighbors
        self._west_neighbors = horizontal_west_link_neighbor(
            self._grid.shape, _horizontal_active_link_ids
        )
        self._east_neighbors = horizontal_east_link_neighbor(
            self._grid.shape, _horizontal_active_link_ids
        )

        # replace bdy condition links
        (ids,) = np.where(self._west_neighbors[horiz_bdy_ids] == -1)
        ids = horiz_bdy_ids[ids]
        self._west_neighbors[ids] = _horizontal_active_link_ids[ids]

        (ids,) = np.where(self._east_neighbors[horiz_bdy_ids] == -1)
        ids = horiz_bdy_ids[ids]
        self._east_neighbors[ids] = _horizontal_active_link_ids[ids]

        (ids,) = np.where(self._north_neighbors[vert_bdy_ids] == -1)
        ids = vert_bdy_ids[ids]
        self._north_neighbors[ids] = _vertical_active_link_ids[ids]

        (ids,) = np.where(self._south_neighbors[vert_bdy_ids] == -1)
        ids = vert_bdy_ids[ids]
        self._south_neighbors[ids] = _vertical_active_link_ids[ids]

        # Once the neighbor arrays are set up, we change the flag to True!
        self._neighbor_flag = True

    def overland_flow(self, dt=None):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.

        Parameters
        ----------
        dt : float, optional
            The duration over which to simulate overland flow. If not provided,
            the duration will be chosen to be the maximum time step that
            ensures stability.

        Returns
        -------
        elapsed : float
            The elapsed time that was simulated.
        """
        if dt is None:
            try:
                duration = self.calc_time_step()
            except NoWaterError as error:
                raise ValueError(
                    "no water on landscape and dt not provided,"
                    " unable to determine time step"
                ) from error
        else:
            duration = dt

        if isinstance(self._mannings_n, str):
            mannings_at_link = self.grid.at_link[self._mannings_n]
        else:
            mannings_at_link = self._mannings_n
        h_at_node = self._grid.at_node["surface_water__depth"]
        z_at_node = self._grid.at_node["topographic__elevation"]
        q_at_link = self._grid.at_link["surface_water__discharge"]
        h_at_link = self._h_links
        water_surface_slope = self._grid.at_link["water_surface__gradient"]
        q_mean_at_link = self.grid.empty(at="link")

        core_nodes = self._grid.core_nodes
        active_links = self._grid.active_links

        elapsed = 0.0
        while elapsed < duration:
            remaining = duration - elapsed
            try:
                dt_local = min(self.calc_time_step(), remaining)
            except NoWaterError:
                elapsed = duration
                break

            elapsed += dt_local

            # Per Bates et al., 2010, this solution needs to find difference
            # between the highest water surface in the two cells and the
            # highest bed elevation
            h_at_link = calc_bates_flow_height(
                z_at_node,
                h_at_node,
                nodes_at_link=self.grid.nodes_at_link,
                where=active_links,
                out=h_at_link,
            )

            h_at_link = np.clip(h_at_link, self._h_init, None, out=h_at_link)

            # Now we calculate the slope of the water surface elevation at
            # active links
            water_surface_slope = calc_grad_at_link(
                h_at_node + z_at_node,
                length_of_link=self.grid.length_of_link,
                nodes_at_link=self.grid.nodes_at_link,
                where=active_links,
                out=water_surface_slope,
            )

            q_at_link = zero_out_dry_links(h_at_link, where=active_links, out=q_at_link)

            q_mean_at_link.fill(0.0)
            q_mean_at_link = calc_weighted_mean_of_parallel_links(
                q_at_link,
                parallel_links_at_link=self.grid.parallel_links_at_link,
                status_at_link=self.grid.status_at_link,
                theta=self._theta,
                where=self._grid.active_links,
                out=q_mean_at_link,
            )

            q_at_link = calc_discharge_at_links(
                q_at_link,
                q_mean_at_link,
                h_at_link,
                water_surface_slope,
                mannings_at_link,
                g=self._g,
                dt=dt_local,
                where=self._grid.active_links,
                out=q_at_link,
            )

            if self._steep_slopes is True:
                # To prevent water from draining too fast for our time steps...
                # Our two limiting factors, the froude number and courant
                # number.

                # IDs where water discharge is positive...
                (positive_q,) = np.where(q_at_link > 0)

                # ... and negative.
                (negative_q,) = np.where(q_at_link < 0)

                # Where does our calculated q exceed the Froude number? If q
                # does exceed the Froude number, we are getting supercritical
                # flow and discharge needs to be reduced to maintain stability.
                q_at_link = adjust_supercritical_discharge(
                    q_at_link,
                    h_at_link,
                    g=self._g,
                    froude=1.0,
                    where=active_links,
                    out=q_at_link,
                )

                # Where does our calculated q exceed the Courant number and
                # water depth divided amongst 4 links? If the calculated q
                # exceeds the Courant number and is greater than the water
                # depth divided by 4 links, we reduce discharge to maintain
                # stability.
                q_at_link = adjust_unstable_discharge(
                    q_at_link,
                    h_at_link,
                    dx=self._grid.dx,
                    dt=dt_local,
                    where=active_links,
                    out=q_at_link,
                )

            # Once stability has been restored, we calculate the change in
            # water depths on all core nodes by finding the difference between
            # inputs (rainfall) and the inputs/outputs (flux divergence of
            # discharge)
            dhdt = self._rainfall_intensity - self._grid.calc_flux_div_at_node(
                q_at_link
            )

            # Updating our water depths...
            h_at_node[core_nodes] = h_at_node[core_nodes] + dhdt[core_nodes] * dt_local

            # To prevent divide by zero errors, a minimum threshold water depth
            # must be maintained. To reduce mass imbalances, this is set to
            # find locations where water depth is smaller than h_init (default
            # is 0.001) and the new value is self._h_init * 10^-3. This was set
            # as it showed the smallest amount of mass creation in the grid
            # during testing.
            if self._steep_slopes is True:
                h_at_node[h_at_node < self._h_init] = self._h_init * 10.0**-3

            # And reset our field values with the newest water depth and
            # discharge.
            self._grid.at_link["surface_water__discharge"][:] = q_at_link
            #
            #
            #  self._helper_q = self._grid.map_upwind_node_link_max_to_node(self._q)
            #  self._helper_s = self._grid.map_upwind_node_link_max_to_node(
            #      self._water_surface_slope)
            #
            #  self._helper_q = self._grid.map_max_of_link_nodes_to_link(self._helper_q)
            #  self._helper_s = self._grid.map_max_of_link_nodes_to_link(self._helper_s)
            #
            #  self._grid['link']['surface_water__discharge'][
            #     self._active_links_at_open_bdy] = self._helper_q[
            #     self._active_links_at_open_bdy]
            #
            #  self._grid['link']['water_surface__gradient'][
            #     self._active_links_at_open_bdy] = self._helper_s[
            #     self._active_links_at_open_bdy]
            # Update nodes near boundary locations - nodes adjacent to
            # boundaries may have discharge and water surface slopes
            # artifically reduced due to boundary effects. This step removes
            # those errors.

        return elapsed

    def run_one_step(self, dt=None):
        """Generate overland flow across a grid.

        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.

        Parameters
        ----------
        dt : float, optional
            The duration over which to simulate overland flow. If not provided,
            the duration will be chosen to be the maximum time step that
            ensures stability.

        Returns
        -------
        elapsed : float
            The elapsed time that was simulated.
        """
        return self.overland_flow(dt=dt)

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
