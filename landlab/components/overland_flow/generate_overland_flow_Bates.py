"""generate_overland_flow.py.

This component simulates overland flow using
the 2-D numerical model of shallow-water flow
over topography using the Bates et al. (2010)
algorithm for storage-cell inundation modeling.

Written by Jordan Adams, based on code written by Greg Tucker.

Last updated: April 21, 2016
"""

import numpy as np
import scipy.constants

from landlab import Component


class OverlandFlowBates(Component):
    """Simulate overland flow using Bates et al. (2010).

    Landlab component that simulates overland flow using the Bates et al.,
    (2010) approximations of the 1D shallow water equations to be used for 2D
    flood inundation modeling.

    This component calculates discharge, depth and shear stress after some
    precipitation event across any raster grid. Default input file is named
    "overland_flow_input.txt' and is contained in the
    landlab.components.overland_flow folder.

    Parameters
    ----------
    grid : RasterGridModel
        A grid.
    input_file : str
        Contains necessary and optional inputs. If not given, default input
        file is used.

        -  Manning's n is *required*.
        -  Storm duration is needed *if* rainfall_duration is not passed in the
           initialization
        -  Rainfall intensity is needed *if* rainfall_intensity is not passed
           in the initialization
        -  Model run time can be provided in initialization. If not it is set
           to the storm duration

    h_init : float, optional
        Some initial depth in the channels. Default = 0.001 m
    g : float, optional
        Gravitational acceleration, :math:`m / s^2`
    alpha : float, optional
        Non-dimensional time step factor from Bates et al., (2010)
    rho : integer, optional
        Density of water, :math:`kg / m^3`
    ten_thirds : float, optional
        Precalculated value of :math:`10 / 3` which is used in the
        implicit shallow water equation.

    Examples
    --------
    >>> DEM_name = "DEM_name.asc"
    >>> (rg, z) = read_esri_ascii(DEM_name)  # doctest: +SKIP
    >>> of = OverlandFlowBates(rg)  # doctest: +SKIP

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Bates, P., Horritt, M., Fewtrell, T. (2010). A simple inertial formulation
    of the shallow water equations for efficient two-dimensional flood
    inundation modelling Journal of Hydrology  387(1-2), 33-45.
    https://dx.doi.org/10.1016/j.jhydrol.2010.03.027

    """

    _name = "OverlandFlowBates"

    _unit_agnostic = False

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
    }

    def __init__(
        self,
        grid,
        h_init=0.00001,
        alpha=0.7,
        mannings_n=0.03,
        g=scipy.constants.g,
        rainfall_intensity=0.0,
    ):
        super().__init__(grid)

        # First we copy our grid

        self._h_init = h_init
        self._alpha = alpha
        self._mannings_n = mannings_n
        self._g = g
        self._rainfall_intensity = rainfall_intensity

        # Now setting up fields at the links...
        # For water discharge

        self._surface_water__discharge = grid.add_zeros(
            "surface_water__discharge",
            at="link",
            units=self._info["surface_water__discharge"]["units"],
        )

        # Pre-calculated values included for speed.
        self._ten_thirds = 10.0 / 3.0
        self._mannings_n_squared = self._mannings_n * self._mannings_n

        # Start time of simulation is at 1.0 s
        self._elapsed_time = 1.0

        # Assigning a class variable to the water depth field and adding the
        # initial thin water depth
        self._h = self._grid["node"]["surface_water__depth"] = (
            self._grid["node"]["surface_water__depth"] + self._h_init
        )

        # Assigning a class variable to the water discharge field.
        self._q = self._grid["link"]["surface_water__discharge"]

        # Assiging a class variable to the elevation field.
        self._z = self._grid.at_node["topographic__elevation"]

    @property
    def surface_water__discharge(self):
        """The discharge of water on active links."""
        return self._surface_water__discharge

    @property
    def h(self):
        """The depth of water at each node."""
        return self._h

    @property
    def dt(self):
        """dt: component timestep."""
        return self._dt

    @dt.setter
    def dt(self, dt):
        assert dt > 0
        self._dt = dt

    def calc_time_step(self):
        # Adaptive time stepper from Bates et al., 2010 and de Almeida et al.,
        # 2012
        self._dt = (
            self._alpha
            * self._grid.dx
            / np.sqrt(self._g * np.amax(self._grid.at_node["surface_water__depth"]))
        )

        return self._dt

    def overland_flow(self, dt=None):
        """For one time step, this generates 'overland flow' across a given
        grid by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every
        node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.


        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        dt : float, optional
            Time step. Either set when called or the component will do it for
            you.
        """

        # If no dt is provided, one will be calculated using
        # self._gear_time_step()
        if dt is None:
            self.calc_time_step()

        # In case another component has added data to the fields, we just reset
        # our water depths, topographic elevations and water discharge
        # variables to the fields.
        # self._h = self._grid['node']['surface_water__depth']
        self._z = self._grid["node"]["topographic__elevation"]
        self._q = self._grid["link"]["surface_water__discharge"]

        # Here we identify the core nodes and active link ids for later use.
        self._core_nodes = self._grid.core_nodes
        self._active_links = self._grid.active_links

        # Per Bates et al., 2010, this solution needs to find the difference
        # between the highest water surface in the two cells and the highest
        # bed elevation
        zmax = self._grid.map_max_of_link_nodes_to_link(self._z)
        w = self._h + self._z
        wmax = self._grid.map_max_of_link_nodes_to_link(w)
        hflow = wmax[self._grid.active_links] - zmax[self._grid.active_links]

        # Now we calculate the slope of the water surface elevation at active
        # links
        water_surface_slope = self._grid.calc_grad_at_link(w)[self._grid.active_links]

        # Here we calculate discharge at all active links using Eq. 11 from
        # Bates et al., 2010
        self._q[self._active_links] = (
            self._q[self._active_links]
            - self._g * hflow * self._dt * water_surface_slope
        ) / (
            1.0
            + self._g
            * hflow
            * self._dt
            * self._mannings_n_squared
            * abs(self._q[self._active_links])
            / hflow**self._ten_thirds
        )

        # Update our water depths
        dhdt = self._rainfall_intensity - self._grid.calc_flux_div_at_node(self._q)

        self._h[self._core_nodes] = (
            self._h[self._core_nodes] + dhdt[self._core_nodes] * self._dt
        )

        # And reset our field values with the newest water depth and discharge.
        self._grid.at_node["surface_water__depth"] = self._h
        self._grid.at_link["surface_water__discharge"] = self._q
