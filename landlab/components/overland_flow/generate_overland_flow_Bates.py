""" generate_overland_flow.py

 This component simulates overland flow using
 the 2-D numerical model of shallow-water flow
 over topography using the Bates et al. (2010)
 algorithm for storage-cell inundation modeling.

Written by Jordan Adams, based on code written by Greg Tucker.

Last updated: April 21, 2016

"""
from landlab import Component, ModelParameterDictionary
import numpy as np
import os

class OverlandFlowBates(Component):
    u"""Simulate overland flow using Base et al. (2010).

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
    >>> DEM_name = 'DEM_name.asc'
    >>> (rg, z) = read_esri_ascii(DEM_name) # doctest: +SKIP
    >>> of = OverlandFlowBates(rg) # doctest: +SKIP
    """
    _name = 'OverlandFlowBates'

    _input_var_names = ('water__depth', 'topographic__elevation')

    _output_var_names = ('water__depth',
     'water__discharge',
     'water_surface__gradient')

    _var_units = {'water__depth': 'm',
     'water__discharge': 'm3/s',
     'water_surface__gradient': 'm/m',
     'topographic__elevation': 'm'}

    _var_mapping = {'water__depth': 'node',
     'topographic__elevtation': 'node',
     'water__discharge': 'active_link',
     'water_surface__gradient': 'node'}

    _var_mapping = {'water__depth': 'The depth of water at each node.',
     'topographic__elevtation': 'The land surface elevation.',
     'water__discharge': 'The discharge of water on active links.',
     'water_surface__gradient': 'The slope of the water surface at each node.'}

    def __init__(self, grid, h_init=0.00001, alpha=0.7,
                 mannings_n=0.03, g=9.81, rainfall_intensity=0.0,
                 **kwds):

        super(OverlandFlowBates, self).__init__(grid, **kwds)

        # First we copy our grid
        self._grid = grid

        self.h_init = h_init
        self.alpha = alpha
        self.mannings_n = mannings_n
        self.g = g
        self.rainfall_intensity = rainfall_intensity

        # Now setting up fields at the links...
        # For water discharge
        self.water__discharge = grid.add_zeros('link',
                'water__discharge', units=self._var_units['water__discharge'])

        # Pre-calculated values included for speed.
        self.ten_thirds = 10.0 / 3.0
        self.mannings_n_squared = self.mannings_n * self.mannings_n

        # Start time of simulation is at 1.0 s
        self.elapsed_time = 1.0

        # Assigning a class variable to the water depth field and adding the initial thin water depth
        self.h = self._grid['node']['water__depth'] = (
            self._grid['node']['water__depth'] + self.h_init)

        # Assigning a class variable to the water discharge field.
        self.q = self._grid['link']['water__discharge']

        # Assiging a class variable to the elevation field.
        self.z = self._grid.at_node['topographic__elevation']

    def calc_time_step(self):

        # Adaptive time stepper from Bates et al., 2010 and de Almeida et al., 2012
        self.dt = self.alpha * self._grid.dx / np.sqrt(self.g * np.amax(
            self._grid.at_node['water__depth']))

        return self.dt

    def overland_flow(self, dt = None, **kwds):
        """
        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every node.

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

        # If no dt is provided, one will be calculated using self.gear_time_step()
        if dt is None:
            self.calc_time_step()

        # In case another component has added data to the fields, we just reset our
        # water depths, topographic elevations and water discharge variables to the fields.         self.h = self._grid['node']['water__depth']
        self.z = self._grid['node']['topographic__elevation']
        self.q = self._grid['link']['water__discharge']

        # Here we identify the core nodes and active link ids for later use.
        self.core_nodes = self._grid.core_nodes
        self.active_links = self._grid.active_links

        # Per Bates et al., 2010, this solution needs to find the difference between the highest
        # water surface in the two cells and the highest bed elevation
        zmax = self._grid.map_max_of_link_nodes_to_link(self.z)
        w = self.h + self.z
        wmax = self._grid.map_max_of_link_nodes_to_link(w)
        hflow = wmax[self._grid.active_links] - zmax[self._grid.active_links]

        # Now we calculate the slope of the water surface elevation at active links
        water_surface_slope = (
            self._grid.calc_grad_at_link(w)[self._grid.active_links])

        # Here we calculate discharge at all active links using Eq. 11 from Bates et al., 2010
        self.q[self.active_links] = ((self.q[self.active_links] - self.g *
            hflow * self.dt * water_surface_slope) / (1.0 + self.g * hflow *
            self.dt * self.mannings_n_squared * abs(self.q[self.active_links])
            / hflow ** self.ten_thirds))

        # Update our water depths
        dhdt = (self.rainfall_intensity - self._grid.calc_flux_div_at_node(
            self.q))

        self.h[self.core_nodes] = (self.h[self.core_nodes] +
            dhdt[self.core_nodes] * self.dt)

        # And reset our field values with the newest water depth and discharge.
        self._grid.at_node['water__depth'] = self.h
        self._grid.at_link['water__discharge'] = self.q

    @property
    def input_var_names(self):
        return self._input_var_names

    @property
    def output_var_names(self):
        return self._output_var_names

    @property
    def var_units(self):
        return self._var_units

    @property
    def var_mapping(self):
        return self._var_mapping
