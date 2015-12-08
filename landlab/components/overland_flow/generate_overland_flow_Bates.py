#Embedded file name: /Users/Jordan/Documents/landlab/landlab/components/overland_flow/generate_overland_flow_Bates.py
""" generate_overland_flow.py

 This component simulates overland flow using
 the 2-D numerical model of shallow-water flow
 over topography using the Bates et al. (2010)
 algorithm for storage-cell inundation modeling.

Written by Jordan Adams, based on code written by Greg Tucker.
Last updated: July 17, 2015

"""
from landlab import Component, ModelParameterDictionary
import numpy as np
import os

class OverlandFlowBates(Component):
    """  Landlab component that simulates overland flow using the Bates et al., (2010) approximations
    of the 1D shallow water equations to be used for 2D flood inundation modeling.

    This component calculates discharge, depth and shear stress after some precipitation event across
    any raster grid. Default input file is named "overland_flow_input.txt' and is contained in the
    landlab.components.overland_flow folder.

        Inputs
        ------
        grid : Requires a RasterGridModel instance

        input_file : Contains necessary and optional inputs. If not given, default input file is used.
            - Manning's n is REQUIRED.
            - Storm duration is needed IF rainfall_duration is not passed in the initialization
            - Rainfall intensity is needed IF rainfall_intensity is not passed in the initialization
            - Model run time can be provided in initialization. If not it is set to the storm duration

        Constants
        ---------
        h_init : float
            Some initial depth in the channels. Default = 0.001 m
        g : float
            Gravitational acceleration, \x0crac{m}{s^2}
        alpha : float
            Non-dimensional time step factor from Bates et al., (2010)
        rho : integer
            Density of water, \x0crac{kg}{m^3}
        ten_thirds : float
            Precalculated value of \x0crac{10}{3} which is used in the implicit shallow water equation.


        >>> DEM_name = 'DEM_name.asc'
        >>> (rg, z) = read_esri_ascii(DEM_name) # doctest: +SKIP
        >>> of = OverlandFlowBates(rg) # doctest: +SKIP

    """
    _name = 'OverlandFlowBates'

    _input_var_names = set(['water_depth', 'topographic__elevation'])

    _output_var_names = set(['water_depth',
     'water_discharge',
     'shear_stress',
     'water_discharge_at_nodes',
     'water_surface_slope_at_nodes'])

    _var_units = {'water_depth': 'm',
     'water_discharge': 'm3/s',
     'shear_stress': 'Pa',
     'water_discharge_at_nodes': 'm3/s',
     'water_surface_slope_at_nodes': 'm/m',
     'topographic__elevation': 'm'}

    _var_mapping = {'water_depth': 'node',
     'topographic__elevtation': 'node',
     'water_discharge': 'active_link',
     'shear_stress': 'node',
     'water_discharge_at_nodes': 'node',
     'water_surface_slope_at_nodes': 'node'}

    _var_mapping = {'water_depth': 'The depth of water at each node.',
     'topographic__elevtation': 'The land surface elevation.',
     'water_discharge': 'The discharge of water on active links.',
     'shear_stress': 'The calculated shear stress at each node.',
     'water_discharge_at_nodes': 'The water discharge from surrounding links mapped onto nodes.',
     'water_surface_slope_at_nodes': 'The slope of the water surface at each node.'}

    def __init__(self, grid, input_file = None, **kwds):

        super(OverlandFlowBates, self).__init__(grid, **kwds)

        # First we copy our grid
        self._grid = grid

        # Then, we look for a input file...
        if input_file is not None:
            inputs = ModelParameterDictionary(input_file)
        else:
            print("No input file provided! Default file and default values will be used")
            _DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__), 'overland_flow_input.txt')
            input_file = _DEFAULT_INPUT_FILE
            inputs = ModelParameterDictionary(input_file)

        # And here we look to see what parameters are set in the input file.
        # If a parameter is not found, default values are set below.

        #This is an initial thin layer of water to prevent divide by zero errors
        try:
            self.h_init = inputs.read_float('h_init')
        except:
            self.h_init = 0.001

        # This is the time step coeffcient, described in Bates et al., 2010 and
        # de Almeida et al., 2012
        try:
            self.alpha = inputs.read_float('alpha')
        except:
            self.alpha = 0.7

        # Manning's roughness coefficient or Manning's n
        try:
            self.mannings_n = inputs.read_float('Mannings_n')
        except:
            self.mannings_n = 0.01

        # Rainfall intensity
        try:
            self.rainfall_intensity = inputs.read_float('rainfall_intensity')
        except:
            self.rainfall_intensity = 0.0

        try:
            self.g = inputs.read_float('g')
        except:
            self.g = 9.8

        # Setting up all fields found at nodes.
        for name in self._input_var_names:
            if name not in self._grid.at_node:
                self._grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self._grid.at_node:
                self._grid.add_zeros('node', name, units=self._var_units[name])

        # Now setting up fields at the links...
        # For water discharge
        self.water_discharge = grid.add_zeros('link', 'water_discharge', units=self._var_units['water_discharge'])

        # Pre-calculated values included for speed.
        self.ten_thirds = 10.0 / 3.0
        self.mannings_n_squared = self.mannings_n * self.mannings_n

        # Start time of simulation is at 1.0 s
        self.elapsed_time = 1.0

        # Assigning a class variable to the water depth field and adding the initial thin water depth
        self.h = self._grid['node']['water_depth'] = self._grid['node']['water_depth'] + self.h_init

        # Assigning a class variable to the water discharge field.
        self.q = self._grid['link']['water_discharge']

        # Assiging a class variable to the elevation field.
        self.z = self._grid.at_node['topographic__elevation']

    def gear_time_step(self, grid):

        # Adaptive time stepper from Bates et al., 2010 and de Almeida et al., 2012
        dt = self.alpha * self._grid.dx / np.sqrt(self.g * np.amax(self._grid.at_node['water_depth']))

        return dt

    def overland_flow(self, grid, dt = None, **kwds):
        """
        For one time step, this generates 'overland flow' across a given grid
        by calculating discharge at each node.

        Using the depth slope product, shear stress is calculated at every node.

        Outputs water depth, discharge and shear stress values through time at
        every point in the input grid.


        Inputs
        ------
        grid : Requires a RasterGridModel instance

        dt : either set when called or the fxn will do it for you.

        """

        # If no dt is provided, one will be calculated using self.gear_time_step()
        if dt is None:
            dt = self.gear_time_step(grid)

        # In case another component has added data to the fields, we just reset our
        # water depths, topographic elevations and water discharge variables to the fields.         self.h = self._grid['node']['water_depth']
        self.z = self._grid['node']['topographic__elevation']
        self.q = self._grid['link']['water_discharge']

        # Here we identify the core nodes and active link ids for later use.
        self.core_nodes = self._grid.core_nodes
        self.active_links = self._grid.active_links

        # Per Bates et al., 2010, this solution needs to find the difference between the highest
        # water surface in the two cells and the highest bed elevation
        zmax = self._grid.max_of_link_end_node_values(self.z)
        w = self.h + self.z
        wmax = self._grid.max_of_link_end_node_values(w)
        hflow = wmax - zmax

        # Now we calculate the slope of the water surface elevation at active links
        water_surface_slope = self._grid.calculate_gradients_at_active_links(w)

        # Here we calculate discharge at all active links using Eq. 11 from Bates et al., 2010
        self.q[self.active_links] = (self.q[self.active_links] - self.g * hflow * dt * water_surface_slope) / (1.0 + self.g * hflow * dt * self.mannings_n_squared * abs(self.q[self.active_links]) / hflow ** self.ten_thirds)

        # Update our water depths
        dhdt = self.rainfall_intensity - self._grid.calculate_flux_divergence_at_nodes(self.q[self.active_links])
        self.h[self.core_nodes] = self.h[self.core_nodes] + dhdt[self.core_nodes] * dt

        # And reset our field values with the newest water depth and discharge.
        self._grid.at_node['water_depth'] = self.h
        self._grid.at_link['water_discharge'] = self.q

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
