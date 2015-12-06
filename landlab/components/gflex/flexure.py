# -*- coding: utf-8 -*-
"""
This is a Landlab wrapper for A Wickert's gFlex flexure model (Wickert et al.,
submitted to Geoscientific Model Development). The most up-to-date version of
his code can be found at github.com/awickert/gFlex.

This Landlab wrapper will use a snapshot of that code, which YOU need to
install on your own machine.
A stable snapshot of gFlex is hosted on PyPI, which is the recommended version
to install.
If you have pip (the Python package install tool), simply run
'pip install gFlex' from a command prompt.
Alternatively, you can download and unpack the code (from github, or with PyPI,
pypi.python.org/pypi/gFlex/), then run 'python setup.py install'.

Created on Thu Feb 19 18:47:11 2015

@author: daniel.hobley (SiccarPoint @Github)

...following AW's run_in_script_2D.py.
"""
from __future__ import print_function

import numpy as np
import inspect
from landlab import RasterModelGrid, Component
from landlab import ModelParameterDictionary
from landlab import FieldError

try:
    import gflex
except ImportError:
    import warnings, sys

    warnings.warn("gFlex not installed.")
    print("""
To use the gFlex component you must have gFlex installed on your machine.
For installation instructions see gFlex on GitHub:

  https://github.com/awickert/gFlex
          """.strip(), file=sys.stderr)


class gFlex(Component):
    """
    This is the Landlab component controlling the gFlex code.

    Note that gFlex maintains its own internal version if the grid, but this
    should not affect performance.
    """
    _name = 'gFlex'

    _input_var_names = set(['surface_load__stress',
                        ])

    _output_var_names = set(['lithosphere__vertical_displacement',
                             'topographic__elevation',
                        ])

    _var_units = {'earth_material_load__magnitude_of_stress' : 'Pa',
                  'earth_material_load__x_positions' : 'm',
                  'earth_material_load__y_positions' : 'm',
                  'earth_material_load__force' : 'N',
                  'lithosphere__elastic_thickness' : 'm',
                  'lithosphere__vertical_displacement' : 'm',
                  }

    _var_mapping = {'earth_material_load__magnitude_of_stress' : 'node',
                  'earth_material_load__x_positions' : 'node',
                  'earth_material_load__y_positions' : 'node',
                  'earth_material_load__force' : 'node',
                  'lithosphere__elastic_thickness' : 'node',
                  'lithosphere__vertical_displacement' : 'node',
                  }

    _var_doc = {'earth_material_load__magnitude_of_stress' : 'Magnitude of stress exerted by surface load',
                  'earth_material_load__x_positions' : 'x position of any surface load',
                  'earth_material_load__y_positions' : 'y position of any surface load',
                  'earth_material_load__force' : 'Force exerted by surface load',
                  'lithosphere__elastic_thickness' : 'Elastic thickness of the lithosphere',
                  'lithosphere__vertical_displacement' : 'Vertical deflection of the surface and of the lithospheric plate',
                  }

    def __init__(self, grid, params):
        self.initialize(grid, params)

    def initialize(self, grid, params):
        """
        grid must be a landlab RasterModelGrid to use gFlex.

        params may be a dictionary of input parameters, or a standard Landlab
        input text file.

        gFlex requires:
        * E, Young's modulus
        * nu, Poisson's ratio
        * rho_m, the mantle density
        * rho_fill, the infill material density
        * Te, the elastic thickness. Can be scalar or an nnodes-long array.
            If you want an array, supply a grid field name where the data can
            be found.
        * BC_W, BC_E, BC_S, BC_N, strings describing the boundary conditions.
            Choose from ('Dirichlet0', '0Moment0Shear', 'Periodic').


        gFlex can take as options:
        * g, the gravitational acceleration (defaults to 9.81)

        gFlex takes as input fields:
        * surface_load__stress

        gFlex modifies/returns:
        * topographic__elevation (if it exists already)
        * lithosphere__vertical_displacement


        """
        assert RasterModelGrid in inspect.getmro(grid.__class__)
        BC_options = ('Dirichlet0', '0Moment0Shear', '0Slope0Shear', 'Periodic')

        if type(params) == str:
            input_dict = ModelParameterDictionary(params)
        else:
            assert type(params) == dict
            input_dict = params

        #instantiate the module:
        self.flex = gflex.F2D()
        flex = self.flex

        #set up the grid variables:
        self._grid = grid
        flex.dx = grid.dx
        flex.dy = grid.dy

        #we assume these properties are fixed in this relatively
        #straightforward implementation, but they can still be set if you want:
        try:
            flex.Method = input_dict['Method']
        except KeyError:
            flex.Method = 'FD'
        try:
            flex.PlateSolutionType = input_dict['PlateSolutionType']
        except KeyError:
            flex.PlateSolutionType = 'vWC1994'
        try:
            flex.Solver = input_dict['Solver']
        except KeyError:
            flex.Solver = 'direct'
        try:
            quiet = input_dict['Quiet']
        except KeyError:
            flex.Quiet = True
        else:
            flex.Quiet = bool(quiet)

        flex.E = float(input_dict['E'])
        flex.nu = float(input_dict['nu'])
        flex.rho_m = float(input_dict['rho_m'])
        flex.rho_fill = float(input_dict['rho_fill'])
        try:
            flex.g = input_dict['g']
        except KeyError:
            flex.g = 9.81
        flex.BC_W = input_dict['BC_W']
        flex.BC_E = input_dict['BC_E']
        flex.BC_S = input_dict['BC_S']
        flex.BC_N = input_dict['BC_N']
        for i in (flex.BC_E, flex.BC_W, flex.BC_N, flex.BC_S):
            assert i in BC_options

        Te_in = input_dict['Te']
        try:
            flex.Te = float(Te_in)
        except ValueError:
            flex.Te = grid.at_node[Te_in].view().reshape((grid.number_of_node_rows, grid.number_of_node_columns))
            self._input_var_names.add(Te_in)
            self._output_var_names.add(Te_in)

        #set up the link between surface load stresses in the gFlex component and the LL grid field:
        flex.qs = grid.at_node['surface_load__stress'].view().reshape((grid.number_of_node_rows, grid.number_of_node_columns))

        #create a holder for the "pre-flexure" state of the grid, to allow updating of elevs:
        self.pre_flex = np.zeros(grid.number_of_nodes, dtype=float)



    def flex_lithosphere(self, **kwds):
        """
        Executes (& finalizes, from the perspective of gFlex) the core method
        of gFlex. Note that flexure of the lithosphere proceeds to steady state
        in a single timestep.
        """
        #note kwds is redundant at the moment, but could be used subsequently for dynamic control over params
        self.flex.qs = self._grid.at_node['surface_load__stress'].view().reshape((self._grid.number_of_node_rows, self._grid.number_of_node_columns))
        self.flex.initialize()
        self.flex.run()
        self.flex.finalize()

        self._grid.at_node['lithosphere__vertical_displacement'] = self.flex.w.view().ravel()

        try:
            self._grid.at_node['topographic__elevation']
            #a topo exists...
        except FieldError:
            pass
        else:
            topo_diff = self._grid.at_node['lithosphere__vertical_displacement'] - self.pre_flex
            self._grid.at_node['topographic__elevation'] += topo_diff
            self.pre_flex += topo_diff
