#!/usr/bin/env python

import numpy as np

from landlab import Component
from .funcs import get_flexure_parameter


_VALID_METHODS = set(['airy', 'flexure'])


def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)


class FlexureComponent(Component):
    """
    Landlab component that implements a 1 and 2D lithospheric flexure
    model.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flexure import FlexureComponent
    >>> grid = RasterModelGrid(5, 4, 1.e4)
    >>> flex = FlexureComponent(grid)
    >>> flex.name
    'Flexure'
    >>> sorted(flex.input_var_names)
    ['lithosphere__elevation', 'lithosphere__overlying_pressure', 'planet_surface_sediment__deposition_increment']
    >>> sorted(flex.output_var_names)
    ['lithosphere__elevation', 'lithosphere__elevation_increment']
    >>> for var in sorted(flex.units): flex.units[var]
    'm'
    'm'
    'Pa'
    'm'

    >>> flex.grid.number_of_node_rows
    5
    >>> flex.grid.number_of_node_columns
    4
    >>> flex.grid is grid
    True

    >>> np.all(grid.at_node['lithosphere__elevation'] == 0.)
    True

    >>> np.all(grid.at_node['lithosphere__elevation'] == 0.)
    True
    >>> np.all(grid.at_node['lithosphere__overlying_pressure'] == 0.)
    True
    >>> flex.update()
    >>> np.all(grid.at_node['lithosphere__elevation_increment'] == 0.)
    True

    >>> load = grid.at_node['lithosphere__overlying_pressure']
    >>> load[4] = 1e9
    >>> dz = grid.at_node['lithosphere__elevation_increment']
    >>> np.all(dz == 0.)
    True

    >>> flex.update()
    >>> np.all(grid.at_node['lithosphere__elevation_increment'] == 0.)
    False
    """

    ###############THIS IS THE STANDARD DECLARATIONS EXAMPLE.
    #YOU NEED TO DECLARE ALL THE FOLLOWING LIKE THIS HERE IN YOUR COMPONENT!!
    _name = 'Flexure'
    #...give the component a name

    _input_var_names = set([
        'lithosphere__overlying_pressure',
        'lithosphere__elevation',
        'planet_surface_sediment__deposition_increment',
    ])
    #...the component requires these values to do its calculation

    _output_var_names = set([
        'lithosphere__elevation_increment',
        'lithosphere__elevation',
    ])
    #...the component modifies these values

    # the next three dictionaries must each have the union of the two var_name
    # sets above as their keys
    _var_units = {
        'lithosphere__overlying_pressure': 'Pa',
        'lithosphere__elevation': 'm',
        'lithosphere__elevation_increment': 'm',
        'planet_surface_sediment__deposition_increment': 'm',
    }
    #...the units for each field. In the future, there may be unit casting
    # required, as there's nothing stopping you passing in the wrong units!!!

    _var_mapping = {
        'lithosphere__overlying_pressure': 'node',
        'lithosphere__elevation': 'node',
        'lithosphere__elevation_increment': 'node',
        'planet_surface_sediment__deposition_increment': 'node',
    }
    #...the grid centering of each name. Note you CANNOT and SHOULD NOT TRY
    # TO define the same field name but with different centering within one
    # component; you'll need to distinguish the names.

    _var_doc = {
        'lithosphere__overlying_pressure': 'The pressure at the base of the lithosphere',
        'lithosphere__elevation': 'The elevation of the top of the lithosphere, i.e., the land surface',
        'lithosphere__elevation_increment': 'The change in elevation of the top of the lithosphere (the land surface) in one timestep',
        'planet_surface_sediment__deposition_increment': 'The amount of sediment deposited at the land surface in one timestep',
    }
    #...give the names a short description. [We still need to work out how
    # we deal with same name, but different description between components]

    ###lastly, ---> 1. did you remember to import and inherit from Component?
    #do the import up top ->
    #from landlab import Component
    #& in the class declaration ->
    #   class FlexureComponent(Component):
    #        ...
    ### ---> 2. Make sure that self._grid is an alias for the grid after
    # initialization, below
    ################################

    def __init__(self, grid, **kwds):
        self._eet = kwds.pop('eet', 65000.)
        self._youngs = kwds.pop('youngs', 7e10)
        self._method = kwds.pop('method', 'airy')
        self._grid = grid

        assert_method_is_valid(self._method)

        super(FlexureComponent, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        self._last_load = self.grid.field_values('node', 'lithosphere__overlying_pressure').copy()

        self._nodal_values = self.grid['node']

        self._r = self._set_kei_func_grid()

    def _set_kei_func_grid(self):
        from scipy.special import kei

        alpha = get_flexure_parameter(self._eet, self._youngs, 2)
        dx, dy = np.meshgrid(
            np.arange(self._grid.number_of_node_columns) * self._grid.dx,
            np.arange(self._grid.number_of_node_rows) * self._grid.dy)

        return kei(np.sqrt(dx ** 2 + dy ** 2) / alpha)

    def update(self, n_procs=1):
        elevation = self._nodal_values['lithosphere__elevation']
        load = self._nodal_values['lithosphere__overlying_pressure']
        deflection = self._nodal_values['lithosphere__elevation_increment']
        deposition = self._nodal_values['planet_surface_sediment__deposition_increment']

        new_load = ((load - self._last_load) +
                    (deposition * 2650. * 9.81).flat)

        self._last_load = load.copy()

        deflection.fill(0.)

        if self._method == 'airy':
            deflection[:] = new_load / (3300. * 9.81)
        else:
            self.subside_loads(new_load, deflection=deflection,
                               n_procs=n_procs)
        np.subtract(elevation, deflection, out=elevation)

    def subside_loads(self, loads, deflection=None, n_procs=1):
        if deflection is None:
            deflection = np.empty(self.shape, dtype=np.float)

        from .cfuncs import subside_grid_in_parallel

        w = deflection.reshape(self._grid.shape)
        load = loads.reshape(self._grid.shape)
        alpha = get_flexure_parameter(self._eet, self._youngs, 2)

        subside_grid_in_parallel(w, load, self._r, alpha, n_procs)

        return deflection


if __name__ == "__main__":
    import doctest
    doctest.testmod()
