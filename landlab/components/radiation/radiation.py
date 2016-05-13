"""Landlab component that simulates relative incidence shortwave radiation
on sloped surface.

Landlab component that computes 1D and 2D total incident shortwave
radiation. This code also computes relative incidence shortwave radiation
compared to a flat surface. Ref: Flores-Cervantes et al., 2012

.. codeauthor:: Sai Nudurupati & Erkan Istanbulluoglu

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components import Radiation

Create a grid on which to calculate incident shortwave radiation

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> Radiation.input_var_names
('topographic__elevation',)

Check the units for the fields.

>>> Radiation.var_units('topographic__elevation')
'm'

Create the input fields.

>>> grid['node']['topographic__elevation'] = np.array([
...       0., 0., 0., 0.,
...       1., 1., 1., 1.,
...       2., 2., 2., 2.,
...       3., 4., 4., 3.,
...       4., 4., 4., 4.])

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> Radiation.var_help('topographic__elevation')
name: topographic__elevation
description:
  elevation of the ground surface relative to some datum
units: m
at: node
intent: in

Check the output variable names

>>> sorted(Radiation.output_var_names) # doctest: +NORMALIZE_WHITESPACE
['radiation__incoming_shortwave',
 'radiation__net_shortwave',
 'radiation__ratio_to_flat_surface']

Instantiate the 'Radiation' component to work on this grid, and run it.

>>> rad = Radiation(grid)

Run the *update* method to update output variables with current time
>>> current_time = 0.5
>>> rad.update(current_time)

>>> rad.grid.at_cell['radiation__ratio_to_flat_surface']
array([ 0.38488566,  0.38488566,  0.33309785,  0.33309785,  0.37381705,
        0.37381705])

>>> rad.grid.at_cell['radiation__incoming_shortwave']
array([ 398.33664988,  398.33664988,  344.73895668,  344.73895668,
        386.88120966,  386.88120966])
"""

from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import numpy as np


_VALID_METHODS = set(['Grid'])


def _assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)


class Radiation(Component):

    """Compute 1D and 2D total incident shortwave radiation.

    Landlab component that computes 1D and 2D total incident shortwave
    radiation. This code also computes relative incidence shortwave radiation
    compared to a flat surface.

    Construction::

        Radiation(grid, method='Grid', cloudiness=0.2, latitude=34., 
                  albedo=0.2, solar_constant=1366.67, 
                  clearsky_turbidity=2., opt_airmass=0.)

    Parameters
    ----------
    grid: RasterModelGrid
        A grid.
    method: {'Grid'}, optional
        Currently, only default is available.
    cloudiness: float, optional
        Cloudiness.
    latitude: float, optional
        Latitude (radians).
    albedo: float, optional
        Albedo.
    solar_constant: float, optional
        Solar Constant (W/m^2).
    clearsky_turbidity: float, optional
        Clear sky turbidity.
    opt_airmass: float, optional
        Optical air mass.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import Radiation
    >>> import numpy as np

    >>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    >>> rad = Radiation(grid)
    >>> rad.name
    'Radiation'
    >>> rad.input_var_names
    ('topographic__elevation',)
    >>> sorted(rad.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['radiation__incoming_shortwave',
     'radiation__net_shortwave',
     'radiation__ratio_to_flat_surface']
    >>> sorted(rad.units) # doctest: +NORMALIZE_WHITESPACE
    [('radiation__incoming_shortwave', 'W/m^2'),
     ('radiation__net_shortwave', 'W/m^2'),
     ('radiation__ratio_to_flat_surface', 'None'),
     ('topographic__elevation', 'm')]

    >>> rad.grid.number_of_node_rows
    5
    >>> rad.grid.number_of_node_columns
    4
    >>> rad.grid is grid
    True
    >>> np.all(grid.at_cell['radiation__ratio_to_flat_surface'] == 0.)
    True
    >>> np.all(grid.at_node['topographic__elevation'] == 0.)
    True

    >>> grid['node']['topographic__elevation'] = np.array([
    ...       0., 0., 0., 0.,
    ...       1., 1., 1., 1.,
    ...       2., 2., 2., 2.,
    ...       3., 4., 4., 3.,
    ...       4., 4., 4., 4.])
    >>> current_time = 0.5
    >>> rad.update(current_time)
    >>> np.all(grid.at_cell['radiation__ratio_to_flat_surface'] == 0.)
    False
    """

    _name = 'Radiation'

    _input_var_names = (
        'topographic__elevation',
    )

    _output_var_names = (
        'radiation__incoming_shortwave',
        'radiation__ratio_to_flat_surface',
        'radiation__net_shortwave',
    )

    _var_units = {
        'topographic__elevation': 'm',
        'radiation__incoming_shortwave': 'W/m^2',
        'radiation__ratio_to_flat_surface': 'None',
        'radiation__net_shortwave': 'W/m^2',
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'radiation__incoming_shortwave': 'cell',
        'radiation__ratio_to_flat_surface': 'cell',
        'radiation__net_shortwave': 'cell',
    }

    _var_doc = {
        'topographic__elevation':
            'elevation of the ground surface relative to some datum',
        'radiation__incoming_shortwave':
            'total incident shortwave radiation over the time step',
        'radiation__ratio_to_flat_surface':
            'ratio of total incident shortwave radiation on sloped surface \
             to flat surface',
        'radiation__net_shortwave':
            'net incident shortwave radiation over the time step',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, method='Grid', cloudiness=0.2,
                 latitude=34., albedo=0.2, solar_constant=1366.67,
                 clearsky_turbidity=2., opt_airmass=0., **kwds):
        """
        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        method : {'Grid'}, optional
            Currently, only default is available.
        cloudiness: float, optional
            Cloudiness.
        latitude: float, optional
            Latitude (Radians).
        albedo: float, optional
            Albedo.
        solar_constant: float, optional
            Solar Constant (W/m^2).
        clearsky_turbidity: float, optional
            Clear sky turbidity.
        opt_airmass: float, optional
            Optical air mass.
        """

        self._method = method
        self._N = cloudiness
        self._latitude = latitude
        self._A = albedo
        self._Io = solar_constant
        self._n = clearsky_turbidity
        self._m = opt_airmass

        _assert_method_is_valid(self._method)

        super(Radiation, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros(name, at='node', units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros(name, at='cell', units=self._var_units[name])

        if 'Slope' not in self.grid.at_cell:
            self.grid.add_zeros('Slope', at='cell', units='radians')

        if 'Aspect' not in self.grid.at_cell:
            self.grid.add_zeros('Aspect', at='cell', units='radians')

        self._nodal_values = self.grid['node']
        self._cell_values = self.grid['cell']
        self._slope, self._aspect = \
            grid.calculate_slope_aspect_at_nodes_burrough(
                vals='topographic__elevation')
#        self._slope = grid.calc_slope_of_node( \
#                                elevs = 'topographic__elevation')
#        self._aspect =
        self._cell_values['Slope'] = self._slope
        self._cell_values['Aspect'] = self._aspect

    def update(self, current_time, hour=12., **kwds):
        """Update fields with current loading conditions.

        Parameters
        ----------
        current_time: float
              Current time (years).
        hour: float, optional
              Hour of the day.
        """
        self._t = hour
        self._radf = self._cell_values['radiation__ratio_to_flat_surface']
        self._Rs = self._cell_values['radiation__incoming_shortwave']
        self._Rnet = self._cell_values['radiation__net_shortwave']

        self._julian = np.floor(
            (current_time - np.floor(current_time)) * 365.25) # Julian day

        self._phi = np.radians(self._latitude) # Latitude in Radians

        self._delta = 23.45 * np.radians(
            np.cos(2*np.pi / 365 * (172 - self._julian))) # Declination angle

        self._tau = (self._t + 12.0) * np.pi / 12.0 # Hour angle

        self._alpha = np.arcsin(np.sin(self._delta) * np.sin(self._phi) +
                                np.cos(self._delta) * np.cos(self._phi) *
                                np.cos(self._tau)) # Solar Altitude

        if self._alpha <= 0.25 * np.pi / 180.0: # If altitude is -ve,
            self._alpha = 0.25 * np.pi / 180.0 # sun is beyond the horizon

        self._Rgl = (self._Io * np.exp((-1) * self._n * (
            0.128 - 0.054 * np.log10(1. / np.sin(self._alpha)))*(
                1. / np.sin(self._alpha))))
        # Counting for Albedo, Cloudiness and Atmospheric turbidity

        self._phisun = (np.arctan(- np.sin(self._tau) / (np.tan(self._delta) *
                        np.cos(self._phi) - np.sin(self._phi) *
                        np.cos(self._tau)))) # Sun's Azhimuth

        if (self._phisun >= 0 and - np.sin(self._tau) <= 0):
            self._phisun = self._phisun + np.pi

        elif (self._phisun <= 0 and - np.sin(self._tau) >= 0):
            self._phisun = self._phisun + np.pi

        self._flat = (np.cos(np.arctan(0)) * np.sin(self._alpha) +
                      np.sin(np.arctan(0)) * np.cos(self._alpha) *
                      np.cos(self._phisun - 0)) # flat surface reference

        self._Rsflat = self._Rgl * self._flat
        # flat surface total incoming shortwave radiation

        self._Rnetflat = ((1 - self._A) * (1 - 0.65 * (self._N ** 2)) *
                          self._Rsflat)
        # flat surface Net incoming shortwave radiation

        self._sloped = (np.cos(self._slope) * np.sin(self._alpha) +
                        np.sin(self._slope) * np.cos(self._alpha) *
                        np.cos(self._phisun - self._aspect))

        self._radf = self._sloped / self._flat

        self._radf[self._radf <= 0.] = 0.
        self._radf[self._radf > 6.] = 6.

        self._Rs = self._Rsflat * self._radf
        # Sloped surface Toatl Incoming Shortwave Radn
        self._Rnet = self._Rnetflat * self._radf

        self._cell_values['radiation__ratio_to_flat_surface'] = self._radf
        self._cell_values['radiation__incoming_shortwave'] = self._Rs
        self._cell_values['radiation__net_shortwave'] = self._Rnet
