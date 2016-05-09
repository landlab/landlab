
#################################################################
##
##  'Field' concept is implemented for Radiation component.
##
##  Sai Nudurupati and Erkan Istanbulluoglu - 14May2014
#################################################################

from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import numpy as np

_VALID_METHODS = set(['Grid'])

def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)

class Radiation(Component):
    """
    Landlab component that computes 1D and 2D total incident shortwave
    radiation. This code also computes relative incidence shortwave radiation
    compared to a flat surface.

    Construction::
        Radiation(grid, method='Grid', cloudiness=0.2, latitude=34., \
                  albedo=0.2, solar_constant=1366.67, \
                  clearsky_turbidity=2., opt_airmass=0.)

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

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.radiation import Radiation
    >>> import numpy as np
    >>> grid = RasterModelGrid( 5, 4, 0.2 )
    >>> grid['node']['topographic__elevation'] = np.random.rand( grid.number_of_nodes ) * 1000
    >>> rad = Radiation( grid )
    >>> rad.name
    'Radiation'
    >>> current_time = 0.5
    >>> rad.update( current_time )

    >>> x = grid['cell']['topographic__total_shortwave_radiation']
    >>> isinstance(x, np.ndarray)
    True
    >>> x.shape == (6, )
    True

    >>> x = grid['cell']['topographic__radiation_factor']
    >>> isinstance(x, np.ndarray)
    True
    >>> x.shape == (6, )
    True
    """

    _name = 'Radiation'

    _input_var_names = set([
        'topographic__elevation',
    ])

    _output_var_names = set([
        'topographic__total_shortwave_radiation',
        'topographic__radiation_factor',
        'topographic__net_shortwave_radiation',
    ])

    _var_units = {
        'topographic__elevation' : 'm',
        'topographic__total_shortwave_radiation' : 'W/m^2',
        'topographic__radiation_factor' : 'None',
        'topographic__net_shortwave_radiation' : 'W/m^2',
    }
    
    _var_mapping = {
        'topographic__elevation' : 'node',
        'topographic__total_shortwave_radiation' : 'cell',
        'topographic__radiation_factor' : 'cell',
        'topographic__net_shortwave_radiation' : 'cell',
    }

    _var_doc = {
        'topographic__elevation' : 
            'elevation of the ground surface relative to some datum',
        'topographic__total_shortwave_radiation' : 
            'total incident shortwave radiation over the time step',
        'topographic__radiation_factor' : 
            'ratio of total incident shortwave radiation on sloped surface \
             to flat surface',
        'topographic__net_shortwave_radiation' : 
            'net incident shortwave radiation over the time step',
    }

    @use_file_name_or_kwds
    def __init__( self, grid, method='Grid', cloudiness=0.2, \
                  latitude=34., albedo=0.2, solar_constant=1366.67, \
                  clearsky_turbidity=2., opt_airmass=0., **kwds ):
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

        assert_method_is_valid(self._method)

        super(Radiation, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if not name in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        if not 'Slope' in self.grid.at_cell:
            self.grid.add_zeros('cell', 'Slope', units='radians' )

        if not 'Aspect' in self.grid.at_cell:
            self.grid.add_zeros('cell', 'Aspect', units='radians' )

        self._nodal_values = self.grid['node']
        self._cell_values = self.grid['cell']
        self._slope, self._aspect = \
            grid.calculate_slope_aspect_at_nodes_burrough( \
                vals='topographic__elevation')
#        self._slope = grid.calc_slope_of_node( \
#                                elevs = 'topographic__elevation')
#        self._aspect = 
        self._cell_values['Slope'] = self._slope
        self._cell_values['Aspect'] = self._aspect

    def update( self, current_time, **kwds ):
        """
        Update fields with current loading conditions.
        Parameters
        ----------
        current_time: float
              Current time (years).
        hour: float, optional
              Hour of the day. 
        """
        self._t = kwds.pop('Hour', 12.)
        self._radf = self._cell_values['topographic__radiation_factor']
        self._Rs = self._cell_values['topographic__total_shortwave_radiation']
        self._Rnet = self._cell_values['topographic__net_shortwave_radiation']

        self._julian = np.floor( ( current_time - np.floor( current_time ) )  \
                                  * 365.25 )                     # Julian day

        self._phi = np.radians(self._latitude)           # Latitude in Radians

        self._delta = 23.45 * np.radians(np.cos(2*np.pi/365                 \
                        * (172 - self._julian)))           # Declination angle


        self._tau = (self._t + 12.0) * np.pi/12.0                # Hour angle

        self._alpha = np.arcsin(np.sin(self._delta) * np.sin(self._phi)       \
                        + np.cos(self._delta) * np.cos(self._phi)             \
                        * np.cos(self._tau))                  # Solar Altitude

        if self._alpha <= 0.25 * np.pi/180.0:         # If altitude is -ve,
            self._alpha = 0.25 * np.pi/180.0       # sun is beyond the horizon

        self._Rgl = (self._Io*np.exp((-1) * self._n * (0.128 - 0.054 *        \
                    np.log10(1/np.sin(self._alpha)))*(1/np.sin(self._alpha))))
                    # Counting for Albedo, Cloudiness and Atmospheric turbidity

        self._phisun = np.arctan(-np.sin(self._tau)/(np.tan(self._delta)      \
                    *np.cos(self._phi) - np.sin(self._phi)                    \
                        *np.cos(self._tau)))
                                                              # Sun's Azhimuth

        if ( self._phisun >= 0 and -np.sin(self._tau) <= 0 ):
            self._phisun = self._phisun + np.pi

        elif ( self._phisun <= 0 and -np.sin(self._tau) >= 0 ):
            self._phisun = self._phisun + np.pi

        self._flat = np.cos(np.arctan(0)) * np.sin(self._alpha) +             \
                      np.sin(np.arctan(0)) * np.cos(self._alpha)              \
                       * np.cos(self._phisun - 0)
                                                       # flat surface reference

        self._Rsflat = self._Rgl * self._flat
                            # flat surface total incoming shortwave radiation

        self._Rnetflat = (1 - self._A) * (1 - 0.65 * (self._N**2)) * \
                                self._Rsflat
                            # flat surface Net incoming shortwave radiation

        self._sloped = (np.cos(self._slope) *                        \
                    np.sin(self._alpha) + np.sin(self._slope) *    \
                    np.cos(self._alpha) *                          \
                    np.cos(self._phisun - self._aspect))

        self._radf = self._sloped/self._flat

        self._radf[self._radf<=0.] = 0.
        self._radf[self._radf>6.] = 6.

        self._Rs = self._Rsflat * self._radf
                    # Sloped surface Toatl Incoming Shortwave Radn
        self._Rnet = self._Rnetflat * self._radf

        self._cell_values['topographic__radiation_factor'] = self._radf
        self._cell_values['topographic__total_shortwave_radiation'] = self._Rs
        self._cell_values['topographic__net_shortwave_radiation'] = self._Rnet
