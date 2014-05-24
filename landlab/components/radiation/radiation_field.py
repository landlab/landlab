#################################################################
##
##  'Field' concept is implemented for Radiation component.
##
##  Sai Nudurupati and Erkan Istanbulluoglu - 14May2014
#################################################################

from landlab import Component

import numpy as np

_VALID_METHODS = set(['Grid'])

def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)

class Radiation( Component ):
    """
    Landlab component that computes 1D and 2D total incident shortwave
    radiation. This code also computes relative incidence shortwave radiation
    compared to a flat surface.
    
    Radiation(grid, **kwds)
        Adds two 'cellular' fields on grid 'TotalShortWaveRadiation' and 
        'RadiationFactor'
        
    Parameters:
        grid : RasterModelGrid (might work on other grids but not tested yet)
        
      Optional (**kwds):
        method : Currently, only default is available
        CLOUDINESS: set cloudiness value. default value is 0.5
        LATITUDE: set Latitude. default value is 34.0
        ALBEDO: set albedo. default value is 0.8
        SOLARCONSTANT: default value is 1353.0
        CLRSKYTURBIDITY: set clear sky turbidity. default value is 2.
        OPTAIRMASS: set optical air mass. default value is 0.0
        
    >>> from landlab import RasterModelGrid        
    >>> from landlab.components.radiation.radiation_field import Radiation        
    >>> import numpy as np        
    >>> grid = RasterModelGrid( 5, 4, 0.2 )
    >>> grid['node']['Elevation'] = np.random.rand( grid.number_of_nodes ) * 1000
    >>> rad = Radiation( grid )
    >>> rad.name
        Radiation
    >>> current_time = 0.5
    >>> rad.update( current_time )
    
    >>> grid['cell']['TotalShortWaveRadiation']
        Out[11]: 
        array([ 1.90762961,  1.90701536,  1.90732026,  1.90617417,  1.90762696,
        1.90718441])
    >>> grid['cell']['RadiationFactor']
        Out[12]: 
        array([ 0.01091503,  0.01091151,  0.01091326,  0.0109067 ,  0.01091501,
        0.01091248])
    """
    
    _name = 'Radiation'
    
    _input_var_names = set([
        'Elevation',
    ])
    
    _output_var_names = set([
        'TotalShortWaveRadiation',
        'RadiationFactor',
    ])
    
    _var_units = {
        'Elevation' : 'm',
        'TotalShortWaveRadiation' : 'W/m^2',
        'RadiationFactor' : 'None',
    }
    
    def __init__( self, grid, **kwds ):
        self._method = kwds.pop('method', 'Grid')
        self._N = kwds.pop('CLOUDINESS', 0.5)
        self._latitude = kwds.pop('LATITUDE', 34.0)
        self._A = kwds.pop('ALBEDO', 0.8)
        self._Io = kwds.pop('SOLARCONSTANT', 1353.0)
        self._n = kwds.pop('CLRSKYTURBIDITY', 2.0)
        self._m = kwds.pop('OPTAIRMASS', 0.0)

        assert_method_is_valid(self._method)
        
        super(Radiation, self).__init__(grid, **kwds)
        
        for name in self._input_var_names:
            if not name in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])
                
        for name in self._output_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])
                
        self._nodal_values = self.grid['node']

        if not 'Slope' in self.grid.at_cell:
            self.grid.add_zeros('cell', 'Slope', units='None' )

        self._cell_values = self.grid['cell']       
             
    def update( self, current_time ):        
        
        self._elev = self._nodal_values['Elevation']
        self._slope = self._cell_values['Slope']
        self._radf = self._cell_values['RadiationFactor']
        self._Rs = self._cell_values['TotalShortWaveRadiation']
      
        self._julian = np.floor( ( current_time - np.floor( current_time ) )  \
                                  * 365 )                          # Julian day

        self._t = 12                                            # Assuming noon
                                                                                                                                
        self._phi = np.pi/180.0 * self._latitude        # Latitude in Radians
        
        self._delta = 23.45 * np.pi/180.0 * np.cos(2*np.pi/365                 \
                        * (172 - self._julian))             # Declination angle

        
        self._tau = (self._t + 12.0) * np.pi/12.0                # Hour angle

        self._alpha = np.arcsin(np.sin(self._delta) * np.sin(self._phi)        \
                        + np.cos(self._delta) * np.cos(self._phi)              \
                        * np.cos(self._tau))                  # Solar Altitude        

        if self._alpha <= 0.25 * np.pi/180.0:         # If altitude is -ve, 
            self._alpha = 0.25 * np.pi/180.0       # sun is beyond the horizon
            
        self._Rgl = (1 - self._A) * (1 - 0.65 * (self._N**2)) *                \
                 (self._Io*np.exp((-1) * self._n * (0.128 - 0.054 *            \
                    np.log10(1/np.sin(self._alpha)))*(1/np.sin(self._alpha))))
                    # Counting for Albedo, Cloudiness and Atmospheric turbidity

        self._phisun = np.arctan(-np.sin(self._tau)/(np.tan(self._delta)       \
                    *np.cos(self._phi) - np.sin(self._phi)                     \
                        *np.cos(self._tau))) 
                                                                # Sun's Azhimuth

        if ( self._phisun >= 0 and -np.sin(self._tau) <= 0 ):
            self._phisun = self._phisun + np.pi
            
        elif ( self._phisun <= 0 and -np.sin(self._tau) >= 0 ):
            self._phisun = self._phisun + np.pi

        self._flat = np.cos(np.arctan(0)) * np.sin(self._alpha) +              \
                      np.sin(np.arctan(0)) * np.cos(self._alpha)               \
                       * np.cos(self._phisun - 0)
                                                       # flat surface reference
        
        
        
        for i in range( 0, self.grid.number_of_cells ):
            
                       
            Slope, Aspect = self.grid.calculate_slope_aspect_at_nodes_bestFitPlane([self.grid.node_index_at_cells[i]],
                                            self._elev)
            self._slope[i] = Slope[0]               
            aspect = Aspect[0] * np.pi/180.
            
            slope_in_radians = np.arctan(self._slope[i])
            
            self._radf[i] = np.cos(slope_in_radians) *             \
                        np.sin(self._alpha) + np.sin(np.arctan(0)) *        \
                        np.cos(self._alpha)                                 \
                       * np.cos(self._phisun - aspect)/self._flat
                          
            if self._radf[i] <= 0:
                self._radf[i] = 0.1
            if self._radf[i] > 6.0:
                self._radf[i] = 6.0       
        
            self._Rs[i] = self._Rgl * self._radf[i]    # Incoming Shortwave Radn        