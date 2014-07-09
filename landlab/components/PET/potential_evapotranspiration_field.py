#################################################################
##
##  'Field' concept is implemented for Potential Evapotranspiration component.
##
##  Sai Nudurupati and Erkan Istanbulluoglu - 16May2014
## Test - Branching! 09Jul14
#################################################################

from landlab import Component

import numpy as np

_VALID_METHODS = set(['Constant', 'PriestlyTaylor'])

def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)
        

class PotentialEvapotranspiration( Component ):
    """
    Landlab component that calculates Potential Evapotranspiration.
    
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid(5, 4, 1.e4)
    >>> PET = PotentialEvapotranspiration(grid)
    >>> PET.name
    'Potential Evapotranspiration'
    """
    _name = 'Potential Evapotranspiration'
    
    _input_var_names = set([
        'RadiationFactor',        
    ])
    
    _output_var_names = set([
        'PotentialEvapotranspiration',
    ])
    
    _var_units = {
        'PotentialEvapotranspiration' : 'mm',
        'RadiationFactor' : 'None',        
    }
    
    def __init__(self, grid, **kwds):
        self._method = kwds.pop('method', 'Constant')
        # For Priestly Taylor
        self._alpha = kwds.pop('PriestlyTaylorConstant', 1.26)
        self._a = kwds.pop('Albedo', 0.2)
        self._pwhv = kwds.pop('LatentHeatofVaporization', 28.34)
        self._y = kwds.pop('PsychometricConstant', 0.066)
        self._sigma = kwds.pop('StefanBoltzmannConstant', 0.0000000567)
        self._Gsc = kwds.pop('SolarConstant', 1366.67)
        self._lat = kwds.pop('Latitude', 34.0)
        self._z = kwds.pop('ElevationofMeasurement', 300)
        self._Krs = kwds.pop('AdjustmentCoefficient', 0.18)
                
        assert_method_is_valid(self._method)
        
        super(PotentialEvapotranspiration, self).__init__(grid, **kwds)
        
        for name in self._input_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])
                
        for name in self._output_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])
                
        self._cell_values = self.grid['cell']
        
    def update(self, current_time, **kwds):        
                
        if self._method == 'Constant':
            PET_value = kwds.pop('ConstantPotentialEvapotranspiration', 6.)
        elif self._method == 'PriestlyTaylor':
            PET_value = self.PriestlyTaylor( current_time )
        self._PET = PET_value * self._cell_values['RadiationFactor']
        self._cell_values['PotentialEvapotranspiration'] = self._PET
        
    def PriestlyTaylor(self, current_time):
        
        """
            Julian Day - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (25)
        """        
        self._J = np.floor( (current_time - np.floor( current_time)) * 365 ) 

        """
            Solar Declination Angle - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (51)
        """
        self._sdecl = 0.409*np.sin((((2.0*3.14)/365.0)*self._J)-1.39)

        """
            Inverse Relative Distance Factor - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (50)
        """
        self._dr = 1 + (0.033*np.cos((2.0*3.14/365.0)*self._J))

        """
            To calculate ws - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (61)
        """
        self._x = 1.0-((((np.tan(self._phi))**2.0))*((np.tan(self._sdecl)**2.0)))  
        if self._x <= 0:
            self._x = 0.00001;
            
        """
            Sunset Hour Angle - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (60)
        """
        self._ws = (3.14/2.0)-np.atan((-1*np.tan(self._phi)*np.tan(self._sdecl))/(self._x**2.0))

        """
            Extraterrestrial radmodel.docx - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (48)
        """
        self._Ra = 11.57*(24.0/3.14)*4.92*self._dr*((self._ws*np.sin(self._phi)*np.sin(self._sdecl))+ \
                        (np.cos(self._phi)*np.cos(self._sdecl)*(np.sin(self._ws))))

        """
            Clear-sky Solar Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (47)
        """
        self._Rs = (0.75+((2.0*(10**-5.0))*self._z))*self._Ra

        """
            Net Short Wave Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (43)
        """
        self._Rns = self._Rs*(1-self._a)

        """
            Relative Cloudiness - ASCE-EWRI Task Committee Report, Jan-2005 - Page 35
        """
        self._u = 0.3
        
        """
            Cloudiness Function - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (45)
        """
        self._fcd = (1.35*self._u)-0.35

        """
            Net Radiation - ASCE-EWRI Task Committee Report, Jan-2005 - Eqn (42)
        """                                         
        self._Rn = self._Rns 

        """
            Priestly Taylor Evapotranspiration in mm/day - Priestly et. al. 1972
        """
        self._ETp = max(self._alpha*(self._delta/(self._delta+self._y) \
                        )*(self._Rn/self._pwhv), 0)
        
        return( self._ETp )