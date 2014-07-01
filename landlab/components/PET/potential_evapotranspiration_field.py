#################################################################
##
##  'Field' concept is implemented for Potential Evapotranspiration component.
##
##  Sai Nudurupati and Erkan Istanbulluoglu - 16May2014
#################################################################

from landlab import Component

#import numpy as np

_VALID_METHODS = set(['Grid'])

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
        self._method = kwds.pop('method', 'Grid')
                
        assert_method_is_valid(self._method)
        
        super(PotentialEvapotranspiration, self).__init__(grid, **kwds)
        
        for name in self._input_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])
                
        for name in self._output_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])
                
        self._cell_values = self.grid['cell']
        
    def update(self, **kwds):
        
        PET_constant = kwds.pop('ConstantPotentialEvapotranspiration', 6.)        
        self._PET = PET_constant * self._cell_values['RadiationFactor']
        self._cell_values['PotentialEvapotranspiration'] = self._PET
