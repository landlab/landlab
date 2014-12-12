
#################################################################
##
##  Cellular Automaton component design for Vegetation establishment/mortality
##
##  Sai Nudurupati and Erkan Istanbulluoglu - 26 Nov 2014
#################################################################

from landlab import Component
import numpy as np

_VALID_METHODS = set(['Grid'])

def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)


class VegCA( Component ):
    """
    Landlab component that implements 1D and 2D vegetation dynamics
    model.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid(5, 4, 1.e4)
    >>> veg_ca = VegCA(grid)
    >>> veg_ca.name
    'VegCA'
    """
    _name = 'VegCA'

    _input_var_names = set([
        'CumulativeWaterStress',
        'VegetationType',
    ])

    _output_var_names = set([
        
    ])

    _var_units = {
        'CumulativeWaterStress' : 'Pa',
        'VegetationType'  : 'None',
    }

    def __init__(self, grid, **kwds):
        self._method = kwds.pop('method', 'Grid') 
        self._Pemaxg = kwds.pop('Pemaxg', 0.35)   # Pe-max-grass - max probability 
        self._Pemaxsh = kwds.pop('Pemaxsh', 0.2)    # Pe-max-shrub
        self._Pemaxtr = kwds.pop('Pemaxtr', 0.25)   # Pe-max-tree
        self._INg = kwds.pop('ING', 2)  # Allelopathic effect on grass from creosotebush
        
        assert_method_is_valid(self._method)

        super(VegCA, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        for name in self._output_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        self._cell_values = self.grid['cell']


    def update(self, **kwds): 
        # 0 - Bare Soil ; 1 - Shrub Seedling ; 2 - Tree Seedling ;
        # 3 - Grass ; 4 - Shrub ; 5 - Tree
        
        self._VegType = self._cell_values['VegetationType']
        self._CumWS   = self._cell_values['CumulativeWaterStress']
        
        bare_cells = np.where(self._VegType == 0)[0]
        n_bare = len(bare_cells)
        first_ring = self.grid.get_looped_cell_neighbor_list(bare_cells)
        second_ring =                                                   \
            self.grid.get_second_ring_looped_cell_neighbor_list(bare_cells)
        veg_type_fr = self._VegType[first_ring]
        veg_type_sr = self._VegType[second_ring]
        Sh_WS_fr = WS_PFT( veg_type_fr, 4, self._CumWS[first_ring] )
        Tr_WS_fr = WS_PFT( veg_type_fr, 5, self._CumWS[first_ring] )
        Tr_WS_sr = WS_PFT( veg_type_sr, 5, self._CumWS[second_ring] )
                
        n = count(veg_type_fr, 4)
        Phi_sh = Sh_WS_fr/8.
        Phi_tr = (Tr_WS_fr + Tr_WS_sr/2.)/8.
        Phi_g = np.mean(self._CumWS[np.where(self._VegType == 3)])
        Pemaxg = self._Pemaxg * np.ones(n_bare)
        Pemaxsh = self._Pemaxsh * np.ones(n_bare)
        Pemaxtr = self._Pemaxtr * np.ones(n_bare)
        Peg = np.amin(np.vstack((Phi_g/(n*self._INg),Pemaxg)),axis = 0)
        Pesh = np.amin(np.vstack((Phi_sh, Pemaxsh)), axis = 0)
        Petr = np.amin(np.vstack((Phi_tr, Pemaxtr)), axis = 0)
        Select_PFT_E = np.random.randint(0,3,n_bare) # Grass - 0; Shrub - 1; Tree - 2
        Pest = np.choose(Select_PFT_E, [Peg, Pesh, Petr])   # Probability of establishment
        R_Est = np.random.rand(n_bare)  # Random number for comparison to establish
        Establish = np.int32(np.where(np.greater_equal(Pest, R_Est)==True)[0])
        self._VegType[bare_cells[Establish]] = Select_PFT_E[Establish] + 3
        
        
        
        
def count( Arr, value ):
    Res = np.zeros(Arr.shape[0],dtype = int)
    x,y = Arr.shape
    for i in range(0,x):
        for j in range(0,y):
            if Arr[i][j] == value:
                Res[i] += 1
    return Res

def WS_PFT( VegType, PlantType, WS ):
    Phi = np.zeros(WS.shape[0])
    x,y = WS.shape
    for i in range(0,x):
        for j in range(0,y):
            if VegType[i][j] == PlantType:
                Phi[i] += WS[i][j]
    return Phi

                
