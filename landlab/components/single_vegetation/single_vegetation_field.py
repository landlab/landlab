#################################################################
##
##  'Field' concept is implemented for Single Species Vegetation.
##
##  Sai Nudurupati and Erkan Istanbulluoglu - 04Mar2014
#################################################################

from landlab import Component

import numpy as np

_VALID_METHODS = set(['Grid'])

def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)


class Vegetation( Component ):
    """
    Landlab component that implements 1D and 2D vegetation dynamics
    model.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.single_vegetation import Vegetation
    >>> grid = RasterModelGrid(5, 4, 1.e4)
    >>> veg = Vegetation(grid)
    >>> veg.name
    'Vegetation'
    """
    _name = 'Vegetation'

    _input_var_names = set([
        'ActualEvapotranspiration',
        'WaterStress',
        'PotentialEvapotranspiration',
    ])

    _output_var_names = set([
        'LiveLeafAreaIndex',
        'DeadLeafAreaIndex',
        'VegetationCover',
        'LiveBiomass',
        'DeadBiomass',
    ])

    _var_units = {
        'LiveLeafAreaIndex' : 'None',
        'DeadLeafAreaIndex' : 'None',
        'VegetationCover'  : 'None',
        'ActualEvapotranspiration': 'mm',
        'PotentialEvapotranspiration': 'mm',
        'WaterStress': 'Pa',
        'LiveBiomass' : 'g DM m^-2 d^-1',
        'DeadBiomass' : 'g DM m^-2 d^-1',
    }

    def __init__(self, grid, **kwds):
        self._method = kwds.pop('method', 'Grid')
        self._WUE = kwds.pop('WUE', 0.01)
        self._LAI_max = kwds.pop('LAI_MAX', 2.)
        self._cb = kwds.pop('CB', 0.0047)
        self._cd = kwds.pop('CD', 0.009)
        self._ksg = kwds.pop('KSG', 0.012)
        self._kdd = kwds.pop('KDD', 0.013)
        self._kws = kwds.pop('KWS', 0.02)
        self._Blive_init = kwds.pop('BLIVE_INI', 102.)
        self._Bdead_init = kwds.pop('BDEAD_INI', 450.)
        self._ETthresholdup = kwds.pop('ETTup', 3.8)
        self._ETthresholddown = kwds.pop('ETTdwn', 6.5)
        self._ETdmax = kwds.pop('ETdmax', 10)

        assert_method_is_valid(self._method)

        super(Vegetation, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        self._cell_values = self.grid['cell']

        self._Blive_ini = self._Blive_init * np.ones(self.grid.number_of_cells)
        self._Bdead_ini = self._Bdead_init * np.ones(self.grid.number_of_cells)

    def update(self, **kwds):

        PETthreshold_ = kwds.pop('PotentialEvapotranspirationThreshold', 0)
        Tb = kwds.pop('Tb', 24.)
        Tr = kwds.pop('Tr', 0.01)
        PET = self._cell_values['PotentialEvapotranspiration']
        ActualET = self._cell_values['ActualEvapotranspiration']
        Water_stress = self._cell_values['WaterStress']

        self._LAIlive = self._cell_values['LiveLeafAreaIndex']
        self._LAIdead = self._cell_values['DeadLeafAreaIndex']
        self._Blive = self._cell_values['LiveBiomass']
        self._Bdead = self._cell_values['DeadBiomass']
        self._VegCov = self._cell_values['VegetationCover']

        if PETthreshold_ == 1:
            PETthreshold = self._ETthresholdup
        else:
            PETthreshold = self._ETthresholddown

        for cell in range(0, self.grid.number_of_cells):

            LAIlive =  min (self._cb * self._Blive_ini[cell], self._LAI_max)
            LAIdead =  min (self._cd * self._Bdead_ini[cell], (self._LAI_max   \
                            - LAIlive))

            if PET[cell] > PETthreshold:  # Growing Season

                NPP = max((ActualET[cell]/(Tb+Tr))*                            \
                            self._WUE*24.*0.55*1000., 0.001)
                Bmax = (self._LAI_max - LAIdead)/self._cb
                Yconst = (1./((1./Bmax)+(((self._kws*Water_stress[cell])+      \
                            self._ksg)/NPP)))
                Blive = (self._Blive_ini[cell] - Yconst) *np.exp(-(NPP/Yconst)*\
                            ((Tb+Tr)/24.)) + Yconst
                Bdead = (self._Bdead_ini[cell] + (Blive - max(Blive *          \
                            np.exp(-self._ksg * Tb/24),0.00001)))*             \
                            np.exp(-self._kdd * min( PET[cell]/10., 1. ) *     \
                            Tb/24.)

            else:                                 # Senescense

                Blive = max(self._Blive_ini[cell] * np.exp((-2) * self._ksg *  \
                            Tb/24.), 1. )
                Bdead = max( (self._Bdead_ini[cell] + ( self._Blive_ini[cell]  \
                            -max( self._Blive_ini[cell]*np.exp( (-2) *         \
                            self._ksg * Tb/24.), 0.000001)))*                  \
                            np.exp( (-1)*self._kdd *                           \
                            min( PET[cell]/10., 1. ) * Tb/24.), 0. )

            LAIlive =  min( self._cb * Blive, self._LAI_max )
            LAIdead =  min( self._cd * Bdead, ( self._LAI_max - LAIlive ) )
            Vt = 1. - np.exp(-0.75 * (LAIlive + LAIdead))

            self._LAIlive[cell] = LAIlive
            self._LAIdead[cell] = LAIdead
            self._VegCov[cell] = Vt
            self._Blive[cell] = Blive
            self._Bdead[cell] = Bdead

        self._Blive_ini = self._Blive
        self._Bdead_ini = self._Bdead
