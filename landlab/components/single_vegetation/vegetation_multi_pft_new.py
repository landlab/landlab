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


class Vegetation(Component):
    """1D and 2D vegetation dynamics.

    Landlab component that implements 1D and 2D vegetation dynamics model.
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

    def __init__(self, grid, data, **kwds):
        self._method = kwds.pop('method', 'Grid')

        assert_method_is_valid(self._method)

        super(Vegetation, self).__init__(grid)

        self.initialize( data, VEGTYPE = grid['cell']['VegetationType'], \
                            **kwds )

        for name in self._input_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        self._cell_values = self.grid['cell']

        self._Blive_ini = self._Blive_init * np.ones(self.grid.number_of_cells)
        self._Bdead_ini = self._Bdead_init * np.ones(self.grid.number_of_cells)

    def initialize( self, data, **kwds ):
        # GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
        # SHRUBSEEDLING = 4; TREESEEDLING = 5
        self._vegtype = \
          kwds.pop('VEGTYPE', np.zeros(self.grid.number_of_cells,dtype = int))
        self._WUE = np.choose(self._vegtype, kwds.pop('WUE',
                [ data['WUE_grass'], data['WUE_shrub'], data['WUE_tree'],
                  data['WUE_bare'], data['WUE_shrub'], data['WUE_tree'] ]))      # Water Use Efficiency  KgCO2kg-1H2O
        self._LAI_max = np.choose( self._vegtype, kwds.pop('LAI_MAX',
                [ data['LAI_MAX_grass'], data['LAI_MAX_shrub'],
                   data['LAI_MAX_tree'], data['LAI_MAX_bare'],
                   data['LAI_MAX_shrub'], data['LAI_MAX_tree'] ]))               # Maximum leaf area index (m2/m2)
        self._cb = np.choose( self._vegtype, kwds.pop('CB',
                [ data['CB_grass'], data['CB_shrub'], data['CB_tree'],
                  data['CB_bare'], data['CB_shrub'], data['CB_tree'] ]))         # Specific leaf area for green/live biomass (m2 leaf g-1 DM)
        self._cd = np.choose( self._vegtype, kwds.pop('CD',
                [ data['CD_grass'], data['CD_shrub'], data['CD_tree'],
                  data['CD_bare'], data['CD_shrub'], data['CD_tree'] ]))         # Specific leaf area for dead biomass (m2 leaf g-1 DM)
        self._ksg = np.choose( self._vegtype, kwds.pop('KSG',
                [ data['KSG_grass'], data['KSG_shrub'], data['KSG_tree'],
                  data['KSG_bare'], data['KSG_shrub'], data['KSG_tree'] ]))      # Senescence coefficient of green/live biomass (d-1)
        self._kdd = np.choose( self._vegtype, kwds.pop('KDD',
                [ data['KDD_grass'], data['KDD_shrub'], data['KDD_tree'],
                  data['KDD_bare'], data['KDD_shrub'], data['KDD_tree'] ]))      # Decay coefficient of aboveground dead biomass (d-1)
        self._kws = np.choose( self._vegtype, kwds.pop('KWS',
                [ data['KWS_grass'], data['KWS_shrub'], data['KWS_tree'],
                  data['KWS_bare'], data['KWS_shrub'], data['KWS_tree'] ]))      # Maximum drought induced foliage loss rates (d-1)
        self._Blive_init = kwds.pop('BLIVE_INI', data['BLIVE_INI'])
        self._Bdead_init = kwds.pop('BDEAD_INI', data['BDEAD_INI'])
        self._ETthresholdup = kwds.pop('ETTup', data['ETTup'])                   # Growth threshold (mm/d)
        self._ETthresholddown = kwds.pop('ETTdwn', data['ETTdwn'])               # Dormancy threshold (mm/d)
        self._Tdmax = kwds.pop('Tdmax', data['Tdmax'])                           # Constant for dead biomass loss adjustment
        self._w = kwds.pop('w', data['w'])                                       # Conversion factor of CO2 to dry biomass
        #self._ETdmax = np.choose( self._vegtype, kwds.pop('ETdmax',
        #                        [ 10, 10, 10, 10, 10, 10 ]))                    # Constant for dead biomass loss adjustment (mm/d)

        self._cell_values = self.grid['cell']

    def update(self, **kwds):

        PETthreshold_ = kwds.pop('PotentialEvapotranspirationThreshold', 0)
        Tb = kwds.pop('Tb', 24.)
        Tr = kwds.pop('Tr', 0.01)
        PET = self._cell_values['PotentialEvapotranspiration']
        PET30_ = self._cell_values['PotentialEvapotranspiration30']
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

            WUE = self._WUE[cell]
            LAImax = self._LAI_max[cell]
            cb = self._cb[cell]
            cd = self._cd[cell]
            ksg = self._ksg[cell]
            kdd = self._kdd[cell]
            kws = self._kws[cell]
            #ETdmax = self._ETdmax[cell]

            LAIlive =  min (cb * self._Blive_ini[cell], LAImax)
            LAIdead =  min (cd * self._Bdead_ini[cell], (LAImax                \
                            - LAIlive))
            NPP = max((ActualET[cell]/(Tb+Tr))*                                \
                                WUE*24.*self._w*1000, 0.001)

            if self._vegtype[cell] == 0:

                if PET30_[cell] > PETthreshold:          # Growing Season

                    Bmax = (LAImax - LAIdead)/cb
                    Yconst = (1/((1/Bmax)+(((kws*Water_stress[cell])+          \
                                ksg)/NPP)))
                    Blive = (self._Blive_ini[cell] - Yconst) *                 \
                                np.exp(-(NPP/Yconst) * ((Tb+Tr)/24.)) + Yconst
                    Bdead = (self._Bdead_ini[cell] + (Blive - max(Blive *      \
                                np.exp(-ksg * Tb/24.),0.00001)))*              \
                                np.exp(-kdd * min( PET[cell]/self._Tdmax, 1. )*\
                                Tb/24.)

                else:                                 # Senescense

                    Blive = max(self._Blive_ini[cell] * np.exp((-2) * ksg *    \
                                Tb/24.), 1 )
                    Bdead = max( (self._Bdead_ini[cell] +                      \
                                ( self._Blive_ini[cell] -                      \
                                max( self._Blive_ini[cell]*np.exp( (-2) *      \
                                ksg * Tb/24.), 0.000001)))*                    \
                                np.exp( (-1) * kdd *                           \
                                min( PET[cell]/self._Tdmax, 1. ) * Tb/24.), 0. )

            elif self._vegtype[cell] == 3:

                Blive = 0.
                Bdead = 0.

            else:

                Bmax = LAImax/cb
                Yconst = (1/((1/Bmax)+(((kws*Water_stress[cell]) +             \
                            ksg)/NPP)))
                Blive = (self._Blive_ini[cell] - Yconst) *                     \
                            np.exp(-(NPP/Yconst) * ((Tb+Tr)/24.)) + Yconst
                Bdead = (self._Bdead_ini[cell] + (Blive - max(Blive *          \
                            np.exp(-ksg * Tb/24.),0.00001))) *                 \
                            np.exp(-kdd * min( PET[cell]/self._Tdmax, 1. ) *   \
                            Tb/24.)

            LAIlive =  min( cb * (Blive + self._Blive_ini[cell])/2., LAImax )
            LAIdead =  min( cd * (Bdead + self._Bdead_ini[cell])/2.,
                                                       ( LAImax - LAIlive ) )
            if self._vegtype[cell] == 0:
                Vt = 1. - np.exp(-0.75 * (LAIlive + LAIdead))
            else:
                #Vt = 1 - np.exp(-0.75 * LAIlive)
                Vt = 1.

            self._LAIlive[cell] = LAIlive
            self._LAIdead[cell] = LAIdead
            self._VegCov[cell] = Vt
            self._Blive[cell] = Blive
            self._Bdead[cell] = Bdead

        self._Blive_ini = self._Blive
        self._Bdead_ini = self._Bdead
