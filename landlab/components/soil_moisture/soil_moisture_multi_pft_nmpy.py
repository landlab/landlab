#################################################################
##
##  Modification of soil_moisture_field.py to accomodate
##  multiple Plant Functional Types (PFTs)
##
##  Sai Nudurupati and Erkan Istanbulluoglu - 31Oct2014
#################################################################
GRASS = 0
SHRUB = 1
TREE = 2
BARE = 3
SHRUBSEEDLING = 4
TREESEEDLING = 5


from landlab import Component

import numpy as np

_VALID_METHODS = set(['Grid', 'Multi'])

def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)

class SoilMoistureNMPY( Component ):

    _name = 'Soil Moisture'

    _input_var_names = set([
        'VegetationCover',
        'LiveLeafAreaIndex',
        'PotentialEvapotranspiraton',
    ])

    _output_var_names = set([
        'WaterStress',
        'SaturationFraction',
        'Drainage',
        'Runoff',
        'ActualEvapotranspiration',
    ])

    _var_units = {
        'VegetationCover' : 'None',
        'LiveLeafAreaIndex': 'None',
        'PotentialEvapotranspiraton' : 'mm',
        'WaterStress' : 'Pa',
        'SaturationFraction' : 'None',
        'Drainage' : 'mm',
        'Runoff' : 'mm',
        'ActualEvapotranspiration' : 'mm',
    }


    def __init__( self, grid, data, **kwds ):

        self._method = kwds.pop('method', 'Grid')

        assert_method_is_valid(self._method)

        super(SoilMoistureNMPY, self).__init__(grid)

        self.initialize( data, VEGTYPE = grid['cell']['VegetationType'], \
                            **kwds )

        for name in self._input_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        for name in self._output_var_names:
            if not name in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        self._nodal_values = self.grid['node']

        if not 'InitialSaturationFraction' in self.grid.at_cell:
            self.grid.add_zeros('cell', 'InitialSaturationFraction',
                                    units='None' )

        self._cell_values = self.grid['cell']


    def initialize( self, data, **kwds ):
        # GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
        # SHRUBSEEDLING = 4; TREESEEDLING = 5
        self._vegtype = \
          kwds.pop('VEGTYPE', self.grid['cell']['VegetationType'])
        self._runon = kwds.pop('RUNON', 0.)
        self._fbare = kwds.pop('F_BARE', data['F_BARE'])

        self._interception_cap = \
                np.choose(self._vegtype, kwds.pop('INTERCEPT_CAP',
                [ data['INTERCEPT_CAP_grass'], data['INTERCEPT_CAP_shrub'],
                  data['INTERCEPT_CAP_tree'], data['INTERCEPT_CAP_bare'],
                  data['INTERCEPT_CAP_shrub'], data['INTERCEPT_CAP_tree'] ]))   # Full canopy interception (mm)
        self._zr = np.choose(self._vegtype, kwds.pop('ZR',
                [ data['ZR_grass'], data['ZR_shrub'], data['ZR_tree'],
                  data['ZR_bare'], data['ZR_shrub'], data['ZR_tree'] ]))        # Root depth (m)
        self._soil_Ib = np.choose(self._vegtype, kwds.pop('I_B',
                [ data['I_B_grass'], data['I_B_shrub'], data['I_B_tree'],
                  data['I_B_bare'], data['I_B_shrub'], data['I_B_tree'] ]))     # Infiltration capacity of bare soil (mm/h)
        self._soil_Iv = np.choose(self._vegtype, kwds.pop('I_V',
                [ data['I_V_grass'], data['I_V_shrub'], data['I_V_tree'],
                  data['I_V_bare'], data['I_V_shrub'], data['I_V_tree'] ]))     # Infiltration capacity of vegetated soil (mm/h)
        self._soil_Ew = kwds.pop('EW', 0.1)
        self._soil_pc = np.choose(self._vegtype, kwds.pop('PC',
                [ data['PC_grass'], data['PC_shrub'], data['PC_tree'],
                  data['PC_bare'], data['PC_shrub'], data['PC_tree'] ]))        # Soil porosity
        self._soil_fc = np.choose(self._vegtype, kwds.pop('FC',
                [ data['FC_grass'], data['FC_shrub'], data['FC_tree'],
                  data['FC_bare'], data['FC_shrub'], data['FC_tree'] ]))        # Saturation degree at soil field capacity
        self._soil_sc = np.choose(self._vegtype, kwds.pop('SC',
                [ data['SC_grass'], data['SC_shrub'], data['SC_tree'],
                  data['SC_bare'], data['SC_shrub'], data['SC_tree'] ]))        # Saturation degree at soil stomatal closure
        self._soil_wp = np.choose(self._vegtype, kwds.pop('WP',
                [ data['WP_grass'], data['WP_shrub'], data['WP_tree'],
                  data['WP_bare'], data['WP_shrub'], data['WP_tree'] ]))        # Saturation degree at soil wilting point
        self._soil_hgw = np.choose(self._vegtype, kwds.pop('HGW',
                [ data['HGW_grass'], data['HGW_shrub'], data['HGW_tree'],
                  data['HGW_bare'], data['HGW_shrub'], data['HGW_tree'] ]))     # Saturation degree at soil hygroscopic point
        self._soil_beta = np.choose(self._vegtype, kwds.pop('BETA',
                [ data['BETA_grass'], data['BETA_shrub'], data['BETA_tree'],
                  data['BETA_bare'], data['BETA_shrub'], data['BETA_tree'] ]))  # Deep percolation constant
        self._LAI_max = np.choose( self._vegtype, kwds.pop('LAI_MAX',
                [ data['LAI_MAX_grass'], data['LAI_MAX_shrub'],
                   data['LAI_MAX_tree'], data['LAI_MAX_bare'],
                   data['LAI_MAX_shrub'], data['LAI_MAX_tree'] ]))              # Maximum leaf area index (m2/m2)
        self._LAIR_max = np.choose( self._vegtype, kwds.pop('LAIR_MAX',
                [ data['LAIR_MAX_grass'], data['LAIR_MAX_shrub'],
                   data['LAIR_MAX_tree'], data['LAIR_MAX_bare'],
                   data['LAIR_MAX_shrub'], data['LAIR_MAX_tree'] ]))            # Reference leaf area index (m2/m2)


    def update( self, current_time, **kwds ):
        #DEBUGG = 0

        P = kwds.pop('P', 5.)*np.ones(self.grid.number_of_cells)
        Tb = kwds.pop('Tb', 24.)*np.ones(self.grid.number_of_cells)
        Tr = kwds.pop('Tr', 0.0)*np.ones(self.grid.number_of_cells)
        self._PET = self._cell_values['PotentialEvapotranspiration']
        self._SO = self._cell_values['InitialSaturationFraction']
        self._vegcover = self._cell_values['VegetationCover']
        self._water_stress = self._cell_values['WaterStress']
        self._S = self._cell_values['SaturationFraction']
        self._D = self._cell_values['Drainage']
        self._ETA = self._cell_values['ActualEvapotranspiration']
        self._fr = self._cell_values['LiveLeafAreaIndex']/self._LAIR_max
        self._fr[self._fr > 1.] = 1.
        self._Sini = np.zeros(self._SO.shape)
        self._ETmax = np.zeros(self._SO.shape) # record ETmax - Eq 5 - Zhou et al.

#
#        for cell in range(0,self.grid.number_of_cells):
#
#            #print cell
        s = self._SO

        fbare = self._fbare
        ZR = self._zr
        pc = self._soil_pc
        fc = self._soil_fc
        sc = self._soil_sc
        wp = self._soil_wp
        hgw = self._soil_hgw
        beta = self._soil_beta
        gr = np.where(self._vegtype == GRASS)
        sc[gr] = sc[gr]*self._fr[gr]+(1-self._fr[gr])*fc[gr]


        Inf_cap = self._soil_Ib *(1-self._vegcover) +         \
                                self._soil_Iv*self._vegcover
                                                    # Infiltration capacity
        Int_cap = np.minimum(self._vegcover*self._interception_cap,
                        P)                # Interception capacity
        Peff = np.maximum(P-Int_cap, np.zeros(self.grid.number_of_cells))
                                          # Effective precipitation depth
        mu = (Inf_cap/1000.0)/(pc*ZR*(np.exp(beta*(1.-fc))-1.))
        Ep = np.maximum((self._PET*self._fr
                +fbare*self._PET*(1.-self._fr)) - Int_cap, 0.0001*
                 np.ones(self.grid.number_of_cells))  # mm/d
        self._ETmax = Ep
        nu = ((Ep/24.)/1000.)/(pc*ZR) # Loss function parameter
        nuw = ((self._soil_Ew/24.)/1000.)/(pc*ZR) # Loss function parameter
        sini = self._SO + ((Peff+self._runon)/(pc*ZR*1000.))

        self._runoff = np.zeros(self.grid.number_of_cells)
        a1 = np.where(sini>1.)
        self._runoff[a1] = (sini[a1]-1.)*pc[a1]*ZR[a1]*1000.
        sini[a1] = np.ones(a1[0].shape,dtype = float)


        tfc = np.zeros(self.grid.number_of_cells,dtype = float)
        tsc = np.zeros(self.grid.number_of_cells,dtype = float)
        twp = np.zeros(self.grid.number_of_cells,dtype = float)

        c1 = np.where(sini>=fc)
        #if sini>=fc:
        tfc[c1] = (1./(beta[c1]*(mu[c1]-nu[c1])))*(beta[c1]*(fc[c1]-sini[c1])+ \
                    np.log((nu[c1]-mu[c1]+mu[c1]*np.exp(beta[c1]*
                    (sini[c1]-fc[c1])))/nu[c1]))
        tsc[c1] = ((fc[c1]-sc[c1])/nu[c1])+tfc[c1]
        twp[c1] = ((sc[c1]-wp[c1])/(nu[c1]-nuw[c1]))*       \
                    np.log(nu[c1]/nuw[c1])+tsc[c1]

        c2 = np.where(Tb[c1]<tfc[c1])
        #if Tb<tfc:
        s[c1][c2] = abs(sini[c1][c2]-(1./beta[c1][c2])*np.log(((nu[c1][c2]-    \
                        mu[c1][c2]+mu[c1][c2]* np.exp(beta[c1][c2]*
                        (sini[c1][c2]-fc[c1][c2])))*np.exp(beta[c1][c2]*
                        (nu[c1][c2]-mu[c1][c2])*Tb[c1][c2])-mu[c1][c2]*
                        np.exp(beta[c1][c2]*(sini[c1][c2]-fc[c1][c2])))/
                        (nu[c1][c2]-mu[c1][c2])))

        self._D[c1][c2] = ((pc[c1][c2]*ZR[c1][c2]*1000.)*
                            (sini[c1][c2]-s[c1][c2]))-(Tb[c1][c2]*
                            (Ep[c1][c2]/24.))
        self._ETA[c1][c2] = (Tb[c1][c2]*(Ep[c1][c2]/24.))

        c3 = np.where(np.logical_and(Tb[c1]>=tfc[c1],Tb[c1]<tsc[c1]))
        #elif Tb>=tfc and Tb<tsc:
        s[c1][c3] = fc[c1][c3]-(nu[c1][c3]*(Tb[c1][c3]-tfc[c1][c3]))
        self._D[c1][c3] = ((pc[c1][c3]*ZR[c1][c3]*1000.)*(sini[c1][c3]-        \
                            fc[c1][c3]))-((tfc[c1][c3])*(Ep[c1][c3]/24.))
        self._ETA[c1][c3] = (Tb[c1][c3]*(Ep[c1][c3]/24.))

        c4 = np.where(np.logical_and(Tb[c1]>=tsc[c1],Tb[c1]<twp[c1]))
        #elif Tb>=tsc and Tb<twp:
        s[c1][c4] = wp[c1][c4]+(sc[c1][c4]-wp[c1][c4])*((nu[c1][c4]/
                     (nu[c1][c4]-nuw[c1][c4]))*np.exp((-1)*((nu[c1][c4]-       \
                     nuw[c1][c4])/(sc[c1][c4]-wp[c1][c4]))*(Tb[c1][c4]-        \
                     tsc[c1][c4]))-(nuw[c1][c4]/(nu[c1][c4]-nuw[c1][c4])))
        self._D[c1][c4] = ((pc[c1][c4]*ZR[c1][c4]*1000.)*(sini[c1][c4]-        \
                            fc[c1][c4]))-(tfc[c1][c4]*Ep[c1][c4]/24.)
        self._ETA[c1][c4] = (1000.*ZR[c1][c4]*pc[c1][c4]*(sini[c1][c4]-        \
                             s[c1][c4]))-self._D[c1][c4]

        c5 = np.where(Tb[c1]>=twp[c1])

        #else:
        s[c1][c5] = hgw[c1][c5]+(wp[c1][c5]-hgw[c1][c5])*np.exp((-1)*
                     (nuw[c1][c5]/(wp[c1][c5]-hgw[c1][c5]))*
                      np.maximum(Tb[c1][c5]-twp[c1][c5],np.zeros(c5[0].shape)))
        self._D[c1][c5] = ((pc[c1][c5]*ZR[c1][c5]*1000.)*(sini[c1][c5]-        \
                            fc[c1][c5]))-(tfc[c1][c5]*Ep[c1][c5]/24.)
        self._ETA[c1][c5] = (1000.*ZR[c1][c5]*pc[c1][c5]*(sini[c1][c5]-        \
                             s[c1][c5]))-self._D[c1][c5]

        d1 = np.where(np.logical_and(sini<fc,sini>=sc))
        #elif sini<fc and sini>=sc:
        tfc[d1] = np.zeros(d1[0].shape)
        tsc[d1] = (sini[d1]-sc[d1])/nu[d1]
        twp[d1] = ((sc[d1]-wp[d1])/(nu[d1]-nuw[d1]))*np.log(nu[d1]/nuw[d1])    \
                    +tsc[d1]

        d2 = np.where(Tb[d1]<tsc[d1])
        #if Tb<tsc:
        s[d1][d2] = sini[d1][d2] - nu[d1][d2]*Tb[d1][d2]
        self._D[d1][d2] = np.zeros(d2[0].shape)
        self._ETA[d1][d2] = 1000.*ZR[d1][d2]*pc[d1][d2]*(sini[d1][d2]-s[d1][d2])

        d3 = np.where(np.logical_and(Tb[d1]>=tsc[d1],Tb[d1]<twp[d1]))
        #elif Tb>=tsc and Tb<twp:
        s[d1][d2] = wp[d1][d2]+(sc[d1][d2]-wp[d1][d2])*((nu[d1][d2]/
                     (nu[d1][d2]-nuw[d1][d2]))*np.exp((-1*np.ones(d3[0].shape))
                      *((nu[d1][d2]-nuw[d1][d2])/(sc[d1][d2]-wp[d1][d2]))*
                      (Tb[d1][d2]-tsc[d1][d2]))-(nuw[d1][d2]/(nu[d1][d2]-      \
                       nuw[d1][d2])))
        self._D[d1][d2] = np.zeros(d3[0].shape)
        self._ETA[d1][d2] = (1000.*ZR[d1][d2]*pc[d1][d2]*(sini[d1][d2]-        \
                             s[d1][d2]))

        d4 = np.where(Tb[d1]>=twp[d1])
        #else:
        s[d1][d4] = hgw[d1][d4]+(wp[d1][d4]-hgw[d1][d4])*np.exp((-1)*
                     (nuw[d1][d4]/(wp[d1][d4]-hgw[d1][d4]))*
                      (Tb[d1][d4]-twp[d1][d4]))
        self._D[d1][d4] = np.zeros(d4[0].shape)
        self._ETA[d1][d4] = (1000.*ZR[d1][d4]*pc[d1][d4]*(sini[d1][d4]-        \
                             s[d1][d4]))

        e1 = np.where(np.logical_and(sini<sc, sini>=wp))
        #elif sini<sc and sini>=wp:
        tfc[e1] = np.zeros(e1[0].shape)
        tsc[e1] = np.zeros(e1[0].shape)
        twp[e1] = ((sc[e1]-wp[e1])/(nu[e1]-nuw[e1]))*np.log(1+(nu[e1]-nuw[e1])*
                    (sini[e1]-wp[e1])/(nuw[e1]*(sc[e1]-wp[e1])))

        e2 = np.where(Tb[e1]<twp[e1])
        #if Tb<twp:
        s[e1][e2] = wp[e1][e2]+((sc[e1][e2]-wp[e1][e2])/(nu[e1][e2]-           \
                     nuw[e1][e2]))*((np.exp((-1)*((nu[e1][e2]-nuw[e1][e2])
                     /(sc[e1][e2]-wp[e1][e2]))*Tb[e1][e2]))*(nuw[e1][e2]+      \
                     ((nu[e1][e2]-nuw[e1][e2])/(sc[e1][e2]-wp[e1][e2]))*
                     (sini[e1][e2]-wp[e1][e2]))-nuw[e1][e2])
        self._D[e1][e2] = np.zeros(e2[0].shape)
        self._ETA[e1][e2] = (1000.*ZR[e1][e2]*pc[e1][e2]*(sini[e1][e2]-        \
                              s[e1][e2]))

        e3 = np.where(Tb[e1]>=twp[e1])
        #else:
        s[e1][e3] = hgw[e1][e3]+(wp[e1][e3]-hgw[e1][e3])*np.exp((-1)*
                     (nuw[e1][e3]/(wp[e1][e3]-hgw[e1][e3]))*
                     (Tb[e1][e3]-twp[e1][e3]))
        self._D[e1][e3] = np.zeros(e3[0].shape)
        self._ETA[e1][e3] = (1000.*ZR[e1][e3]*pc[e1][e3]*
                              (sini[e1][e3]-s[e1][e3]))

        f1 = np.where(sini<wp)
        #else:
        tfc[f1] = np.zeros(f1[0].shape)
        tsc[f1] = np.zeros(f1[0].shape)
        twp[f1] = np.zeros(f1[0].shape)

        s[f1] = hgw[f1]+(sini[f1]-hgw[f1])*np.exp((-1)*(nuw[f1]/
                  (wp[f1]-hgw[f1]))*Tb[f1])
        self._D[f1] = np.zeros(f1[0].shape)
        self._ETA[f1] = (1000.*ZR[f1]*pc[f1]*(sini[f1]-s[f1]))

        self._water_stress = np.minimum(((np.maximum(((sc - (s+sini)/2.)
                              /(sc - wp)),np.zeros(self.grid.number_of_cells))) \
                               **4.),np.ones(self.grid.number_of_cells))
        self._S = s
        self._SO = s
        self._Sini = sini

        current_time += (Tb+Tr)/(24.*365.25)
        return( current_time )
