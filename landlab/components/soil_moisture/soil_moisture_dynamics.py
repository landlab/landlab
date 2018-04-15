from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import numpy as np

_VALID_METHODS = set(['Grid', 'Multi'])


def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)


class SoilMoisture(Component):
    """
    Landlab component that simulates root-zone average soil moisture at each
    cell using inputs of potential evapotranspiration, live leaf area index,
    and vegetation cover.

    This component uses a single soil moisture layer and models soil moisture
    loss through transpiration by plants, evaporation by bare soil, and
    leakage. The solution of water balance is based on Laio et. al 2001. The
    component requires fields of initial soil moisture, rainfall input (if
    any), time to the next storm and potential transpiration.

    Ref: Laio, F., Porporato, A., Ridolfi, L., & Rodriguez-Iturbe, I. (2001).
    Plants in water-controlled ecosystems: active role in hydrologic processes
    and response to water stress: II. Probabilistic soil moisture dynamics.
    Advances in Water Resources, 24(7), 707-723.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.soil_moisture import SoilMoisture
    >>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    >>> SoilMoisture.name
    'Soil Moisture'
    >>> sorted(SoilMoisture.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['soil_moisture__root_zone_leakage',
     'soil_moisture__saturation_fraction',
     'surface__evapotranspiration',
     'surface__runoff',
     'vegetation__water_stress']
    >>> sorted(SoilMoisture.units) # doctest: +NORMALIZE_WHITESPACE
    [('rainfall__daily_depth', 'mm'),
     ('soil_moisture__initial_saturation_fraction', 'None'),
     ('soil_moisture__root_zone_leakage', 'mm'),
     ('soil_moisture__saturation_fraction', 'None'),
     ('surface__evapotranspiration', 'mm'),
     ('surface__potential_evapotranspiration_rate', 'mm'),
     ('surface__runoff', 'mm'),
     ('vegetation__cover_fraction', 'None'),
     ('vegetation__live_leaf_area_index', 'None'),
     ('vegetation__plant_functional_type', 'None'),
     ('vegetation__water_stress', 'None')]
    >>> grid['cell']['vegetation__plant_functional_type']= (
    ...            np.zeros(grid.number_of_cells, dtype=int))
    >>> SM = SoilMoisture(grid)
    >>> SM.grid.number_of_cell_rows
    3
    >>> SM.grid.number_of_cell_columns
    2
    >>> SM.grid is grid
    True
    >>> import numpy as np
    >>> np.allclose(grid.at_cell['soil_moisture__saturation_fraction'], 0.)
    True
    >>> grid['cell']['surface__potential_evapotranspiration_rate']= np.array([
    ...            0.2554777, 0.2554777 , 0.22110221, 0.22110221,
    ...            0.24813062, 0.24813062])
    >>> grid['cell']['soil_moisture__initial_saturation_fraction']= (
    ...        0.75 * np.ones(grid.number_of_cells))
    >>> grid['cell']['vegetation__live_leaf_area_index']= (
    ...        2. * np.ones(grid.number_of_cells))
    >>> grid['cell']['vegetation__cover_fraction']= (
    ...        np.ones(grid.number_of_cells))
    >>> current_time = 0.5
    >>> grid['cell']['rainfall__daily_depth'] = (
    ...        25. * np.ones(grid.number_of_cells))
    >>> current_time = SM.update(current_time)
    >>> np.allclose(grid.at_cell['soil_moisture__saturation_fraction'], 0.)
    False
    """
    _name = 'Soil Moisture'

    _input_var_names = (
        'vegetation__cover_fraction',
        'vegetation__live_leaf_area_index',
        'surface__potential_evapotranspiration_rate',
        'soil_moisture__initial_saturation_fraction',
        'vegetation__plant_functional_type',
        'rainfall__daily_depth',
    )

    _output_var_names = (
        'vegetation__water_stress',
        'soil_moisture__saturation_fraction',
        'soil_moisture__root_zone_leakage',
        'surface__runoff',
        'surface__evapotranspiration',
    )

    _var_units = {
        'vegetation__cover_fraction': 'None',
        'vegetation__live_leaf_area_index': 'None',
        'surface__potential_evapotranspiration_rate': 'mm',
        'vegetation__plant_functional_type': 'None',
        'vegetation__water_stress': 'None',
        'soil_moisture__saturation_fraction': 'None',
        'soil_moisture__initial_saturation_fraction': 'None',
        'soil_moisture__root_zone_leakage': 'mm',
        'surface__runoff': 'mm',
        'surface__evapotranspiration': 'mm',
        'rainfall__daily_depth': 'mm',
    }

    _var_mapping = {
        'vegetation__cover_fraction': 'cell',
        'vegetation__live_leaf_area_index': 'cell',
        'surface__potential_evapotranspiration_rate': 'cell',
        'vegetation__plant_functional_type': 'cell',
        'vegetation__water_stress': 'cell',
        'soil_moisture__saturation_fraction': 'cell',
        'soil_moisture__initial_saturation_fraction': 'cell',
        'soil_moisture__root_zone_leakage': 'cell',
        'surface__runoff': 'cell',
        'surface__evapotranspiration': 'cell',
        'rainfall__daily_depth': 'cell',
    }

    _var_doc = {
        'vegetation__cover_fraction':
            'fraction of land covered by vegetation',
        'vegetation__live_leaf_area_index':
            'one-sided green leaf area per unit ground surface area',
        'surface__potential_evapotranspiration_rate':
            'potential sum of evaporation and plant transpiration',
        'vegetation__plant_functional_type':
            'classification of plants (int), grass=0, shrub=1, tree=2, \
             bare=3, shrub_seedling=4, tree_seedling=5',
        'vegetation__water_stress':
            'parameter that represents nonlinear effects of water deficit \
             on plants',
        'soil_moisture__saturation_fraction':
            'relative volumetric water content (theta) - limits=[0,1]',
        'soil_moisture__initial_saturation_fraction':
            'initial soil_moisture__saturation_fraction',
        'soil_moisture__root_zone_leakage':
            'leakage of water into deeper portions of the soil not accessible \
             to the plant',
        'surface__runoff':
            'runoff from ground surface',
        'surface__evapotranspiration':
            'actual sum of evaporation and plant transpiration',
        'rainfall__daily_depth':
            'Rain in (mm) as a field, allowing spatio-temporal soil moisture \
             saturation analysis.',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, runon=0., f_bare=0.7, soil_ew=0.1,
                 intercept_cap_grass=1., zr_grass=0.3, I_B_grass=20.,
                 I_V_grass=24., pc_grass=0.43, fc_grass=0.56, sc_grass=0.33,
                 wp_grass=0.13, hgw_grass=0.1, beta_grass=13.8,
                 LAI_max_grass=2., LAIR_max_grass=2.88,
                 intercept_cap_shrub=1.5, zr_shrub=0.5, I_B_shrub=20.,
                 I_V_shrub=40., pc_shrub=0.43, fc_shrub=0.56, sc_shrub=0.24,
                 wp_shrub=0.13, hgw_shrub=0.1, beta_shrub=13.8,
                 LAI_max_shrub=2., LAIR_max_shrub=2.,
                 intercept_cap_tree=2., zr_tree=1.3, I_B_tree=20.,
                 I_V_tree=40., pc_tree=0.43, fc_tree=0.56, sc_tree=0.22,
                 wp_tree=0.15, hgw_tree=0.1, beta_tree=13.8,
                 LAI_max_tree=4., LAIR_max_tree=4.,
                 intercept_cap_bare=1., zr_bare=0.15, I_B_bare=20.,
                 I_V_bare=20., pc_bare=0.43, fc_bare=0.56, sc_bare=0.33,
                 wp_bare=0.13, hgw_bare=0.1, beta_bare=13.8,
                 LAI_max_bare=0.01, LAIR_max_bare=0.01, **kwds):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        runon: float, optional
            Runon from higher elevation (mm).
        f_bare: float, optional
            Fraction to partition PET for bare soil (None).
        soil_ew: float, optional
            Residual Evaporation after wilting (mm/day).
        intercept_cap: float, optional
            Plant Functional Type (PFT) specific full canopy interception
            capacity.
        zr: float, optional
            Root depth (m).
        I_B: float, optional
            Infiltration capacity of bare soil (mm/h).
        I_V: float, optional
            Infiltration capacity of vegetated soil (mm/h).
        pc: float, optional
            Soil porosity (None).
        fc: float, optional
            Soil saturation degree at field capacity (None).
        sc: float, optional
            Soil saturation degree at stomatal closure (None).
        wp: float, optional
            Soil saturation degree at wilting point (None).
        hgw: float, optional
            Soil saturation degree at hygroscopic point (None).
        beta: float, optional
            Deep percolation constant = 2*b+3 where b is
            water retention (None).
        LAI_max: float, optional
            Maximum leaf area index (m^2/m^2).
        LAIR_max: float, optional
            Reference leaf area index (m^2/m^2).
        """
        self._method = kwds.pop('method', 'Grid')

        assert_method_is_valid(self._method)

        super(SoilMoisture, self).__init__(grid)

        self.initialize(runon=runon, f_bare=f_bare, soil_ew=soil_ew,
                        intercept_cap_grass=intercept_cap_grass,
                        zr_grass=zr_grass, I_B_grass=I_B_grass,
                        I_V_grass=I_V_grass, pc_grass=pc_grass,
                        fc_grass=fc_grass, sc_grass=sc_grass,
                        wp_grass=wp_grass, hgw_grass=hgw_grass,
                        beta_grass=beta_grass, LAI_max_grass=LAI_max_grass,
                        LAIR_max_grass=LAIR_max_grass,
                        intercept_cap_shrub=intercept_cap_shrub,
                        zr_shrub=zr_shrub, I_B_shrub=I_B_shrub,
                        I_V_shrub=I_V_shrub, pc_shrub=pc_shrub,
                        fc_shrub=fc_shrub, sc_shrub=sc_shrub,
                        wp_shrub=wp_shrub, hgw_shrub=hgw_shrub,
                        beta_shrub=beta_shrub, LAI_max_shrub=LAI_max_shrub,
                        LAIR_max_shrub=LAIR_max_shrub,
                        intercept_cap_tree=intercept_cap_tree, zr_tree=zr_tree,
                        I_B_tree=I_B_tree, I_V_tree=I_V_tree, pc_tree=pc_tree,
                        fc_tree=fc_tree, sc_tree=sc_tree, wp_tree=wp_tree,
                        hgw_tree=hgw_tree, beta_tree=beta_tree,
                        LAI_max_tree=LAI_max_tree, LAIR_max_tree=LAIR_max_tree,
                        intercept_cap_bare=intercept_cap_bare, zr_bare=zr_bare,
                        I_B_bare=I_B_bare, I_V_bare=I_V_bare, pc_bare=pc_bare,
                        fc_bare=fc_bare, sc_bare=sc_bare, wp_bare=wp_bare,
                        hgw_bare=hgw_bare, beta_bare=beta_bare,
                        LAI_max_bare=LAI_max_bare,
                        LAIR_max_bare=LAIR_max_bare, **kwds)

        for name in self._input_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_cell:
                self.grid.add_zeros('cell', name, units=self._var_units[name])

        self._nodal_values = self.grid['node']

        self._cell_values = self.grid['cell']

    def initialize(self, runon=0., f_bare=0.7, soil_ew=0.1,
                   intercept_cap_grass=1., zr_grass=0.3, I_B_grass=20.,
                   I_V_grass=24., pc_grass=0.43, fc_grass=0.56, sc_grass=0.33,
                   wp_grass=0.13, hgw_grass=0.1, beta_grass=13.8,
                   LAI_max_grass=2., LAIR_max_grass=2.88,
                   intercept_cap_shrub=1.5, zr_shrub=0.5, I_B_shrub=20.,
                   I_V_shrub=40., pc_shrub=0.43, fc_shrub=0.56, sc_shrub=0.24,
                   wp_shrub=0.13, hgw_shrub=0.1, beta_shrub=13.8,
                   LAI_max_shrub=2., LAIR_max_shrub=2.,
                   intercept_cap_tree=2., zr_tree=1.3, I_B_tree=20.,
                   I_V_tree=40., pc_tree=0.43, fc_tree=0.56, sc_tree=0.22,
                   wp_tree=0.15, hgw_tree=0.1, beta_tree=13.8,
                   LAI_max_tree=4., LAIR_max_tree=4.,
                   intercept_cap_bare=1., zr_bare=0.15, I_B_bare=20.,
                   I_V_bare=20., pc_bare=0.43, fc_bare=0.56, sc_bare=0.33,
                   wp_bare=0.13, hgw_bare=0.1, beta_bare=13.8,
                   LAI_max_bare=0.01, LAIR_max_bare=0.01, **kwds):
        # GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
        # SHRUBSEEDLING = 4; TREESEEDLING = 5
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        runon: float, optional
            Runon from higher elevation (mm).
        f_bare: float, optional
            Fraction to partition PET for bare soil (None).
        soil_ew: float, optional
            Residual Evaporation after wilting (mm/day).
        intercept_cap: float, optional
            Plant Functional Type (PFT) specific full canopy interception
            capacity.
        zr: float, optional
            Root depth (m).
        I_B: float, optional
            Infiltration capacity of bare soil (mm/h).
        I_V: float, optional
            Infiltration capacity of vegetated soil (mm/h).
        pc: float, optional
            Soil porosity (None).
        fc: float, optional
            Soil saturation degree at field capacity (None).
        sc: float, optional
            Soil saturation degree at stomatal closure (None).
        wp: float, optional
            Soil saturation degree at wilting point (None).
        hgw: float, optional
            Soil saturation degree at hygroscopic point (None).
        beta: float, optional
            Deep percolation constant = 2*b+3 where b is
            water retention (None).
        parameter (None)
        LAI_max: float, optional
            Maximum leaf area index (m^2/m^2).
        LAIR_max: float, optional
            Reference leaf area index (m^2/m^2).
        """

        self._vegtype = self.grid['cell']['vegetation__plant_functional_type']
        self._runon = runon
        self._fbare = f_bare
        self._interception_cap = np.choose(self._vegtype, [
             intercept_cap_grass, intercept_cap_shrub, intercept_cap_tree,
             intercept_cap_bare, intercept_cap_shrub, intercept_cap_tree])

        self._zr = np.choose(self._vegtype, [
             zr_grass, zr_shrub, zr_tree, zr_bare, zr_shrub, zr_tree])

        self._soil_Ib = np.choose(self._vegtype, [
             I_B_grass, I_B_shrub, I_B_tree, I_B_bare, I_B_shrub, I_B_tree])

        self._soil_Iv = np.choose(self._vegtype, [
             I_V_grass, I_V_shrub, I_V_tree, I_V_bare, I_V_shrub, I_V_tree])

        self._soil_Ew = soil_ew
        self._soil_pc = np.choose(self._vegtype, [
             pc_grass, pc_shrub, pc_tree, pc_bare, pc_shrub, pc_tree])

        self._soil_fc = np.choose(self._vegtype, [
             fc_grass, fc_shrub, fc_tree, fc_bare, fc_shrub, fc_tree])

        self._soil_sc = np.choose(self._vegtype, [
             sc_grass, sc_shrub, sc_tree, sc_bare, sc_shrub, sc_tree])

        self._soil_wp = np.choose(self._vegtype, [
             wp_grass, wp_shrub, wp_tree, wp_bare, wp_shrub, wp_tree])

        self._soil_hgw = np.choose(self._vegtype, [
             hgw_grass, hgw_shrub, hgw_tree, hgw_bare, hgw_shrub, hgw_tree])

        self._soil_beta = np.choose(self._vegtype, [
             beta_grass, beta_shrub, beta_tree,
             beta_bare, beta_shrub, beta_tree])

        self._LAI_max = np.choose(self._vegtype, [
             LAI_max_grass, LAI_max_shrub, LAI_max_tree,
             LAI_max_bare, LAI_max_shrub, LAI_max_tree])

        self._LAIR_max = np.choose(self._vegtype, [
             LAIR_max_grass, LAIR_max_shrub, LAIR_max_tree,
             LAIR_max_bare, LAIR_max_shrub, LAIR_max_tree])

    def update(self, current_time, Tb=24., Tr=0., **kwds):
        """
        Update fields with current loading conditions.

        Parameters
        ----------
        current_time: float
            Current time (years).
        Tr: float, optional
            Storm duration (hours).
        Tb: float, optional
            Inter-storm duration (hours).
        """
        P_ = self._cell_values['rainfall__daily_depth']
        self._PET = \
            self._cell_values['surface__potential_evapotranspiration_rate']
        self._SO = \
            self._cell_values['soil_moisture__initial_saturation_fraction']
        self._vegcover = self._cell_values['vegetation__cover_fraction']
        self._water_stress = self._cell_values['vegetation__water_stress']
        self._S = self._cell_values['soil_moisture__saturation_fraction']
        self._D = self._cell_values['soil_moisture__root_zone_leakage']
        self._ETA = self._cell_values['surface__evapotranspiration']
        self._fr = (self._cell_values['vegetation__live_leaf_area_index'] /
                    self._LAIR_max)
        self._runoff = self._cell_values['surface__runoff']
        # LAIl = self._cell_values['vegetation__live_leaf_area_index']
        # LAIt = LAIl+self._cell_values['DeadLeafAreaIndex']
        # if LAIt.all() == 0.:
        #     self._fr = np.zeros(self.grid.number_of_cells)
        # else:
        #     self._fr = (self._vegcover[0]*LAIl/LAIt)
        self._fr[self._fr > 1.] = 1.
        self._Sini = np.zeros(self._SO.shape)
        self._ETmax = np.zeros(self._SO.shape)

        for cell in range(0, self.grid.number_of_cells):
            P = P_[cell]
            # print cell
            s = self._SO[cell]
            fbare = self._fbare
            ZR = self._zr[cell]
            pc = self._soil_pc[cell]
            fc = self._soil_fc[cell]
            scc = self._soil_sc[cell]
            wp = self._soil_wp[cell]
            hgw = self._soil_hgw[cell]
            beta = self._soil_beta[cell]
            if self._vegtype[cell] == 0:   # 0 - GRASS
                sc = scc*self._fr[cell]+(1-self._fr[cell])*fc
            else:
                sc = scc

            Inf_cap = (self._soil_Ib[cell]*(1-self._vegcover[cell]) +
                       self._soil_Iv[cell]*self._vegcover[cell])
            # Infiltration capacity
            Int_cap = min(self._vegcover[cell]*self._interception_cap[cell],
                          P)
            # Interception capacity
            Peff = max(P-Int_cap, 0.)         # Effective precipitation depth
            mu = (Inf_cap/1000.0)/(pc*ZR*(np.exp(beta*(1.-fc))-1.))
            Ep = max((self._PET[cell]*self._fr[cell] +
                     fbare*self._PET[cell]*(1.-self._fr[cell])) -
                     Int_cap, 0.0001)  # mm/d
            self._ETmax[cell] = Ep
            nu = ((Ep / 24.) / 1000.) / (pc*ZR)   # Loss function parameter
            nuw = ((self._soil_Ew/24.)/1000.)/(pc*ZR)
            # Loss function parameter
            sini = self._SO[cell] + ((Peff+self._runon)/(pc*ZR*1000.))

            if sini > 1.:
                self._runoff[cell] = (sini-1.)*pc*ZR*1000.
                # print 'Runoff =', self._runoff
                sini = 1.
            else:
                self._runoff[cell] = 0.

            if sini >= fc:
                tfc = (1./(beta*(mu-nu)))*(beta*(fc-sini) + np.log((
                       nu-mu+mu*np.exp(beta*(sini-fc)))/nu))
                tsc = ((fc-sc)/nu)+tfc
                twp = ((sc-wp)/(nu-nuw))*np.log(nu/nuw)+tsc

                if Tb < tfc:
                    s = abs(sini-(1./beta)*np.log(((nu-mu+mu *
                            np.exp(beta*(sini-fc)))*np.exp(beta*(nu-mu)*Tb) -
                            mu*np.exp(beta*(sini-fc)))/(nu-mu)))

                    self._D[cell] = ((pc*ZR*1000.)*(sini-s))-(Tb*(Ep/24.))
                    self._ETA[cell] = (Tb*(Ep/24.))

                elif Tb >= tfc and Tb < tsc:
                    s = fc-(nu*(Tb-tfc))
                    self._D[cell] = ((pc*ZR*1000.)*(sini-fc))-((tfc)*(Ep/24.))
                    self._ETA[cell] = (Tb*(Ep/24.))

                elif Tb >= tsc and Tb < twp:
                    s = (wp+(sc-wp)*((nu/(nu-nuw))*np.exp((-1)*((nu-nuw) /
                         (sc-wp))*(Tb-tsc))-(nuw/(nu-nuw))))
                    self._D[cell] = ((pc*ZR*1000.)*(sini-fc))-(tfc*Ep/24.)
                    self._ETA[cell] = (1000.*ZR*pc*(sini-s))-self._D[cell]

                else:
                    s = (hgw+(wp-hgw)*np.exp((-1)*(nuw/(wp-hgw)) *
                         max(Tb-twp, 0.)))
                    self._D[cell] = ((pc*ZR*1000.)*(sini-fc))-(tfc*Ep/24.)
                    self._ETA[cell] = (1000.*ZR*pc*(sini-s))-self._D[cell]

            elif sini < fc and sini >= sc:
                tfc = 0.
                tsc = (sini-sc)/nu
                twp = ((sc-wp)/(nu-nuw))*np.log(nu/nuw)+tsc

                if Tb < tsc:
                    s = sini - nu*Tb
                    self._D[cell] = 0.
                    self._ETA[cell] = 1000.*ZR*pc*(sini-s)

                elif Tb >= tsc and Tb < twp:
                    s = (wp+(sc-wp)*((nu/(nu-nuw))*np.exp((-1) *
                         ((nu-nuw)/(sc-wp))*(Tb-tsc))-(nuw/(nu-nuw))))
                    self._D[cell] = 0
                    self._ETA[cell] = (1000.*ZR*pc*(sini-s))

                else:
                    s = hgw+(wp-hgw)*np.exp((-1)*(nuw/(wp-hgw))*(Tb-twp))
                    self._D[cell] = 0.
                    self._ETA[cell] = (1000.*ZR*pc*(sini-s))

            elif sini < sc and sini >= wp:
                tfc = 0
                tsc = 0
                twp = (((sc-wp)/(nu-nuw))*np.log(1+(nu-nuw)*(sini-wp) /
                       (nuw*(sc-wp))))

                if Tb < twp:
                    s = (wp+((sc-wp)/(nu-nuw))*((np.exp((-1)*((nu-nuw) /
                         (sc-wp))*Tb))*(nuw+((nu-nuw)/(sc-wp))*(sini-wp))-nuw))
                    self._D[cell] = 0.
                    self._ETA[cell] = (1000.*ZR*pc*(sini-s))

                else:
                    s = hgw+(wp-hgw)*np.exp((-1)*(nuw/(wp-hgw))*(Tb-twp))
                    self._D[cell] = 0.
                    self._ETA[cell] = (1000.*ZR*pc*(sini-s))

            else:
                tfc = 0.
                tsc = 0.
                twp = 0.

                s = hgw+(sini-hgw)*np.exp((-1)*(nuw/(wp-hgw))*Tb)
                self._D[cell] = 0.
                self._ETA[cell] = (1000.*ZR*pc*(sini-s))

            self._water_stress[cell] = min(((max(((sc - (s+sini)/2.) /
                                           (sc - wp)), 0.))**4.), 1.0)
            self._S[cell] = s
            self._SO[cell] = s
            self._Sini[cell] = sini

        current_time += (Tb+Tr)/(24.*365.25)
        return current_time
