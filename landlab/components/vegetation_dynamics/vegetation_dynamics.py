import numpy as np

from landlab import Component

_VALID_METHODS = {"Grid"}


def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class Vegetation(Component):
    """Landlab component that simulates net primary productivity, biomass and
    leaf area index at each cell based on inputs of root-zone average soil
    moisture.

    Ref: Zhou, X., Istanbulluoglu, E., & Vivoni, E. R. (2013). Modeling the
    ecohydrological role of aspect controlled radiation on tree grass shrub
    coexistence in a semiarid climate. Water Resources Research,
    49(5), 2872-2895.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import Vegetation

    Create a grid on which to simulate vegetation dynamics.

    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))

    The grid will need some input data. To check the names of the fields
    that provide the input to this component, use the *input_var_names*
    class property.

    >>> sorted(Vegetation.input_var_names)
    ['surface__evapotranspiration',
     'surface__potential_evapotranspiration_30day_mean',
     'surface__potential_evapotranspiration_rate',
     'vegetation__plant_functional_type',
     'vegetation__water_stress']

    >>> sorted(Vegetation.units)
    [('surface__evapotranspiration', 'mm'),
     ('surface__potential_evapotranspiration_30day_mean', 'mm'),
     ('surface__potential_evapotranspiration_rate', 'mm'),
     ('vegetation__cover_fraction', 'None'),
     ('vegetation__dead_biomass', 'g m^-2 d^-1'),
     ('vegetation__dead_leaf_area_index', 'None'),
     ('vegetation__live_biomass', 'g m^-2 d^-1'),
     ('vegetation__live_leaf_area_index', 'None'),
     ('vegetation__plant_functional_type', 'None'),
     ('vegetation__water_stress', 'None')]

    Provide all the input fields.

    >>> grid["cell"]["vegetation__plant_functional_type"] = np.zeros(
    ...     grid.number_of_cells, dtype=int
    ... )
    >>> grid["cell"]["surface__evapotranspiration"] = 0.2 * np.ones(
    ...     grid.number_of_cells
    ... )
    >>> grid["cell"]["surface__potential_evapotranspiration_rate"] = np.array(
    ...     [0.25547770, 0.25547770, 0.22110221, 0.22110221, 0.24813062, 0.24813062]
    ... )
    >>> grid["cell"]["surface__potential_evapotranspiration_30day_mean"] = np.array(
    ...     [0.25547770, 0.25547770, 0.22110221, 0.22110221, 0.24813062, 0.24813062]
    ... )
    >>> grid["cell"]["vegetation__water_stress"] = 0.01 * np.ones(grid.number_of_cells)

    Instantiate the 'Vegetation' component.

    >>> Veg = Vegetation(grid)

    >>> Veg.grid.number_of_cell_rows
    3
    >>> Veg.grid.number_of_cell_columns
    2
    >>> Veg.grid is grid
    True
    >>> import numpy as np
    >>> sorted(Vegetation.output_var_names)
    ['vegetation__cover_fraction',
     'vegetation__dead_biomass',
     'vegetation__dead_leaf_area_index',
     'vegetation__live_biomass',
     'vegetation__live_leaf_area_index']

    >>> np.all(grid.at_cell["vegetation__live_leaf_area_index"] == 0.0)
    True

    >>> Veg.update()

    >>> np.all(grid.at_cell["vegetation__live_leaf_area_index"] == 0.0)
    False

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Zhou, X., Istanbulluoglu, E., and Vivoni, E. R.: Modeling the
    ecohydrological role of aspect-controlled radiation on tree-grass-shrub
    coexistence in a semiarid climate, Water Resour. Res., 49, 2872– 2895,
    doi:10.1002/wrcr.20259, 2013.

    """

    _name = "Vegetation"

    _unit_agnostic = False

    _info = {
        "surface__evapotranspiration": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "mm",
            "mapping": "cell",
            "doc": "actual sum of evaporation and plant transpiration",
        },
        "surface__potential_evapotranspiration_30day_mean": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "mm",
            "mapping": "cell",
            "doc": "30 day mean of surface__potential_evapotranspiration",
        },
        "surface__potential_evapotranspiration_rate": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "mm",
            "mapping": "cell",
            "doc": "potential sum of evaporation and potential transpiration",
        },
        "vegetation__cover_fraction": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": "fraction of land covered by vegetation",
        },
        "vegetation__dead_biomass": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "g m^-2 d^-1",
            "mapping": "cell",
            "doc": (
                "weight of dead organic mass per unit area - measured in terms "
                "of dry matter"
            ),
        },
        "vegetation__dead_leaf_area_index": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": "one-sided dead leaf area per unit ground surface area",
        },
        "vegetation__live_biomass": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "g m^-2 d^-1",
            "mapping": "cell",
            "doc": (
                "weight of green organic mass per unit area - measured in terms "
                "of dry matter"
            ),
        },
        "vegetation__live_leaf_area_index": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": "one-sided green leaf area per unit ground surface area",
        },
        "vegetation__plant_functional_type": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": (
                "classification of plants (int), grass=0, shrub=1, tree=2, bare=3, "
                "shrub_seedling=4, tree_seedling=5"
            ),
        },
        "vegetation__water_stress": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": "parameter that represents nonlinear effects of water deficit on plants",
        },
    }

    def __init__(
        self,
        grid,
        Blive_init=102.0,
        Bdead_init=450.0,
        ETthreshold_up=3.8,
        ETthreshold_down=6.8,
        Tdmax=10.0,
        w=0.55,
        WUE_grass=0.01,
        LAI_max_grass=2.0,
        cb_grass=0.0047,
        cd_grass=0.009,
        ksg_grass=0.012,
        kdd_grass=0.013,
        kws_grass=0.02,
        WUE_shrub=0.0025,
        LAI_max_shrub=2.0,
        cb_shrub=0.004,
        cd_shrub=0.01,
        ksg_shrub=0.002,
        kdd_shrub=0.013,
        kws_shrub=0.02,
        WUE_tree=0.0045,
        LAI_max_tree=4.0,
        cb_tree=0.004,
        cd_tree=0.01,
        ksg_tree=0.002,
        kdd_tree=0.013,
        kws_tree=0.01,
        WUE_bare=0.01,
        LAI_max_bare=0.01,
        cb_bare=0.0047,
        cd_bare=0.009,
        ksg_bare=0.012,
        kdd_bare=0.013,
        kws_bare=0.02,
        method="Grid",
        PETthreshold_switch=0,
        Tb=24.0,
        Tr=0.01,
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        Blive_init: float, optional
            Initial value for vegetation__live_biomass. Converted to field.
        Bdead_init: float, optional
            Initial value for vegetation__dead_biomass. Coverted to field.
        ETthreshold_up: float, optional
            Potential Evapotranspiration (PET) threshold for
            growing season (mm/d).
        ETthreshold_down: float, optional
            PET threshold for dormant season (mm/d).
        Tdmax: float, optional
            Constant for dead biomass loss adjustment (mm/d).
        w: float, optional
            Conversion factor of CO2 to dry biomass (Kg DM/Kg CO2).
        WUE: float, optional
            Water Use Efficiency - ratio of water used in plant water
            lost by the plant through transpiration (KgCO2Kg-1H2O).
        LAI_max: float, optional
            Maximum leaf area index (m2/m2).
        cb: float, optional
            Specific leaf area for green/live biomass (m2 leaf g-1 DM).
        cd: float, optional
            Specific leaf area for dead biomass (m2 leaf g-1 DM).
        ksg: float, optional
            Senescence coefficient of green/live biomass (d-1).
        kdd: float, optional
            Decay coefficient of aboveground dead biomass (d-1).
        kws: float, optional
            Maximum drought induced foliage loss rate (d-1).
        method: str
            Method name.
        Tr: float, optional
            Storm duration (hours).
        Tb: float, optional
            Inter-storm duration (hours).
        PETthreshold_switch: int, optional
            Flag to indiate the PET threshold. This controls whether the
            threshold is for growth (1) or dormancy (any other value).
        """
        super().__init__(grid)

        self.Tb = Tb
        self.Tr = Tr
        self.PETthreshold_switch = PETthreshold_switch

        self._method = method

        assert_method_is_valid(self._method)

        self.initialize(
            Blive_init=Blive_init,
            Bdead_init=Bdead_init,
            ETthreshold_up=ETthreshold_up,
            ETthreshold_down=ETthreshold_down,
            Tdmax=Tdmax,
            w=w,
            WUE_grass=WUE_grass,
            LAI_max_grass=LAI_max_grass,
            cb_grass=cb_grass,
            cd_grass=cd_grass,
            ksg_grass=ksg_grass,
            kdd_grass=kdd_grass,
            kws_grass=kws_grass,
            WUE_shrub=WUE_shrub,
            LAI_max_shrub=LAI_max_shrub,
            cb_shrub=cb_shrub,
            cd_shrub=cd_shrub,
            ksg_shrub=ksg_shrub,
            kdd_shrub=kdd_shrub,
            kws_shrub=kws_shrub,
            WUE_tree=WUE_tree,
            LAI_max_tree=LAI_max_tree,
            cb_tree=cb_tree,
            cd_tree=cd_tree,
            ksg_tree=ksg_tree,
            kdd_tree=kdd_tree,
            kws_tree=kws_tree,
            WUE_bare=WUE_bare,
            LAI_max_bare=LAI_max_bare,
            cb_bare=cb_bare,
            cd_bare=cd_bare,
            ksg_bare=ksg_bare,
            kdd_bare=kdd_bare,
            kws_bare=kws_bare,
        )

        self.initialize_output_fields()

        self._cell_values = self._grid["cell"]

        self._Blive_ini = self._Blive_init * np.ones(self._grid.number_of_cells)
        self._Bdead_ini = self._Bdead_init * np.ones(self._grid.number_of_cells)

    @property
    def Tb(self):
        """Storm duration (hours)."""
        return self._Tb

    @Tb.setter
    def Tb(self, Tb):
        if Tb < 0.0:
            raise ValueError("Tb must be non-negative")
        self._Tb = Tb

    @property
    def Tr(self):
        """Inter-storm duration (hours)."""
        return self._Tr

    @Tr.setter
    def Tr(self, Tr):
        assert Tr >= 0
        self._Tr = Tr

    @property
    def PETthreshold_switch(self):
        """Flag to indiate the PET threshold.

        This controls whether the threshold is for growth (1) or
        dormancy (any other value).
        """
        return self._PETthreshold_switch

    @PETthreshold_switch.setter
    def PETthreshold_switch(self, PETthreshold_switch):
        self._PETthreshold_switch = PETthreshold_switch

    def initialize(
        self,
        Blive_init=102.0,
        Bdead_init=450.0,
        ETthreshold_up=3.8,
        ETthreshold_down=6.8,
        Tdmax=10.0,
        w=0.55,
        WUE_grass=0.01,
        LAI_max_grass=2.0,
        cb_grass=0.0047,
        cd_grass=0.009,
        ksg_grass=0.012,
        kdd_grass=0.013,
        kws_grass=0.02,
        WUE_shrub=0.0025,
        LAI_max_shrub=2.0,
        cb_shrub=0.004,
        cd_shrub=0.01,
        ksg_shrub=0.002,
        kdd_shrub=0.013,
        kws_shrub=0.02,
        WUE_tree=0.0045,
        LAI_max_tree=4.0,
        cb_tree=0.004,
        cd_tree=0.01,
        ksg_tree=0.002,
        kdd_tree=0.013,
        kws_tree=0.01,
        WUE_bare=0.01,
        LAI_max_bare=0.01,
        cb_bare=0.0047,
        cd_bare=0.009,
        ksg_bare=0.012,
        kdd_bare=0.013,
        kws_bare=0.02,
    ):
        # GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
        # SHRUBSEEDLING = 4; TREESEEDLING = 5
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        Blive_init: float, optional
            Initial value for vegetation__live_biomass. Converted to field.
        Bdead_init: float, optional
            Initial value for vegetation__dead_biomass. Coverted to field.
        ETthreshold_up: float, optional
            Potential Evapotranspiration (PET) threshold for
            growing season (mm/d).
        ETthreshold_down: float, optional
            PET threshold for dormant season (mm/d).
        Tdmax: float, optional
            Constant for dead biomass loss adjustment (mm/d).
        w: float, optional
            Conversion factor of CO2 to dry biomass (Kg DM/Kg CO2).
        WUE: float, optional
            Water Use Efficiency - ratio of water used in plant water
            lost by the plant through transpiration (KgCO2Kg-1H2O).
        LAI_max: float, optional
            Maximum leaf area index (m2/m2).
        cb: float, optional
            Specific leaf area for green/live biomass (m2 leaf g-1 DM).
        cd: float, optional
            Specific leaf area for dead biomass (m2 leaf g-1 DM).
        ksg: float, optional
            Senescence coefficient of green/live biomass (d-1).
        kdd: float, optional
            Decay coefficient of aboveground dead biomass (d-1).
        kws: float, optional
            Maximum drought induced foliage loss rate (d-1).
        """
        self._vegtype = self._grid["cell"]["vegetation__plant_functional_type"]
        self._WUE = np.choose(
            self._vegtype,
            [WUE_grass, WUE_shrub, WUE_tree, WUE_bare, WUE_shrub, WUE_tree],
        )
        # Water Use Efficiency  KgCO2kg-1H2O
        self._LAI_max = np.choose(
            self._vegtype,
            [
                LAI_max_grass,
                LAI_max_shrub,
                LAI_max_tree,
                LAI_max_bare,
                LAI_max_shrub,
                LAI_max_tree,
            ],
        )
        # Maximum leaf area index (m2/m2)
        self._cb = np.choose(
            self._vegtype, [cb_grass, cb_shrub, cb_tree, cb_bare, cb_shrub, cb_tree]
        )
        # Specific leaf area for green/live biomass (m2 leaf g-1 DM)
        self._cd = np.choose(
            self._vegtype, [cd_grass, cd_shrub, cd_tree, cd_bare, cd_shrub, cd_tree]
        )
        # Specific leaf area for dead biomass (m2 leaf g-1 DM)
        self._ksg = np.choose(
            self._vegtype,
            [ksg_grass, ksg_shrub, ksg_tree, ksg_bare, ksg_shrub, ksg_tree],
        )
        # Senescence coefficient of green/live biomass (d-1)
        self._kdd = np.choose(
            self._vegtype,
            [kdd_grass, kdd_shrub, kdd_tree, kdd_bare, kdd_shrub, kdd_tree],
        )
        # Decay coefficient of aboveground dead biomass (d-1)
        self._kws = np.choose(
            self._vegtype,
            [kws_grass, kws_shrub, kws_tree, kws_bare, kws_shrub, kws_tree],
        )
        # Maximum drought induced foliage loss rates (d-1)
        self._Blive_init = Blive_init
        self._Bdead_init = Bdead_init
        self._ETthresholdup = ETthreshold_up  # Growth threshold (mm/d)
        self._ETthresholddown = ETthreshold_down  # Dormancy threshold (mm/d)
        self._Tdmax = Tdmax  # Constant for dead biomass loss adjustment
        self._w = w  # Conversion factor of CO2 to dry biomass

        self._Blive_ini = self._Blive_init * np.ones(self._grid.number_of_cells)
        self._Bdead_ini = self._Bdead_init * np.ones(self._grid.number_of_cells)

    def update(self):
        """Update fields with current loading conditions.

        This method looks to the properties ``PETthreshold_switch``,
        ``Tb``, and ``Tr`` and uses their values to calculate the new
        field values.
        """
        PETthreshold_ = self._PETthreshold_switch
        Tb = self._Tb
        Tr = self._Tr

        PET = self._cell_values["surface__potential_evapotranspiration_rate"]
        PET30_ = self._cell_values["surface__potential_evapotranspiration_30day_mean"]
        ActualET = self._cell_values["surface__evapotranspiration"]
        Water_stress = self._cell_values["vegetation__water_stress"]

        self._LAIlive = self._cell_values["vegetation__live_leaf_area_index"]
        self._LAIdead = self._cell_values["vegetation__dead_leaf_area_index"]
        self._Blive = self._cell_values["vegetation__live_biomass"]
        self._Bdead = self._cell_values["vegetation__dead_biomass"]
        self._VegCov = self._cell_values["vegetation__cover_fraction"]

        if PETthreshold_ == 1:
            PETthreshold = self._ETthresholdup
        else:
            PETthreshold = self._ETthresholddown

        for cell in range(0, self._grid.number_of_cells):
            WUE = self._WUE[cell]
            LAImax = self._LAI_max[cell]
            cb = self._cb[cell]
            cd = self._cd[cell]
            ksg = self._ksg[cell]
            kdd = self._kdd[cell]
            kws = self._kws[cell]
            # ETdmax = self._ETdmax[cell]
            LAIlive = min(cb * self._Blive_ini[cell], LAImax)
            LAIdead = min(cd * self._Bdead_ini[cell], (LAImax - LAIlive))
            NPP = max((ActualET[cell] / (Tb + Tr)) * WUE * 24.0 * self._w * 1000, 0.001)

            if self._vegtype[cell] == 0:
                if PET30_[cell] > PETthreshold:
                    # Growing Season
                    Bmax = (LAImax - LAIdead) / cb
                    Yconst = 1 / (
                        (1 / Bmax) + (((kws * Water_stress[cell]) + ksg) / NPP)
                    )
                    Blive = (self._Blive_ini[cell] - Yconst) * np.exp(
                        -(NPP / Yconst) * ((Tb + Tr) / 24.0)
                    ) + Yconst
                    Bdead = (
                        self._Bdead_ini[cell]
                        + (Blive - max(Blive * np.exp(-1 * ksg * Tb / 24.0), 0.00001))
                    ) * np.exp(-1 * kdd * min(PET[cell] / self._Tdmax, 1.0) * Tb / 24.0)
                else:  # Senescense
                    Blive = max(
                        self._Blive_ini[cell] * np.exp((-2) * ksg * Tb / 24.0), 1
                    )
                    Bdead = max(
                        (
                            self._Bdead_ini[cell]
                            + (
                                self._Blive_ini[cell]
                                - (
                                    max(
                                        self._Blive_ini[cell]
                                        * np.exp((-2) * ksg * Tb / 24.0),
                                        0.000001,
                                    )
                                )
                            )
                            * np.exp(
                                (-1)
                                * kdd
                                * min(PET[cell] / self._Tdmax, 1.0)
                                * Tb
                                / 24.0
                            ),
                            0.0,
                        )
                    )

            elif self._vegtype[cell] == 3:
                Blive = 0.0
                Bdead = 0.0

            else:
                Bmax = LAImax / cb
                Yconst = 1.0 / (
                    (1.0 / Bmax) + (((kws * Water_stress[cell]) + ksg) / NPP)
                )
                Blive = (self._Blive_ini[cell] - Yconst) * np.exp(
                    -(NPP / Yconst) * ((Tb + Tr) / 24.0)
                ) + Yconst
                Bdead = (
                    self._Bdead_ini[cell]
                    + (Blive - max(Blive * np.exp(-ksg * Tb / 24.0), 0.00001))
                ) * np.exp(-kdd * min(PET[cell] / self._Tdmax, 1.0) * Tb / 24.0)

            LAIlive = min(cb * (Blive + self._Blive_ini[cell]) / 2.0, LAImax)
            LAIdead = min(
                cd * (Bdead + self._Bdead_ini[cell]) / 2.0, (LAImax - LAIlive)
            )
            if self._vegtype[cell] == 0:
                Vt = 1.0 - np.exp(-0.75 * (LAIlive + LAIdead))
            else:
                # Vt = 1 - np.exp(-0.75 * LAIlive)
                Vt = 1.0

            self._LAIlive[cell] = LAIlive
            self._LAIdead[cell] = LAIdead
            self._VegCov[cell] = Vt
            self._Blive[cell] = Blive
            self._Bdead[cell] = Bdead

        self._Blive_ini = self._Blive
        self._Bdead_ini = self._Bdead
