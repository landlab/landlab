import numpy as np

from landlab import Component

_VALID_METHODS = {"Grid"}
GRASS = 0
SHRUB = 1
TREE = 2
BARE = 3
SHRUBSEEDLING = 4
TREESEEDLING = 5


def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError("%s: Invalid method name" % method)


class VegCA(Component):
    """Landlab component that simulates inter-species plant competition using a
    2D cellular automata model.

    This code is based on Cellular Automata Tree Grass Shrub Simulator (CATGraSS).
    It simulates spatial competition of multiple plant functional types through
    establishment and mortality. In the current code, tree, grass and
    shrubs are used.

    Ref: Zhou, X., Istanbulluoglu, E., & Vivoni, E. R. (2013). Modeling the
    ecohydrological role of aspect controlled radiation on tree grass shrub
    coexistence in a semiarid climate. Water Resources Research,
    49(5), 2872-2895.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import VegCA
    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))
    >>> VegCA.name
    'Cellular Automata Plant Competition'
    >>> sorted(VegCA.output_var_names)
    ['plant__age', 'plant__live_index']
    >>> sorted(VegCA.units)
    [('plant__age', 'Years'),
     ('plant__live_index', 'None'),
     ('vegetation__cumulative_water_stress', 'None'),
     ('vegetation__plant_functional_type', 'None')]
    >>> grid["cell"]["vegetation__plant_functional_type"] = np.arange(
    ...     0, grid.number_of_cells, dtype=int
    ... )
    >>> grid["cell"]["vegetation__cumulative_water_stress"] = np.ones(
    ...     grid.number_of_cells
    ... )
    >>> ca_veg = VegCA(grid)
    >>> ca_veg.grid.number_of_cell_rows
    3
    >>> ca_veg.grid.number_of_cell_columns
    2
    >>> ca_veg.grid is grid
    True
    >>> import numpy as np
    >>> A = np.copy(grid["cell"]["plant__age"])
    >>> ca_veg.update()
    >>> np.alltrue(grid["cell"]["plant__age"] == A)
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

    _name = "Cellular Automata Plant Competition"

    _unit_agnostic = False

    _info = {
        "plant__age": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "Years",
            "mapping": "cell",
            "doc": "Age of plant",
        },
        "plant__live_index": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": "1 - vegetation__cumulative_water_stress",
        },
        "vegetation__cumulative_water_stress": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": "cumulative vegetation__water_stress over the growing season",
        },
        "vegetation__plant_functional_type": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "None",
            "mapping": "cell",
            "doc": (
                "classification of plants (int), grass=0, shrub=1, tree=2, "
                "bare=3, shrub_seedling=4, tree_seedling=5"
            ),
        },
    }

    def __init__(
        self,
        grid,
        Pemaxg=0.35,
        ING=2.0,
        ThetaGrass=0.62,
        PmbGrass=0.05,
        Pemaxsh=0.2,
        ThetaShrub=0.8,
        PmbShrub=0.01,
        tpmaxShrub=600,
        Pemaxtr=0.25,
        ThetaTree=0.72,
        PmbTree=0.01,
        tpmaxTree=350,
        ThetaShrubSeedling=0.64,
        PmbShrubSeedling=0.03,
        tpmaxShrubSeedling=18,
        ThetaTreeSeedling=0.64,
        PmbTreeSeedling=0.03,
        tpmaxTreeSeedling=18,
        method="Grid",
        Edit_VegCov=True,
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        Pemaxg: float, optional
            Maximal establishment probability of grass.
        ING: float, optional
            Parameter to define allelopathic effect of creosote on grass.
        ThetaGrass: float, optional
            Drought resistance threshold of grass.
        PmbGrass: float, optional
            Background mortality probability of grass.
        Pemaxsh: float, optional
            Maximal establishment probability of shrub.
        ThetaShrub: float, optional
            Drought resistance threshold of shrub.
        PmbShrub: float, optional
            Background mortality probability of shrub.
        tpmaxShrub: float, optional
            Maximum age of shrub (years).
        Pemaxtr: float, optional
            Maximal establishment probability of tree.
        Thetatree: float, optional
            Drought resistance threshold of tree.
        PmbTree: float, optional
            Background mortality probability of tree.
        tpmaxTree: float, optional
            Maximum age of tree (years).
        ThetaShrubSeedling: float, optional
            Drought resistance threshold of shrub seedling.
        PmbShrubSeedling: float, optional
            Background mortality probability of shrub seedling.
        tpmaxShrubSeedling: float, optional
            Maximum age of shrub seedling (years).
        ThetaTreeSeedling: float, optional
            Drought resistance threshold of tree seedling.
        PmbTreeSeedling: float, optional
            Background mortality probability of tree seedling.
        tpmaxTreeSeedling: float, optional
            Maximum age of tree seedling (years).
        method: str, optional
            Method used.
        Edit_VegCov: bool, optional
            If Edit_VegCov=True, an optional field
            'vegetation__boolean_vegetated' will be output, (i.e.) if a cell is
            vegetated the corresponding cell of the field will be 1, otherwise
            it will be 0.

        """
        super().__init__(grid)

        self.Edit_VegCov = Edit_VegCov

        self._Pemaxg = Pemaxg  # Pe-max-grass - max probability
        self._Pemaxsh = Pemaxsh  # Pe-max-shrub
        self._Pemaxtr = Pemaxtr  # Pe-max-tree
        self._INg = ING  # Allelopathic effect on grass from creosotebush
        self._th_g = ThetaGrass  # grass
        self._th_sh = ThetaShrub  # shrub - Creosote
        self._th_tr = ThetaTree  # Juniper pine
        self._th_sh_s = ThetaShrubSeedling  # shrub seedling
        self._th_tr_s = ThetaTreeSeedling  # Juniper pine seedling
        self._Pmb_g = PmbGrass  # Background mortality probability - grass
        self._Pmb_sh = PmbShrub  # shrub
        self._Pmb_tr = PmbTree  # tree
        self._Pmb_sh_s = PmbShrubSeedling  # shrub seedling
        self._Pmb_tr_s = PmbTreeSeedling  # tree seedling
        self._tpmax_sh = tpmaxShrub  # Maximum age - shrub
        self._tpmax_tr = tpmaxTree  # Maximum age - tree
        self._tpmax_sh_s = tpmaxShrubSeedling  # Maximum age - shrub seedling
        self._tpmax_tr_s = tpmaxTreeSeedling  # Maximum age - tree seedling

        self._method = method

        assert_method_is_valid(self._method)

        # if "vegetation__plant_functional_type" not in self._grid.at_cell:
        #     grid["cell"]["vegetation__plant_functional_type"] = np.random.randint(
        #         0, 6, grid.number_of_cells
        #     )

        self.initialize_output_fields()

        self._cell_values = self._grid["cell"]

        VegType = grid["cell"]["vegetation__plant_functional_type"]

        tp = np.zeros(grid.number_of_cells, dtype=int)
        tp[VegType == TREE] = np.random.randint(
            0, self._tpmax_tr, np.where(VegType == TREE)[0].shape
        )
        tp[VegType == SHRUB] = np.random.randint(
            0, self._tpmax_sh, np.where(VegType == SHRUB)[0].shape
        )
        locs_trees = np.where(VegType == TREE)[0]
        locs_shrubs = np.where(VegType == SHRUB)[0]
        VegType[locs_trees[tp[locs_trees] < self._tpmax_tr_s]] = TREESEEDLING
        VegType[locs_shrubs[tp[locs_shrubs] < self._tpmax_sh_s]] = SHRUBSEEDLING
        grid["cell"]["plant__age"] = tp.astype(float)

    @property
    def Edit_VegCov(self):
        """Flag to indicate whether an optional field is created.

        If Edit_VegCov=True, an optional field
        'vegetation__boolean_vegetated' will be output, (i.e.) if a cell
        is vegetated the corresponding cell of the field will be 1,
        otherwise it will be 0.
        """
        return self._Edit_VegCov

    @Edit_VegCov.setter
    def Edit_VegCov(self, Edit_VegCov):
        assert isinstance(Edit_VegCov, bool)
        self._Edit_VegCov = Edit_VegCov

    def update(self, dt=1):
        """Update fields with current loading conditions.

        Parameters
        ----------
        dt: int, optional
            Time elapsed - time step (years).
        """
        self._VegType = self._cell_values["vegetation__plant_functional_type"]
        self._CumWS = self._cell_values["vegetation__cumulative_water_stress"]
        self._live_index = self._cell_values["plant__live_index"]
        self._tp = self._cell_values["plant__age"] + dt

        # Check if shrub and tree seedlings have matured
        shrub_seedlings = np.where(self._VegType == SHRUBSEEDLING)[0]
        tree_seedlings = np.where(self._VegType == TREESEEDLING)[0]
        matured_shrubs = np.where(self._tp[shrub_seedlings] > self._tpmax_sh_s)[0]
        matured_trees = np.where(self._tp[tree_seedlings] > self._tpmax_tr_s)[0]
        self._VegType[shrub_seedlings[matured_shrubs]] = SHRUB
        self._VegType[tree_seedlings[matured_trees]] = TREE
        self._tp[shrub_seedlings[matured_shrubs]] = 0
        self._tp[tree_seedlings[matured_trees]] = 0

        # Establishment
        self._live_index = 1 - self._CumWS  # Plant live index = 1 - WS
        bare_cells = np.where(self._VegType == BARE)[0]
        n_bare = len(bare_cells)
        first_ring = self._grid.looped_neighbors_at_cell[bare_cells]
        second_ring = self._grid.second_ring_looped_neighbors_at_cell[bare_cells]
        veg_type_fr = self._VegType[first_ring]
        veg_type_sr = self._VegType[second_ring]
        Sh_WS_fr = WS_PFT(veg_type_fr, SHRUB, self._live_index[first_ring])
        Tr_WS_fr = WS_PFT(veg_type_fr, TREE, self._live_index[first_ring])
        Tr_WS_sr = WS_PFT(veg_type_sr, TREE, self._live_index[second_ring])

        n = count(veg_type_fr, SHRUB)
        Phi_sh = Sh_WS_fr / 8.0
        Phi_tr = (Tr_WS_fr + Tr_WS_sr / 2.0) / 8.0
        Phi_g = np.mean(self._live_index[np.where(self._VegType == GRASS)])
        Pemaxg = self._Pemaxg * np.ones(n_bare)
        Pemaxsh = self._Pemaxsh * np.ones(n_bare)
        Pemaxtr = self._Pemaxtr * np.ones(n_bare)
        Peg = np.amin(np.vstack((Phi_g / (n * self._INg), Pemaxg)), axis=0)
        Pesh = np.amin(np.vstack((Phi_sh, Pemaxsh)), axis=0)
        Petr = np.amin(np.vstack((Phi_tr, Pemaxtr)), axis=0)
        Select_PFT_E = np.random.choice([GRASS, SHRUBSEEDLING, TREESEEDLING], n_bare)
        # Grass - 0; Shrub Seedling - 4; Tree Seedling - 5
        Pest = np.choose(Select_PFT_E, [Peg, 0, 0, 0, Pesh, Petr])
        # Probability of establishment
        R_Est = np.random.rand(n_bare)
        # Random number for comparison to establish
        Establish = np.int32(np.where(np.greater_equal(Pest, R_Est))[0])
        self._VegType[bare_cells[Establish]] = Select_PFT_E[Establish]
        self._tp[bare_cells[Establish]] = 0

        # Mortality
        plant_cells = np.where(self._VegType != BARE)[0]
        n_plant = len(plant_cells)
        Theta = np.choose(
            self._VegType[plant_cells],
            [self._th_g, self._th_sh, self._th_tr, 0, self._th_sh_s, self._th_tr_s],
        )
        PMd = self._CumWS[plant_cells] - Theta
        PMd[PMd < 0.0] = 0.0
        tpmax = np.choose(
            self._VegType[plant_cells],
            [
                200000,
                self._tpmax_sh,
                self._tpmax_tr,
                0,
                self._tpmax_sh_s,
                self._tpmax_tr_s,
            ],
        )
        PMa = np.zeros(n_plant)
        tp_plant = self._tp[plant_cells]
        tp_greater = np.where(tp_plant > 0.5 * tpmax)[0]
        PMa[tp_greater] = (
            (tp_plant[tp_greater] - 0.5 * tpmax[tp_greater]) / (0.5 * tpmax[tp_greater])
        ) - 1
        PMb = np.choose(
            self._VegType[plant_cells],
            [
                self._Pmb_g,
                self._Pmb_sh,
                self._Pmb_tr,
                0,
                self._Pmb_sh_s,
                self._Pmb_tr_s,
            ],
        )
        PM = PMd + PMa + PMb
        PM[PM > 1.0] = 1.0
        R_Mor = np.random.rand(n_plant)  # Random number for comparison to kill
        Mortality = np.int32(np.where(np.greater_equal(PM, R_Mor))[0])
        self._VegType[plant_cells[Mortality]] = BARE
        self._tp[plant_cells[Mortality]] = 0

        self._cell_values["plant__age"] = self._tp

        if self._Edit_VegCov:
            self._grid["cell"]["vegetation__boolean_vegetated"] = np.zeros(
                self._grid.number_of_cells, dtype=int
            )
            self._grid["cell"]["vegetation__boolean_vegetated"][
                self._VegType != BARE
            ] = 1

        # For debugging purposes
        self._bare_cells = bare_cells
        self._Established = bare_cells[Establish]
        self._plant_cells = plant_cells
        self._Mortified = plant_cells[Mortality]


def count(Arr, value):
    Res = np.zeros(Arr.shape[0], dtype=int)
    x, y = Arr.shape
    for i in range(0, x):
        for j in range(0, y):
            if Arr[i][j] == value:
                Res[i] += 1
    return Res


def WS_PFT(VegType, PlantType, WS):
    Phi = np.zeros(WS.shape[0])
    x, y = WS.shape
    for i in range(0, x):
        for j in range(0, y):
            if VegType[i][j] == PlantType:
                Phi[i] += WS[i][j]
    return Phi
