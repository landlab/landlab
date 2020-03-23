# Import Packages
import numpy as np
from landlab import FieldError, Component
from ...utils.decorators import use_file_name_or_kwds
from .funcs import convert_phy_pft_to_distr_pft, convert_distr_pft_to_phy_pft

_VALID_SCHEMES = set(["ravi_et_al_2009", "zhou_et_al_2013"])


def _assert_pft_scheme_is_valid(scheme):
    if scheme not in _VALID_SCHEMES:
        raise ValueError("%s: Invalid PFT scheme" % scheme)


# Declare Global Variables (If any)
# 'ravi_et_al_2009' pft_scheme - used internally
BARE = 0
GRASS = 1
SHRUB = 2
BURNTGRASS = 3
BURNTSHRUB = 4
TREE = 5
BURNTTREE = 6
SHRUBSEED = 7
TREESEED = 8


class SpatialDisturbance(Component):
    """
    Landlab component that implements that implements
    spatial disturbances, such as wildfires and
    grazing, conceptualized by Ravi and D'Odorico (2009).
    These disturbances modify the vegetation occupying
    the cells (e.g. convert SHRUB to BURNTSHRUB) depending on
    their spatial reach.

    This component is designed to work with two different
    plant functional type (PFT) schemes: 'ravi_et_al_2009',
    and 'zhou_et_al_2013'.
    'ravi_et_al_2009': Each cell in the RasterModelGrid object
    can take an integer from 0 through 8 to
    represent the cell states for bare soil [0],
    grass [1], shrub [2], burnt grass [3], burnt shrub [4],
    tree [5], burnt tree [6], shrub seed [7], and tree seed [8].
    This scheme is an extended version of the PFTs used in
    Ravi and D'Odorico (2009).
    'zhou_et_al_2013' (default): Each cell in the
    RasterModelGrid object can take an integer
    from 0 through 5 to represent the cell states
    for grass [0], shrub [1], tree [2], bare soil [3],
    shrub seedling [4], and tree seedling [5].
    To use the 'zhou_et_al_2013' scheme,
    select this scheme while instantiating the
    component and pass the input vegetation field
    "vegetation__plant_functional_type". To use
    the other scheme, pass the vegetation
    field as input to the methods (i.e.,
    SpatialDisturbance.initiate_fires(V=veg_field)).
    NOTE: The methods of this component internally
    convert the input PFTs in 'zhou_et_al_2013' scheme
    to 'ravi_et_al_2009', implement the algorithm,
    and convert the PFTs back to 'zhou_et_al_2013'
    when the component is instantiated with
    'zhou_et_al_2013' scheme.

    There are two key process-representing methods in this
    component: graze(), and initiate_fires().

    References:
    Ravi, S., & D’Odorico, P. (2009). Post-fire
    resource redistribution and fertility island dynamics
    in shrub encroached desert grasslands: a modeling approach.
    Landscape Ecology, 24(3), 325-335.
    Zhou, X., Istanbulluoglu, E., & Vivoni, E. R. (2013). Modeling the
    ecohydrological role of aspect controlled radiation on tree grass shrub
    coexistence in a semiarid climate. Water Resources Research,
    49(5), 2872-2895.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(0)
    >>> from landlab import RasterModelGrid as rmg
    >>> from landlab.components import SpatialDisturbance
    >>> grid = rmg((10, 10), xy_spacing=(0.2, 0.2))
    >>> SpatialDisturbance.name
    'Spatial Disturbance'
    >>> sorted(SpatialDisturbance.output_var_names)
    ['vegetation__plant_functional_type']
    >>> sorted(SpatialDisturbance.units)
    [('vegetation__plant_functional_type', 'None')]

    Let us look at an example using the default
    'zhou_et_al_2013' PFT scheme.

    >>> grid.at_cell["vegetation__plant_functional_type"] = (
    ...     np.random.randint(0, 4, size=grid.number_of_cells))
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"] == 0)[0].shape
    (15,)
    >>> sd = SpatialDisturbance(grid)
    >>> sd._grid.number_of_cell_rows
    8
    >>> sd._grid.number_of_cell_columns
    8
    >>> sd._grid is grid
    True
    >>> (V, grazed_cells) = sd.graze(grazing_pressure=0.5)
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==0)[0].shape
    (9,)
    >>> grid.at_cell["vegetation__plant_functional_type"] = (
    ...     np.random.randint(0, 3, size=grid.number_of_cells))
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==0)[0].shape
    (21,)
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==1)[0].shape
    (18,)
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==2)[0].shape
    (25,)
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==3)[0].shape
    (0,)
    >>> (V, burnt_locs, ignition_cells) = sd.initiate_fires(
    ...        n_fires=10,
    ...        fire_area_mean=0.0625,
    ...        sh_susc=0.8,
    ...        gr_susc=1.,
    ...        tr_susc=0.,
    ... )
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==0)[0].shape
    (11,)
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==1)[0].shape
    (15,)
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==2)[0].shape
    (25,)
    >>> np.where(
    ...     grid.at_cell["vegetation__plant_functional_type"]==3)[0].shape
    (13,)

    Let us look at an example using the
    'ravi_et_al_2009' PFT scheme.

    >>> new_grid = rmg((10, 10), xy_spacing=(0.2, 0.2))
    >>> V = (np.random.randint(1, 3, size=new_grid.number_of_cells))
    >>> np.where(V == GRASS)[0].shape
    (25,)
    >>> np.where(V == BARE)[0].shape
    (0,)
    >>> new_sd = SpatialDisturbance(new_grid, pft_scheme="ravi_et_al_2009")
    >>> (V, grazed_cells) = new_sd.graze(V=V, grazing_pressure=0.5)
    >>> np.where(V == GRASS)[0].shape
    (16,)
    >>> np.where(V == BARE)[0].shape
    (9,)
    >>> V = (np.random.randint(1, 3, size=new_grid.number_of_cells))
    >>> np.where(V == GRASS)[0].shape
    (27,)
    >>> np.where(V == SHRUB)[0].shape
    (37,)
    >>> np.where(V == BARE)[0].shape
    (0,)
    >>> np.where(V == BURNTGRASS)[0].shape
    (0,)
    >>> np.where(V == BURNTSHRUB)[0].shape
    (0,)
    >>> (V, burnt_locs, ignition_cells) = new_sd.initiate_fires(
    ...        V=V,
    ...        n_fires=10,
    ...        fire_area_mean=0.0625,
    ...        sh_susc=0.8,
    ...        gr_susc=1.,
    ... )
    >>> np.where(V == GRASS)[0].shape
    (17,)
    >>> np.where(V == SHRUB)[0].shape
    (16,)
    >>> np.where(V == BARE)[0].shape
    (0,)
    >>> np.where(V == BURNTGRASS)[0].shape
    (10,)
    >>> np.where(V == BURNTSHRUB)[0].shape
    (21,)

    Note that the initiate_fires method converted
    grass and shrubs, to burnt_grass and burnt_shrubs
    respectively, when we used 'ravi_et_al_2009' PFT scheme
    but grass and shrubs were converted to bare cells
    when we used 'zhou_et_al_2013' PFT scheme.
    """

    _name = "Spatial Disturbance"

    _input_var_names = ("vegetation__plant_functional_type",)

    _output_var_names = ("vegetation__plant_functional_type",)

    _var_units = {
        "vegetation__plant_functional_type": "None",
    }

    _var_mapping = {
        "vegetation__plant_functional_type": "cell",
    }

    _var_doc = {
        "vegetation__plant_functional_type": "classification of plant type - zhou_et_al_2013 (int)"
        + "grass=0, shrub=1, tree=2, bare=3,"
        + "shrub_seedling=4, tree_seedling=5",
    }

    @use_file_name_or_kwds
    def __init__(self, grid, pft_scheme="zhou_et_al_2013", **kwds):
        """
        Parameters:
        ----------
        grid: RasterModelGrid
            grid, Landlab's RasterModelGrid object
        pft_scheme: str
            Vegetation Plant Functional Type (PFT) scheme;
            Either 'ravi_et_al_2009' or 'zhou_et_al_2013'.
        """
        self._pft_scheme = pft_scheme
        _assert_pft_scheme_is_valid(self._pft_scheme)
        self._grid = grid

        if self._pft_scheme == "zhou_et_al_2013":
            if "vegetation__plant_functional_type" not in self._grid.at_cell:
                raise FieldError(
                    "Cellular field of 'Plant Functional Type'" + " is required!"
                )

    def graze(self, V=None, grazing_pressure=0.01):
        """
        removes grass from  a cell where grazing_pressure at
        grass-occupied cell is greater than a
        uniformly distributed random number U[0, 1]
        generated at each grass cell.

        Parameters:
        ----------
        V: numpy array, shape = [grid.number_of_cells],
           compulsory if pft_scheme = 'ravi_et_al_2009'.
            Vegetation Plant Functional Type;
            BARE = 0; GRASS = 1; SHRUB = 2; BURNTGRASS = 3; BURNTSHRUB = 4;
            TREE = 5; BURNTTREE = 6; SHRUBSEED = 7; TREESEED = 8
        grazing_pressure: float , optional
            Probability of grazing at GRASS cells [0, 1].

        Returns:
        -------
        A tuple: (V, grazed_cells);
        V: numpy array, shape = [grid.number_of_cells],
           Useful only if pft_scheme = 'ravi_et_al_2009'.
            Vegetation Plant Functional Type;
            BARE = 0; GRASS = 1; SHRUB = 2; BURNTGRASS = 3; BURNTSHRUB = 4;
            TREE = 5; BURNTTREE = 6; SHRUBSEED = 7; TREESEED = 8
        grazed_cells: numpy array of ints, (-)
            cell_ids of cells where GRASS cells have
            been converted to BARE due to
            grazing. (purpose - debugging)
        """
        if self._pft_scheme == "zhou_et_al_2013":
            vegtype = self._grid.at_cell["vegetation__plant_functional_type"]
            V = convert_phy_pft_to_distr_pft(self._grid, vegtype)
        elif self._pft_scheme == "ravi_et_al_2009":
            if V is None:
                raise ValueError(
                    "Cellular field of 'Plant Functional Type'" + " should be provided!"
                )
        grz_prob = (
            0.6 * grazing_pressure
            + 2 * 0.4 * grazing_pressure * np.random.random_sample()
        )
        grass_cells = np.where(V == 1)[0]
        compute_ = np.random.random(grass_cells.shape)
        grazed_cells = grass_cells[compute_ < grz_prob]
        V[grazed_cells] = 0
        if self._pft_scheme == "zhou_et_al_2013":
            vegtype = convert_distr_pft_to_phy_pft(self._grid, V)
            self._grid.at_cell["vegetation__plant_functional_type"] = vegtype
        return (V, grazed_cells)

    def initiate_fires(
        self,
        V=None,
        n_fires=2,
        fire_area_mean=0.0625,
        fire_area_dev=0.01,
        sh_susc=1.0,
        tr_susc=1.0,
        gr_susc=1.0,
        sh_seed_susc=1.0,
        tr_seed_susc=1.0,
        gr_vuln=1.0,
        sh_vuln=0.0,
        sh_seed_vuln=0.0,
        tr_vuln=0.0,
        tr_seed_vuln=0.0,
    ):
        """
        This method simulates a desired number of lightning strikes
        by randomly placing them across the grid. Fire starts depenging
        on the vegetation occupying the cell struck by lightning and
        it's vulnerability. For example, fire starts with a probability
        of 90% if lightning strikes a cell occupied by SHRUB
        when sh_vuln is 0.9. BARE cells cannot start a fire.
        Fire spreads to the neighbors based on 'susceptibility'
        of the vegetated neighbor. For example, fire spreads
        with a probability of 75% to a TREE neighbor with
        tr_susc of 0.75. Fire does not spread to a BARE cell.
        The size of the fire is limited
        by fire_area_mean (i.e., the fire stops
        spreading if the burnt area exceeds fire_area_mean).
        This method updaes the field "vegetation__plant_functional_type".

        Parameters:
        ----------
        grid: RasterModelGrid
            grid, Landlab's RasterModelGrid object
        V: numpy array, shape = [grid.number_of_cells],
           compulsory if pft_scheme = 'ravi_et_al_2009'.
            Vegetation Plant Functional Type;
            BARE = 0; GRASS = 1; SHRUB = 2; BURNTGRASS = 3; BURNTSHRUB = 4;
            TREE = 5; BURNTTREE = 6; SHRUBSEED = 7; TREESEED = 8
        n_fires: int, optional
            Number of lightning strikes to be created
        fire_area_mean: float, optional
            mean area of uniform distribution to sample fire size
        fire_area_dev: float, optional
            standard deviation of uniform distribution to sample fire size
        sh_susc: float, optional
            susceptibility of SHRUB to fire
        tr_susc: float, optional
            susceptibility of TREE to fire
        gr_susc: float, optional
            susceptibility of GRASS to fire
        sh_seed_susc: float, optional
            susceptibility of SHRUBSEED to fire
        tr_seed_susc: float, optional
            susceptibility of TREESEED to fire
        gr_vuln: float, optional
            probability of GRASS cell to catch fire due to
            lightning
        sh_vuln: float, optional
            probability of SHRUB cell to catch fire due to
            lightning
        sh_seed_vuln: float, optional
            probability of SHRUBSEED cell to catch fire due to
            lightning
        tr_vuln: float, optional
            probability of TREE cell to catch fire due to
            lightning
        tr_seed_vuln: float, optional
            probability of TREESEED cell to catch fire due to
            lightning

        Returns:
        -------
        A tuple: (V, burnt_locs, ignition_cells);
        V: numpy array, shape = [grid.number_of_cells],
           Useful only if pft_scheme = 'ravi_et_al_2009'.
            Vegetation Plant Functional Type;
            BARE = 0; GRASS = 1; SHRUB = 2; BURNTGRASS = 3; BURNTSHRUB = 4;
            TREE = 5; BURNTTREE = 6; SHRUBSEED = 7; TREESEED = 8
        burnt_locs: numpy array of ints, (-)
            cell_ids of cells where fire has spread. (purpose - debugging)
        ignition_cells: numpy array of ints, (-)
            cell_ids of cells where lightning struck. (purpose - debugging)
        """
        if self._pft_scheme == "zhou_et_al_2013":
            vegtype = self._grid.at_cell["vegetation__plant_functional_type"]
            V = convert_phy_pft_to_distr_pft(self._grid, vegtype)
        elif self._pft_scheme == "ravi_et_al_2009":
            if V is None:
                raise ValueError(
                    "Cellular field of 'Plant Functional Type'" + " should be provided!"
                )
        susc = self._set_susceptibility(
            V,
            sh_susc=sh_susc,
            tr_susc=tr_susc,
            gr_susc=gr_susc,
            sh_seed_susc=sh_seed_susc,
            tr_seed_susc=tr_seed_susc,
        )
        ignition_cells = []
        burnt_locs = []  # Total burnt locations for all fires
        for i in range(0, n_fires):
            ignition_cell = np.random.choice(self._grid.number_of_cells, 1)
            if self._is_cell_ignitable(
                V,
                ignition_cell,
                gr_vuln=gr_vuln,
                sh_vuln=sh_vuln,
                sh_seed_vuln=sh_seed_vuln,
                tr_vuln=tr_vuln,
                tr_seed_vuln=tr_seed_vuln,
            ):
                (fire_locs, V) = self._spread_fire(
                    V,
                    ignition_cell,
                    fire_area_mean=fire_area_mean,
                    fire_area_dev=fire_area_dev,
                    susc=susc,
                )
            else:
                fire_locs = []
            burnt_locs += fire_locs
            ignition_cells += list(ignition_cell)

        if self._pft_scheme == "zhou_et_al_2013":
            vegtype = convert_distr_pft_to_phy_pft(self._grid, V)
            self._grid.at_cell["vegetation__plant_functional_type"] = vegtype
        return (V, burnt_locs, ignition_cells)

    def _spread_fire(
        self, V, ignition_cell, fire_area_mean=0.0625, fire_area_dev=0.01, susc=None
    ):
        """
        This private method implements the fire algorithm.
        """
        if susc is None:
            susc = np.ones(self._grid.number_of_cells)
        fire_burnt = 0  # To check how many cells are being burnt
        grass_cells = np.where(V == GRASS)[0]
        if int(grass_cells.shape[0]) == 1:
            return [], V, []
        fire_locs = []  # record all the cell ids where fire has spread
        fire_locs += list(ignition_cell)
        burning_cells = [ignition_cell]
        V = self._burn_veg(V, burning_cells)
        fire_burnt += 1
        alr_cntd = []
        # loop to propagate fires one ring at a time
        while burning_cells != []:
            newly_burnt = []  # Cells to be burnt in the sub-loop
            for cell in burning_cells:
                neigh_ = self._grid.looped_neighbors_at_cell[cell]
                veg_neighbors = neigh_[np.where(V[neigh_] != BARE)]
                unique_neigh = np.setdiff1d(veg_neighbors, alr_cntd)
                alr_cntd += list(unique_neigh)
                susc_neigh = self._check_susc(unique_neigh, susc[unique_neigh])
                newly_burnt += susc_neigh
            if newly_burnt == []:
                break
            burning_cells = np.unique(np.array(newly_burnt))
            fire_locs += list(burning_cells)
            V = self._burn_veg(V, burning_cells)
            fire_burnt += int(burning_cells.shape[0])
            fire_area_sample = self._fetch_uniform_random_fire_area(
                fire_area_mean, fire_area_dev
            )
            if fire_burnt > fire_area_sample * self._grid.number_of_cells:
                break
        return (fire_locs, V)

    def _fetch_uniform_random_fire_area(self, fire_area_mean, fire_area_dev):
        a = fire_area_mean - fire_area_dev
        return a + 2 * fire_area_dev * np.random.random_sample()

    def _burn_veg(self, V, newly_burnt):
        newly_burnt = np.array(newly_burnt, dtype=int)
        burnt_grass = newly_burnt[np.where(V[newly_burnt] == GRASS)[0]]
        burnt_shrub = newly_burnt[np.where(V[newly_burnt] == SHRUB)[0]]
        burnt_tree = newly_burnt[np.where(V[newly_burnt] == TREE)[0]]
        burnt_shrub_seed = newly_burnt[np.where(V[newly_burnt] == SHRUBSEED)[0]]
        burnt_tree_seed = newly_burnt[np.where(V[newly_burnt] == TREESEED)[0]]
        V[burnt_grass] = BURNTGRASS
        V[burnt_shrub] = BURNTSHRUB
        V[burnt_tree] = BURNTTREE
        V[burnt_shrub_seed] = BURNTSHRUB
        V[burnt_tree_seed] = BURNTTREE
        return V

    def _check_susc(self, some_neighbors, susc):
        if some_neighbors.shape[0] == 0:
            susc_neighbors = []
        else:
            rand_val = np.random.rand(some_neighbors.shape[0])
            susc_neighbors = some_neighbors[rand_val < susc]
        return list(susc_neighbors)

    def _set_susceptibility(
        self,
        V=None,
        sh_susc=1.0,
        tr_susc=1.0,
        gr_susc=1.0,
        sh_seed_susc=1.0,
        tr_seed_susc=1.0,
    ):
        susc = np.zeros(self._grid.number_of_cells)
        susc[V == SHRUB] = sh_susc
        susc[V == TREE] = tr_susc
        susc[V == GRASS] = gr_susc
        susc[V == SHRUBSEED] = sh_seed_susc
        susc[V == TREESEED] = tr_seed_susc
        return susc

    def _is_cell_ignitable(
        self,
        V,
        ignition_cell,
        gr_vuln=1.0,
        sh_vuln=0.0,
        sh_seed_vuln=0.0,
        tr_vuln=0.0,
        tr_seed_vuln=0.0,
    ):
        vulnerability = np.choose(
            V[ignition_cell],
            [0.0, gr_vuln, sh_vuln, 0.0, 0.0, tr_vuln, 0.0, sh_seed_vuln, tr_seed_vuln],
        )
        rand_val = np.random.rand()
        if rand_val < vulnerability:
            return True
        else:
            return False
