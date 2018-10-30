
"""
21 May 2017
Authors: Sai Nudurupati & Erkan Istanbulluoglu

This is to replicate create_fires function from
Sujith Ravi's model as implemented in resource_redistribution_funcs.py.
In this fire creation method, we will also consider trees. Also, we
introduce fire suscessibility thresholds for shrubs and trees.

0: Bare; 1: Grass; 2: Shrub; 3: Burnt Grass; 4: Burnt Shrub; 5: Trees
6: Burnt Trees; 7: Shrub Seed; 8: Tree Seed
"""

#%% Import Packages
import numpy as np
from landlab import FieldError, Component
from ...utils.decorators import use_file_name_or_kwds
from funcs import (convert_phy_pft_to_distr_pft,
                   convert_distr_pft_to_phy_pft)
#from __future__ import print_function


_VALID_SCHEMES = set(['zhou_et_al_2013', 'ravi_et_al_2009'])


def _assert_pft_scheme_is_valid(scheme):
    if scheme not in _VALID_SCHEMES:
        raise ValueError('%s: Invalid PFT scheme' % scheme)


# %% Declare Global Variables (If any)
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

# %%


class SpatialDisturbance(Component):
    """
    What is this component (brief description)?

    A little more information about this component!

    What does this component output? Are there multiple
    processes? What are the key methods?

    Construction::
        SpatialDisturbance()

    Parameters:
    ----------
    grid: RasterModelGrid
        grid, Landlab's RasterModelGrid object
    V: array_like
        Vegetation Plant Functional Type; shape = [grid.number_of_cells]
        BARE = 0; GRASS = 1; SHRUB = 2; BURNTGRASS = 3; BURNTSHRUB = 4;
        TREE = 5; BURNTTREE = 6; SHRUBSEED = 7; TREESEED = 8
    n_fires: int, optional
        Number of fires to be created
    fire_area_mean: float, optional
        mean area of uniform distribution to sample fire size
    fire_area_dev: float, optional
        standard deviation of uniform distribution to sample fire size
    sh_susc: float, optional
        susceptibility of shrubs to burn
    tr_susc: float, optional
        susceptibility of trees to burn
    gr_susc: float, optional
        susceptibility of grass to burn
    sh_seed_susc: float, optional
        susceptibility of shrub seedlings to burn
    tr_seed_susc: float, optional
        susceptibility of tree seedlings to burn

    Examples
    --------
    - Add an example usage of this class.
    - Make sure that this example is run as if
    you are running commands from 'python command window'
    or 'ipython command window'.
    - The goal here is to familiarize user with different
    *input_fields*, and *output_fields*. More detailed
    the examples are, the better!
    Note: there might be issues with formatting this area.
    Deal with them!
    """

    _name = 'SpatialDisturbance'

    __version__ = 'landlab version in which this component is created in'

    _input_var_names = (
            'vegetation__plant_functional_type',
            )

    _output_var_names = (
            'vegetation__plant_functional_type',
            )

    _var_units = {
            'vegetation__plant_functional_type': 'None',
            }

    _var_mapping = {
            'vegetation__plant_functional_type': 'cell',
            }

    _var_doc = {
            'vegetation__plant_functional_type':
            'classification of plant type - zhou_et_al_2013 (int)' +
            'grass=0, shrub=1, tree=2, bare=3,' +
            'shrub_seedling=4, tree_seedling=5',
            }

    @use_file_name_or_kwds
    def __init__(self, grid, pft_scheme='zhou_et_al_2013',
                 **kwds):
        self._pft_scheme = pft_scheme
        _assert_pft_scheme_is_valid(self._pft_scheme)
        super(SpatialDisturbance, self).__init__(grid, **kwds)

        if self._pft_scheme == 'zhou_et_al_2013':
            if 'vegetation__plant_functional_type' not in self.grid.at_cell:
                raise FieldError("Cellular field of 'Plant Functional Type'" +
                                 " is required!")

    def graze(self, V=None, grazing_pressure=0.01):
        """
        Function to implement grazing
        """
        if self._pft_scheme == 'zhou_et_al_2013':
            vegtype = self._grid.at_cell['vegetation__plant_functional_type']
            V = convert_phy_pft_to_distr_pft(self._grid, vegtype)
        elif self._pft_scheme == 'ravi_et_al_2009':
            if V is None:
                raise ValueError("Cellular field of 'Plant Functional Type'" +
                                 " should be provided!")
        grz_prob = (0.6 * grazing_pressure +
                    2 * 0.4 * grazing_pressure * np.random.random_sample())
        grass_cells = np.where(V == 1)[0]
        compute_ = np.random.random(grass_cells.shape)
        grazed_cells = grass_cells[compute_ < grz_prob]
        V[grazed_cells] = 0
        if self._pft_scheme == 'zhou_et_al_2013':
            vegtype = convert_distr_pft_to_phy_pft(self._grid, V)
            self._grid.at_cell['vegetation__plant_functional_type'] = vegtype
        return (V, grazed_cells)

    def initiate_fires(self, V=None, n_fires=2, fire_area_mean=0.0625,
                     fire_area_dev=0.01, sh_susc=1., tr_susc=1.,
                     gr_susc=1., sh_seed_susc=1., tr_seed_susc=1.):
        """
        - Add description to this method. If this is the main method or
        one of the main methods, i.e. if this method performs a
        process for which this component is written, make sure you
        mention it.
        - Search for BasicModelInterface (BMI) on CSDMS website.
        We try to follow this interface to enable easier coupling
        with other models.

        Parameters:
        ----------
        grid: RasterModelGrid
            grid, Landlab's RasterModelGrid object
        V: array_like
            Vegetation Plant Functional Type; shape = [grid.number_of_cells]
            BARE = 0; GRASS = 1; SHRUB = 2; BURNTGRASS = 3; BURNTSHRUB = 4;
            TREE = 5; BURNTTREE = 6; SHRUBSEED = 7; TREESEED = 8
        n_fires: int, optional
            Number of fires to be created
        fire_area_mean: float, optional
            mean area of uniform distribution to sample fire size
        fire_area_dev: float, optional
            standard deviation of uniform distribution to sample fire size
        sh_susc: float, optional
            susceptibility of shrubs to burn
        tr_susc: float, optional
            susceptibility of trees to burn
        gr_susc: float, optional
            susceptibility of grass to burn
        sh_seed_susc: float, optional
            susceptibility of shrub seedlings to burn
        tr_seed_susc: float, optional
            susceptibility of tree seedlings to burn
        """
        if self._pft_scheme == 'zhou_et_al_2013':
            vegtype = self._grid.at_cell['vegetation__plant_functional_type']
            V = convert_phy_pft_to_distr_pft(self._grid, vegtype)
        elif self._pft_scheme == 'ravi_et_al_2009':
            if V is None:
                raise ValueError("Cellular field of 'Plant Functional Type'" +
                                 " should be provided!")
        susc = self._set_susceptibility(V, sh_susc=sh_susc, 
                                        tr_susc=tr_susc,
                                        gr_susc=gr_susc,
                                        sh_seed_susc=sh_seed_susc,
                                        tr_seed_susc=tr_seed_susc)
        ignition_cells = []
        burnt_locs = []  # Total burnt locations for all fires
        for i in range(0, n_fires):
            ignition_cell = np.random.choice(self._grid.number_of_cells, 1)
            if V[ignition_cell] == GRASS:
                (fire_locs, V) = self._spread_fire(
                                      V, ignition_cell,
                                      fire_area_mean=fire_area_mean,
                                      fire_area_dev=fire_area_dev,
                                      susc=susc)
            else:
                fire_locs = []
            burnt_locs += fire_locs
            ignition_cells += list(ignition_cell)

        if self._pft_scheme == 'zhou_et_al_2013':
            vegtype = convert_distr_pft_to_phy_pft(self._grid, V)
            self._grid.at_cell['vegetation__plant_functional_type'] = vegtype
        return (V, burnt_locs, ignition_cells)

    def _spread_fire(self, V, ignition_cell, fire_area_mean=0.0625,
                     fire_area_dev=0.01,
                     susc=None):
        """
        - An underscore in the front signals the user to stay away
        from using this method (is intended for internal use).
        - Works just like a regular method but implies a hidden method.
        - You should still document these methods though.

        Parameters:
        ----------
        grid: RasterModelGrid
            grid, Landlab's RasterModelGrid object
        V: array_like
            Vegetation Plant Functional Type; shape = [grid.number_of_cells]
            BARE = 0; GRASS = 1; SHRUB = 2; BURNTGRASS = 3; BURNTSHRUB = 4;
            TREE = 5; BURNTTREE = 6; SHRUBSEED = 7; TREESEED = 8
        ignition_cell: int
            cell id where the fire starts
        fire_area_mean: float, optional
            mean area of uniform distribution to sample fire size
        fire_area_dev: float, optional
            standard deviation of uniform distribution to sample fire size
        sh_susc: float, optional
            susceptibility of shrubs to burn
        tr_susc: float, optional
            susceptibility of trees to burn
        gr_susc: float, optional
            susceptibility of grass to burn
        sh_seed_susc: float, optional
            susceptibility of shrub seedlings to burn
        tr_seed_susc: float, optional
            susceptibility of tree seedlings to burn
        """
        if susc == None:
            susc = np.ones(self.grid.number_of_cells)
        fire_burnt = 0    # To check how many cells are being burnt
        grass_cells = np.where(V == GRASS)[0]
        if int(grass_cells.shape[0]) == 1:
            return [], V, []
        fire_locs = []       # record all the cell ids where fire has spread
        fire_locs += list(ignition_cell)
        burning_cells = [ignition_cell]
        V = self._burn_veg(V, burning_cells)
        fire_burnt += 1
        alr_cntd = []
        # loop to propagate fires one ring at a time
        while (burning_cells != []):
            newly_burnt = []   # Cells to be burnt in the sub-loop
            for cell in burning_cells:
                neigh_ = self._grid.looped_neighbors_at_cell[cell]
                veg_neighbors = (neigh_[np.where(V[neigh_] != BARE)])
                unique_neigh = np.setdiff1d(veg_neighbors, alr_cntd)
                alr_cntd += list(unique_neigh)
                susc_neigh = self._check_susc(unique_neigh,
                                              susc[unique_neigh])
                newly_burnt += (susc_neigh)
            if newly_burnt == []:
                break
            burning_cells = np.unique(np.array(newly_burnt))
            fire_locs += list(burning_cells)
            V = self._burn_veg(V, burning_cells)
            fire_burnt += int(burning_cells.shape[0])
            fire_area_sample = (self._fetch_uniform_random_fire_area(
                                        fire_area_mean, fire_area_dev))
            if fire_burnt > fire_area_sample*self._grid.number_of_cells:
                break
        return (fire_locs, V)

    def _fetch_uniform_random_fire_area(self, fire_area_mean, fire_area_dev):
        a = fire_area_mean - fire_area_dev
        return (a+2*fire_area_dev*np.random.random_sample())

    def _burn_veg(self, V, newly_burnt):
        newly_burnt = np.array(newly_burnt, dtype=int)
        burnt_grass = newly_burnt[np.where(V[newly_burnt] == GRASS)[0]]
        burnt_shrub = newly_burnt[np.where(V[newly_burnt] == SHRUB)[0]]
        burnt_tree = newly_burnt[np.where(V[newly_burnt] == TREE)[0]]
        burnt_shrub_seed = newly_burnt[np.where(
                                        V[newly_burnt] == SHRUBSEED)[0]]
        burnt_tree_seed = newly_burnt[np.where(V[newly_burnt] == TREESEED)[0]]
        V[burnt_grass] = BURNTGRASS
        V[burnt_shrub] = BURNTSHRUB
        V[burnt_tree] = BURNTTREE
        V[burnt_shrub_seed] = BURNTSHRUB
        V[burnt_tree_seed] = BURNTTREE
        return (V)

    def _check_susc(self, some_neighbors, susc):
        if some_neighbors.shape[0] == 0:
            susc_neighbors = []
        else:
            rand_val = np.random.rand(some_neighbors.shape[0])
            susc_neighbors = some_neighbors[rand_val < susc]
        return (list(susc_neighbors))

    def _set_susceptibility(self, V=None, sh_susc=1.,
                            tr_susc=1., gr_susc=1.,
                            sh_seed_susc=1., tr_seed_susc=1.):
        susc = np.zeros(self.grid.number_of_cells)
        susc[V==SHRUB] = sh_susc
        susc[V==TREE] = tr_susc
        susc[V==GRASS] = gr_susc
        susc[V==SHRUBSEED] = sh_seed_susc
        susc[V==TREESEED] = tr_seed_susc
        return susc
