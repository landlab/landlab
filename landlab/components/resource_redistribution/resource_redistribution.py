# -*- coding: utf-8 -*-
"""
06 Jul 2018
Authors: Sai Nudurupati & Erkan Istanbulluoglu

Resource Redistribution component for landlab. This component
has erode, deposit, establishment and mortality routines from
Sujith Ravi's 2009 Landscape Ecology paper.
0: Bare; 1: Grass; 2: Shrub; 3: Burnt Grass; 4: Burnt Shrub;
"""

#%% Import Packages
import numpy as np
from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
#from __future__ import print_function


#%% Declare Global Variables (If any)
# 'ravi_et_al_2009' pft_scheme - used internally
BARE = 0
GRASS = 1
SHRUB = 2
BURNTGRASS = 3
BURNTSHRUB = 4

#%%
class ResourceRedistribution(Component):
    """
    What is this component (brief description)?
    
    A little more information about this component!

    What does this component output? Are there multiple
    processes? What are the key methods?

    Construction::
        ResourceRedistribution()

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
            'soil__resources',
            )

    _output_var_names = (
            'vegetation__plant_functional_type',
            'soil__resources',
            )

    _var_units = {
            'vegetation__plant_functional_type': None,
            'soil__resources': None,
            }

    _var_mapping = {
            'vegetation__plant_functional_type': 'cell',
            'soil__resources': 'cell',
            }

    _var_doc = {
            'vegetation__plant_functional_type':
            'classification of plant type - bare=0, grass=1,' +
            'shrub=2, burntgrass=3, burntshrub=4',
            'soil__resources': 'level of soil resources',
            }

    @use_file_name_or_kwds
    def __init__(self, grid, e=0.1, R_low_threshold=-2.,
                 R_threshold=2., R_dep_threshold=1.,
                 Rth_gr=0.4, Rth_sh=0.8, P_gr_regrwth=0.25,
                 P_sh_regrwth=0.25, Pen=0.05, Pgrz=0.01,
                 P_gr=0.5, sh_max_age=600, sh_seedling_max_age=18,
                 sh_seedling_mor_dis=0., sh_mor_dis_low_thresh_age=300,
                 sh_mor_dis_low_slp=0.01, sh_mor_dis_high_slp=0.99,
                 P_sh_fire_mor=0.75, P_gr_fire_mor=1.,
                 sh_mor_ws_thresh=0.01, gr_mor_ws_thresh=0.08,
                 **kwds):
        super(ResourceRedistribution, self).__init__(grid, **kwds)

        name = 'vegetation__plant_functional_type'
        if name not in self.grid.at_cell:
            print('Since a cellular field of PFTs ' +
                  'is not provided, the field ' +
                  name + ' is initialized to 0 - all bare cells!')
            self.grid.add_zeros('cell', name, units=self._var_units[name])
        name = 'soil__resources'
        if name not in self.grid.at_cell:
            print('Since a cellular field for Resources ' +
                  'is not provided, the field ' +
                  name + ' is initialized to 1!')
            self.grid.add_ones('cell', name, units=self._var_units[name])

        self.initialize(e=e, R_low_threshold=R_low_threshold,
                        R_threshold=R_threshold,
                        R_dep_threshold=R_dep_threshold, Rth_gr=Rth_gr,
                        Rth_sh=Rth_sh, P_gr_regrwth=P_gr_regrwth,
                        P_sh_regrwth=P_sh_regrwth, Pen=Pen, Pgrz=Pgrz,
                        P_gr=P_gr, sh_max_age=sh_max_age,
                        sh_seedling_max_age=sh_seedling_max_age,
                        sh_seedling_mor_dis=sh_seedling_mor_dis,
                        sh_mor_dis_low_thresh_age=sh_mor_dis_low_thresh_age,
                        sh_mor_dis_low_slp=sh_mor_dis_low_slp,
                        sh_mor_dis_high_slp=sh_mor_dis_high_slp,
                        P_sh_fire_mor=P_sh_fire_mor,
                        P_gr_fire_mor=P_gr_fire_mor,
                        sh_mor_ws_thresh=sh_mor_ws_thresh,
                        gr_mor_ws_thresh=gr_mor_ws_thresh,
                        **kwds)

    def initialize(self, e=0.1, R_low_threshold=-2.,
                   R_threshold=2., R_dep_threshold=1.,
                   Rth_gr=0.4, Rth_sh=0.8, P_gr_regrwth=0.25,
                   P_sh_regrwth=0.25, Pen=0.05, Pgrz=0.01,
                   P_gr=0.5, sh_max_age=600, sh_seedling_max_age=18,
                   sh_seedling_mor_dis=0., sh_mor_dis_low_thresh_age=300,
                   sh_mor_dis_low_slp=0.01, sh_mor_dis_high_slp=0.99,
                   P_sh_fire_mor=0.75, P_gr_fire_mor=1.,
                   sh_mor_ws_thresh=0.01, gr_mor_ws_thresh=0.08,
                   **kwds):
        self._e = e
        self._R_low_threshold = R_low_threshold
        self._R_threshold = R_threshold
        self._R_dep_threshold = R_dep_threshold
        self._Rth_gr = Rth_gr
        self._Rth_sh = Rth_sh
        self._P_gr_regrwth = P_gr_regrwth
        self._P_sh_regrwth = P_sh_regrwth
        self._Pen = Pen
        self._Pgrz = Pgrz
        self._P_gr = P_gr
        self._sh_max_age = sh_max_age
        self._sh_seedling_max_age = sh_seedling_max_age
        self._sh_seedling_mor_dis = sh_seedling_mor_dis
        self._sh_mor_dis_low_thresh_age = sh_mor_dis_low_thresh_age
        self._sh_mor_dis_low_slp = sh_mor_dis_low_slp
        self._sh_mor_dis_high_slp = sh_mor_dis_high_slp
        self._P_sh_fire_mor = P_sh_fire_mor
        self._P_gr_fire_mor = P_gr_fire_mor
        self._sh_mor_ws_thresh = sh_mor_ws_thresh
        self._gr_mor_ws_thresh = gr_mor_ws_thresh

    def erode(self):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        R = self.grid.at_cell['soil__resources']
        Elig_R = np.where(R > self._R_low_threshold)[0]
        # Deal with shrubs first
        burnt_shrubs = Elig_R[V[Elig_R] == BURNTSHRUB]
        eroded_soil_shrub = (np.int(burnt_shrubs.shape[0]) * 4 * self._e)
        R[burnt_shrubs] -= (4 * self._e)
        # Deal with erosion of other cells
        burnt_grass = Elig_R[V[Elig_R]==BURNTGRASS]
        bare_cells = Elig_R[V[Elig_R]==BARE]
        R[burnt_grass] -= (0.2 * self._e)
        R[bare_cells] -= self._e
        eroded_soil = (np.int(burnt_grass.shape[0]) * 0.2 * self._e +
                       np.int(bare_cells.shape[0]) * self._e)
        return (eroded_soil, eroded_soil_shrub, burnt_shrubs,
                burnt_grass, bare_cells)

    def deposit(self, eroded_soil, eroded_soil_shrub):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        R = self.grid.at_cell['soil__resources']
        exclusive = np.arange(0, self.grid.number_of_cells)
        burnt_shrubs = np.where(V == BURNTSHRUB)[0]
        if int(burnt_shrubs.shape[0]) > 0:
            burnt_shrubs_neigh = np.unique(
                    self.grid.looped_neighbors_at_cell[burnt_shrubs])
            R[burnt_shrubs_neigh] += (
                    eroded_soil_shrub / float(burnt_shrubs_neigh.shape[0]))
            exclusive = (
                    exclusive[np.in1d(exclusive, burnt_shrubs, invert=True)])
            exclusive = (
                exclusive[np.in1d(exclusive, burnt_shrubs_neigh, invert=True)])
        else:
            burnt_shrubs_neigh = []
        shrub_exclusive = exclusive[np.where(V[exclusive] == SHRUB)[0]]
        grass_exclusive = exclusive[np.where(V[exclusive] == GRASS)[0]]
        bare_exclusive = exclusive[np.where(V[exclusive] == BARE)[0]]
        ## Calculating how much each bare cell will get from eroded_soil 'Eb'
        ## Shrubs will receive thrice Es = 3*Eb && Eg = 2*Eb
        weighted_parts = (1 * int(bare_exclusive.shape[0]) +
                          3 * int(shrub_exclusive.shape[0]) +
                          2 * int(grass_exclusive.shape[0]))
        eroded_soil_part = (eroded_soil / float(weighted_parts))
        R[bare_exclusive] += eroded_soil_part
        R[shrub_exclusive] += (3. * eroded_soil_part)
        R[grass_exclusive] += (2. * eroded_soil_part)
        return (burnt_shrubs_neigh, exclusive, shrub_exclusive,
                grass_exclusive, bare_exclusive, eroded_soil_part)

    def re_adjust_resource(self):
        R = self.grid.at_cell['soil__resources']
    ## Resource exceeding R_threshold will be distributed to its neighbors                                
        resource_adjusted = 0.
        eligible_locs_to_adj_neigh = np.array([])
        locs_to_adj = np.where(R > self._R_threshold)[0]
        if int(locs_to_adj.shape[0]) > 0:        
            resource_adjusted = np.sum(R[locs_to_adj] - self._R_threshold)
            locs_to_adj_neigh = np.unique(
                    self.grid.looped_neighbors_at_cell[locs_to_adj])
            eligible_locs_to_adj_neigh = (
                    locs_to_adj_neigh[R[locs_to_adj_neigh] < self._R_dep_threshold])
            if int(eligible_locs_to_adj_neigh.shape[0]) > 0:
                R[eligible_locs_to_adj_neigh] += (
                    resource_adjusted / float(
                            eligible_locs_to_adj_neigh.shape[0])) 
                R[locs_to_adj] = self._R_threshold
    ## Resource below R_low_threshold is raised to R_low_threshold 
        sed_to_borrow = 0.
        Elig_locs = np.where(R < self._R_low_threshold)[0]
        if int(Elig_locs.shape[0]) > 0:
            sed_to_borrow = np.absolute(
                    np.sum(self._R_low_threshold - R[Elig_locs]))
            R[Elig_locs] = self._R_low_threshold        
            locs_to_borrow = np.where(R > 0.)[0]
            R[locs_to_borrow] -= (
                    sed_to_borrow / float(locs_to_borrow.shape[0]))    
        return (resource_adjusted, eligible_locs_to_adj_neigh,
                Elig_locs, sed_to_borrow)

    def establish(self, V_age):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        R = self.grid.at_cell['soil__resources']
        burnt_shrubs = np.where(V == BURNTSHRUB)[0]
        burnt_grass = np.where(V == BURNTGRASS)[0]
        ## Regrowth in burnt area
        shrubs_r_regrwth = burnt_shrubs[np.where(R[burnt_shrubs] >  self._Rth_sh)[0]]
        grass_r_regrwth = burnt_grass[np.where(R[burnt_grass] >  self._Rth_gr)[0]]
        P_check_1 = np.random.random(shrubs_r_regrwth.shape)
        est_1 = shrubs_r_regrwth[np.where(P_check_1 <  self._P_sh_regrwth)[0]]
        V[est_1] = SHRUB
        V_age[est_1] = 0
        P_check_2 = np.random.random(grass_r_regrwth.shape)
        est_2 = grass_r_regrwth[np.where(P_check_2 <  self._P_gr_regrwth)[0]]
        V[est_2] = GRASS
        V_age[est_2] = 0
        ## Regrowth in grazed area
        # shrub encroachment due to neighbors
        bare_cells = np.where(V == BARE)[0]
        shrub_grz_regrwth = bare_cells[np.where(R[bare_cells] >  self._Rth_sh)[0]]
        neigh_sh_grz_regrwth = (
            self.grid.looped_neighbors_at_cell[shrub_grz_regrwth])
        ns_sh = self._np_ndarray_count(neigh_sh_grz_regrwth)
        ns_Pgrz = (ns_sh *  self._Pen)      # ns * P2
        P_check_3 = np.random.random(ns_Pgrz.shape)
        est_3 = shrub_grz_regrwth[np.where(P_check_3 < ns_Pgrz)[0]]
        V[est_3] = SHRUB    # Establish shrubs
        V_age[est_3] = 0
        # shrub encroachment due to grazing
        bare_cells_ = np.where(V == BARE)[0]
        shrub_grz_regrwth_ = bare_cells_[np.where(R[bare_cells_] >  self._Rth_sh)[0]]
        P_check_4 = np.random.random(shrub_grz_regrwth_.shape)
        est_4 = shrub_grz_regrwth_[np.where(P_check_4 <  self._Pgrz)[0]]
        V[est_4] = SHRUB    # Establish grass
        V_age[est_4] = 0
        # grass growth where shrubs haven't encroached - grass seed dispersal
        bare_cells_2 = np.where(V == BARE)[0]
        grass_grz_regrwth = bare_cells_2[np.where(R[bare_cells_2] >  self._Rth_gr)[0]]
        P_check_5 = np.random.random(grass_grz_regrwth.shape)
        est_5 = grass_grz_regrwth[np.where(P_check_5 <  self._P_gr)[0]]
        V[est_5] = GRASS    # Establish grass
        V_age[est_5] = 0
        return (V_age, est_1, est_2, est_3, est_4, est_5)

    def _np_ndarray_count(self, neigh_sh_grz_regrwth):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        ns_sh = np.zeros(neigh_sh_grz_regrwth.shape[0], dtype=int)
        for i in range(0, ns_sh.shape[0]):
            ns_sh[i] = (
                np.int(
                    np.where(V[neigh_sh_grz_regrwth[i,:]] == 2)[0].shape[0]))
        return (ns_sh)

    def _compute_Prob_mortality_age(self, V_age, Pmor_age):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        # Age mortality for Vegetation other than shrub = 0.
        # Shrub seedling
        shrubs_ = np.where(V == SHRUB)[0]
        shrub_seedlings = (
                    shrubs_[np.where(V_age[shrubs_] < self._sh_seedling_max_age)[0]])
        Pmor_age[shrub_seedlings] = self._sh_seedling_mor_dis
        young_shrubs = (
                shrubs_[np.logical_and(V_age[shrubs_] > self._sh_seedling_max_age,
                        V_age[shrubs_] < self._sh_mor_dis_low_thresh_age)])
        Pmor_age[young_shrubs] = self._sh_mor_dis_low_slp  # This is constant
        adult_shrubs = (
                shrubs_[np.where(V_age[shrubs_] > self._sh_mor_dis_low_thresh_age)])
        Pmor_age[adult_shrubs] = (
                 self._sh_mor_dis_low_slp +
                 ((np.asfarray(V_age[adult_shrubs]) - self._sh_mor_dis_low_thresh_age) /
                 (self._sh_max_age - self._sh_mor_dis_low_thresh_age) *
                 (self._sh_mor_dis_high_slp)))
        return (Pmor_age)

    def mortality(self, V_age):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        # Killing burnt vegetation    
        burnt_shrubs = np.where(V == BURNTSHRUB)[0]
        P_check_1 = np.random.random(burnt_shrubs.shape)
        kill_1 = burnt_shrubs[np.where(P_check_1 < self._P_sh_fire_mor)[0]]
        V[kill_1] = BARE
        V_age[kill_1] = 0
        burnt_grass = np.where(V == BURNTGRASS)[0]
        P_check_2 = np.random.random(burnt_grass.shape)
        kill_2 = burnt_grass[np.where(P_check_2 < self._P_gr_fire_mor)[0]]
        V[kill_2] = BARE
        V_age[kill_2] = 0
        
        # Mortality due to Water Stress and age
        Pmor_age_ws = np.zeros(V.shape)
              # Creating a field for Prob_mortality for Water Stress & Age
        Pmor_age = np.zeros(V.shape)
              # Creating a field for Prob_mortality due to age alone 
        Pmor_age_ws[V == GRASS] += (0.6 * self._gr_mor_ws_thresh +
                                    2 * 0.4 * self._gr_mor_ws_thresh *
                                    np.random.random_sample())
              # Grass' water stess probability
        Pmor_age_ws[V == SHRUB] += (0.6 * self._sh_mor_ws_thresh +
                                    2 * 0.4 * self._sh_mor_ws_thresh *
                                    np.random.random_sample())
              # Shrubs' water stress probability
        Pmor_age = self._compute_Prob_mortality_age(V_age, Pmor_age)
        Pmor_age_ws += Pmor_age
        P_check_3 = np.random.random(V.shape)
        kill_3 = np.where(P_check_3 < Pmor_age_ws)[0]
        V[kill_3] = BARE
        V_age[kill_3] = 0
        return (V_age, Pmor_age, Pmor_age_ws)

    def initialize_Veg_age(self, V_age):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        V_age[V == SHRUB] = np.random.randint(0, self._sh_max_age,
                                              V_age[V == SHRUB].shape)
        return (V_age)
        
    def update_Veg_age(self, V_age):
        V = self.grid.at_cell['vegetation__plant_functional_type']
        V_age[V != BARE] += 1
        return (V_age)
