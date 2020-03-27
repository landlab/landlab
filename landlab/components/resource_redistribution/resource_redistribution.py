# -*- coding: utf-8 -*-

import numpy as np
from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import warnings


# Declare Global Variables (If any)
# 'ravi_et_al_2009' pft_scheme - used internally
BARE = 0
GRASS = 1
SHRUB = 2
BURNTGRASS = 3
BURNTSHRUB = 4


class ResourceRedistribution(Component):
    """
    Landlab component that implements a two-state cellular automata
    model, that couples categorical vegetation (V) states and soil
    resource level (R), conceptualized by Ravi and D'Odorico (2009).

    Each cell in the RasterModelGrid object can take an integer from
    0 through 4 to represent the cell states for bare soil [0],
    grass [1], shrub [2], burnt grass [3], and burnt shrub [4].
    R is conceptualized to represent a soil resource value relative
    to the mean resource level of a productiive ecosystem (default:
    -2 to 2).

    There are four key process-representing methods in this
    component: establish(), mortality(), erode(), and deposit().
    re_adjust_resource() is a helpful utility function that helps
    limit R within the thresholds.

    Ref: Ravi, S., & D’Odorico, P. (2009). Post-fire
    resource redistribution and fertility island dynamics
    in shrub encroached desert grasslands: a modeling approach.
    Landscape Ecology, 24(3), 325-335.

    .. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(0)
    >>> from landlab import RasterModelGrid as rmg
    >>> from landlab.components import ResourceRedistribution
    >>> grid = rmg((5, 4), xy_spacing=(0.2, 0.2))
    >>> ResourceRedistribution.name
    'Resource Redistribution'
    >>> sorted(ResourceRedistribution.output_var_names)
    ['soil__resources',
     'vegetation__plant_functional_type']
    >>> sorted(ResourceRedistribution.units)
    [('soil__resources', 'None'),
     ('vegetation__plant_functional_type', 'None')]
    >>> grid.at_cell["vegetation__plant_functional_type"] = (
    ...     np.random.randint(0, 5, size=grid.number_of_cells))
    >>> np.allclose(
    ...    grid.at_cell["vegetation__plant_functional_type"],
    ...    np.array([4, 0, 3, 3, 3, 1]))
    True
    >>> grid.at_cell["soil__resources"] = (
    ...     np.ones(grid.number_of_cells, dtype=float))
    >>> rr = ResourceRedistribution(grid)
    >>> rr.grid.number_of_cell_rows
    3
    >>> rr.grid.number_of_cell_columns
    2
    >>> rr.grid is grid
    True
    >>> (eroded_soil,
    ...  eroded_soil_shrub,
    ...  burnt_shrub,
    ...  burnt_grass,
    ...  bare_cells) = rr.erode()
    >>> np.round(eroded_soil, decimals=2) == 0.16
    True
    >>> burnt_shrub.shape == (1,)
    True
    >>> (burnt_shrubs_neigh,
    ...  exclusive,
    ...  shrub_exclusive,
    ...  grass_exclusive,
    ...  bare_exclusive,
    ...  eroded_soil_part) = rr.deposit(eroded_soil, eroded_soil_shrub)
    >>> np.allclose(burnt_shrubs_neigh, np.array([1, 2, 3, 4, 5]))
    True
    >>> eroded_soil_part == 0
    True
    >>> (resource_adjusted,
    ...  eligible_locs_to_adj_neigh,
    ...  elig_locs,
    ...  sed_to_borrow) = rr.re_adjust_resource()
    >>> resource_adjusted == 0.
    True
    >>> V_age = np.zeros(rr.grid.number_of_cells, dtype=int)
    >>> V_age = rr.initialize_Veg_age(V_age=V_age)
    >>> np.allclose(V_age, np.zeros(rr.grid.number_of_cells, dtype=int))
    True
    >>> (V_age, est_1, est_2, est_3, est_4, est_5) = rr.establish(V_age)
    >>> np.allclose(grid.at_cell["vegetation__plant_functional_type"],
    ...             np.array([4, 1, 3, 3, 3, 1]))
    True
    >>> (V_age, Pmor_age, Pmor_age_ws) = rr.mortality(V_age)
    >>> np.allclose(grid.at_cell["vegetation__plant_functional_type"],
    ...             np.array([4, 1, 0, 0, 0, 1]))
    True
    >>> V_age = rr.update_Veg_age(V_age)
    >>> np.allclose(V_age, np.array([1, 1, 0, 0, 0, 1]))
    True
    """

    _name = "Resource Redistribution"

    _input_var_names = (
        "vegetation__plant_functional_type",
        "soil__resources",
    )

    _output_var_names = (
        "vegetation__plant_functional_type",
        "soil__resources",
    )

    _var_units = {
        "vegetation__plant_functional_type": "None",
        "soil__resources": "None",
    }

    _var_mapping = {
        "vegetation__plant_functional_type": "cell",
        "soil__resources": "cell",
    }

    _var_doc = {
        "vegetation__plant_functional_type": "classification of plant type - bare=0, grass=1,"
        + "shrub=2, burntgrass=3, burntshrub=4",
        "soil__resources": "level of soil resources",
    }

    @use_file_name_or_kwds
    def __init__(
        self,
        grid,
        e=0.1,
        R_low_threshold=-2.0,
        R_threshold=2.0,
        R_dep_threshold=1.0,
        Rth_gr=0.4,
        Rth_sh=0.8,
        P_gr_regrwth=0.25,
        P_sh_regrwth=0.25,
        Pen=0.05,
        Pgrz=0.01,
        P_gr=0.5,
        sh_max_age=600,
        sh_seedling_max_age=18,
        sh_seedling_mor_dis=0.0,
        sh_mor_dis_low_thresh_age=300,
        sh_mor_dis_low_slp=0.01,
        sh_mor_dis_high_slp=0.99,
        P_sh_fire_mor=0.75,
        P_gr_fire_mor=1.0,
        sh_mor_ws_thresh=0.01,
        gr_mor_ws_thresh=0.08,
        **kwds,
    ):
        """
        Parameters:
        ----------
        grid: RasterModelGrid
            grid, Landlab's RasterModelGrid object
        e: float, optional
            soil erosion at bare soil cell (-)
        R_low_threshold: float, optional
            Minimum value of R (-)
        R_threshold: float, optional
            Maximum value of R (-)
        R_dep_threshold: float, optional
            threshold for deposition during adjustment (-)
        Rth_gr: float, optional
            Minimum R for establishment of Grass (-)
        Rth_sh: float, optional
            Minimum R for establishment of Shrub (-)
        P_gr_regrwth: float, optional
            Probability of grass regrowth in burnt grass patches (-)
        P_sh_regrwth: float, optional
            Probability of shrub regrowth in burnt shrub patches (-)
        Pen: float, optional
            Probability of shrub establishment (neighborhood) (-)
        Pgrz: float, optional
            Probability of shrub establishment
            (seed dispersal due to herbivores) (-)
        P_gr: float, optional
            Probability of grass establishment (-)
        sh_max_age: int, optional
            Maximum age of shrubs (yr)
        sh_seedling_max_age: int, optional
            Maximum age of shrub seedlings (yr)
        sh_seedling_mor_dis: float, optional
            Probability of mortality of shrub seedlings due to disease (-)
        sh_mor_dis_low_thresh_age: int, optional
            Age at which shrubs start experiencing disease (yr)
        sh_mor_dis_low_slp: float, optional
            Probability of mortality of shrubs due to disease
            at  sh_mor_dis_low_thresh_age (-)
        sh_mor_dis_high_slp: float, optional
            Probability of mortality of shrubs due to disease
            at sh_max_age (-)
        P_sh_fire_mor: float, optional
            Probability of shrub mortality due to fire (-)
        P_gr_fire_mor: float, optional
            Probability of grass mortality due to fire (-)
        sh_mor_ws_thresh: float, optional
            Background mortality probability of shrubs due
            to waterstress (P_mor_ws) (-)
        gr_mor_ws_thresh: float, optional
            Background mortality probability of grass due
            to waterstress (-)
        """

        super(ResourceRedistribution, self).__init__(grid, **kwds)

        name = "vegetation__plant_functional_type"
        if name not in self.grid.at_cell:
            warnings.warn(
                "Since a cellular field of PFTs "
                + "is not provided, the field "
                + name
                + " is initialized to 0 - all bare cells!"
            )
            self.grid.add_zeros("cell", name, units=self._var_units[name])
        name = "soil__resources"
        if name not in self.grid.at_cell:
            warnings.warn(
                "Since a cellular field for Resources "
                + "is not provided, the field "
                + name
                + " is initialized to 1!"
            )
            self.grid.add_ones("cell", name, units=self._var_units[name])

        self.initialize(
            e=e,
            R_low_threshold=R_low_threshold,
            R_threshold=R_threshold,
            R_dep_threshold=R_dep_threshold,
            Rth_gr=Rth_gr,
            Rth_sh=Rth_sh,
            P_gr_regrwth=P_gr_regrwth,
            P_sh_regrwth=P_sh_regrwth,
            Pen=Pen,
            Pgrz=Pgrz,
            P_gr=P_gr,
            sh_max_age=sh_max_age,
            sh_seedling_max_age=sh_seedling_max_age,
            sh_seedling_mor_dis=sh_seedling_mor_dis,
            sh_mor_dis_low_thresh_age=sh_mor_dis_low_thresh_age,
            sh_mor_dis_low_slp=sh_mor_dis_low_slp,
            sh_mor_dis_high_slp=sh_mor_dis_high_slp,
            P_sh_fire_mor=P_sh_fire_mor,
            P_gr_fire_mor=P_gr_fire_mor,
            sh_mor_ws_thresh=sh_mor_ws_thresh,
            gr_mor_ws_thresh=gr_mor_ws_thresh,
            **kwds,
        )

    def initialize(
        self,
        e=0.1,
        R_low_threshold=-2.0,
        R_threshold=2.0,
        R_dep_threshold=1.0,
        Rth_gr=0.4,
        Rth_sh=0.8,
        P_gr_regrwth=0.25,
        P_sh_regrwth=0.25,
        Pen=0.05,
        Pgrz=0.01,
        P_gr=0.5,
        sh_max_age=600,
        sh_seedling_max_age=18,
        sh_seedling_mor_dis=0.0,
        sh_mor_dis_low_thresh_age=300,
        sh_mor_dis_low_slp=0.01,
        sh_mor_dis_high_slp=0.99,
        P_sh_fire_mor=0.75,
        P_gr_fire_mor=1.0,
        sh_mor_ws_thresh=0.01,
        gr_mor_ws_thresh=0.08,
        **kwds,
    ):
        """
        Parameters:
        ----------
        grid: RasterModelGrid
            grid, Landlab's RasterModelGrid object
        e: float, optional
            soil erosion at bare soil cell (-)
        R_low_threshold: float, optional
            Minimum value of R (-)
        R_threshold: float, optional
            Maximum value of R (-)
        R_dep_threshold: float, optional
            threshold for deposition during adjustment (-)
        Rth_gr: float, optional
            Minimum R for establishment of Grass (-)
        Rth_sh: float, optional
            Minimum R for establishment of Shrub (-)
        P_gr_regrwth: float, optional
            Probability of grass regrowth in burnt grass patches (-)
        P_sh_regrwth: float, optional
            Probability of shrub regrowth in burnt shrub patches (-)
        Pen: float, optional
            Probability of shrub establishment (neighborhood) (-)
        Pgrz: float, optional
            Probability of shrub establishment
            (seed dispersal due to herbivores) (-)
        P_gr: float, optional
            Probability of grass establishment (-)
        sh_max_age: int, optional
            Maximum age of shrubs (yr)
        sh_seedling_max_age: int, optional
            Maximum age of shrub seedlings (yr)
        sh_seedling_mor_dis: float, optional
            Probability of mortality of shrub seedlings due to disease (-)
        sh_mor_dis_low_thresh_age: int, optional
            Age at which shrubs start experiencing disease (yr)
        sh_mor_dis_low_slp: float, optional
            Probability of mortality of shrubs due to disease
            at  sh_mor_dis_low_thresh_age (-)
        sh_mor_dis_high_slp: float, optional
            Probability of mortality of shrubs due to disease
            at sh_max_age (-)
        P_sh_fire_mor: float, optional
            Probability of shrub mortality due to fire (-)
        P_gr_fire_mor: float, optional
            Probability of grass mortality due to fire (-)
        sh_mor_ws_thresh: float, optional
            Background mortality probability of shrubs due
            to waterstress (P_mor_ws) (-)
        gr_mor_ws_thresh: float, optional
            Background mortality probability of grass due
            to waterstress (-)
        """
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
        """
        This method removes resource R from a cell depending on the
        vegetation (PFT) occupying each cell in discrete units of e.
        """
        V = self.grid.at_cell["vegetation__plant_functional_type"]
        R = self.grid.at_cell["soil__resources"]
        Elig_R = np.where(R > self._R_low_threshold)[0]
        # Deal with shrubs first
        burnt_shrubs = Elig_R[V[Elig_R] == BURNTSHRUB]
        eroded_soil_shrub = np.int(burnt_shrubs.shape[0]) * 4 * self._e
        R[burnt_shrubs] -= 4 * self._e
        # Deal with erosion of other cells
        burnt_grass = Elig_R[V[Elig_R] == BURNTGRASS]
        bare_cells = Elig_R[V[Elig_R] == BARE]
        R[burnt_grass] -= 0.2 * self._e
        R[bare_cells] -= self._e
        eroded_soil = (
            np.int(burnt_grass.shape[0]) * 0.2 * self._e
            + np.int(bare_cells.shape[0]) * self._e
        )
        return (eroded_soil, eroded_soil_shrub, burnt_shrubs, burnt_grass, bare_cells)

    def deposit(self, eroded_soil, eroded_soil_shrub):
        """
        This method deposits resource R (generally eroded using
        the erode method) in discrete units of e.

        Parameters
        ----------
        eroded_soil: float, (-)
            R eroded from cells with burnt grass and bare cells.
            This R is deposited based on the algorithm described
            in Ravi et. al 2009.
        eroded_soil_shrub: float, (-)
            R eroded from cells with burnt shrubs. This R
            is deposited in the cells neighboring burnt shrubs.

        Returns
        -------
        A tuple: (burnt_shrubs_neigh,
                  exclusive,
                  shrub_exclusive,
                  grass_exclusive,
                  bare_exclusive,
                  eroded_soil_part);
        burnt_shrubs_neigh: numpy array of ints, (-)
            cell_ids of looped shrub neighbors (neighbors of cells
            occupied by BURNTSHRUB).
        exclusive: numpy array of ints, (-)
            cell_ids of cells except those occupied by BURNTSHRUB and
            burnt_shrubs_neigh. (purpose - debugging)
        shrub_exclusive: numpy array of ints, (-)
            cell_ids of cells within exclusive that are occupied
            by SHRUB. (purpose - debugging)
        grass_exclusive: numpy array of ints, (-)
            cell_ids of cells within exclusive that are occupied
            by GRASS. (purpose - debugging)
        bare_exclusive: numpy array of ints, (-)
            cell_ids of cells within exclusive that are occupied
            by BARE. (purpose - debugging)
        """
        V = self.grid.at_cell["vegetation__plant_functional_type"]
        R = self.grid.at_cell["soil__resources"]
        exclusive = np.arange(0, self.grid.number_of_cells)
        burnt_shrubs = np.where(V == BURNTSHRUB)[0]
        if int(burnt_shrubs.shape[0]) > 0:
            burnt_shrubs_neigh = np.unique(
                self.grid.looped_neighbors_at_cell[burnt_shrubs]
            )
            R[burnt_shrubs_neigh] += eroded_soil_shrub / float(
                burnt_shrubs_neigh.shape[0]
            )
            exclusive = exclusive[np.in1d(exclusive, burnt_shrubs, invert=True)]
            exclusive = exclusive[np.in1d(exclusive, burnt_shrubs_neigh, invert=True)]
        else:
            burnt_shrubs_neigh = []
        shrub_exclusive = exclusive[np.where(V[exclusive] == SHRUB)[0]]
        grass_exclusive = exclusive[np.where(V[exclusive] == GRASS)[0]]
        bare_exclusive = exclusive[np.where(V[exclusive] == BARE)[0]]
        # Calculating how much each bare cell will get from eroded_soil 'Eb'
        # Shrubs will receive thrice Es = 3*Eb && Eg = 2*Eb
        weighted_parts = (
            1 * int(bare_exclusive.shape[0])
            + 3 * int(shrub_exclusive.shape[0])
            + 2 * int(grass_exclusive.shape[0])
        )
        if weighted_parts != 0:
            eroded_soil_part = eroded_soil / float(weighted_parts)
        else:
            eroded_soil_part = 0
        R[bare_exclusive] += eroded_soil_part
        R[shrub_exclusive] += 3.0 * eroded_soil_part
        R[grass_exclusive] += 2.0 * eroded_soil_part
        return (
            burnt_shrubs_neigh,
            exclusive,
            shrub_exclusive,
            grass_exclusive,
            bare_exclusive,
            eroded_soil_part,
        )

    def re_adjust_resource(self):
        """
        This method re-adjusts resource R based on the
        user-set maximum and minimum thresholds.
        Returns
        -------
        A tuple: (resource_adjusted,
                  eligible_locs_to_adj_neigh,
                  Elig_locs,
                  sed_to_borrow);
        resource_adjusted: float, (-)
            cumulative amount of R in excess of R_threshold.
            (purpose - debugging)
        eligible_locs_to_adj_neigh: numpy array of ints, (-)
            cell_ids of neighbors of cells with
            R > R_threshold. (purpose - debugging)
        Elig_locs: numpy array of ints, (-)
            cell_ids of cells where R < R_threshold.
            R is deposited into these cells. (purpose - debugging)
        sed_to_borrow: float, (-)
            cumulative amount of R < R_threshold that
            is borrowed. (purpose - debugging)
        """
        R = self.grid.at_cell["soil__resources"]
        # Resource exceeding R_threshold will be distributed to its neighbors
        resource_adjusted = 0.0
        eligible_locs_to_adj_neigh = np.array([])
        locs_to_adj = np.where(R > self._R_threshold)[0]
        if int(locs_to_adj.shape[0]) > 0:
            resource_adjusted = np.sum(R[locs_to_adj] - self._R_threshold)
            locs_to_adj_neigh = np.unique(
                self.grid.looped_neighbors_at_cell[locs_to_adj]
            )
            eligible_locs_to_adj_neigh = locs_to_adj_neigh[
                R[locs_to_adj_neigh] < self._R_dep_threshold
            ]
            if int(eligible_locs_to_adj_neigh.shape[0]) > 0:
                R[eligible_locs_to_adj_neigh] += resource_adjusted / float(
                    eligible_locs_to_adj_neigh.shape[0]
                )
                R[locs_to_adj] = self._R_threshold
        # Resource below R_low_threshold is raised to R_low_threshold
        sed_to_borrow = 0.0
        Elig_locs = np.where(R < self._R_low_threshold)[0]
        if int(Elig_locs.shape[0]) > 0:
            sed_to_borrow = np.absolute(np.sum(self._R_low_threshold - R[Elig_locs]))
            R[Elig_locs] = self._R_low_threshold
            locs_to_borrow = np.where(R > 0.0)[0]
            R[locs_to_borrow] -= sed_to_borrow / float(locs_to_borrow.shape[0])
        return (resource_adjusted, eligible_locs_to_adj_neigh, Elig_locs, sed_to_borrow)

    def establish(self, V_age):
        """
        This method implements the establishment algorithm for PFT
        based on the rules outlined in Ravi et al. 2009. This
        method updates the field "vegetation__plant_functional_type".

        Parameters
        ----------
        V_age: numpy array of ints, shape=[number_of_cells] (yrs)
            age of the PFT (V) occupying each cell.

        Returns
        -------
        A tuple: (V_age, est_1, est_2, est_3, est_4, est_5);
        V_age: numpy array of ints, (yrs)
            age of the PFT (V) occupying each cell. This
            array is updated by this method.
        est_1: numpy array of ints, (-)
            cell_ids of BURNTSHRUB cells where SHRUB regrows
            (established). (purpose - debugging)
        est_2: numpy array of ints, (-)
            cell_ids of BURNTGRASS cells where GRASS regrows
            (established). (purpose - debugging)
        est_3: numpy array of ints, (-)
            cell_ids of BARE cells where SHRUB establishes
            due to seed dispersal. (purpose - debugging)
        est_4: numpy array of ints, (-)
            cell_ids of BARE cells where SHRUB establishes
            due to grazing. (purpose - debugging)
        est_5: numpy array of ints, (-)
            cell_ids of BARE cells where GRASS establishes
            due to seed dispersal. (purpose - debugging)
        """
        V = self.grid.at_cell["vegetation__plant_functional_type"]
        R = self.grid.at_cell["soil__resources"]
        burnt_shrubs = np.where(V == BURNTSHRUB)[0]
        burnt_grass = np.where(V == BURNTGRASS)[0]
        # Regrowth in burnt area
        shrubs_r_regrwth = burnt_shrubs[np.where(R[burnt_shrubs] > self._Rth_sh)[0]]
        grass_r_regrwth = burnt_grass[np.where(R[burnt_grass] > self._Rth_gr)[0]]
        P_check_1 = np.random.random(shrubs_r_regrwth.shape)
        est_1 = shrubs_r_regrwth[np.where(P_check_1 < self._P_sh_regrwth)[0]]
        V[est_1] = SHRUB
        V_age[est_1] = 0
        P_check_2 = np.random.random(grass_r_regrwth.shape)
        est_2 = grass_r_regrwth[np.where(P_check_2 < self._P_gr_regrwth)[0]]
        V[est_2] = GRASS
        V_age[est_2] = 0
        # Regrowth in grazed area
        # shrub encroachment due to neighbors
        bare_cells = np.where(V == BARE)[0]
        shrub_grz_regrwth = bare_cells[np.where(R[bare_cells] > self._Rth_sh)[0]]
        neigh_sh_grz_regrwth = self.grid.looped_neighbors_at_cell[shrub_grz_regrwth]
        ns_sh = self._np_ndarray_count(neigh_sh_grz_regrwth)
        ns_Pgrz = ns_sh * self._Pen  # ns * P2
        P_check_3 = np.random.random(ns_Pgrz.shape)
        est_3 = shrub_grz_regrwth[np.where(P_check_3 < ns_Pgrz)[0]]
        V[est_3] = SHRUB  # Establish shrubs
        V_age[est_3] = 0
        # shrub encroachment due to grazing
        bare_cells_ = np.where(V == BARE)[0]
        shrub_grz_regrwth_ = bare_cells_[np.where(R[bare_cells_] > self._Rth_sh)[0]]
        P_check_4 = np.random.random(shrub_grz_regrwth_.shape)
        est_4 = shrub_grz_regrwth_[np.where(P_check_4 < self._Pgrz)[0]]
        V[est_4] = SHRUB  # Establish grass
        V_age[est_4] = 0
        # grass growth where shrubs haven't encroached - grass seed dispersal
        bare_cells_2 = np.where(V == BARE)[0]
        grass_grz_regrwth = bare_cells_2[np.where(R[bare_cells_2] > self._Rth_gr)[0]]
        P_check_5 = np.random.random(grass_grz_regrwth.shape)
        est_5 = grass_grz_regrwth[np.where(P_check_5 < self._P_gr)[0]]
        V[est_5] = GRASS  # Establish grass
        V_age[est_5] = 0
        return (V_age, est_1, est_2, est_3, est_4, est_5)

    def _np_ndarray_count(self, neigh_sh_grz_regrwth):
        V = self.grid.at_cell["vegetation__plant_functional_type"]
        ns_sh = np.zeros(neigh_sh_grz_regrwth.shape[0], dtype=int)
        for i in range(0, ns_sh.shape[0]):
            ns_sh[i] = np.int(np.where(V[neigh_sh_grz_regrwth[i, :]] == 2)[0].shape[0])
        return ns_sh

    def _compute_Prob_mortality_age(self, V_age, Pmor_age):
        V = self.grid.at_cell["vegetation__plant_functional_type"]
        # Age mortality for Vegetation other than shrub = 0.
        # Shrub seedling
        shrubs_ = np.where(V == SHRUB)[0]
        shrub_seedlings = shrubs_[
            np.where(V_age[shrubs_] < self._sh_seedling_max_age)[0]
        ]
        Pmor_age[shrub_seedlings] = self._sh_seedling_mor_dis
        young_shrubs = shrubs_[
            np.logical_and(
                V_age[shrubs_] > self._sh_seedling_max_age,
                V_age[shrubs_] < self._sh_mor_dis_low_thresh_age,
            )
        ]
        Pmor_age[young_shrubs] = self._sh_mor_dis_low_slp  # This is constant
        adult_shrubs = shrubs_[
            np.where(V_age[shrubs_] > self._sh_mor_dis_low_thresh_age)
        ]
        Pmor_age[adult_shrubs] = self._sh_mor_dis_low_slp + (
            (np.asfarray(V_age[adult_shrubs]) - self._sh_mor_dis_low_thresh_age)
            / (self._sh_max_age - self._sh_mor_dis_low_thresh_age)
            * (self._sh_mor_dis_high_slp)
        )
        return Pmor_age

    def mortality(self, V_age):
        """
        This method implements the mortality algorithm for PFT
        based on the rules outlined in Ravi et al. 2009. This
        method updates the field "vegetation__plant_functional_type".

        Parameters
        ----------
        V_age: numpy array of ints, shape=[number_of_cells] (yrs)
            Age of the PFT (V) occupying each cell. This
            array is updated by this method.

        Returns
        -------
        A tuple: (V_age, Pmor_age, Pmor_age_ws);
        V_age: numpy array of ints, (yrs)
            age of the PFT (V) occupying each cell. This
            array is updated by this method.
        Pmor_age: numpy array of floats, (-)
            probability of mortality due to age
            at each cell. Note that only SHRUB are
            susceptible to this probability.
            (purpose - debugging)
        Pmor_age_ws: numpy array of floats, (-)
            probability of mortality due to age and
            water stress at each cell. Note that
            probability of mortality due to waterstress
            can be obtained by Pmor_age_ws - Pmor_age.
            (purpose - debugging)
        """
        V = self.grid.at_cell["vegetation__plant_functional_type"]
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
        Pmor_age_ws[V == GRASS] += (
            0.6 * self._gr_mor_ws_thresh
            + 2 * 0.4 * self._gr_mor_ws_thresh * np.random.random_sample()
        )
        # Grass' water stess probability
        Pmor_age_ws[V == SHRUB] += (
            0.6 * self._sh_mor_ws_thresh
            + 2 * 0.4 * self._sh_mor_ws_thresh * np.random.random_sample()
        )
        # Shrubs' water stress probability
        Pmor_age = self._compute_Prob_mortality_age(V_age, Pmor_age)
        Pmor_age_ws += Pmor_age
        P_check_3 = np.random.random(V.shape)
        kill_3 = np.where(P_check_3 < Pmor_age_ws)[0]
        V[kill_3] = BARE
        V_age[kill_3] = 0
        return (V_age, Pmor_age, Pmor_age_ws)

    def initialize_Veg_age(self, V_age):
        """
        This method initializes random age for shrubs. The input
        array, is used for shape. TODO: Use self.grid.number_of_cells
        for shape.

        Parameters
        ----------
        V_age: numpy array of ints, shape=[number_of_cells] (yrs)
            Age of the PFT (V) occupying each cell.

        Returns
        -------
        V_age: numpy array of ints, (yrs)
            age of the PFT (V) occupying each cell. This
            array is updated by this method.
        """
        V = self.grid.at_cell["vegetation__plant_functional_type"]
        V_age[V == SHRUB] = np.random.randint(
            0, self._sh_max_age, V_age[V == SHRUB].shape
        )
        return V_age

    def update_Veg_age(self, V_age):
        """
        This method updates V_age by one year at
        each vegetated cell (V != BARE)

        Parameters
        ----------
        V_age: numpy array of ints, shape=[number_of_cells] (yrs)
            Age of the PFT (V) occupying each cell.

        Returns
        -------
        V_age: numpy array of ints, (yrs)
            age of the PFT (V) occupying each cell. This
            array is updated by this method.
        """
        V = self.grid.at_cell["vegetation__plant_functional_type"]
        V_age[V != BARE] += 1
        return V_age
