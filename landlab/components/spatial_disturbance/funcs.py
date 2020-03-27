
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 22:23:12 2018

@author: saisiddu
"""
import numpy as np


def convert_phy_pft_to_distr_pft(grid, vegtype):
    """
    Convert PFT in Physical Ecohyd model to PFT in Create_fires/resource_
    redistribution codes.
    """
    BARE_distr = 0
    GRASS_distr = 1
    SHRUB_distr = 2
    TREE_distr = 5
    SHRUBSEED_distr = 7
    TREESEED_distr = 8
#    BARE_phy = 3
    GRASS_phy = 0
    SHRUB_phy = 1
    TREE_phy = 2
    SHRUBSEED_phy = 4
    TREESEED_phy = 5
    V = (BARE_distr * np.ones(grid.number_of_cells, dtype=int))
    grass_phy_cells = np.where(vegtype == GRASS_phy)[0]
    shrub_phy_cells = np.where(vegtype == SHRUB_phy)[0]
    tree_phy_cells = np.where(vegtype == TREE_phy)[0]
    shrubseed_phy_cells = np.where(vegtype == SHRUBSEED_phy)[0]
    treeseed_phy_cells = np.where(vegtype == TREESEED_phy)[0]
    V[grass_phy_cells] = GRASS_distr
    V[shrub_phy_cells] = SHRUB_distr
    V[tree_phy_cells] = TREE_distr
    V[shrubseed_phy_cells] = SHRUBSEED_distr
    V[treeseed_phy_cells] = TREESEED_distr
    return V


def convert_distr_pft_to_phy_pft(grid, V):
    """
    Convert PFT in Create_fires/resource_redistribution codes to
    PFT in Physical Ecohyd model.
    """
    #    BARE_distr = 0
    GRASS_distr = 1
    SHRUB_distr = 2
    TREE_distr = 5
    SHRUBSEED_distr = 7
    TREESEED_distr = 8
    BARE_phy = 3
    GRASS_phy = 0
    SHRUB_phy = 1
    TREE_phy = 2
    SHRUBSEED_phy = 4
    TREESEED_phy = 5
    vegtype = (BARE_phy * np.ones(grid.number_of_cells, dtype=int))
    grass_distr_cells = np.where(V == GRASS_distr)[0]
    shrub_distr_cells = np.where(V == SHRUB_distr)[0]
    tree_distr_cells = np.where(V == TREE_distr)[0]
    shrubseed_distr_cells = np.where(V == SHRUBSEED_distr)[0]
    treeseed_distr_cells = np.where(V == TREESEED_distr)[0]
    vegtype[grass_distr_cells] = GRASS_phy
    vegtype[shrub_distr_cells] = SHRUB_phy
    vegtype[tree_distr_cells] = TREE_phy
    vegtype[shrubseed_distr_cells] = SHRUBSEED_phy
    vegtype[treeseed_distr_cells] = TREESEED_phy
    return vegtype
