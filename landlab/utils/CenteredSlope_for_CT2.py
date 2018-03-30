#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 10:14:52 2018

@author: margaux
"""

import numpy as np
import gc
from landlab.item_collection import ItemCollection


class ClastCollection(ItemCollection):
    """
    """
    def __init__(self,
                 grid,
                 clast_x=[],
                 clast_y=[],
                 clast_elev=[],
                 clast_radius=[]):
        """
        """

        # Save a reference to the grid
        self._grid = grid
#
#        # Store data
#        self._clast__x = clast_x
#        self._clast__y = clast_y
#        self._clast__elev = clast_elev
#        self._clast__radius = clast_radius
#

        # Determine reference (closest) node for each clast
        self._clast__node = self._grid.find_nearest_node((self._clast__x, self._clast__y))
        # Clast set size:
        self._nb_of_clast = len(self._clast__node)

        # Store the input information and others in a dictionary:
        
        clast_data = {'clast__x' : clast_x,
                      'clast__y' : clast_y,
                      'clast__elev' : clast_elev,
                      'clast__node' : self._clast__node,
                      'lambda_0' : np.zeros(self._nb_of_clast),
                      'lambda_mean' : np.zeros(self._nb_of_clast),
                      'slope__WE' : np.zeros(self._nb_of_clast),
                      'slope__SN' : np.zeros(self._nb_of_clast),
                      'slope__steepest_azimuth' : np.full(self._nb_of_clast, np.NaN),
                      'slope__steepest_dip' : np.zeros(self._nb_of_clast),
                      'distance__to_exit' : np.full(self._nb_of_clast, np.NaN),
                      'target_node' : -np.ones(self._nb_of_clast, dtype=int),
                      'change_x' : np.zeros(self._nb_of_clast),
                      'change_y' : np.zeros(self._nb_of_clast),
                      'hop_length' : np.zeros(self._nb_of_clast),
                      'total_travelled_dist' : np.zeros(self._nb_of_clast)
                      }

        # Build ItemCollection containing clast data:
        ItemCollection.__init__(self,               # necessary?
                                self._grid,
                                data=clast_data,
                                grid_element='node',
                                element_id=self._clast__node)


    def _neighborhood(self, clast):
        self._clast__node = self._grid.find_nearest_node((self._clast__x, self._clast__y)) # needs this to update after clast has moved
        _node = self._clast__node[clast]
        _grid = self._grid
        ### Calculation of slopes: W to E and S to N
        we_slope=(_grid.at_node['topographic__elevation'][west_node]-
                  _grid.at_node['topographic__elevation'][east_node])/(2*rmg.dx)
        
        sn_slope=(self._grid.at_node['topographic__elevation'][south_node]-
                  self._grid.at_node['topographic__elevation'][north_node])/(2*self._grid.dy)
        
        
        ### COORDINATES OF CORNERS OF CELL
        x_of_corners_at_node=np.zeros((self._grid.number_of_nodes,4))
        y_of_corners_at_node=np.zeros((self._grid.number_of_nodes,4))
        
        for i in range(self._grid.number_of_nodes):
            x_of_corners_at_node[i,:]=[self._grid.node_x[i]+self._grid.dx/2, self._grid.node_x[i]-self._grid.dx/2, self._grid.node_x[i]-self._grid.dx/2, self._grid.node_x[i]+self._grid.dx/2]
            y_of_corners_at_node[i,:]=[self._grid.node_y[i]+self._grid.dx/2, self._grid.node_x[i]+self._grid.dx/2, self._grid.node_x[i]-self._grid.dx/2, self._grid.node_x[i]-self._grid.dx/2]
        
        ### Adjacent row and col nodes
        _row_col_adjacent_nodes_at_node = self._grid.neighbors_at_node[_node]
        ### Adjacent diagonal nodes
        _diagonal_adjacent_nodes_at_node = self._grid.diagonal_adjacent_nodes_at_node[_node]
    
    
    ### Determination of direction and value of steepest slope (ss)
    def _move_to(self, clast):
        _node = self.at[clast, 'clast__node']
        _z_node = self._grid.at_node['topographic__elevation'][_node]
        if we_slope == 0:
            if sn_slope == 0:
                # centered slope is null, clast does not move:
                ss_dip = 0
                ss_azimuth = None
                dist_to_exit = None
                target_node = -1
    #            # OPTION2 (to develop):
    #            # clast moves to random direction, according to lambda_0:

            else: # SN slope is not 0, ss direction is S or N
                if sn_slope < 0: # ss direction is South
                    ss_dip = np.abs(sn_slope)
                    ss_azimuth = np.radians(270) # South
                    dist_to_exit = self._clast_y - (self.node_y[_node]-(dy/2))
                    target_node = _row_col_adjacent_nodes_at_node[3]
                else: # ss direction is North
                    ss_dip = np.abs(sn_slope)
                    ss_azimuth = np.radians(90) # North
                    dist_to_exit = (self.node_y[_node]+(dy/2)) - self._clast_y
                    target_node = _row_col_adjacent_nodes_at_node[1]
        else: # we_slope is not 0
            if sn_slope == 0:
                if we_slope < 0: # ss direction is West
                    ss_dip = np.abs(we_slope)
                    ss_azimuth = np.radians(180) # West
                    dist_to_exit = self._clast_x - (self.node_x[_node]-(dx/2))
                    target_node = _row_col_adjacent_nodes_at_node[2]
                else: # ss direction is East
                    ss_dip = np.abs(we_slope)
                    ss_azimuth = 0 # East
                    dist_to_exit = (self.node_y[_node]+(dx/2)) - self._clast_x
                    target_node = _row_col_adjacent_nodes_at_node[0]
            else: # sn_slope is not 0
                ss_dip = np.sqrt(np.power(sn_slope, 2) + np.power(we_slope, 2))
                if sn_slope > 0 and we_slope > 0: # Quarter = NE
                    ss_azimuth = np.arctan(np.abs(sn_slope / we_slope))
                    corner = 0
                    clast-corner_azimuth = np.arctan(
                            np.abs(self._clast_y-y_of_corners_at_node[node, corner])/
                                   np.abs(self._clast_x-x_of_corners_at_node[node, corner]))
                    if ss_azimuth < clast-corner_azimuth: # Eigth = NE-row
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(x_of_corners_at_node[node, corner]-self._clast_x)) / (np.cos(ss_azimuth)))
                        target_node = _row_col_adjacent_nodes_at_node[0]
                        # Coordinates of vector clast-to-border:
                        [change_x, change_y] = [dist_to_exit / np.cos(ss_azimuth), dist_to_exit / np.sin(ss_azimuth)]
                    elif ss_azimuth > clast-corner_azimuth: # Eigth = NE-col
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(y_of_corners_at_node[node, corner]-self._clast_y)) / (np.cos(np.radians(90)-ss_azimuth)))
                        target_node = _row_col_adjacent_nodes_at_node[1]
                        [change_x, change_y] = [-dist_to_exit / np.sin(np.radians(90)-ss_azimuth), dist_to_exit / np.cos(np.radians(90)-ss_azimuth)]
                    elif ss_azimuth == clast-corner_azimuth: # exit direction is diagonal
                        dist_to_exit = np.sqrt(np.power(np.abs(x_of_corners_at_node[node, corner]-self._clast_x), 2) + np.power(np.abs(y_of_corners_at_node[node, corner]-self._clast_y), 2))
                        target_node = _diagonal_adjacent_nodes_at_node[0]
                        # Coordinates of vector clast-to-border:
                        [change_x, change_y] = [x_of_corners_at_node[node, corner]-self._clast_x, y_of_corners_at_node[node, corner] - self._clast_y]
                elif sn_slope > 0 and we_slope < 0: # Quarter = NW
                    ss_azimuth = np.radians(90) + np.arctan(np.abs(sn_slope / we_slope))
                    corner = 1
                    clast-corner_azimuth = np.radians(180) - np.arctan(
                            np.abs(self._clast_y-y_of_corners_at_node[node, corner])/
                                   np.abs(self._clast_x-x_of_corners_at_node[node, corner]))
                    if ss_azimuth < clast-corner_azimuth: # Eigth = NW-col
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(y_of_corners_at_node[node, corner]-self._clast_y)) / (np.cos(ss_azimuth-np.radians(90))))
                        [change_x, change_y] = [-dist_to_exit / np.sin(ss_azimuth-np.radians(90)), dist_to_exit / np.cos(ss_azimuth-np.radians(90))]
                        target_node = _row_col_adjacent_nodes_at_node[1]
                    elif ss_azimuth > clast-corner_azimuth: # Eigth = NW-row
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(x_of_corners_at_node[node, corner]-self._clast_x)) / (np.cos(np.radians(180)-ss_azimuth)))
                        [change_x, change_y] = [-dist_to_exit / np.cos(np.radians(180)-ss_azimuth), dist_to_exit / np.sin(np.radians(180)-ss_azimuth)]
                        target_node = _row_col_adjacent_nodes_at_node[2]
                    elif ss_azimuth == clast-corner_azimuth: # exit direction is diagonal
                        dist_to_exit = np.sqrt(np.power(np.abs(x_of_corners_at_node[node, corner]-self._clast_x), 2) + np.power(np.abs(y_of_corners_at_node[node, corner]-self._clast_y), 2))
                        [change_x, change_y] = [-(x_of_corners_at_node[node, corner]-self._clast_x), y_of_corners_at_node[node, corner] - self._clast_y]
                        target_node = _diagonal_adjacent_nodes_at_node[1]
    
                elif sn_slope < 0 and we_slope < 0: # Quarter = SW
                    ss_azimuth = np.radians(180) + np.arctan(np.abs(sn_slope / we_slope))
                    corner = 2
                    clast-corner_azimuth = np.radians(180) + np.arctan(
                            np.abs(self._clast_y-y_of_corners_at_node[node, corner])/
                                   np.abs(self._clast_x-x_of_corners_at_node[node, corner]))
                    if ss_azimuth < clast-corner_azimuth: # Eigth = SW-row
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(x_of_corners_at_node[node, corner]-self._clast_x)) / (np.cos(ss_azimuth-np.radians(180))))
                        [change_x, change_y] = [-dist_to_exit / np.cos(ss_azimuth- np.radians(180)), -dist_to_exit / np.sin(ss_azimuth-np.radians(180))]
                        target_node = _row_col_adjacent_nodes_at_node[2]
                    elif ss_azimuth > clast-corner_azimuth: # Eigth = SW-col
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(y_of_corners_at_node[node, corner]-self._clast_y)) / (np.cos(np.radians(270)-ss_azimuth)))
                        [change_x, change_y] = [-dist_to_exit / np.sin(np.radians(270)-ss_azimuth), -dist_to_exit / np.cos(np.radians(270)-ss_azimuth)]
                        target_node = _row_col_adjacent_nodes_at_node[3]
                    elif ss_azimuth == clast-corner_azimuth: # exit direction is diagonal
                        dist_to_exit = np.sqrt(np.power(np.abs(x_of_corners_at_node[node, corner]-self._clast_x), 2) + np.power(np.abs(y_of_corners_at_node[node, corner]-self._clast_y), 2))
                        [change_x, change_y] = [-(x_of_corners_at_node[node, corner]-self._clast_x), -(y_of_corners_at_node[node, corner] - self._clast_y)]
                        target_node = _diagonal_adjacent_nodes_at_node[2]
    
                elif sn_slope < 0 and we_slope < 0: # Quarter = SE
                    ss_azimuth = np.radians(270) + np.arctan(np.abs(sn_slope / we_slope))
                    corner = 3
                    clast-corner_azimuth = np.radians(360) - np.arctan(
                            np.abs(self._clast_y-y_of_corners_at_node[node, corner])/
                                   np.abs(self._clast_x-x_of_corners_at_node[node, corner]))
                    if ss_azimuth < clast-corner_azimuth: # Eigth = SE-col
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(y_of_corners_at_node[node, corner]-self._clast_y)) / (np.cos(ss_azimuth-np.radians(270))))
                        [change_x, change_y] = [dist_to_exit / np.sin(ss_azimuth-np.radians(270)), -dist_to_exit / np.cos(ss_azimuth-np.radians(270))]
                        target_node = _row_col_adjacent_nodes_at_node[3]
                    elif ss_azimuth > clast-corner_azimuth: # Eigth = SE-row
                        dist_to_exit = (1 / cos(ss_dip)) * ((np.abs(x_of_corners_at_node[node, corner]-self._clast_yx) / (np.cos(np.radians(360)-ss_azimuth))))
                        [change_x, change_y] = [dist_to_exit / np.cos(np.radians(360)-ss_azimuth), -dist_to_exit / np.sin(np.radians(360)-ss_azimuth)]
                        target_node = _row_col_adjacent_nodes_at_node[0]
                    elif ss_azimuth == clast-corner_azimuth: # exit direction is diagonal
                        dist_to_exit = np.sqrt(np.power(np.abs(x_of_corners_at_node[node, corner]-self._clast_x), 2) + np.power(np.abs(y_of_corners_at_node[node, corner]-self._clast_y), 2))
                        [change_x, change_y] = [(x_of_corners_at_node[node, corner]-self._clast_x), -(y_of_corners_at_node[node, corner] - self._clast_y)]
                        target_node = _diagonal_adjacent_nodes_at_node[3]
        # Calculate lambda_mean:
        _lambda_mean = self._lambda_0[i] * (Si + S) / (Si - S)
        self.at[clast, 'lambda_mean']


        # Save values to dataframe:
        self.at[clast, 'slope__steepest_azimuth'] = ss_azimuth
        self.at[clast, 'slope__steepest_dip'] = ss_dip
        self.at[clast, 'distance__to_exit'] = dist_to_exit
        self.at[clast, 'target_node'] = target_node
        self.at[clast, 'change_x'] = change_x
        self.at[clast, 'change_y'] = change_y






    #def _test_leave_cell(self, clast):
    #    _node = self._clast__node[clast]
    #    _z_node = self._grid.at_node['topographic__elevation'][_node]
    #    S=ss_dip[clast]
    #    Si = self._Si
    #    _azimuth=azimuth
    #    lambda_mean = self._lambda_0[clast] * (Si + S) / (Si - S)
    #    R = np.random.rand(1)
    #
    #    if S >= Si:
    #        proba_leave_cell = 1
    #    else:
    #        proba_leave_cell = np.exp(((dist_to_exit) / np.cos(np.arctan(S)) / lambda_mean))
    #
    #    _move = np.zeros(1, dtype=bool)
    #
    #    if proba_leave_cell >= R:
    #        # Clast leaves node:
    #        _move = True
    #    else:
    #        _move = False
    #
    #    return _move
    
    
    #
    #
    #def _move_clast_change_cell(self, clast):
    #        self._clast_x[clast]+=change_x
    #        self._clast_y[clast]+=change_y
    #        hop_length+=dist_to_exit[clast]
    #        
    #        
    #        
    #def _move_clast_in_cell(self, clast):
    #    
    
    def _change_cell_proba(self, clast):
        lambda_mean = self.at[clast, 'lambda_mean']
        dist_to_exit = self.at[clast, 'dist_to_exit']

        # Draw a random sample in the probability distribution of travel distances:
        self.rand_length = np.random.exponential(scale=lambda_mean, 1)

        if self.rand_length < dist_to_exit: # clast stays in cell
            _change_cell = False
        else: # self.rand_length >= dist_to_exit: clast leaves cell
            _change_cell = True

        return _change_cell

    def _move_in_cell(self, clast):
        # clast stays in cell, move of distance rand_length along slope
        ss_azimuth = self.at[clast, 'slope__steepest_azimuth']
        ss_dip = self.at[clast, 'slope__steepest_dip']
        x_horizontal = self.rand_length * np.cos(ss_dip)
        if ss_azimuth <= np.radians(90):
            [change_x, change_y] = [x_horizontal * np.cos(ss_azimuth), x_horizontal * np.sin(ss_azimuth)]
        elif ss_azimuth <= np.radians(180):
            [change_x, change_y] = [-x_horizontal * np.cos(np.radians(180)-ss_azimuth), x_horizontal * np.sin(np.radians(180)-ss_azimuth)]
        elif ss_azimuth <= np.radians(270):
            [change_x, change_y] = [-x_horizontal * np.sin(np.radians(270)-ss_azimuth), -x_horizontal * np.cos(np.radians(270)-ss_azimuth)]
        else: # ss_azimuth <= np.radians(360)
            [change_x, change_y] = [x_horizontal * np.cos(np.radians(360)-ss_azimuth), -x_horizontal * np.sin(np.radians(360)-ss_azimuth)]
        # Update clast coordinates:
        self.at[clast, 'clast__x'] += change_x
        self.at[clast, 'clast__y'] += change_y
        self.at[clast, 'change_x'] = change_x
        self.at[clast, 'change_y'] = change_y


        self.at[clast, 'hop_length'] += self.rand_length

    def _move_out_of_cell(self, clast):
        # clast leaves cell, move of distance dist_to_exit along slope
        self.at[clast, 'clast__x'] += change_x
        self.at[clast, 'clast__y'] += change_y
        self.at[clast, 'hop_length'] += self.at[clast, 'dist_to_exit']
        self.at[clast, 'clast__node'] = self.at[clast, 'target_node']
        # Save values to dataframe:

########JUST FOR TESTING PURPOSE##############################################
        if self.at[clast, 'clast__node'] != self.at[clast, 'target_node']:
            print('Error: moved to wrong node')
##############################################################################
            
            
    def clast_solver_Exponential(self, dt=1., Si=1.2, kappa=0.0001, tau_0=1, uplift=None): # lambda_0=1, 
    
            # Store various values that will be used
            self._kappa = kappa
            self._dt = dt
            self._erosion__depth = self._erosion_rate * self._dt
            self._deposition__thickness = self._deposition_rate * self._dt
            self._Si = Si #  slope above which particle motion continues indefinetly (not equal to critical slope of diffusion, see Furbish and Haff, 2010)
            hop_length_i = 0
#            self._lambda_0 = lambda_0
    
            self._lambda_0=np.zeros(self._nb_of_clast)
            self._lambda_0[:] = (self._dt * self._kappa * self._grid.dx) / (self._Si * 2 * np.array(self._clast__radius))
#            self._lambda_0[:] += lambda_0

            for i in range(self._nb_of_clast):   # Treat clasts one after the other
                # Get values from dataframe:
                _x_clast_init = self.at[i, 'clast__x']
                _y_clast_init = self.at[i, 'clast__y']
                _z_clast_init = self.at[i, 'clast__y']
                # Test if clast in grid core (=not phantom)
                if ClastCollection.phantom(self, i) == False:
                    #print('not phantom')
                    # Clast is in grid core
                    # Test if clast is detached:
                    if ClastCollection.clast_detach_proba(self, i) == True:
                        #print('detached')
                        # Clast is detached -> Test if moves (leaves node):
                        ClastCollection._neighborhood(self, i)
                        ClastCollection._move_to(self, i)
                        
                        self.at[clast, 'hop_length'] =0
                        self.rand_length = 0.
                        while ClastCollection._change_cell_proba(self, i) == True:
                            ClastCollection._move_out_of_cell(self,i)
                            ClastCollection._neighborhood(self, i)
                            ClastCollection._move_to(self, i)
                        else:
                            ClastCollection._move_in_cell(self,i)
                            
                            
    for i in range(self._set_size):
        
        ClastCollection._neighborhood(i)   # updates neighborhood info
        ClastCollection._move_to(i)        # dertermines direction of mvt
        if np.isnan(ss_azimuth) is False: # if centered slope is not null (not on a flat)
    
            ClastCollection.move(i)
    
    
    
    
    
    
    
    
    

# add delta a little extra travelled distance  (delta) when displacing 
# the clast so that it does change cell, doesn't stay at boundary?
# but change_x, change_y calculated from dist_to_exit which is always longer
# than the horizontal distance so should be ok
# delta = np.power(10, -3)