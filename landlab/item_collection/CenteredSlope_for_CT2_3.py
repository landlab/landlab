#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import gc
from landlab.item_collection import ItemCollection
from matplotlib.pyplot import figure, show, plot, xlabel, ylabel, title, legend

from math import pi

from scipy.stats import truncnorm

class ClastCollection(ItemCollection):
    """
    ClastCollection is a datastructure inheriting from ItemCollection and
    includes methods applying to the data.

    It creates a collection of clasts that live on grid cells and controls their
    movement on the grid.

    ClastCollection does not change the landscape. It must be used together
    with another Landlab component that calculates erosion and deposition on
    every grid node. Currently, ClastCollection can be used with:
        - TransportLengthHillslopeDiffuser
        - SPACE
        (-ErosionDeposition) (TO DO)

    A flow director must also be used (provides topographic__steepest_slope).

    ClastsCollection uses the erosion calculated on each node to determine
    whether a clast is mobile: a clast can be moved if it is entirely within
    the eroded layer.
    If a clast is mobile, it is assigned a potential travel distance, drawn
    from an exponential probability distribution of travel distances as
    described by Furbish and Haff (2010) and Furbish and Roering (2013).
    The distribution of travel distances depends on the slope, the clast size
    and is conditional (takes into account the fact that the clast has already
    travelled), it is thus calculated for each clast, and on each cell it
    crosses during its travel.
    If the assigned travel distance is less than the distance to exit the
    current cell, the clast is moved within the cell according to the assigned
    travel distance.
    If the assigned travel distance is greater than the distance to exit the
    cell, the clast is moved to the next downstream cell and a new potential
    travel distance is calculated and assigned to the clast.

    Clast collection currently only works on RasterModelGrid.

    Methods
    -------
    clast_detach_proba
    clast_solver_Exponential
    phantom

    """
    def __init__(self,
                 grid,
                 clast_x=[],
                 clast_y=[],
                 clast_elev=[],
                 clast_radius=[]):
        """
        Parameters
        ----------
        grid : RasterModelGrid
        clast_x : number-of-clasts long array
            x coordinates of the initial clast positions
        clast_y : number-of-clasts long array
            y coordinates of the initial clast positions
        clast_elev : number-of-clasts long array
            Initial elevation of the clasts
        clast_radius : number-of-clasts long array
            Initial radius of the clasts

        Examples
        --------




        """

        # Save a reference to the grid:
        self._grid = grid

        ### Coordinates of cell corners surrounding each node:
        self.x_of_corners_at_node=np.zeros((self._grid.number_of_nodes,4))
        self.y_of_corners_at_node=np.zeros((self._grid.number_of_nodes,4))

        for i in range(self._grid.number_of_nodes):
            self.x_of_corners_at_node[i,:]=[
                    self._grid.node_x[i]+self._grid.dx/2,
                    self._grid.node_x[i]-self._grid.dx/2,
                    self._grid.node_x[i]-self._grid.dx/2,
                    self._grid.node_x[i]+self._grid.dx/2]
            self.y_of_corners_at_node[i,:]=[
                    self._grid.node_y[i]+self._grid.dx/2,
                    self._grid.node_x[i]+self._grid.dx/2,
                    self._grid.node_x[i]-self._grid.dx/2,
                    self._grid.node_x[i]-self._grid.dx/2]


        # Determine reference (closest) node for each clast:
        _clast__node=[]
        _clast__node[:] = self._grid.find_nearest_node((clast_x[:], clast_y[:]))

        # Determine current cell for each clast:
        _clast__cell=[]
        _clast__cell[:] = self._grid.cell_at_node[_clast__node]

        # Clast collection size:
        self.nb_of_clast = len(_clast__node)

        # Store the input information in a dictionary and create the other
        # fields that will be populated later:
        clast_data = {'clast__x' : clast_x,
                      'clast__y' : clast_y,
                      'clast__elev' : clast_elev,
                      'clast__node' : _clast__node,
                      'clast__radius' : clast_radius,
                      'lambda_0' : np.zeros(self._nb_of_clast),
                      'lambda_mean' : np.zeros(self._nb_of_clast),
                      'slope__WE' : np.zeros(self._nb_of_clast),
                      'slope__SN' : np.zeros(self._nb_of_clast),
                      'slope__steepest_azimuth' : np.full(
                              self._nb_of_clast, np.NaN),
                      'slope__steepest_dip' : np.zeros(self._nb_of_clast),
                      'distance__to_exit' : np.full(self._nb_of_clast, np.NaN),
                      'target_node' : -np.ones(self._nb_of_clast, dtype=int),
                      'target_node_flag': -np.ones(
                              self._nb_of_clast, dtype=int),
                      'change_x' : np.zeros(self._nb_of_clast),
                      'change_y' : np.zeros(self._nb_of_clast),
                      'hop_length' : np.zeros(self._nb_of_clast),
                      'total_travelled_dist' : np.zeros(self._nb_of_clast)}

        # Build ItemCollection containing clast data:
        super(ClastCollection, self).__init__(
                                self._grid,
                                data=clast_data,
                                grid_element='cell',
                                element_id=_clast__cell)

#        ItemCollection.__init__(self,...)

#        super().__init__(
#                self,
#                self._grid,
#                data=clast_data,
#                grid_element='node',
#                element_id=_clast__node)

    def _neighborhood(self, clast):
        """ Determine the neighboring nodes of each clast, the node towards
        whichthe clast will be moved (target node) and the slope used to
        move the clast.

        If the clast is on a hillslope celle and none of the neighbor nodes
        are boundary nodes, calculate the (signed) centered slopes around
        clast's reference node: west-east slope is elevation difference
        between west neighbor and east neighbor, and north-south slope is
        elevation difference between north neighbor and south neighbor.

        If the clast is on a river cell and/or if at least one of the neighbor
        node is a boundary node, the slope used to move the clast is provided
        by Flow director: the slope is either the provided
        topographic__steepest_slope or 0 if clast's node is sink.
        """
        df = self.DataFrame

        # If clast is phantom, slope is 0:
        if ClastCollection.phantom(self, clast) == True:
            df.at[clast, 'slope__WE'] = 0
            df.at[clast, 'slope__SN'] = 0

        else:
            # Determine reference (closest) node for each clast:
            _grid = self._grid
            df.at[clast, 'clast__node'] = _grid.find_nearest_node(
                    (df.at[clast, 'clast__x'],
                         df.at[clast, 'clast__y']))
            _node = df.at[clast, 'clast__node']

            # Adjacent row and col nodes:
            _row_col_adjacent_nodes_at_node = (
                    _grid.neighbors_at_node[_node])
            east_node = _row_col_adjacent_nodes_at_node[0]
            north_node = _row_col_adjacent_nodes_at_node[1]
            west_node = _row_col_adjacent_nodes_at_node[2]
            south_node = _row_col_adjacent_nodes_at_node[3]

            # Diagonal nodes:
            _diagonal_adjacent_nodes_at_node = (
                    _grid.diagonal_adjacent_nodes_at_node[_node])

            # Full neighborhood: E, N, W, S, NE, NW, SW, SE
            _neighbor_nodes = np.concatenate((
                    _row_col_adjacent_nodes_at_node,
                        _diagonal_adjacent_nodes_at_node),
                            axis=0)

            # Case where one of the neighbor is boundary:
            if any(i in _row_col_adjacent_nodes_at_node for\
                       i in _grid.boundary_nodes) == True or \
                           ClastCollection._cell_is_hillslope(self, clast) == \
                               False:
                print('close to boundary and/or river')
                if _grid.at_node['flow__receiver_node'][_node] != _node:
                    # if FlowDirector designates a receiver node other than
                    # the node itself, clast is moving toward receiver node:
                    df.at[clast, 'target_node'] = (
                            _grid.at_node['flow__receiver_node'][_node])
                    df.at[clast, 'target_node_flag'] = (
                            np.where(_neighbor_nodes ==
                                     _grid.at_node['flow__receiver_node']\
                                         [_node])[0])
                    df.at[clast, 'slope__WE'] = np.NaN
                    df.at[clast, 'slope__SN'] = np.NaN
                else: # flow receiver = node itself
                    df.at[clast, 'target_node'] = _node
                    df.at[clast, 'target_node_flag'] = -1
                    df.at[clast, 'slope__WE'] = 0.
                    df.at[clast, 'slope__SN'] = 0.

            else: # if not close to boundary

                # Calculation of slopes: W to E and S to N, units=m/m
                df.at[clast, 'slope__WE'] = ((
                        _grid.at_node['topographic__elevation'][west_node]-(
                                _grid.at_node['topographic__elevation']\
                                    [east_node]))/(2*_grid.dx))

                df.at[clast, 'slope__SN'] = \
                    (_grid.at_node['topographic__elevation'][south_node]-
                          _grid.at_node['topographic__elevation']\
                              [north_node])/(2*_grid.dy)



    def _move_to(self, clast):
        """ Determine the direction and value of slope where the clast will
        travel.

        Depending on the sign of north-south and west-east slopes

        """

        _grid = self._grid
        df = self.DataFrame
        _node = df.at[clast, 'clast__node']
        _node_x = _grid.node_x[_node]
        _node_y = _grid.node_y[_node]
        _clast_x = df.at[clast, 'clast__x']
        _clast_y = df.at[clast, 'clast__y']
        _node_z = _grid.at_node['topographic__elevation'][_node]
        we_slope = df.at[clast, 'slope__WE']
        sn_slope = df.at[clast, 'slope__SN']
        _radius = df.at[clast, 'clast__radius']
        # Adjacent row and col nodes:
        _row_col_adjacent_nodes_at_node = _grid.neighbors_at_node[_node]
        # Adjacent diagonal nodes:
        _diagonal_adjacent_nodes_at_node = _grid.diagonal_adjacent_nodes_at_node[_node]

        # If clast is phantom, it will not move, set everything to 0:
        if ClastCollection.phantom(self, clast) == True:
            target_node = df.at[clast, 'clast__node']
            ss_azimuth = np.NaN
            ss_dip = 0.
            dist_to_exit = np.NaN
            [change_x, change_y] = [0., 0.]

        # Clast is not phantom:
        # Cases where one of the neighbor is boundary node or clast is on a
        # river cell: use the topographis__steepest_slope and
        # flow__receiver_node provided by flow director:
        elif (np.isnan(df.at[clast, 'slope__WE']) == True and np.isnan(df.at[clast, 'slope__SN']) == True) or ClastCollection._cell_is_hillslope(self, clast) == False:
            print('clast is next to boundary or on river cell: calculating slope to receiver')
            df.at[clast, 'target_node'] = self._grid.at_node['flow__receiver_node'][_node]
            target_node = df.at[clast, 'target_node']
            print('receiver is : %s' %target_node)
            target_node_flag = df.at[clast, 'target_node_flag']


            if target_node_flag == 0: # East
                we_slope = (_node_z - self._grid.at_node['topographic__elevation'][df.at[clast, 'target_node']]) / _grid.dx
                sn_slope = 0.
                ss_azimuth = 0. # east

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * ((_node_x + (self._grid.dx/2)) - _clast_x)
                [change_x, change_y] = [dist_to_exit + _radius, 0.0]
################## FOR TESTING PURPOSE ONLY ##################
                if df.at[clast, 'slope__WE'] < 0:
                    print('error')
                else:
                    pass
###############################################################
            elif target_node_flag == 1: # North
                we_slope = 0.
                sn_slope = (_node_z - self._grid.at_node['topographic__elevation'][df.at[clast, 'target_node']]) / _grid.dy
                ss_azimuth = np.radians(90) # north

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * ((_node_y + (self._grid.dy/2)) - _clast_y)
                [change_x, change_y] = [0.0, dist_to_exit + _radius]
            elif target_node_flag == 2: # West
                we_slope = (self._grid.at_node['topographic__elevation'][df.at[clast, 'target_node']] - _node_z) / _grid.dx
                sn_slope = 0.
                ss_azimuth = np.radians(180) # west

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * (_clast_x - (_node_x - (self._grid.dx/2)))
                [change_x, change_y] = [-dist_to_exit - _radius, 0.0]
            elif target_node_flag == 3: # South
                we_slope = 0.
                sn_slope = (self._grid.at_node['topographic__elevation'][df.at[clast, 'target_node']] - _node_z) / _grid.dy
                ss_azimuth = np.radians(270) # south

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * (_clast_y - (_node_y - (self._grid.dy/2)))
                [change_x, change_y] = [0.0, -dist_to_exit - _radius]
            else: # Diagonals
                we_slope = np.NaN
                sn_slope = np.NaN
                _target_node_x = _grid.node_x[df.at[clast, 'target_node']]
                _target_node_y = _grid.node_y[df.at[clast, 'target_node']]
                _target_node_z = _grid.at_node['topographic__elevation'][df.at[clast, 'target_node']]

                ss_dip = np.arctan(_node_z - _target_node_z) / np.sqrt(np.power(_grid.node_x[target_node]-_grid.node_x[_node],2)+np.power(_grid.node_y[target_node]-_grid.node_y[_node],2))
                dist_to_exit = (1 / np.cos(ss_dip)) * np.sqrt(np.power((_target_node_x - _clast_x), 2) + np.power((_target_node_y - _clast_y), 2))

                if target_node_flag == 4: # NE
                    corner = 0
                    ss_azimuth = np.arctan((self.y_of_corners_at_node[_node, corner] - _clast_y) / (self.x_of_corners_at_node[_node, corner] - _clast_x))
                    [change_x, change_y] = [self.x_of_corners_at_node[_node, corner] - _clast_x +  _radius, self.y_of_corners_at_node[_node, corner] - _clast_y + _radius]

                elif target_node_flag == 5: # NW
                    corner = 1
                    ss_azimuth = np.radians(90) + np.arctan((_clast_x - self.x_of_corners_at_node[_node, corner]) / (self.y_of_corners_at_node[_node, corner] - _clast_y))
                    [change_x, change_y] = [-(_clast_x - self.x_of_corners_at_node[_node, corner]) - _radius, self.y_of_corners_at_node[_node, corner] - _clast_y +  _radius]

                elif target_node_flag == 6: # SW
                    corner = 2
                    ss_azimuth = np.radians(180) + np.arctan((_clast_y - self.y_of_corners_at_node[_node, corner]) / (_clast_x - self.x_of_corners_at_node[_node, corner]))
                    [change_x, change_y] = [-(_clast_x - self.x_of_corners_at_node[_node, corner]) - _radius, -(_clast_y - self.y_of_corners_at_node[_node, corner]) - _radius]

                elif target_node_flag == 7: # SE
                    corner = 3
                    ss_azimuth = np.radians(270) + np.arctan((self.x_of_corners_at_node[_node, corner] - _clast_x) / (_clast_y - self.y_of_corners_at_node[_node, corner]))
                    [change_x, change_y] = [(self.x_of_corners_at_node[_node, corner] - _clast_x) + _radius, -(_clast_y - self.y_of_corners_at_node[_node, corner]) - _radius]
                else:
                    print('error target node flag?')

#            print('we=%s' %we_slope)
#            print('sn=%s' %sn_slope)
#            print('SS_HORIZ_NORM=%s' %ss_horiz_norm)
#            print('SteepSlope_norm=%s' %ss_norm)
#            print('ss_dip=%s' %ss_dip)
#            print('SN_norm=%s' %sn_norm)
#            print('WE_norm=%s' %we_norm)


            df.at[clast, 'slope__WE'] = we_slope
            df.at[clast, 'slope__SN'] = sn_slope
            df.at[clast, 'slope__steepest_azimuth'] = ss_azimuth
            df.at[clast, 'slope__steepest_dip'] = ss_dip
            df.at[clast, 'distance__to_exit'] = dist_to_exit
            df.at[clast, 'target_node'] = target_node
            #df.at[clast, 'target_node_flag'] = target_node_flag



        else: # clast is not next to boundary
            # norm of steepest slope vector projected on horizontal plane:
            ss_horiz_norm = np.sqrt(np.power(we_slope, 2) + np.power(sn_slope, 2))

            # norms of vectors SN and WE:
            sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
            we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))

            # norm of steepest slope vector = norm of resultant of SN and WE:
            ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))

            if we_slope == 0:
                if sn_slope == 0:
                    print('slope is null')
                    # centered slope is null, clast does not move:
                    ss_dip = 0
                    ss_azimuth = None
                    dist_to_exit = -1
                    target_node = _node
                    [change_x, change_y] = [0.0, 0.0]
        #            # OPTION2 (to develop):
        #            # clast moves to random direction, according to lambda_0:
        #            ss_azimuth = np.random.uniform(0.0, 2*pi, 1)   # pick a direction at random
        #            df.at[clast, 'slope__steepest_azimuth'] = ss_azimuth
        #            df.at[clast, 'slope__flag'] = '
        #            self.rand_length = np.random.exponential(scale=df.at[clast, 'lambda_0'], 1)
        #            dist_to_exit =


                else: # SN slope is not 0, ss direction is S or N
                    # dip of steepest slope:
                    ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                    if sn_slope < 0: # ss direction is South
                        print('going South')
                        if ss_dip != np.arctan(np.abs(sn_slope)):   # dip in radians
                            print('error, dip is %s' %ss_dip)
                            print('should be:')
                            print(np.arctan(np.abs(sn_slope)))

                        ss_azimuth = np.radians(270) # South
                        dist_to_exit = (1 / np.cos(ss_dip)) * (_clast_y - (_node_y - (self._grid.dy/2))) #
                        target_node = _row_col_adjacent_nodes_at_node[3]
                        [change_x, change_y] = [0.0, -dist_to_exit - _radius]
                    else: # ss direction is North
                        print('going North')
                        if ss_dip != np.arctan(np.abs(sn_slope)):
                            print('error, dip is %s' %ss_dip)

                        ss_azimuth = np.radians(90) # North
                        dist_to_exit = (1 / np.cos(ss_dip)) * ((_node_y + (self._grid.dy/2)) - _clast_y) #
                        target_node = _row_col_adjacent_nodes_at_node[1]
                        [change_x, change_y] = [0.0, dist_to_exit + _radius]

            else: # we_slope is not 0
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                if sn_slope == 0:  # ss is W or E
                    if we_slope < 0: # ss direction is West
                        print('West')
                        ss_dip = np.arctan(np.abs(we_slope))
                        ss_azimuth = np.radians(180) # West
                        dist_to_exit = (1 / np.cos(ss_dip)) * (_clast_x - (_node_x - (self._grid.dx/2))) #
                        target_node = _row_col_adjacent_nodes_at_node[2]
                        [change_x, change_y] = [-dist_to_exit - _radius, 0.0]
                    else: # ss direction is East
                        print('East')
                        ss_dip = np.arctan(np.abs(we_slope))
                        ss_azimuth = 0 # East
                        dist_to_exit = (1 / np.cos(ss_dip)) * ((_node_x + (self._grid.dx/2)) - _clast_x) #
                        target_node = _row_col_adjacent_nodes_at_node[0]
                        [change_x, change_y] = [dist_to_exit + _radius, 0.0]

                else: # sn_slope is not 0
                    if sn_slope > 0 and we_slope > 0: # Quarter = NE
                        ss_azimuth = np.arctan(np.abs(sn_slope / we_slope))
                        corner = 0
                        clast_to_corner_azimuth = np.arctan(
                                np.abs(_clast_y - self.y_of_corners_at_node[_node, corner])/
                                       np.abs(_clast_x - self.x_of_corners_at_node[_node, corner]))
                        if ss_azimuth < clast_to_corner_azimuth: # Eigth = NE-row
                            print('NE-R')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x)) / (np.cos(ss_azimuth))) #
                            target_node = _row_col_adjacent_nodes_at_node[0]
                            [change_x, change_y] = [(dist_to_exit * abs(we_slope) / ss_horiz_norm) + _radius, (dist_to_exit * abs(sn_slope) / ss_horiz_norm) + _radius]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = NE-col
                            print('NE-C')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(np.radians(90) - ss_azimuth))) #
                            target_node = _row_col_adjacent_nodes_at_node[1]
                            [change_x, change_y] = [(dist_to_exit * abs(we_slope) / ss_horiz_norm) + _radius, (dist_to_exit * abs(sn_slope) / ss_horiz_norm) + _radius]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('NE-corner')
                            dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            target_node = _diagonal_adjacent_nodes_at_node[0]
                            [change_x, change_y] = [(dist_to_exit * abs(we_slope) / ss_horiz_norm) + _radius, (dist_to_exit * abs(sn_slope) / ss_horiz_norm) + _radius]
                    elif sn_slope > 0 and we_slope < 0: # Quarter = NW
                        ss_azimuth = np.radians(90) + np.arctan(np.abs(we_slope / sn_slope))
                        corner = 1
                        clast_to_corner_azimuth = np.radians(180) - np.arctan(
                                np.abs(_clast_y - self.y_of_corners_at_node[_node, corner])/
                                       np.abs(_clast_x - self.x_of_corners_at_node[_node, corner]))
                        if ss_azimuth < clast_to_corner_azimuth: # Eigth = NW-col
                            print('NW-C')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(ss_azimuth - np.radians(90)))) #
                            [change_x, change_y] = [-(dist_to_exit * abs(we_slope) / ss_horiz_norm) - _radius, (dist_to_exit * abs(sn_slope) / ss_horiz_norm) + _radius]
                            target_node = _row_col_adjacent_nodes_at_node[1]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = NW-row
                            print('NW-R')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x)) / (np.cos(np.radians(180) - ss_azimuth))) #
                            [change_x, change_y] = [-(dist_to_exit * abs(we_slope) / ss_horiz_norm) - _radius, (dist_to_exit * abs(sn_slope) / ss_horiz_norm) + _radius]
                            target_node = _row_col_adjacent_nodes_at_node[2]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('NW-corner')
                            dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            [change_x, change_y] = [-(dist_to_exit * abs(we_slope) / ss_horiz_norm) - _radius, (dist_to_exit * abs(sn_slope) / ss_horiz_norm) + _radius]
                            target_node = _diagonal_adjacent_nodes_at_node[1]

                    elif sn_slope < 0 and we_slope < 0: # Quarter = SW
                        ss_azimuth = np.radians(180) + np.arctan(np.abs(sn_slope / we_slope))
                        corner = 2
                        clast_to_corner_azimuth = np.radians(180) + np.arctan(
                                np.abs(_clast_y - self.y_of_corners_at_node[_node, corner])/
                                       np.abs(_clast_x - self.x_of_corners_at_node[_node, corner]))
                        if ss_azimuth < clast_to_corner_azimuth: # Eigth = SW-row
                            print('SW-R')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x)) / (np.cos(ss_azimuth - np.radians(180)))) #
                            [change_x, change_y] = [-(dist_to_exit * abs(we_slope) / ss_horiz_norm) - _radius, -(dist_to_exit * abs(sn_slope) / ss_horiz_norm) - _radius]
                            target_node = _row_col_adjacent_nodes_at_node[2]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = SW-col
                            print('SW-C')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(np.radians(270) - ss_azimuth)))  #
                            [change_x, change_y] = [-(dist_to_exit * abs(we_slope) / ss_horiz_norm) - _radius, -(dist_to_exit * abs(sn_slope) / ss_horiz_norm) - _radius]
                            target_node = _row_col_adjacent_nodes_at_node[3]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('SW-corner')
                            dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            [change_x, change_y] = [-(dist_to_exit * abs(we_slope) / ss_horiz_norm) - _radius, -(dist_to_exit * abs(sn_slope) / ss_horiz_norm) - _radius]
                            target_node = _diagonal_adjacent_nodes_at_node[2]

                    elif sn_slope < 0 and we_slope > 0: # Quarter = SE
                        ss_azimuth = np.radians(270) + np.arctan(np.abs(we_slope / sn_slope))
                        corner = 3
                        clast_to_corner_azimuth = np.radians(360) - np.arctan(
                                np.abs(_clast_y - self.y_of_corners_at_node[_node, corner])/
                                       np.abs(_clast_x - self.x_of_corners_at_node[_node, corner]))
                        if ss_azimuth < clast_to_corner_azimuth: # Eigth = SE-col
                            print('SE-C')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(ss_azimuth - np.radians(270)))) #
                            [change_x, change_y] = [(dist_to_exit * abs(we_slope) / ss_horiz_norm) + _radius, -(dist_to_exit * abs(sn_slope) / ss_horiz_norm) - _radius]
                            target_node = _row_col_adjacent_nodes_at_node[3]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = SE-row
                            print('SE-R')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x) / (np.cos(np.radians(360) - ss_azimuth))))  #
                            [change_x, change_y] = [(dist_to_exit * abs(we_slope) / ss_horiz_norm) + _radius, -(dist_to_exit * abs(sn_slope) / ss_horiz_norm) - _radius]
                            target_node = _row_col_adjacent_nodes_at_node[0]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('SE-corner')
                            dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            [change_x, change_y] = [(dist_to_exit * abs(we_slope) / ss_horiz_norm) + _radius, -(dist_to_exit * abs(sn_slope) / ss_horiz_norm) - _radius]
                            target_node = _diagonal_adjacent_nodes_at_node[3]

        # Save values to dataframe:
        df.at[clast, 'slope__steepest_azimuth'] = ss_azimuth
        df.at[clast, 'slope__steepest_dip'] = ss_dip
        df.at[clast, 'distance__to_exit'] = dist_to_exit
        df.at[clast, 'target_node'] = target_node


        # Calculate lambda_mean:
        df.at[clast, 'lambda_0'] = ((self._dt * self._kappa * _grid.dx) / (self._Si *  2 * df.at[clast,'clast__radius']))
        #((self._dt * self._kappa * _grid.dx) / (np.tan(df.at[clast, 'slope__steepest_dip']) *  2 * df.at[clast,'clast__radius'])) #* (np.tan(df.at[clast, 'slope__steepest_dip']) * ((self._Si-np.tan(df.at[clast, 'slope__steepest_dip']))/(self._Si+np.tan(df.at[clast, 'slope__steepest_dip']))))  # self._Si * 20 *
        #df.at[clast, 'lambda_0'] = _grid.at_node['sediment__deposition_coeff'][_node] * _grid.dx
        #print(df.at[clast, 'lambda_0'])

        if ClastCollection._cell_is_hillslope(self, clast) == True:   # Cell is hillslope
            if np.tan(df.at[clast, 'slope__steepest_dip']) >= self._Si:
                lambda_0 = df.at[clast, 'lambda_0']
                lambda_mean = np.power(10, 10)
            else:
                lambda_0  = df.at[clast,'lambda_0']  # 1
                lambda_mean = lambda_0 * (self._Si + np.tan(df.at[clast, 'slope__steepest_dip'])) / (self._Si - np.tan(df.at[clast, 'slope__steepest_dip']))
        else:   # Cell is river
            if np.tan(df.at[clast, 'slope__steepest_dip'] >= self._Si):
                lambda_0 = df.at[clast, 'lambda_0']
                lambda_mean = np.power(10, 10)
            else:
                lambda_0 = 10 * df.at[clast,'lambda_0'] # 100
                lambda_mean = lambda_0 * (self._Si + np.tan(df.at[clast, 'slope__steepest_dip'])) / (self._Si - np.tan(df.at[clast, 'slope__steepest_dip']))


        df.at[clast, 'lambda_mean'] = lambda_mean
        df.at[clast, 'lambda_0'] = lambda_0


        print('change_x = %s' %change_x)
        print('change_y = %s' %change_y)

        x=truncnorm(-(df.at[clast, 'lambda_0']/2), (df.at[clast, 'lambda_0']/2), loc=0, scale=(df.at[clast, 'lambda_0']/2)/4)
        [dx,dy] = x.rvs(2)/10
        df.at[clast, 'change_x'] = change_x #+ dx  #np.random.normal(loc=0.0, scale=(df.at[clast, 'lambda_mean']/_grid.dx), size=1)   # adds a randomness for lateral spreading
        df.at[clast, 'change_y'] = change_y #+ dy  #np.random.normal(loc=0.0, scale=(df.at[clast, 'lambda_mean']/_grid.dx), size=1)   # adds a randomness for lateral spreading



    def _change_cell_proba(self, clast):
        lambda_mean = self.df.at[clast, 'lambda_mean']
        dist_to_exit = self.df.at[clast, 'distance__to_exit']

        if dist_to_exit == -1: # case where slope is null
            _change_cell = False
        else:
            # Draw a random sample in the probability distribution of travel distances:
            self.rand_length = np.random.exponential(scale=lambda_mean, size=1)

            if self.rand_length < dist_to_exit: # clast stays in cell
                _change_cell = False
            else: # self.rand_length >= dist_to_exit: clast leaves cell
                _change_cell = True

        print('rand_length = %s' % self.rand_length)
        print('dist_to_exit= %s' % dist_to_exit)

        return _change_cell


    def _move_in_cell(self, clast):
        # clast stays in cell, move of distance rand_length along slope
        ss_azimuth = self.df.at[clast, 'slope__steepest_azimuth']
        ss_dip = self.df.at[clast, 'slope__steepest_dip']
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
        self.df.at[clast, 'clast__x'] += change_x
        self.df.at[clast, 'clast__y'] += change_y
        #self.df.at[clast, 'change_x'] = change_x
        #self.df.at[clast, 'change_y'] = change_y


        self.df.at[clast, 'hop_length'] += self.rand_length

    def _move_out_of_cell(self, clast):
        # clast leaves cell, move of distance dist_to_exit along slope
        self.df.at[clast, 'clast__x'] += self.df.at[clast, 'change_x']
        self.df.at[clast, 'clast__y'] += self.df.at[clast, 'change_y']
        self.df.at[clast, 'hop_length'] += (np.sqrt(np.power(self.df.at[clast, 'change_x'], 2)+ np.power(self.df.at[clast, 'change_y'],2))) / np.cos(self.df.at[clast, 'slope__steepest_dip']) # self.df.at[clast, 'distance__to_exit']
        self.df.at[clast, 'clast__node'] = self.df.at[clast, 'target_node']
        self.df.at[clast, 'element_id'] = self._grid.cell_at_node[self.df.at[clast, 'clast__node']]


    def phantom(self, clast):
        """Check if the clast has exited the grid.

        When a clast reaches a boundary node, it exits the grid and is thus
        flagged as phantom. It is not erased from the dataframe but is not
        moved anymore.

        Parameters
        ----------
        clast : int
            Clast ID.

        """
        # To add: Also phantom when totally dissovled (radius = 0)
        clast_node = self.df.at[clast, 'clast__node']
        clast_radius = self.df.at[clast, 'clast__radius']
        boundary_nodes = self._grid.boundary_nodes
        _phantom = np.zeros(1, dtype=bool)

        if clast_node in boundary_nodes:
            _phantom = True
        elif clast_radius == 0.:
            _phantom = True
        else:
            _phantom = False

        return _phantom

    def clast_detach_proba(self, clast):
        """ Test if clast is detached (mobile): clast is detached if erosion
        is sufficient to expose its base.

        Parameters
        ----------
        clast : int
            Clast ID

        """

        clast__node = self.df.at[clast, 'clast__node']
        clast__elev = self.df.at[clast, 'clast__elev']
        topo__elev = self._grid.at_node['topographic__elevation'][clast__node]
        erosion = self._erosion__depth[clast__node]

        _detach = np.zeros(1, dtype=bool)

        if erosion >= topo__elev - clast__elev:
            _detach = True
        else:
            _detach = False

        return _detach

    def _cell_is_hillslope(self, clast):
        """Test if cell currently occupied by clast is hillslope (or river).

        The area/slope ratio can be used to determine if a cell belongs to the
        hillslope or river domain.

        See Tucker and Bras, 1998, Water Resources Research, for complete
        justification.

        Parameters
        ----------
        clast : int
            Clast ID

        """

        # Should threshold be an input parameter (user-defined)?

        _node = self.df.at[clast, 'clast__node']
        _grid = self._grid
        threshold = 1e4

        if np.isnan(self.df.at[clast, 'slope__steepest_azimuth']) == True:
            area_over_slope = 1   # clast next to boundary
        else:
            area_over_slope = _grid.at_node['drainage_area'][_node] / _grid.at_node['topographic__steepest_slope'][_node]

        if area_over_slope < threshold:
            _cell_is_hillslope = True
        else:
            _cell_is_hillslope = False

        return _cell_is_hillslope

    def clast_solver_Exponential(self, dt=1., Si=1.2, kappa=0.0001, uplift=None, erosion_method='TLDiff'): # lambda_0=1,

        self.df=self.DataFrame
        # Method loop: Depending on the method used to evolve the landscape,
        # get sediment influx, erosion and deposition
        self.erosion_method = erosion_method
        if self.erosion_method == 'TLDiff':
            self._erosion_rate = self._grid.at_node['sediment__erosion_rate']
            self._deposition_rate = self._grid.at_node['sediment__deposition_rate']
            self._sediment__flux_in = self._grid.at_node['sediment__flux_in']
        elif self.erosion_method == 'Space':
            from landlab.components import Space
            for obj in gc.get_objects():
                if isinstance(obj, Space):    # look for instance of Space
                    self._erosion_rate = obj.Es + obj.Er   # Works if only one instance of Space was made
                    self._deposition_rate = obj.depo_rate
                    self._sediment__flux_in = obj.qs_in
            # self.space = space_name CAN'T CALL INSTANCE FROM INPUT STRING NAME

        # Future version: multiple compo -> add fluxes?
        # Store various values that will be used
        self._kappa = kappa
        self._dt = dt
        self._erosion__depth = self._erosion_rate * self._dt
        self._deposition__thickness = self._deposition_rate * self._dt
        self._Si = Si #  slope above which particle motion continues indefinetly (not equal to critical slope of diffusion, see Furbish and Haff, 2010)
#            self._lambda_0 = lambda_0

#        self._lambda_0=np.zeros(self._nb_of_clast)
#        for i in range(self._nb_of_clast-1):
#            self._lambda_0[i] = (self._dt * self._kappa * max(self._grid.dx, self._grid.dy)) / (self._Si * self.df.at[i, 'clast__radius'])
##            self._lambda_0[:] += lambda_0

        # Uplift:
        if uplift is not None:
            if type(uplift) is str:
                self._uplift = self._grid.at_node[uplift]
            elif type(uplift) in (float, int):
                self._uplift = np.ones(self._grid.number_of_nodes) * (
                        float(uplift))
            elif len(uplift) == self._grid.number_of_nodes:
                self._uplift = np.array(uplift)
            else:
                raise TypeError('Supplied type of uplift is not recognized')

        for clast in range(self._nb_of_clast):   # Treat clasts one after the other
            print('CLAST %s' %clast)
            # Get values from dataframe:
#                _x_clast_init = self.df.at[i, 'clast__x']
#                _y_clast_init = self.df.at[i, 'clast__y']
#                _z_clast_init = self.df.at[i, 'clast__y']
            # Test if clast in grid core (=not phantom)
            if ClastCollection.phantom(self, clast) == False:
                print('not phantom')
                # Clast is in grid core
                # Test if clast is detached:
                if ClastCollection.clast_detach_proba(self, clast) == True:
                    print('detached')
                    # Clast is detached -> update neighborhood info:
                    ClastCollection._neighborhood(self, clast)
                    print('neighborhood')
                    ClastCollection._move_to(self, clast)
                    print('move_to')

                    self.df.at[clast, 'hop_length'] =0
                    self.rand_length = 0.
                    #Test if moves (leaves node):
                    if np.isnan(self.df.at[clast,'slope__steepest_azimuth']) == False: # if centered slope is not null (not on a flat)
                        print('not on flat')
                        while ClastCollection._change_cell_proba(self, clast) == True:
                            if ClastCollection.phantom(self, clast) == False:
                                print('HERE')
                                print(self.df)
                                print('must change cell')
                                ClastCollection._move_out_of_cell(self,clast)
                                print('moved out')
                                print('clastx= %s' %self.DataFrame.at[clast, 'clast__x'])
                                print('clasty= %s' %self.DataFrame.at[clast, 'clast__y'])
                                figure(1)
                                plot(self.DataFrame.at[clast, 'clast__x'], self.DataFrame.at[clast, 'clast__y'], 'o', color='gray', markersize=1)
                                ########JUST FOR TESTING PURPOSE##################################
                                if self.df.at[clast, 'clast__node'] != self.df.at[clast, 'target_node']:
                                    print('Error: moved to wrong node')
                                ##################################################################
                                self.df.at[clast, 'target_node_flag'] = -1
                                ClastCollection._neighborhood(self, clast)
                                print('done new neighborhood')
                                ClastCollection._move_to(self, clast)
                                print('done new move_to, test change_cell_proba')
                            else:
                                print('clast has gone out of grid')
                                break


                        else:
                            if ClastCollection.phantom(self, clast) == False:
                                if self.df.at[clast,'distance__to_exit'] == -1:# case where slope is null
                                    pass
                                else:
                                    print(self.df)
                                    print('must move in cell')
                                    ClastCollection._move_in_cell(self,clast)
                                    print('moved in cell')
                            else:
                                break
                    else: # if centered slope is null
                        pass   # go to next clast

                # Update elevation:
                self.df.at[clast, 'clast__elev'] = self._grid.at_node['topographic__elevation'][self.df.at[clast, 'clast__node']] - (self._deposition__thickness[self.df.at[clast, 'clast__node']] * np.random.rand(1))
                print('elevation updated')
                # Update total travelled distance:
                self.df.at[clast, 'total_travelled_dist'] += self.df.at[clast, 'hop_length']
                if hasattr(self, '_uplift') is True:   # uplift clast if necessary
                    self.df.at[clast, 'clast__elev'] += self._uplift[self.df.at[clast, 'clast__node']] * self._dt
                else:
                    print('not detached')
                    pass
            else: # clast is phantom, go to next clast
                # for display purpose: phantom clast has a radius of 0
                # self._clast__radius(i) = 0.
                pass








# add delta a little extra travelled distance  (delta) when displacing
# the clast so that it does change cell, doesn't stay at boundary?
# but change_x, change_y calculated from dist_to_exit which is always longer
# than the horizontal distance so should be ok
# delta = np.power(10, -3)




            #                 distances = _grid.all_node_distances_map
