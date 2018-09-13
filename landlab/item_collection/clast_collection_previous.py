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
    including methods to move the clasts.

    It creates a collection of clasts that live on grid cells and controls
    their movement on the grid.

    ClastCollection does not change the landscape. It must be used together
    with another Landlab component that computes erosion and deposition on
    every grid node. Currently, ClastCollection can be used with:
        - TransportLengthHillslopeDiffuser
        - SPACE
        (-ErosionDeposition) (TO DO)

    A flow director must also be used (to provide flow__receiver_node,
    and topographic__steepest_slope).

    ClastCollection uses the erosion calculated on each node to determine
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
    travel distance, along the steepest slope direction.
    If the assigned travel distance is greater than the distance to exit the
    cell, the clast is moved to the next downstream cell and a new potential
    travel distance is calculated and assigned to the clast.

    Clast collection currently only works on RasterModelGrid.

    As per ItemCollection, initializing ClastCollection creates a Pandas
    object and clast data is stored as a Dataframe. In addition to entries
    created through the ItemCollection inheritance (i.e. grid_element and
    element_id), the ClastCollection dataframe contains:
         clast__x = x coordinate of clast
         clast__y = y coordinate of clast
         clast__elev = clast elevation
         clast__node = clast reference (closest) node
         clast__radius = radius (m)
         lambda_0 = characteristic length scale
         lambda_mean = mean travel distance
         slope__WE = west-east slope at the current position of the clast
         slope__SN = south-north slope at the current position of the clast
         slope__steepest_azimuth = azimuth (anticlockwise, from East) of
                                                           the steepest slope
         slope__steepest_dip = angle (rads, positive) of the steepest slope
         distance__to_exit = distance from the clast's position to the cell
                   face or corner through which the clast will exit the cell,
                   along the steepest slope (m)
         target_node = node of the cell towards which the clast is moving
         target_node_flag = digit representing the position of the target node
               relative to the current reference node : 0=E, 1=N, 2=W, 3=S,
               4=NE, 5=NW, 6=SW, 7=SE
        # REMOVE NEXT TWO?? not used ##########################################
        # change_x = change in the x coordinate of the clast due to travel step
        # change_y = change in the y coordinate of the clast due to travel step
        #######################################################################
         hop_length = distance travelled in the current time step
         total_travelled_dist = total distance that the clast has travelled

    Methods
    -------
    clast_detach_proba
    clast_solver_Exponential
    phantom

    """

# TO DO #######################################################################
    _name = 'ClastCollection'

    _cite_as = """Insert a BibTeX formatted reference here"""

    _input_var_names = (
        'name',
    )

    _output_var_names = (
        'name',
    )

    _var_units = {
        'name' : 'unit',
    }

    _var_mapping = {
        'name' : 'location',
    }

    _var_doc = {
        'name':
            'description',
    }
###############################################################################


    def __init__(self,
                 grid,
                 clast_x=[],
                 clast_y=[],
                 clast_elev=[],
                 clast_radius=[]):
        """ Initialize ClastCollection.

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

        >>> import numpy as np
        >>> from landlab.item_collection import ClastCollection
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowRouter,
        ...                                TransportLengthHillslopeDiffuser

        Define grid and initial topography:
        >>> mg = RasterModelGrid(5,5)
        >>> z = mg.node_y*0.2
        >>> _ = mg.add_field('node', 'topographic__elevation', z)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
        >>> mg.set_watershed_boundary_condition_outlet_id(2,
        ...                                               mg['node']['topographic__elevation'],
        ...                                               -9999.)

        Instantiate flow router and hillslope diffuser:
        >>> fr = FlowRouter(mg, method='D4')
        >>> kappa = 0.001
        >>> tl_diff = TransportLengthHillslopeDiffuser(mg,
        ...                                            erodibility=kappa,
        ...                                            slope_crit=0.6)

        Instantiate clast collection, one clast, starting position near the
        middle of the grid, its base near the surface:
        >>> init_position = [1.8,
        ...                  2.2,
        ...                  mg.at_node['topographic__elevation'][12]-0.02]
        >>> CC = ClastCollection(mg,
        ...                 clast_x=[init_position[0]],
        ...                 clast_y=[init_position[1]],
        ...                 clast_elev=[init_position[2]],
        ...                 clast_radius=np.ones(1)*0.5)

        Set seed for randomness in clast collection methods
        >>> np.random.seed(seed = 5000)

        Run for one timestep of 50y
        >>> fr.run_one_step()
        >>> tl_diff.run_one_step(50)
        >>> CC.run_one_step(dt=50,
        ...                 Si=1.2,
        ...                 kappa=kappa,
        ...                 uplift=None,
        ...                 erosion_method='TLDiff',
        ...                 hillslope_river__threshold=1e10,
        ...                 lateral_spreading='off')

        Clast did not move (because not entirely in eroded sediment layer)
        Check that position has not changed:
        >>> np.array_equal(
        ...     [CC.DataFrame.at[0,'clast__x'],
        ...      CC.DataFrame.at[0,'clast__y'],
        ...      CC.DataFrame.at[0,'clast__elev']],
        ...      init_position)
        True

        Run again for one timestep, clast will move:
        >>> fr.run_one_step()
        >>> tl_diff.run_one_step(50)
        >>> CC.run_one_step(dt=50,
        ...                 Si=1.2,
        ...                 kappa=kappa,
        ...                 uplift=None,
        ...                 erosion_method='TLDiff',
        ...                 hillslope_river__threshold=1e10,
        ...                 lateral_spreading='off')

        Clast has moved within same cell:
        >>> np.array_equal(
        ...     [CC.DataFrame.at[0,'clast__x'],
        ...      CC.DataFrame.at[0,'clast__y'],
        ...      CC.DataFrame.at[0,'clast__elev']],
        ...      init_position)
        False
        >>> CC.DataFrame.at[0,'element_id']
        4
        >>> np.greater(CC.DataFrame.at[0,'total_travelled_dist'], 0.0)
        True


        Run again for longer time, with lateral spreading,
        clast will move to next cell:

        >>> for t in range(0, 10001, 50):
        ...    np.random.seed(seed = 5000)
        ...    fr.run_one_step()
        ...    tl_diff.run_one_step(50)
        ...    CC.run_one_step(dt=50,
        ...                    Si=1.2,
        ...                    kappa=kappa,
        ...                    uplift=None,
        ...                    erosion_method='TLDiff',
        ...                    hillslope_river__threshold=1e10,
        ...                    lateral_spreading='on')

        Clast has moved to next cell:
        >>> CC.DataFrame.at[0,'element_id']
        1

        Run again, will move to next cell which is the outlet:

        >>> mg.at_node['topographic__elevation'][2]=-1.

        >>> for t in range(0, 501, 50):
        ...    fr.run_one_step()
        ...    tl_diff.run_one_step(50)
        ...    CC.run_one_step(dt=50,
        ...                    Si=1.2,
        ...                    kappa=kappa,
        ...                    uplift=None,
        ...                    erosion_method='TLDiff',
        ...                    hillslope_river__threshold=1e10,
        ...                    lateral_spreading='off')

        Clast has moved to next cell (reference node is now boundary node 2):
        >>> [CC.DataFrame.at[0,'element_id'], CC.DataFrame.at[0,'clast__node']]
        [-1, 2]
        Clast is at the outlet, has gone out of grid so it is considered
        as 'phantom'
        >>> CC.phantom(0)
        True


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
        self._nb_of_clast = len(_clast__node)

        # Store the input information in a dictionary and create the other
        # fields that will be populated later:

        # MOVED TO DOCUMENTATION ##############################################
        # clast__x = x coordinate of clast
        # clast__y = y coordinate of clast
        # clast__elev = clast elevation
        # clast__node = clast reference (closest) node
        # clast__radius = radius (m)
        # lambda_0 = characteristic length scale
        # lambda_mean = mean travel distance
        # slope__WE = west-east slope at the current position of the clast
        # slope__SN = south-north slope at the current position of the clast
        # slope__steepest_azimuth = direction (anticlockwise, from East) of
        #                                                   the steepest slope
        # slope__steepest_dip = angle (rads, positive) of the steepest slope
        # distance__to_exit = distance from the clast's position to the cell
        #           face or corner through which the clast will exit the cell,
        #           along the steepest slope (m)
        # target_node = node of the cell towards which the clast is moving
        # target_node_flag = digit representing the position of the target node
        #       relative to the current reference node : 0=E, 1=N, 2=W, 3=S,
        #       4=NE, 5=NW, 6=SW, 7=SE
        # REMOVE NEXT TWO?? not used ##########################################
        # change_x = change in the x coordinate of the clast due to travel step
        # change_y = change in the y coordinate of the clast due to travel step
        #######################################################################
        # hop_length = distance travelled in the current time step
        # total_travelled_dist = total distance that the clast has travelled
        #######################################################################

        clast_data = {'clast__x' : clast_x,
                      'clast__y' : clast_y,
                      'clast__elev' : clast_elev,
                      'clast__node' : _clast__node,
                      'clast__initial_radius' : clast_radius,
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
                      'change_x' : np.zeros(self._nb_of_clast), #remove?
                      'change_y' : np.zeros(self._nb_of_clast), #remove?
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
        which the clast will be moved (target node) and the slope used to
        move the clast.

        If the clast is on a hillslop cell and none of the neighbor nodes
        are boundary nodes, calculate the (signed) centered slopes around
        clast's reference node: west-east slope is elevation difference
        between west neighbor and east neighbor, and north-south slope is
        elevation difference between north neighbor and south neighbor.

        If the clast is on a river cell and/or if at least one of the neighbor
        node is a boundary node, the slope used to move the clast is provided
        by Flow director: the slope is either the provided
        topographic__steepest_slope or 0 if clast's node is sink.

        Parameters
        ----------
        clast : int
            Clast ID
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
                    # slopes will be calculated by _move_to depending on the
                    # target node, they are set as NaN for now:
                    df.at[clast, 'slope__WE'] = np.NaN
                    df.at[clast, 'slope__SN'] = np.NaN
                else: # flow receiver = node itself
                    df.at[clast, 'target_node'] = _node
                    df.at[clast, 'target_node_flag'] = -1
                    df.at[clast, 'slope__WE'] = np.NaN
                    df.at[clast, 'slope__SN'] = np.NaN

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
        travel, and calculate the travel distance for the clast to exit the
        current cell.

        Parameters
        ----------
        clast : int
            Clast ID
        """

        # Saving references:
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
        _diagonal_adjacent_nodes_at_node = (
                _grid.diagonal_adjacent_nodes_at_node[_node])

        # If clast is phantom, it will not move, set everything to 0:
        if ClastCollection.phantom(self, clast) == True:
            target_node = df.at[clast, 'clast__node']
            ss_azimuth = np.NaN
            ss_dip = 0.
            dist_to_exit = np.NaN
            [change_x, change_y] = [0., 0.]

        # Clast is not phantom:
        # Cases where one of the neighbor is boundary node or clast is on a
        # river cell: use the topographic__steepest_slope and
        # flow__receiver_node provided by flow director:
        # TO DO: use Dinf? calculate ss_azimuth within the steepest triangle
        elif (np.isnan(df.at[clast, 'slope__WE']) == True and\
              np.isnan(df.at[clast, 'slope__SN']) == True) or \
              ClastCollection._cell_is_hillslope(self, clast) == False:
            print('clast is next to boundary or on river cell')
            df.at[clast, 'target_node'] = (
                    self._grid.at_node['flow__receiver_node'][_node])
            target_node = df.at[clast, 'target_node']
            print('receiver is : %s' %target_node)
            target_node_flag = df.at[clast, 'target_node_flag']

            if target_node_flag == -1: # node is sink, clast does not move
                we_slope = np.NaN
                sn_slope =np.NaN
                ss_azimuth = np.NaN
                ss_dip = 0.
                dist_to_exit = np.inf
                [change_x, change_y] = [0.,0.]
                target_node = _node

            elif target_node_flag == 0: # East
                we_slope = (
                        _node_z - self._grid.at_node[
                                'topographic__elevation'][
                                        df.at[clast, 'target_node']]) / (
                                _grid.dx) # = topographic__steepest_slope?
                sn_slope = 0.
                ss_azimuth = 0. # east

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(
                        np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector
                # = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * (
                        (_node_x + (self._grid.dx/2)) - _clast_x)
                [change_x, change_y] = [dist_to_exit + _radius, 0.0]
################## FOR TESTING ONLY ##################
                if df.at[clast, 'slope__WE'] < 0:
                    print('error')
                else:
                    pass
###############################################################
            elif target_node_flag == 1: # North
                we_slope = 0.
                sn_slope = (_node_z - (
                        self._grid.at_node['topographic__elevation'][
                                df.at[clast, 'target_node']])) / (
                                _grid.dy)
                ss_azimuth = np.radians(90) # north

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(
                        np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector
                # = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * (
                        (_node_y + (self._grid.dy/2)) - _clast_y)
                [change_x, change_y] = [0.0, dist_to_exit + _radius]
            elif target_node_flag == 2: # West
                we_slope = (self._grid.at_node['topographic__elevation'][
                        df.at[clast, 'target_node']] - _node_z) / _grid.dx
                sn_slope = 0.
                ss_azimuth = np.radians(180) # west

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(
                        np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector
                # = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * (
                        _clast_x - (_node_x - (self._grid.dx/2)))
                [change_x, change_y] = [-dist_to_exit - _radius, 0.0]
            elif target_node_flag == 3: # South
                we_slope = 0.
                sn_slope = (self._grid.at_node['topographic__elevation'][
                        df.at[clast, 'target_node']] - _node_z) / _grid.dy
                ss_azimuth = np.radians(270) # south

                # norm of steepest slope vector projected on horizontal plane:
                ss_horiz_norm = np.sqrt(
                        np.power(we_slope, 2) + np.power(sn_slope, 2))
                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
                # norm of steepest slope vector
                # = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
                # dip of steepest slope:
                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                dist_to_exit = (1 / np.cos(ss_dip)) * (
                        _clast_y - (_node_y - (self._grid.dy/2)))
                [change_x, change_y] = [0.0, -dist_to_exit - _radius]

            else: # Diagonals
                we_slope = np.NaN
                sn_slope = np.NaN
                _target_node_x = _grid.node_x[df.at[clast, 'target_node']]
                _target_node_y = _grid.node_y[df.at[clast, 'target_node']]
                _target_node_z = _grid.at_node['topographic__elevation'][
                        df.at[clast, 'target_node']]

## TO DO : Format to Pep8
                ss_dip = np.arctan(_node_z - _target_node_z) / (np.sqrt(np.power(_grid.node_x[target_node]-_grid.node_x[_node],2)+np.power(_grid.node_y[target_node]-_grid.node_y[_node],2)))
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

                else: #FOR TESTING ONLY
                    print('error target node flag')

#FOR TESTING ONLY
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
            ss_horiz_norm = np.sqrt(
                    np.power(we_slope, 2) + np.power(sn_slope, 2))

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
                    ss_azimuth = np.NaN
                    dist_to_exit = np.inf
                    target_node = _node
                    [change_x, change_y] = [0.0, 0.0]
        #            # OPTION2 (to do?):
        #            # clast moves to random direction, according to lambda_0:
        #            ss_azimuth = np.random.uniform(0.0, 2*pi, 1)
        #            df.at[clast, 'slope__steepest_azimuth'] = ss_azimuth
        #            df.at[clast, 'slope__flag'] = '
        #            self.rand_length = np.random.exponential(
        #                                   scale=df.at[clast, 'lambda_0'], 1)
        #            dist_to_exit =


                else: # SN slope is not 0, ss direction is S or N
                    # dip of steepest slope:
                    ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                    if sn_slope < 0: # ss direction is South
                        print('going South')
#                        if ss_dip != np.arctan(np.abs(sn_slope)):   # radians
#                            print('error, dip is %s' %ss_dip)
#                            print('should be:')
#                            print(np.arctan(np.abs(sn_slope)))

                        ss_azimuth = np.radians(270) # South
                        dist_to_exit = (1 / np.cos(ss_dip)) * (
                                _clast_y - (_node_y - (self._grid.dy/2)))
                        target_node = _row_col_adjacent_nodes_at_node[3]
                        [change_x, change_y] = [0.0, -dist_to_exit - _radius]
                    else: # ss direction is North
                        print('going North')
                        if ss_dip != np.arctan(np.abs(sn_slope)):
                            print('error, dip is %s' %ss_dip)

                        ss_azimuth = np.radians(90) # North
                        dist_to_exit = (1 / np.cos(ss_dip)) * (
                                (_node_y + (self._grid.dy/2)) - _clast_y)
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
                        dist_to_exit = (1 / np.cos(ss_dip)) * (
                                _clast_x - (_node_x - (self._grid.dx/2)))
                        target_node = _row_col_adjacent_nodes_at_node[2]
                        [change_x, change_y] = [-dist_to_exit - _radius, 0.0]
                    else: # ss direction is East
                        print('East')
                        ss_dip = np.arctan(np.abs(we_slope))
                        ss_azimuth = 0 # East
                        dist_to_exit = (1 / np.cos(ss_dip)) * (
                                (_node_x + (self._grid.dx/2)) - _clast_x)
                        target_node = _row_col_adjacent_nodes_at_node[0]
                        [change_x, change_y] = [dist_to_exit + _radius, 0.0]

                else: # sn_slope is not 0
                    if sn_slope > 0 and we_slope > 0: # Quarter = NE
                        ss_azimuth = np.arctan(np.abs(sn_slope / we_slope))
                        corner = 0
                        clast_to_corner_azimuth = np.arctan(
                                np.abs(_clast_y - self.y_of_corners_at_node[
                                        _node, corner])/
                                       np.abs(_clast_x - (
                                               self.x_of_corners_at_node[
                                                       _node, corner])))
                        if ss_azimuth < clast_to_corner_azimuth: # Eigth = NE-row
                            print('NE-R')
                            dist_to_exit = (1 / np.cos(ss_dip)) * (
                                    (np.abs(self.x_of_corners_at_node[
                                            _node, corner] - _clast_x)) / (
                                                        np.cos(ss_azimuth)))
                            target_node = _row_col_adjacent_nodes_at_node[0]
# TO DO: Format to Pep8
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


        # Calculate lambda_0 and lambda_mean :
        df.at[clast, 'lambda_0'] = ((self._dt * self._kappa * _grid.dx) / (2 * df.at[clast,'clast__radius'])) * (np.tan(df.at[clast, 'slope__steepest_dip']) * ((self._Si-np.tan(df.at[clast, 'slope__steepest_dip']))/(self._Si+np.tan(df.at[clast, 'slope__steepest_dip']))))
                #(self._dt * self._kappa * _grid.dx) / (
#                        2 * df.at[clast,'clast__radius'])) #self._Si *
# self._Si *
#
        if ClastCollection._cell_is_hillslope(self, clast) == True:
            # Cell is hillslope:
            if np.tan(df.at[clast, 'slope__steepest_dip']) >= self._Si:
                lambda_0 = df.at[clast, 'lambda_0']
                lambda_mean = np.power(10, 10)
            else:
                lambda_0  = df.at[clast,'lambda_0']
                lambda_mean = lambda_0 * (
                        self._Si + np.tan(
                                df.at[clast, 'slope__steepest_dip'])) / (
                                self._Si - np.tan(
                                        df.at[clast, 'slope__steepest_dip']))
        else:   # Cell is river:
            if np.tan(df.at[clast, 'slope__steepest_dip'] >= self._Si):
                lambda_0 = df.at[clast, 'lambda_0']
                lambda_mean = np.power(10, 10)
            else:
                lambda_0 = 10 * df.at[clast,'lambda_0']
                lambda_mean = lambda_0 * (
                        self._Si + np.tan(
                                df.at[clast, 'slope__steepest_dip'])) / (
                                self._Si - np.tan(
                                        df.at[clast, 'slope__steepest_dip']))

        df.at[clast, 'lambda_mean'] = lambda_mean
        df.at[clast, 'lambda_0'] = lambda_0

# For testing only: #####################
        print('change_x = %s' %change_x)
        print('change_y = %s' %change_y)
#########################################

        # Add randomness for lateral spreading (if option on):
        if self._lateral_spreading == 'on':
            if df.at[clast, 'target_node'] in _grid.boundary_nodes:
                [delta_x,delta_y] = [0,0] # to avoid clast going out of the grid
            else:
                # Create truncated normal distribution
                # (centered in 0; max, min = +/- lambda0/2)
                x=truncnorm(-
                            (df.at[clast, 'lambda_0']/2), (
                                    df.at[clast, 'lambda_0']/2), loc=0, scale=(
                                            df.at[clast, 'lambda_0']/4))
                [delta_x,delta_y] = x.rvs(2)/10

                df.at[clast, 'change_x'] = change_x + delta_x
                df.at[clast, 'change_y'] = change_y + delta_y
                # formerly: + np.random.normal(loc=0.0, scale=(df.at[clast, 'lambda_mean']/_grid.dx), size=1)   # adds a randomness for lateral spreading
        else:
            df.at[clast, 'change_x'] = change_x
            df.at[clast, 'change_y'] = change_y


    def _change_cell_proba(self, clast):
        """ Determine the probability for a clast to exit from its current
        cell. A random travel distance is drawn from a probability distribution
        calculated for each clast, on its current cell. If this travel distance
        is higher than the distance for the clast to exit the cell, then the
        clast will change cell.

        Parameters
        ----------
        clast : int
            Clast ID
        """
        lambda_mean = self.df.at[clast, 'lambda_mean']
        dist_to_exit = self.df.at[clast, 'distance__to_exit']

        if np.isinf(dist_to_exit): # cases where slope is null or sink
            _change_cell = False
        else:
            # Draw a random sample in the probability distribution
            # of travel distances:
            self.rand_length = np.random.exponential(scale=lambda_mean, size=1)

            if self.rand_length < dist_to_exit: # clast stays in cell
                _change_cell = False
            else: # self.rand_length >= dist_to_exit: clast leaves cell
                _change_cell = True

# For testing only:##################################
        print('rand_length = %s' % self.rand_length)
        print('dist_to_exit= %s' % dist_to_exit)
#####################################################

        return _change_cell


    def _move_in_cell(self, clast):
        """ Move the clast within the cell, along the path of steepest slope,
        of a distance corresponding to the drawn travel distance if that travel
        distance drawn from the probability distribution is smaller than the
        distance to exit the cell.

        Parameters
        ----------
        clast : int
            Clast ID
        """
        # clast stays in cell, move of distance rand_length along slope
        ss_azimuth = self.df.at[clast, 'slope__steepest_azimuth']
        ss_dip = self.df.at[clast, 'slope__steepest_dip']
        x_horizontal = self.rand_length * np.cos(ss_dip)

        if np.isnan(ss_azimuth): # clast is in sink
            [change_x, change_y] = [0., 0.]
        elif ss_azimuth <= np.radians(90):
            [change_x, change_y] = [(
                    x_horizontal * np.cos(ss_azimuth)), (
                            x_horizontal * np.sin(ss_azimuth))]
        elif ss_azimuth <= np.radians(180):
            [change_x, change_y] = [(
                    -x_horizontal * np.cos(np.radians(180)-ss_azimuth)), (
                            x_horizontal * np.sin(np.radians(180)-ss_azimuth))]
        elif ss_azimuth <= np.radians(270):
            [change_x, change_y] = [(
                    -x_horizontal * np.sin(np.radians(270)-ss_azimuth)), (
                            -x_horizontal * np.cos(
                                    np.radians(270)-ss_azimuth))]
        else: # ss_azimuth <= np.radians(360)
            [change_x, change_y] = [(
                    x_horizontal * np.cos(np.radians(360)-ss_azimuth)), (
                            -x_horizontal * np.sin(
                                    np.radians(360)-ss_azimuth))]

        # Update clast coordinates:
        self.df.at[clast, 'clast__x'] += change_x
        self.df.at[clast, 'clast__y'] += change_y
        # Update hop length:
        self.df.at[clast, 'hop_length'] += self.rand_length

    def _move_out_of_cell(self, clast):
        """ Move the clast out of the cell, along the path of steepest slope.
        The clast gets placed fully in the target cell, just along the face it
        just crossed. This happens if the travel distance drawn from the
        probability distribution is larger than the distance to exit the cell.

        Parameters
        ----------
        clast : int
            Clast ID
        """
        # clast leaves cell, move of distance dist_to_exit along slope
        self.df.at[clast, 'clast__x'] += self.df.at[clast, 'change_x']
        self.df.at[clast, 'clast__y'] += self.df.at[clast, 'change_y']
        self.df.at[clast, 'hop_length'] += (
                np.sqrt(np.power(
                        self.df.at[clast, 'change_x'], 2)+ np.power(
                                self.df.at[clast, 'change_y'],2))) / (
                        np.cos(self.df.at[clast, 'slope__steepest_dip']))

        self.df.at[clast, 'clast__node'] = self.df.at[clast, 'target_node']
        self.df.at[clast, 'element_id'] = self._grid.cell_at_node[self.df.at[clast, 'clast__node']]


    def phantom(self, clast):
        """Check if the clast has exited the grid.

        When a clast reaches a boundary node, it exits the grid and is thus
        flagged as phantom. Also, when a clast's radius is reduced to zero by
        abrasion, it is considered as phantom. In both cases, the clast is not
        erased from the dataframe but is not moved anymore.

        Parameters
        ----------
        clast : int
            Clast ID.

        """

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

        The area/slope ratio is used to determine if a cell belongs to the
        hillslope or river domain.

        See Tucker and Bras, 1998, Water Resources Research, for complete
        justification.

        Parameters
        ----------
        clast : int
            Clast ID

        """

        df = self.DataFrame
        _node = df.at[clast, 'clast__node']
        _grid = self._grid
        threshold = self._threshold

        # If next to boundary, set A/S to 1:
        if np.isnan(df.at[clast, 'slope__steepest_azimuth']) == True:
            area_over_slope = 1   # clast next to boundary
        else:
            area_over_slope = (_grid.at_node['drainage_area'][_node] /
                               _grid.at_node[
                                       'topographic__steepest_slope'][_node])

        if area_over_slope < threshold:
            _cell_is_hillslope = True
        else:
            _cell_is_hillslope = False

        return _cell_is_hillslope

    def clast_solver_Exponential(self,
                                 dt=1.,
                                 Si=1.2,
                                 kappa=0.0001,
                                 uplift=None,
                                 erosion_method='TLDiff',
                                 hillslope_river__threshold=1e4,
                                 lateral_spreading='off'):
        """See :func:`run_one_step`

        """
        self.df=self.DataFrame
        # Method loop: Depending on the method used to evolve the landscape,
        # get sediment influx, erosion and deposition
        self.erosion_method = erosion_method
        if self.erosion_method == 'TLDiff':
            self._erosion_rate = self._grid.at_node['sediment__erosion_rate']
            self._deposition_rate = self._grid.at_node[
                    'sediment__deposition_rate']
            self._sediment__flux_in = self._grid.at_node['sediment__flux_in'] # NOT USED

        elif self.erosion_method == 'Space':
            from landlab.components import Space
            for obj in gc.get_objects():
                if isinstance(obj, Space):    # look for instance of Space
                    self._erosion_rate = obj.Es + obj.Er   # Works if only one instance of Space was made
                    self._deposition_rate = obj.depo_rate
                    self._sediment__flux_in = obj.qs_in  # NOT USED

        else:
            raise ValueError('Erosion method must be "TLDiff" or "Space"')

        # Future version: multiple compo -> add fluxes?
        # Store values that will be used:
        self._kappa = kappa # if using SPACE, kappa=K_sed. TO DO: can be a field or array
        self._dt = dt
        self._erosion__depth = self._erosion_rate * self._dt
        self._deposition__thickness = self._deposition_rate * self._dt
        #  slope above which particle motion continues indefinetly (not equal
        # to critical slope of diffusion, see Furbish and Haff, 2010):
        self._Si = Si
        self._threshold = hillslope_river__threshold
        self._lateral_spreading = lateral_spreading


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

        # Treat clasts one after the other:
        for clast in range(self._nb_of_clast):
            # FOR TESTING ONLY ######
            print('CLAST %s' %clast)
            #########################

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
                                if np.isinf(self.df.at[clast,'distance__to_exit']):# case where slope is null or simk
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


            else: # clast is phantom, go to next clast
                # for display purpose: phantom clast has a radius of 0
                # self._clast__radius(i) = 0.
                pass


    def run_one_step(self, dt=1.,
                     Si=1.2,
                     kappa=0.0001,
                     uplift=None,
                     erosion_method='TLDiff',
                     hillslope_river__threshold=1e4,
                     lateral_spreading='off'):

        """ Run the clast-tracker for one timestep dt: evaluate the type of
        cell each clast is currently in, wheteher each clast is detached,
        whether it is moved, where and how far it is moved and update each
        clast's characteristics.

        Parameters
        ----------
        dt : float
            Timestep
        Si : float
            Magnitude of a critical slope beyond which particles in motion
            are not stopped [L/L]
        kappa : float
            Erodibility (units vary)
        uplift : float, int, field or array
            Uplift rate [L/T]
        erosion_method : str
            Method (component) used to evolve the landscape. Options at present
            include:
                - 'TLDiff' (default): TransportLengthHillslopeDiffuser
                - 'Space': Space
        hillslope_river__threshold : float
            Threshold value of the drainage area/slope ratio over which a cell
            will be considered a river cell [L^2]
        lateral_spreading : str
            Optional, enables the lateral spreading of clasts (random deviation
            from steepest slope direction, scaled to lambda_0) when 'on', no
            deviation if 'off'. Default is 'off'.

        """
        self.clast_solver_Exponential(dt,
                                      Si,
                                      kappa,
                                      uplift,
                                      erosion_method,
                                      hillslope_river__threshold,
                                      lateral_spreading)