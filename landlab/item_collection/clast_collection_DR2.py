#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import gc
from landlab.data_record import DataRecord
from matplotlib.pyplot import figure, show, plot, xlabel, ylabel, title, legend

from math import pi

from scipy.stats import truncnorm

class ClastCollection(DataRecord):
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
        ...                                               mg['node'][
        ...                                         'topographic__elevation'],
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
        # Clast collection size:
        self._nb_of_clast = len(_clast__node)
        # CHeck clast elev (not above node elev):
        _clast__elev=clast_elev
        for clast in range(self._nb_of_clast):
            if _clast__elev[clast] > \
            self._grid.at_node['topographic__elevation'][_clast__node[clast]]:
                _clast__elev[clast] = \
                self._grid.at_node['topographic__elevation'][
                        _clast__node[clast]] # snaps clast to grid if was above
        # Determine current cell for each clast:
        _clast__cell=[]
        _clast__cell[:] = self._grid.cell_at_node[_clast__node]

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
        # change_x = change in the x coordinate of the clast due to travel step
        # change_y = change in the y coordinate of the clast due to travel step
        # hop_length = distance travelled in the current time step
        # total_travelled_dist = total distance that the clast has travelled
        #######################################################################

        clast_data = {'clast__x' : (['item_id', 'time'], clast_x.reshape((self._nb_of_clast,1))),
                      'clast__y' : (['item_id', 'time'], clast_y.reshape((self._nb_of_clast,1))),
                      'clast__elev' : (['item_id', 'time'], clast_elev.reshape((self._nb_of_clast,1))),
                      'clast__node' : (['item_id'], np.array(_clast__node)),
                      'clast__initial_radius' : (['item_id'], clast_radius),
                      'clast__radius' : (['item_id', 'time'], clast_radius.reshape((self._nb_of_clast,1))),
                      'lambda_0' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'lambda_mean' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'slope__WE' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'slope__SN' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'slope__steepest_azimuth' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'slope__steepest_dip' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'distance__to_exit' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'target_node' : (['item_id'], -np.ones(self._nb_of_clast, dtype=int)),
                      'target_node_flag': (['item_id'], -np.ones(self._nb_of_clast, dtype=int)),
                      'distance__to_travel' : (['item_id'], np.full(self._nb_of_clast, np.NaN)),
                      'change_x' : (['item_id'], np.zeros(self._nb_of_clast)), #remove?
                      'change_y' : (['item_id'], np.zeros(self._nb_of_clast)), #remove?
                      'hop_length' : (['item_id', 'time'], np.zeros((self._nb_of_clast,1))),
                      'total_travelled_dist' : (['item_id'], np.zeros(self._nb_of_clast)),
                      'close2boundary' : (['item_id'], np.full(self._nb_of_clast, np.NaN))}

        # Build ItemCollection containing clast data:
        super(ClastCollection, self).__init__(
                                self._grid,
                                time=[0.],
                                items={'grid_element': 'cell',
                                       'element_id': np.array(
                                               (_clast__cell),
                                               dtype=int).reshape((
                                                       self._nb_of_clast,1))},
                                data_vars=clast_data,
                                attrs=None,
                                compat='broadcast_equals')


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

        # Calculate lambda_0 :
        self['lambda_0'][clast] = self.attrs['kappa'] / (2 * self['clast__radius'][clast, -1] * self.attrs['disturbance_fqcy'])


        # If clast is phantom, slope is 0:
        if ClastCollection.phantom(self, clast) == True or self['clast__radius'][clast, -1] == 0:
            self['slope__WE'][clast] = 0
            self['slope__SN'][clast] = 0

        else:
            # Determine reference (closest) node for each clast:
            _grid = self._grid
            _node = self['clast__node'][clast] = \
                _grid.find_nearest_node((self['clast__x'][clast,-1],\
                                       self['clast__y'][clast,-1]))

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
            print(_neighbor_nodes)
            # Case where one of the neighbor is boundary:
            if any(i in _row_col_adjacent_nodes_at_node for\
                       i in _grid.boundary_nodes) == True or \
                           ClastCollection._cell_is_hillslope(self, clast) == \
                               False:
                print('close to boundary and/or river')
                self['close2boundary'][clast].values=True
                if _grid.at_node['flow__receiver_node'][_node] != _node:
                    # if FlowDirector designates a receiver node other than
                    # the node itself, clast is moving toward receiver node:
                    self['target_node'][clast] = (
                            _grid.at_node['flow__receiver_node'][_node])
                    print('target node: ')
                    print(self['target_node'][clast])
                    self['target_node_flag'][clast] = np.where(
                            _neighbor_nodes == _grid.at_node[
                                    'flow__receiver_node'][_node])[0][0]
                    # slopes will be calculated by _move_to depending on the
                    # target node, they are set as NaN for now:
                    self['slope__WE'][clast] = np.NaN
                    self['slope__SN'][clast] = np.NaN
                else: # flow receiver = node itself
                    self['target_node'][clast] = _node
                    self['target_node_flag'][clast] = -1
                    self['slope__WE'][clast] = np.NaN
                    self['slope__SN'][clast] = np.NaN

            else: # if not close to boundary
                # target_node and target_node_flag determined by _move_to
                # Calculation of slopes: W to E and S to N, units=m/m
                self['slope__WE'][clast] = ((
                        _grid.at_node['topographic__elevation'][west_node]-(
                                _grid.at_node['topographic__elevation']\
                                    [east_node]))/(2*_grid.dx))

                self['slope__SN'][clast] = \
                    (_grid.at_node['topographic__elevation'][south_node]-
                          _grid.at_node['topographic__elevation']\
                              [north_node])/(2*_grid.dy)

        ####### ADD LATERAL SPREADING HERE?

                # norm of steepest slope vector projected on horizontal plane:
                we_slope = self['slope__WE'][clast].values.item()
                sn_slope = self['slope__SN'][clast].values.item()
                ss_horiz_norm = np.sqrt(np.power(we_slope, 2) + np.power(sn_slope, 2))

                # norms of vectors SN and WE:
                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))

                # norm of steepest slope vector = norm of resultant of SN and WE:
                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))

                ss_dip = np.arccos(ss_horiz_norm / ss_norm)


                        # Add randomness for lateral spreading (if option on):
                if self.attrs['lateral_spreading'] == 'on':
                    # maximum angle of deviation from steepest slope direction:
                    max_deviation = np.pi * np.exp(ss_dip/(-np.sqrt(self['lambda_0'][clast].values.item()))) * np.cos(ss_dip)
                    print('max_dev=%s' %max_deviation)
                    # draw random angle of deviation in normal distrib:
#                    draw_dev=truncnorm((-max_deviation),
#                                        max_deviation,
#                                        loc=0.0,
#                                        scale=1).rvs(1)

                    draw_dev=np.random.uniform(low=-max_deviation,
                                                  high=max_deviation,
                                                  size=1)
                    print('draw_dev=%s' %draw_dev)

                    slope_WE_with_dev = self['slope__WE'][clast].values.item() * np.cos(draw_dev[0]) - self['slope__SN'][clast].values.item() * np.sin(draw_dev[0])
                    slope_SN_with_dev = self['slope__WE'][clast].values.item() * np.sin(draw_dev[0]) + self['slope__SN'][clast].values.item() * np.cos(draw_dev[0])
                    print('slope_we_with_dev=%s' %slope_WE_with_dev)
                    self['slope__WE'][clast] = slope_WE_with_dev
                    self['slope__SN'][clast] = slope_SN_with_dev




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
        _node = self['clast__node'][clast]
        _node_x = _grid.node_x[_node]
        _node_y = _grid.node_y[_node]
        _node_z = _grid.at_node['topographic__elevation'][_node]
        _clast_x = self['clast__x'][clast, -1].values.item()
        _clast_y = self['clast__y'][clast, -1].values.item()
        we_slope = self['slope__WE'][clast].values.item()
        sn_slope = self['slope__SN'][clast].values.item()
#        _radius = self['clast__radius'][clast][-2]
        # Adjacent row and col nodes:
        _row_col_adjacent_nodes_at_node = _grid.neighbors_at_node[_node]
        # Adjacent diagonal nodes:
        _diagonal_adjacent_nodes_at_node = (
                _grid.diagonal_adjacent_nodes_at_node[_node])

        # If clast is phantom, it will not move, set everything to 0:
        if ClastCollection.phantom(self, clast) == True:
            target_node = self['clast__node'][clast]
            ss_azimuth = np.NaN
            ss_dip = 0.
            dist_to_exit = np.NaN
            [change_x, change_y] = [0., 0.]

        # Clast is not phantom:
        # Cases where one of the neighbor is boundary node or clast is on a
        # river cell: use the topographic__steepest_slope and
        # flow__receiver_node provided by flow director:
        elif (np.isnan(self['slope__WE'][clast]) == True and\
              np.isnan(self['slope__SN'][clast]) == True) or \
              self._cell_is_hillslope(clast) == False:
            print('clast is next to boundary or on river cell')
            target_node = self['target_node'][clast] = (
                    _grid.at_node['flow__receiver_node'][_node])
            print('receiver is : %s' %target_node)
            target_node_flag = self['target_node_flag'][clast]

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
                                        self['target_node'][clast]]) / (
                                _grid.dx) # = topographic__steepest_slope?
                sn_slope = 0.
                ss_azimuth = 0. # east

#                # norm of steepest slope vector projected on horizontal plane:
#                ss_horiz_norm = np.sqrt(
#                        np.power(we_slope, 2) + np.power(sn_slope, 2))
#                # norms of vectors SN and WE:
#                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
#                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
#                # norm of steepest slope vector
#                # = norm of resultant of SN and WE:
#                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
#                # dip of steepest slope:
#                ss_dip = np.arccos(ss_horiz_norm / ss_norm)

                ss_dip = abs(np.arctan(we_slope))
                dist_to_exit = (1 / np.cos(ss_dip)) * abs(
                        (_node_x + (self._grid.dx/2)) - _clast_x)
                [change_x, change_y] = [abs((_node_x + (self._grid.dx/2)) - _clast_x) + 1e-10, 0.0]

            elif target_node_flag == 1: # North
                we_slope = 0.
                sn_slope = (_node_z - (
                        self._grid.at_node['topographic__elevation'][
                                self['target_node'][clast]])) / (
                                _grid.dy)
                ss_azimuth = np.radians(90) # north

#                # norm of steepest slope vector projected on horizontal plane:
#                ss_horiz_norm = np.sqrt(
#                        np.power(we_slope, 2) + np.power(sn_slope, 2))
#                # norms of vectors SN and WE:
#                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
#                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
#                # norm of steepest slope vector
#                # = norm of resultant of SN and WE:
#                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
#                # dip of steepest slope:
#                ss_dip = np.arccos(ss_horiz_norm / ss_norm)
                ss_dip = abs(np.arctan(sn_slope))

                dist_to_exit = (1 / np.cos(ss_dip)) * abs(
                        (_node_y + (self._grid.dy/2)) - _clast_y)
                [change_x, change_y] = [0.0, abs((_node_y + (self._grid.dy/2)) - _clast_y) + 1e-10]
            elif target_node_flag == 2: # West
                we_slope = (self._grid.at_node['topographic__elevation'][
                        self['target_node'][clast]] - _node_z) / _grid.dx
                sn_slope = 0.
                ss_azimuth = np.radians(180) # west

#                # norm of steepest slope vector projected on horizontal plane:
#                ss_horiz_norm = np.sqrt(
#                        np.power(we_slope, 2) + np.power(sn_slope, 2))
#                # norms of vectors SN and WE:
#                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
#                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
#                # norm of steepest slope vector
#                # = norm of resultant of SN and WE:
#                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
#                # dip of steepest slope:
#                ss_dip = np.arccos(ss_horiz_norm / ss_norm)
                ss_dip = abs(np.arctan(we_slope))
                dist_to_exit = (1 / np.cos(ss_dip)) * abs(
                        _clast_x - (_node_x - (self._grid.dx/2)))
                [change_x, change_y] = [-abs(_clast_x - (_node_x - (self._grid.dx/2))) - 1e-10, 0.0]
            elif target_node_flag == 3: # South
                we_slope = 0.
                sn_slope = (self._grid.at_node['topographic__elevation'][
                        self['target_node'][clast]] - _node_z) / (
                    _grid.dy)
                ss_azimuth = np.radians(270) # south
#
#                # norm of steepest slope vector projected on horizontal plane:
#                ss_horiz_norm = np.sqrt(
#                        np.power(we_slope, 2) + np.power(sn_slope, 2))
#                # norms of vectors SN and WE:
#                sn_norm = abs(sn_slope) / np.cos(np.arctan(abs(sn_slope)))
#                we_norm = abs(we_slope) / np.cos(np.arctan(abs(we_slope)))
#                # norm of steepest slope vector
#                # = norm of resultant of SN and WE:
#                ss_norm = np.sqrt(np.power(sn_norm, 2) + np.power(we_norm, 2))
#                # dip of steepest slope:
#                ss_dip = np.arccos(ss_horiz_norm / ss_norm)
                ss_dip = abs(np.arctan(sn_slope))
                print('ssdip= %s' %ss_dip)
                print('clastnode %s' %_node)
                print('clastnodey %s' %_node_y)
                print('clasty %s' %_clast_y)

                dist_to_exit = (1 / np.cos(ss_dip)) * abs(
                        _clast_y - (_node_y - (self._grid.dy/2)))
                print('dist to exit toward SOUTH: %s' %dist_to_exit)
                [change_x, change_y] = [0.0, -abs(_clast_y - (_node_y - (self._grid.dy/2))) - 1e-10]

            else: # Diagonals
                we_slope = np.NaN
                sn_slope = np.NaN
                _target_node_x = _grid.node_x[
                        self['target_node'][clast]]
                _target_node_y = _grid.node_y[
                        self['target_node'][clast]]
                _target_node_z = _grid.at_node['topographic__elevation'][
                        self['target_node'][clast]]

## TO DO : Format to Pep8
                ss_dip = np.arctan(_node_z - _target_node_z) / (np.sqrt(np.power(_grid.node_x[target_node]-_grid.node_x[_node],2)+np.power(_grid.node_y[target_node]-_grid.node_y[_node],2)))
                dist_to_exit = (1 / np.cos(ss_dip)) * np.sqrt(np.power((_target_node_x - _clast_x), 2) + np.power((_target_node_y - _clast_y), 2))

                if target_node_flag == 4: # NE
                    corner = 0
                    ss_azimuth = np.arctan((self.y_of_corners_at_node[_node, corner] - _clast_y) / (self.x_of_corners_at_node[_node, corner] - _clast_x))
                    [change_x, change_y] = [self.x_of_corners_at_node[_node, corner] - _clast_x +  1e-10, self.y_of_corners_at_node[_node, corner] - _clast_y + 1e-10]

                elif target_node_flag == 5: # NW
                    corner = 1
                    ss_azimuth = np.radians(90) + np.arctan((_clast_x - self.x_of_corners_at_node[_node, corner]) / (self.y_of_corners_at_node[_node, corner] - _clast_y))
                    [change_x, change_y] = [-(_clast_x - self.x_of_corners_at_node[_node, corner]) - 1e-10, self.y_of_corners_at_node[_node, corner] - _clast_y +  1e-10]

                elif target_node_flag == 6: # SW
                    corner = 2
                    ss_azimuth = np.radians(180) + np.arctan((_clast_y - self.y_of_corners_at_node[_node, corner]) / (_clast_x - self.x_of_corners_at_node[_node, corner]))
                    [change_x, change_y] = [-(_clast_x - self.x_of_corners_at_node[_node, corner]) - 1e-10, -(_clast_y - self.y_of_corners_at_node[_node, corner]) - 1e-10]

                elif target_node_flag == 7: # SE
                    corner = 3
                    ss_azimuth = np.radians(270) + np.arctan((self.x_of_corners_at_node[_node, corner] - _clast_x) / (_clast_y - self.y_of_corners_at_node[_node, corner]))
                    [change_x, change_y] = [(self.x_of_corners_at_node[_node, corner] - _clast_x) + 1e-10, -(_clast_y - self.y_of_corners_at_node[_node, corner]) - 1e-10]

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


            self['slope__WE'][clast] = we_slope
            self['slope__SN'][clast] = sn_slope
            self['slope__steepest_azimuth'][clast] = ss_azimuth
            self['slope__steepest_dip'][clast] = ss_dip
            self['distance__to_exit'][clast] = dist_to_exit
            self['target_node'][clast] = target_node
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
                        print(sn_slope)

                        ss_azimuth = np.radians(270) # South
                        dist_to_exit = (1 / np.cos(ss_dip)) * (
                                _clast_y - (_node_y - (self._grid.dy/2)))
                        target_node = _row_col_adjacent_nodes_at_node[3]
                        [change_x, change_y] = [0.0, -dist_to_exit - 1e-10]
                    else: # ss direction is North
                        print('going North')
                        if ss_dip != np.arctan(np.abs(sn_slope)):
                            print('error, dip is %s' %ss_dip)

                        ss_azimuth = np.radians(90) # North
                        dist_to_exit = (1 / np.cos(ss_dip)) * (
                                (_node_y + (self._grid.dy/2)) - _clast_y)
                        target_node = _row_col_adjacent_nodes_at_node[1]
                        [change_x, change_y] = [0.0, dist_to_exit + 1e-10]

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
                        [change_x, change_y] = [-dist_to_exit - 1e-10, 0.0]
                    else: # ss direction is East
                        print('East')
                        ss_dip = np.arctan(np.abs(we_slope))
                        ss_azimuth = 0 # East
                        dist_to_exit = (1 / np.cos(ss_dip)) * (
                                (_node_x + (self._grid.dx/2)) - _clast_x)
                        target_node = _row_col_adjacent_nodes_at_node[0]
                        [change_x, change_y] = [dist_to_exit + 1e-10, 0.0]

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
                            h_dist_to_exit = (
                                    (np.abs(self.x_of_corners_at_node[
                                            _node, corner] - _clast_x)) / (
                                                        np.cos(ss_azimuth)))
                            target_node = _row_col_adjacent_nodes_at_node[0]
# TO DO: Format to Pep8
                            [change_x, change_y] = [h_dist_to_exit * np.cos(ss_azimuth) + 1e-10, h_dist_to_exit * np.sin(ss_azimuth) + 1e-10]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = NE-col
                            print('NE-C')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(np.radians(90) - ss_azimuth)))
                            h_dist_to_exit = ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(np.radians(90) - ss_azimuth)))
                            target_node = _row_col_adjacent_nodes_at_node[1]
                            [change_x, change_y] = [h_dist_to_exit * np.sin(np.radians(90) - ss_azimuth) + 1e-10, h_dist_to_exit * np.cos(np.radians(90) - ss_azimuth) + 1e-10]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('NE-corner')
                            dist_to_exit = (1 / np.cos(ss_dip)) * np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            h_dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            target_node = _diagonal_adjacent_nodes_at_node[0]
                            [change_x, change_y] = [np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x) + 1e-10, np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y) + 1e-10]
                    elif sn_slope > 0 and we_slope < 0: # Quarter = NW
                        ss_azimuth = np.radians(90) + np.arctan(np.abs(we_slope / sn_slope))
                        corner = 1
                        clast_to_corner_azimuth = np.radians(180) - np.arctan(
                                np.abs(_clast_y - self.y_of_corners_at_node[_node, corner])/
                                       np.abs(_clast_x - self.x_of_corners_at_node[_node, corner]))
                        if ss_azimuth < clast_to_corner_azimuth: # Eigth = NW-col
                            print('NW-C')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(ss_azimuth - np.radians(90)))) #
                            h_dist_to_exit = ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(ss_azimuth - np.radians(90)))) #
                            [change_x, change_y] = [-h_dist_to_exit * (np.sin(ss_azimuth - np.radians(90))) - 1e-10, h_dist_to_exit * (np.cos(ss_azimuth - np.radians(90))) + 1e-10]
                            target_node = _row_col_adjacent_nodes_at_node[1]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = NW-row
                            print('NW-R')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x)) / (np.cos(np.radians(180) - ss_azimuth))) #
                            h_dist_to_exit = ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x)) / (np.cos(np.radians(180) - ss_azimuth))) #
                            [change_x, change_y] = [-h_dist_to_exit * (np.cos(np.radians(180) - ss_azimuth)) - 1e-10, h_dist_to_exit * (np.sin(np.radians(180) - ss_azimuth))+ 1e-10]
                            target_node = _row_col_adjacent_nodes_at_node[2]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('NW-corner')
                            dist_to_exit = (1 / np.cos(ss_dip)) * np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            h_dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            [change_x, change_y] = [-np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x) - 1e-10, np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y) + 1e-10]
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
                            h_dist_to_exit = ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x)) / (np.cos(ss_azimuth - np.radians(180)))) #
                            [change_x, change_y] = [-h_dist_to_exit * np.cos(ss_azimuth - np.radians(180)) - 1e-10, -h_dist_to_exit * np.sin(ss_azimuth - np.radians(180))- 1e-10]
                            target_node = _row_col_adjacent_nodes_at_node[2]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = SW-col
                            print('SW-C')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(np.radians(270) - ss_azimuth)))  #
                            h_dist_to_exit = ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(np.radians(270) - ss_azimuth)))  #
                            [change_x, change_y] = [-h_dist_to_exit * np.sin(np.radians(270) - ss_azimuth)- 1e-10, -h_dist_to_exit *np.cos(np.radians(270) - ss_azimuth) - 1e-10]
                            target_node = _row_col_adjacent_nodes_at_node[3]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('SW-corner')
                            dist_to_exit = (1 / np.cos(ss_dip)) * np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            h_dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            [change_x, change_y] = [-np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x) - 1e-10, -np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y) - 1e-10]
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
                            h_dist_to_exit = ((np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y)) / (np.cos(ss_azimuth - np.radians(270)))) #
                            [change_x, change_y] = [h_dist_to_exit * np.cos(ss_azimuth - np.radians(270)) + 1e-10, -h_dist_to_exit * np.sin(ss_azimuth - np.radians(270)) - 1e-10]
                            target_node = _row_col_adjacent_nodes_at_node[3]
                        elif ss_azimuth > clast_to_corner_azimuth: # Eigth = SE-row
                            print('SE-R')
                            dist_to_exit = (1 / np.cos(ss_dip)) * ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x) / (np.cos(np.radians(360) - ss_azimuth))))  #
                            h_dist_to_exit = ((np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x) / (np.cos(np.radians(360) - ss_azimuth))))  #
                            [change_x, change_y] = [dist_to_exit * np.cos(np.radians(360) - ss_azimuth) + 1e-10, -dist_to_exit *np.sin(np.radians(360) - ss_azimuth) - 1e-10]
                            target_node = _row_col_adjacent_nodes_at_node[0]
                        elif ss_azimuth == clast_to_corner_azimuth: # exit direction is diagonal
                            print('SE-corner')
                            dist_to_exit = (1 / np.cos(ss_dip)) * np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            h_dist_to_exit = np.sqrt(np.power(np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x), 2) + np.power(np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y), 2))
                            [change_x, change_y] = [np.abs(self.x_of_corners_at_node[_node, corner] - _clast_x) + 1e-10, -np.abs(self.y_of_corners_at_node[_node, corner] - _clast_y) - 1e-10]
                            target_node = _diagonal_adjacent_nodes_at_node[3]

        # Save values to dataframe:
        self['slope__steepest_azimuth'][clast] = ss_azimuth
        self['slope__steepest_dip'][clast] = ss_dip
        self['distance__to_exit'][clast] = dist_to_exit
        self['target_node'][clast] = target_node

        self['change_x'][clast] = change_x
        self['change_y'][clast] = change_y


        # Calculate lambda_0 and lambda_mean :
        self['lambda_0'][clast] = self.attrs['kappa'] / (2 * self['clast__radius'][clast, -1] * self.attrs['disturbance_fqcy'])
        # ((self._dt * self._kappa * _grid.dx) / (2 * df.at[clast,'clast__radius'])) * (np.tan(df.at[clast, 'slope__steepest_dip']) * ((self._Si-np.tan(df.at[clast, 'slope__steepest_dip']))/(self._Si+np.tan(df.at[clast, 'slope__steepest_dip']))))
                #(self._dt * self._kappa * _grid.dx) / (
#                        2 * df.at[clast,'clast__radius'])) #self._Si *


#
#        # Add randomness for lateral spreading (if option on):
#        if self.attrs['lateral_spreading'] == 'on':
#            if self['target_node'][clast] in _grid.boundary_nodes:
#                [delta_x,delta_y] = [0,0] # to avoid clast going out of the grid
#            else:
##                # Create truncated normal distribution
##                # (centered in 0; max, min = +/- lambda0/2)
##                x=truncnorm(-(self['lambda_0'][clast]/2), (
##                                    self['lambda_0'][clast]/2),
##                            loc=0,
##                            scale=(self['lambda_0'][clast]/4))
##                [delta_x,delta_y] = x.rvs(2)/10
##
##            self['change_x'][clast] = change_x + delta_x
##            self['change_y'][clast] = change_y + delta_y
##                # formerly: + np.random.normal(loc=0.0, scale=(df.at[clast, 'lambda_mean']/_grid.dx), size=1)   # adds a randomness for lateral spreading
#                max_deviation = np.pi * np.exp(self['slope__steepest_dip'][clast]/(-np.sqrt(self['lambda_0'][clast]))) * np.cos(self['slope__steepest_dip'][clast])
#                draw_dev=truncnorm((-max_deviation[1]/0.01),
#                                   max_deviation[1]/0.01,
#                                   loc=0.0,
#                                   scale=0.01).rvs(1)
#
#            self['change_x'][clast] = change_x + delta_x
#            self['change_y'][clast] = change_y + delta_y
#        else:
#            self['change_x'][clast] = change_x
#            self['change_y'][clast] = change_y

# For testing only: #####################
        print('change_x = %s' %change_x)
        print('change_y = %s' %change_y)
#########################################

        # Calculate lambda_0 and lambda_mean :
# If river, lambda0 multiplied:
        if ClastCollection._cell_is_hillslope(self, clast) == True:
            # Cell is hillslope:
            if np.tan(self['slope__steepest_dip'][clast]) >= self.attrs['Si']:
                lambda_0 = self['lambda_0'][clast]
                lambda_mean = np.power(10, 10)
            else:
                lambda_0  = self['lambda_0'][clast]
                lambda_mean = lambda_0 * (
                        self.attrs['Si'] + np.tan(
                                self['slope__steepest_dip'][clast])) / (
                                self.attrs['Si'] - np.tan(
                                    self['slope__steepest_dip'][clast]))
        else:   # Cell is river:
            if np.tan(self['slope__steepest_dip'][clast] >= self.attrs['Si']):
                lambda_0 = self['lambda_0'][clast]
                lambda_mean = np.power(10, 10)
            else:
                lambda_0 = 10 * self['lambda_0'][clast]
                lambda_mean = lambda_0 * (
                        self.attrs['Si'] + np.tan(
                                self['slope__steepest_dip'][clast])) / (
                                self.attrs['Si'] - np.tan(
                                    self['slope__steepest_dip'][clast]))

        self['lambda_mean'][clast] = lambda_mean
        self['lambda_0'][clast] = lambda_0




    def _move(self, clast):
            """ Move the clast along the path of steepest slope. First draws a
            travel distance from a probability distribution. If that travel
            distance is smaller than the distance to exit the cell in the
            direction of steepest slope, then move is handled by
            _move_in_cell. Otherwise, the clast is moved to the edge of the
            next downslope cell. The distance to travel, remaining from the
            original travel distance, is scaled as follows:
                new distance_to_travel = (lambda_mean in new cell) * (
                        distance remaining to travel) / (lambda_mean in
                                                    previous cell)
            This new distance to travel is compared to the distance to exit
            new cell (in the steepest slope direction of this new cell) and
            the clast is move accordingly, in cell or to next cell, re-
            iterating the same process until the distance left to travel is 0.

            Parameters
            ----------
            clast : int
                Clast ID
            """
            #self['hop_length'][clast][-1] = 0

            # cases where slope is null or sink:
            if np.isinf(self['distance__to_exit'][clast]):
                # clast does not move:
                pass
            else:
                # Draw a random sample in the probability distribution
                # of travel distances:
                self['distance__to_travel'][clast] = np.random.exponential(
                            scale= self['lambda_mean'][clast].values,
                            size=1)[0]
                print('dist_to_travel:')
                print(self['distance__to_travel'][clast].values)

                while self['distance__to_travel'][clast].values > 0.:
                    print('moving')

                    if ClastCollection.phantom(self, clast) == False:

                        if self['distance__to_travel'][clast] > (
                                self['distance__to_exit'][clast]):
                            print('MOVE OUT OF CELL')
                            print('from node: %s' %self['clast__node'][clast].values)
                            print('current x: %s' %self['clast__x'][clast, -1].values)
                            print('current y: %s' %self['clast__y'][clast, -1].values)
                            print('going to node: %s' %self['target_node'][clast].values)
                            # clast must change cell

                            # clast leaves cell,
                            # move of dist_to_exit along steepest slope:
                            print('applying changex %s' %self['change_x'][clast].values)
                            print('applying changey %s' %self['change_y'][clast].values)
                            self['clast__x'][clast, -1] += \
                                    self['change_x'][clast]
                            self['clast__y'][clast, -1] += \
                                    self['change_y'][clast]
                            print('new x: %s' %self['clast__x'][clast, -1].values)
                            print('new y: %s' %self['clast__y'][clast, -1].values)
                            self['hop_length'][clast, -1] += (
                                    np.sqrt(np.power(
                                            self['change_x'].values[clast], 2)\
                            + np.power(self['change_y'][clast],2)))\
                            / np.cos(self['slope__steepest_dip'][clast])
                            print('hop length: %s' %self['hop_length'][clast, -1])
                            self['clast__node'][clast] = self._grid.find_nearest_node((self['clast__x'][clast,-1],\
                                       self['clast__y'][clast,-1]))
                            print('new node: %s' %self['clast__node'][clast].values)
                            # check if phantom?
                            self['element_id'][clast] = (
                                    self._grid.cell_at_node[
                                    self['clast__node'][clast]])

                            # update distance to travel left after first hop:
                            self['distance__to_travel'][clast] -= (
                                    self['distance__to_exit'][clast])
                            print('distance left to travel %s' %self['distance__to_travel'][clast])
                            # keep memory of lambda_mean:
                            _lambda_mean_old = self['lambda_mean'][clast]

                            # update new neighborhood and travel info:
                            ClastCollection._neighborhood(self, clast)
                            ClastCollection._move_to(self, clast)

                            # Scale distance_to_travel to new lambda mean:
                            self['distance__to_travel'][clast] = (
                                    self['lambda_mean'][clast] * (
                                            self['distance__to_travel']\
                                            [clast])) / (_lambda_mean_old)
                            print('distance left to travel RESCALED %s' %self['distance__to_travel'][clast])


                        else: # move in cell, once
                            print('MOVE IN CELL')
                            ClastCollection._move_in_cell(self, clast)
                            break

                    else: # clast is phantom, does not move
                        print('clast has gone out of grid')
                        self['hop_length'][clast, -1] = 0.
                        break




    def _move_in_cell(self, clast):
            """ Move the clast within the cell, along the path of steepest
            slope, of a distance corresponding to the drawn travel distance if
            that travel distance drawn from the probability distribution is
            smaller than the distance to exit the cell.

            Parameters
            ----------
            clast : int
                Clast ID
            """
            # clast stays in cell, move of distance rand_length along slope
            ss_azimuth = self['slope__steepest_azimuth'][clast]
            ss_dip = self['slope__steepest_dip'][clast]
            x_horizontal = self['distance__to_travel'][clast] * np.cos(ss_dip)

            if np.isnan(ss_azimuth): # clast is in sink
                [change_x, change_y] = [0., 0.]
            elif ss_azimuth <= np.radians(90):
                [change_x, change_y] = [(
                        x_horizontal * np.cos(ss_azimuth)), (
                                x_horizontal * np.sin(ss_azimuth))]
            elif ss_azimuth <= np.radians(180):
                [change_x, change_y] = [(
                        -x_horizontal * np.cos(np.radians(180)-ss_azimuth)), (
                                x_horizontal * np.sin(
                                        np.radians(180)-ss_azimuth))]
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

            self['change_x'][clast] = change_x
            self['change_y'][clast] = change_y
            # Update clast coordinates:
            self['clast__x'][clast, -1] += change_x
            self['clast__y'][clast, -1] += change_y
            # Update hop length:
            self['hop_length'][clast, -1] += self['distance__to_travel'][clast]

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

        clast_node = self['clast__node'][clast].values
        clast_radius = self['clast__radius'][clast, -1]
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

        clast__node = self['clast__node'][clast]
        clast__elev = self['clast__elev'][clast, -2] #-2?
        topo__elev = self._grid.at_node['topographic__elevation'][clast__node]

        _detach = np.zeros(1, dtype=bool)


        _clast_depth = topo__elev - clast__elev


        proba_mobile = (1-np.exp(-self.attrs['disturbance_fqcy'] * self.attrs['dt'])) * (
                np.exp(-_clast_depth / self.attrs['d_star']))

        if proba_mobile >= 0.35:  ### VALUE TO SET #############################
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

        _node = self['clast__node'][clast]
        _grid = self._grid

        # If next to boundary, set A/S to 1:
        if np.isnan(self['slope__steepest_azimuth'][clast]) == True:
            area_over_slope = 1   # clast next to boundary
        else:
            area_over_slope = (_grid.at_node['drainage_area'][_node] /
                               _grid.at_node[
                                       'topographic__steepest_slope'][_node])

        if area_over_slope < self.attrs['threshold']:
            _cell_is_hillslope = True
        else:
            _cell_is_hillslope = False

        return _cell_is_hillslope

    def clast_solver_Exponential(self,
                                 dt=1.,
                                 Si=1.2,
                                 kappa=0.0001,
                                 uplift=None,
                                 #erosion_method='TLDiff',
                                 hillslope_river__threshold=1e4,
                                 lateral_spreading='off',
                                 disturbance_fqcy=1.,
                                 d_star=1.):
        """See :func:`run_one_step`

        """

        # Store values that will be used:
        self.attrs['kappa'] = kappa # if using SPACE, kappa=K_sed. TO DO: can be a field or array
        self.attrs['dt'] = dt
        #self._erosion__depth = self._erosion_rate * self._dt# NOT USED
#        self._deposition__thickness = self._deposition_rate * self._dt
        #  slope above which particle motion continues indefinetly (not equal
        # to critical slope of diffusion, see Furbish and Haff, 2010):
        self.attrs['Si'] = Si
        self.attrs['threshold'] = hillslope_river__threshold
        self.attrs['lateral_spreading'] = lateral_spreading
        self.attrs['disturbance_fqcy'] = disturbance_fqcy
        self.attrs['d_star'] = d_star


        # Uplift: #### TO MODIF (put in attrs?) AND TEST #########
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

        # Create new time coordinate for clast__x,y,z rad and hop_length (just
        # one needed, the others are filled with NaN automatically):
        self.add_record(time=[self.latest_time+self.attrs['dt']],
                        item_id=[0],
                        new_record={'clast__x' : (
                                ['item_id', 'time'], np.full(
                                        1, np.NaN).reshape(1,1))})
        # for now clasts don't change in size:
        # if self.clast_abrasion = False:
        self['clast__radius'] = self['clast__radius'].ffill('time')
        self['clast__x'] = self['clast__x'].ffill('time')
        self['clast__y'] = self['clast__y'].ffill('time')


        # Treat clasts one after the other:
        for clast in range(self._nb_of_clast):
            # Calculate lambda_0 (for display purposes mostly):
            self['lambda_0'][clast] = self.attrs['kappa'] / (2 * self['clast__radius'][clast, -1] * self.attrs['disturbance_fqcy'])

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


                    self['hop_length'][clast,-1] = 0.
                    self['distance__to_travel'][clast] = 0.

                    ClastCollection._move(self, clast)
                    print('moved')

                    # Update elevation:
                    # Clast is buried in the active layer with inv exp proba
                    self['clast__elev'][clast, -1] = self._grid.at_node['topographic__elevation'][self['clast__node'][clast].values] #- np.random.exponential(scale=(self.attrs['d_star']*self.attrs['disturbance_fqcy']/self['clast__radius'][clast].values[-1]), size=1)[0] #(self.attrs['d_star']/100 * (np.random.rand(1)))
                    print('elevation updated')
                    if hasattr(self, '_uplift') is True:   # uplift clast if necessary
                            self['clast__elev'][clast, -1] += \
                            self._uplift[self['clast__node'][clast]] * \
                            self.attrs['dt']
                            print('elevation updated with uplift')

                    # Update total travelled distance:
                    self['total_travelled_dist'][clast] += \
                        self['hop_length'][clast, -1]
                    # To implement: clast_abarasion
                    #if self.clast_abrasion = True:
                    # self['clast__radius'][clast][-1]=self.clast_abrasion(clast)
                else: # clast not detached, don't move:
                    self['hop_length'][:,-1] = 0
                    self['clast__x'] = self['clast__x'].ffill('time')
                    self['clast__y'] = self['clast__y'].ffill('time')
                    self['clast__elev'] = self['clast__elev'].ffill('time')
                    if hasattr(self, '_uplift') is True:   # uplift clast if necessary
                            self['clast__elev'][clast, -1] += \
                            self._uplift[self['clast__node'][clast]] * \
                            self.attrs['dt']
                            print('elevation updated with uplift')

            else: # clast is phantom, go to next clast
                # for display purpose: phantom clast has a radius of 0
                # self._clast__radius(i) = 0.
                self['hop_length'][clast, -1] = 0.
                pass


    def run_one_step(self, dt=1.,
                     Si=1.2,
                     kappa=0.0001,
                     uplift=None,
                     erosion_method='TLDiff',
                     hillslope_river__threshold=1e4,
                     lateral_spreading='off',
                     disturbance_fqcy=1.,
                     d_star=1.):

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
                                      lateral_spreading,
                                      disturbance_fqcy)