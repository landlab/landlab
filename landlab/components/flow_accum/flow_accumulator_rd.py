#!/usr/env/python

"""
flow_accumulator.py: Component to accumulate flow and calculate drainage area.

Provides the FlowAccumulator component which accumulates flow and calculates
drainage area. FlowAccumulator supports multiple methods for calculating flow
direction. Optionally a depression finding component can be specified and flow
directing, depression finding, and flow routing can all be accomplished
together.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component
from landlab import RasterModelGrid # for type tests
from landlab.utils.decorators import use_field_name_or_array
from landlab.core.messages import warning_message

from landlab.components.flow_accum import flow_accum_bw
from landlab.components.flow_accum import flow_accum_to_n

from landlab import BAD_INDEX_VALUE
import six
import numpy as np

import richdem as rp

@use_field_name_or_array('node')
def _return_surface(grid, surface):
    """
    Private function to return the surface to direct flow over.

    This function exists to take advantange of the 'use_field_name_or_array
    decorator which permits providing the surface as a field name or array.
    """
    return surface


class RdFlowAccumulator(Component):

    """
    Component to accumulate flow and calculate drainage area.

    This is accomplished by first finding flow directions by a user-specified
    method and then calculating the drainage area and discharge.

    Optionally, spatially variable runoff can be set either by the model grid
    field 'water__unit_flux_in' or the input variable *runoff_rate**.

    Optionally a depression finding component can be specified and flow
    directing, depression finding, and flow routing can all be accomplished
    together.


    NOTE: The perimeter nodes  NEVER contribute to the accumulating flux, even
    if the  gradients from them point inwards to the main body of the grid.
    This is because under Landlab definitions, perimeter nodes lack cells, so
    cannot accumulate any discharge.

    FlowAccumulator stores as ModelGrid fields:

        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array of all but the first element of the delta data structure:
            *flow__data_structure_delta*. The first element is always zero.
        -  Link array of the D data structure: *flow__data_structure_D*

    The FlowDirector component will add additional ModelGrid fields.
    DirectToOne methods(Steepest/D4 and D8) and DirectToMany(NAMES HERE) use
    different model grid fields.

    DirectToOne Methods (Steeptest/D4 and D8) store the following as ModelGrid
    fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_node'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Node array containing ID of link that leads from each node to its
           receiver, or BAD_INDEX_VALUE if no link:
           *'flow__link_to_receiver_node'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    DirectToMany Methods (MFD) store the following as ModelGrid
    fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_nodes'*. This array is 2D, and is
           of dimension (number of nodes x max number of receivers).
        -  Node array of flow proportions: *'flow__receiver_proportions'*. This
           array is 2D, and is of dimension (number of nodes x max number of
           receivers).
        -  Node array of links carrying flow:  *'flow__links_to_receiver_nodes'*.
           This array is 2D, and is of dimension (number of nodes x max number of
           receivers).
        -  Node array of the steepest downhill receiver. *'flow__receiver_nodes'*
        -  Node array of steepest downhill slope from each receiver:
           *'topographic__steepest_slope'*
        -  Node array containing ID of steepest link that leads from each node to a
           receiver, or BAD_INDEX_VALUE if no link:
           *'flow__link_to_receiver_node'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class is :func:`run_one_step`

    `run_one_step` takes the optional argument update_flow_director (default is
    True) that determines if the flow_director is re-run before flow is
    accumulated.

    Parameters
    ----------
    grid : ModelGrid
        A grid of type Voroni.
    surface : field name at node or array of length node
        The surface to direct flow across.
    fill_surface : should the surface be filled, or two surfaces maintained
# do same things as in DepressionFinder
    route_method : string
        A string of method. Options are:
        D4
        D8
        MFD-with-diagonals
        MFD-no-diagonals
        MFD-with-diagonals-kinematic-wave
        MFD-no-diagonals-kinematic-wave
        ...
        ...
        Dinf
    depression_method : string
        fill_depression
        breach_lake...
        ...
        ...

    runoff_rate : float, optional (m/time)
        If provided, sets the (spatially constant) runoff rate. If a spatially
        variable runoff rate is desired, use the input field
        'water__unit_flux_in'. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field.
        If neither are set, defaults to spatially constant unit input.
    **kwargs : any additional parameters to pass to ... e.g. epsilon,

    Examples
    --------


    """

    _name = 'RdFlowAccumulator'

    _input_var_names = ('topographic__elevation',
                        'water__unit_flux_in'
                        )

    _output_var_names = ('drainage_area',
                         'surface_water__discharge',
                         'flow__upstream_node_order',
                         'flow__nodes_not_in_stack',
                         'flow__data_structure_delta',
                         'flow__data_structure_D'
                         )

    _var_units = {'topographic__elevation': 'm',
                  'flow__receiver_node': 'm',
                  'water__unit_flux_in': 'm/s',
                  'drainage_area': 'm**2',
                  'surface_water__discharge': 'm**3/s',
                  'flow__upstream_node_order': '-',
                  'flow__data_structure_delta': '-',
                  'flow__data_structure_D': '-',
                  'flow__nodes_not_in_stack': '-'
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'flow__receiver_node': 'node',
                    'water__unit_flux_in': 'node',
                    'drainage_area': 'node',
                    'surface_water__discharge': 'node',
                    'flow__upstream_node_order': 'node',
                    'flow__nodes_not_in_stack': 'grid',
                    'flow__data_structure_delta': 'node',
                    'flow__data_structure_D': 'link',
                    }
    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'surface_water__discharge': 'Discharge of water through each node',
        'water__unit_flux_in':
            'External volume water per area per time input to each node '
            '(e.g., rainfall rate)',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'flow__data_structure_delta':
            'Node array containing the elements delta[1:] of the data '
            'structure "delta" used for construction of the downstream-to-'
            'upstream node array',
        'flow__data_structure_D':
            'Link array containing the data structure D used for construction'
            'of the downstream-to-upstream node array',
        'flow__nodes_not_in_stack':
            'Boolean value indicating if there are any nodes that have not yet'
            'been added to the stack stored in flow__upstream_node_order.'
               }

    def __init__(self,
                 grid,
                 surface='topographic__elevation',
                 filled_surface =
                 route_method='D4',
                 depression_method = 'fill_depression'
                 runoff_rate=None,
                 depression_finder=None,
                 **kwargs):
        """
        Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and compatability
        for the flow_director and depression_finder keyword arguments, tests
        the argument of runoff_rate, and initializes new fields.
        """
        super(RdFlowAccumulator, self).__init__(grid)
        # Keep a local reference to the grid
        self._grid = grid

        # Grid type testing
        if not isinstance(self._grid, RasterModelGrid):
            raise ValueError('RdFlowAccumulator only works with RasterModelGrids.')

        self.kwargs = kwargs
        # STEP 1: Testing of input values, supplied either in function call or
        # as part of the grid.
        self._test_water_inputs(grid, runoff_rate)

        # save elevations and node_cell_area to class properites.
        self.surface = surface
        self.surface_values = _return_surface(grid, surface)

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.

        self.node_cell_area = node_cell_area

        # STEP 2:
        # identify Flow Director method, save name, import and initialize the correct
        # flow director component if necessary
        self._add_director(flow_director)
        self._add_depression_finder(depression_finder)

        # check that all KWARGS are used.
        if len(self.kwargs)>0:
            raise ValueError('unused kwargs passed to FlowAccumulator: ',
                             self.kwargs)

        # This component will track of the following variables.
        # Attempt to create each, if they already exist, assign the existing
        # version to the local copy.

        #   - drainage area at each node
        #   - receiver of each node
        #   - D array
        #   - delta array
        #   - missing nodes in stack.
        try:
            self.drainage_area = grid.add_zeros('drainage_area', at='node',
                                                dtype=float)
        except FieldError:
            self.drainage_area = grid.at_node['drainage_area']

        try:
            self.discharges = grid.add_zeros('surface_water__discharge',
                                             at='node', dtype=float)
        except FieldError:
            self.discharges = grid.at_node['surface_water__discharge']

        try:
            self.upstream_ordered_nodes = grid.add_field('flow__upstream_node_order',
                                                         BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                                         at='node', dtype=int)

        except FieldError:
            self.upstream_ordered_nodes = grid.at_node[
                'flow__upstream_node_order']

        try:
            self.delta_structure = grid.add_field('flow__data_structure_delta',
                                                  BAD_INDEX_VALUE*grid.ones(at='node', dtype=int),
                                                  at='node', dtype=int)
        except FieldError:
            self.delta_structure = grid.at_node['flow__data_structure_delta']

        try:

            if self.flow_director.to_n_receivers == 'many' and self._is_raster:
                # needs to be BAD_INDEX_VALUE
                self.D_structure = grid.add_field('flow__data_structure_D',
                                                  BAD_INDEX_VALUE*np.ones((self._grid.number_of_links, 2),
                                                  dtype=int),
                                                  at='link',
                                                  dtype=int,
                                                  noclobber=False)
            else:

                # needs to be BAD_INDEX_VALUE
                self.D_structure = grid.add_field('flow__data_structure_D',
                                                  BAD_INDEX_VALUE*grid.ones(at='link'),
                                                  at='link', dtype=int)
        except FieldError:
            self.D_structure = grid.at_link['flow__data_structure_D']

        self.nodes_not_in_stack = True

    @property
    def node_drainage_area(self):
        """Return the drainage area."""
        return self._grid['node']['drainage_area']

    @property
    def node_water_discharge(self):
        """Return the surface water discharge."""
        return self._grid['node']['surface_water__discharge']

    @property
    def node_order_upstream(self):
        """Return the upstream node order (drainage stack)."""
        return self._grid['node']['flow__upstream_node_order']

    def _test_water_inputs(self, grid, runoff_rate):
        """Test inputs for runoff_rate and water__unit_flux_in."""
        # testing input for runoff rate, can be None, a string associated with
        # a field at node, a single float or int, or an array of size number of
        # nodes.
        if runoff_rate is not None:
            if isinstance(runoff_rate, str):
                runoff_rate = grid.at_node[runoff_rate]
            elif isinstance(runoff_rate, (float, int)):
                pass
            else:
                assert runoff_rate.size == grid.number_of_nodes

        # test for water__unit_flux_in
        try:
            grid.at_node['water__unit_flux_in']
        except FieldError:
            if runoff_rate is None:
                # assume that if runoff rate is not supplied, that the value
                # should be set to one everywhere.
                grid.add_ones('node', 'water__unit_flux_in', dtype=float)
            else:
                if isinstance(runoff_rate, (float, int)):
                    grid.add_empty('node', 'water__unit_flux_in', dtype=float)
                    grid.at_node['water__unit_flux_in'].fill(runoff_rate)
                else:
                    grid.at_node['water__unit_flux_in'] = runoff_rate
        else:
            if runoff_rate is not None:
                print ("RdFlowAccumulator found both the field " +
                       "'water__unit_flux_in' and a provided float or " +
                       "array for the runoff_rate argument. THE FIELD IS " +
                       "BEING OVERWRITTEN WITH THE SUPPLIED RUNOFF_RATE!")
                if isinstance(runoff_rate, (float, int)):
                    grid.at_node['water__unit_flux_in'].fill(runoff_rate)
                else:
                    grid.at_node['water__unit_flux_in'] = runoff_rate


    def accumulate_flow(self, update_flow_director=True):
        """
        Function to make FlowAccumulator calculate drainage area and discharge.

        Running run_one_step() results in the following to occur:
            1. Flow directions are updated (unless update_flow_director is set
            as False).
            2. Intermediate steps that analyse the drainage network topology
            and create datastructures for efficient drainage area and discharge
            calculations.
            3. Calculation of drainage area and discharge.
            4. Depression finding and mapping, which updates drainage area and
            discharge.
        """
        # step 1. Find flow directions by specified method
        if update_flow_director == True:
            self.flow_director.run_one_step()

        # further steps vary depending on how many recievers are present
        # one set of steps is for route to one (D8, Steepest/D4)
        if self.flow_director.to_n_receivers == 'one':

            # step 3. Run depression finder if passed
            # Depression finder reaccumulates flow at the end of its routine.
            if self.depression_finder_provided is not None:

                self.depression_finder.map_depressions()

                a = self._grid['node']['drainage_area']
                q = self._grid['node']['surface_water__discharge']

            else:
                # step 2. Get r
                r = self._grid['node']['flow__receiver_node']

                # step 2. Stack, D, delta construction
                nd = flow_accum_bw._make_number_of_donors_array(r)
                delta = flow_accum_bw._make_delta_array(nd)
                D = flow_accum_bw._make_array_of_donors(r, delta)
                s = flow_accum_bw.make_ordered_node_array(r)

                # put theese in grid so that depression finder can use it.
                # store the generated data in the grid
                self._grid['node']['flow__data_structure_delta'][:] = delta[1:]
                self._grid['link']['flow__data_structure_D'][:len(D)] = D
                self._grid['node']['flow__upstream_node_order'][:] = s

                # step 4. Accumulate (to one or to N depending on direction method. )
                a, q = flow_accum_bw.find_drainage_area_and_discharge(s,
                                                                      r,
                                                                      self.node_cell_area,
                                                                      self._grid.at_node['water__unit_flux_in'])
                self._grid['node']['drainage_area'][:] = a
                self._grid['node']['surface_water__discharge'][:] = q

        else:
            # step 2. Get r and p
            r = self._grid['node']['flow__receiver_nodes']
            p = self._grid['node']['flow__receiver_proportions']

            # step 2. Stack, D, delta construction
            nd = flow_accum_to_n._make_number_of_donors_array_to_n(r, p)
            delta = flow_accum_to_n._make_delta_array_to_n(nd)
            D = flow_accum_to_n._make_array_of_donors_to_n(r, p, delta)
            s = flow_accum_to_n.make_ordered_node_array_to_n(r, p)

            # put theese in grid so that depression finder can use it.
            # store the generated data in the grid
            self._grid['node']['flow__data_structure_delta'][:] = delta[1:]

            if self._is_raster:
                tempD = BAD_INDEX_VALUE * np.ones((self._grid.number_of_links*2))
                tempD[:len(D)] = D
                self._grid['link']['flow__data_structure_D'][:] = tempD.reshape((self._grid.number_of_links, 2))
            else:
                self._grid['link']['flow__data_structure_D'][:len(D)] = D
            self._grid['node']['flow__upstream_node_order'][:] = s

            # step 3. Run depression finder if passed
            # at present this must go at the end.

            # step 4. Accumulate (to one or to N depending on direction method. )
            a, q = flow_accum_to_n.find_drainage_area_and_discharge_to_n(s,
                                                                         r,
                                                                         p,
                                                                         self.node_cell_area,
                                                                         self._grid.at_node['water__unit_flux_in'])
            # store drainage area and discharge.
            self._grid['node']['drainage_area'][:] = a
            self._grid['node']['surface_water__discharge'][:] = q

            # at the moment, this is where the depression finder needs to live.
            if self.depression_finder_provided is not None:
                self.depression_finder.map_depressions()

        return (a, q)

    def run_one_step(self):
        """
        Accumulate flow and save to the model grid.

        run_one_step() checks for updated boundary conditions, calculates
        slopes on links, finds baselevel nodes based on the status at node,
        calculates flow directions, and accumulates flow and saves results to
        the grid.

        An alternative to run_one_step() is accumulate_flow() which does the
        same things but also returns the drainage area and discharge.
        """
        self.accumulate_flow()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
