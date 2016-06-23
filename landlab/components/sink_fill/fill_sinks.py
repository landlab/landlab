# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19.

@author: dejh
"""
from __future__ import print_function

import landlab
from landlab import (ModelParameterDictionary, Component, FieldError,
                     FIXED_VALUE_BOUNDARY)

from landlab.utils.decorators import use_file_name_or_kwds, deprecated
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.components.flow_routing import (DepressionFinderAndRouter,
                                             FlowRouter)
from landlab.grid.base import BAD_INDEX_VALUE
import numpy as np


class SinkFiller(Component):
    """
    This component identifies depressions in a topographic surface, then fills
    them in in the topography.  No attempt is made to conserve sediment mass.
    User may specify whether the holes should be filled to flat, or with a
    gradient downwards towards the depression outlet. The gradient can be
    spatially variable, and is chosen to not reverse any drainage directions
    at the perimeter of each lake.

    The primary method of this class is 'run_one_step'. 'fill_pits' is a
    synonym.

    Constructor assigns a copy of the grid, and calls the initialize
    method.

    Construction::

        SinkFiller(grid, routing='D8', apply_slope=False, fill_slope=1.e-5):

    Parameters
    ----------
    grid : ModelGrid
        A landlab grid.
    routing : {'D8', 'D4'} (optional)
        If grid is a raster type, controls whether fill connectivity can
        occur on diagonals ('D8', default), or only orthogonally ('D4').
        Has no effect if grid is not a raster.
    apply_slope : bool
        If False (default), leave the top of the filled sink flat. If True,
        apply the slope fill_slope to the top surface to allow subsequent flow
        routing. A test is performed to ensure applying this slope will not
        alter the drainage structure at the edge of the filled region
        (i.e., that we are not accidentally reversing the flow direction
        far from the outlet.)
    fill_slope : float (m/m)
        The slope added to the top surface of filled pits to allow flow
        routing across them, if apply_slope.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab import BAD_INDEX_VALUE as XX
    >>> from landlab.components import FlowRouter, SinkFiller
    >>> import numpy as np
    >>> lake1 = np.array([34, 35, 36, 44, 45, 46, 54, 55, 56, 65, 74])
    >>> lake2 = np.array([78, 87, 88])
    >>> guard_nodes = np.array([23, 33, 53, 63, 73, 83])
    >>> lake = np.concatenate((lake1, lake2))
    >>> mg = RasterModelGrid((10, 10), 1.)
    >>> z = np.ones(100, dtype=float)
    >>> z += mg.node_x  # add a slope
    >>> z[guard_nodes] += 0.001  # forces the flow out of a particular node
    >>> z[lake] = 0.
    >>> field = mg.add_field('node', 'topographic__elevation', z,
    ...                      units='-', copy=True)
    >>> fr = FlowRouter(mg)
    >>> fr.run_one_step()
    >>> mg.at_node['flow__sink_flag'][mg.core_nodes].sum()
    14
    >>> hf = SinkFiller(mg, apply_slope=False)
    >>> hf.run_one_step()
    >>> np.allclose(mg.at_node['topographic__elevation'][lake1], 4.)
    True
    >>> np.allclose(mg.at_node['topographic__elevation'][lake2], 7.)
    True

    Now reset and demonstrate the adding of an inclined surface:

    >>> field[:] = z
    >>> hf = SinkFiller(mg, apply_slope=True)
    >>> hf.run_one_step()
    >>> hole1 = np.array([4.00007692, 4.00015385, 4.00023077, 4.00053846,
    ...                   4.00038462, 4.00030769, 4.00069231, 4.00061538,
    ...                   4.00046154, 4.00076923, 4.00084615])
    >>> hole2 = np.array([7.4, 7.2, 7.6])
    >>> np.allclose(mg.at_node['topographic__elevation'][lake1], hole1)
    True
    >>> np.allclose(mg.at_node['topographic__elevation'][lake2], hole2)
    True
    >>> fr.run_one_step()
    >>> mg.at_node['flow__sink_flag'][mg.core_nodes].sum()
    0
    """
    _name = 'SinkFiller'

    _input_var_names = ('topographic__elevation',
                        )

    _output_var_names = ('topographic__elevation',
                         'sediment_fill__depth',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'sediment_fill__depth': 'm',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'sediment_fill__depth': 'node',
                    }

    _var_doc = {'topographic__elevation': 'Surface topographic elevation',
                'sediment_fill__depth': 'Depth of sediment added at each' +
                                        'node',
                }

    @use_file_name_or_kwds
    def __init__(self, grid, routing='D8', apply_slope=False,
                 fill_slope=1.e-5, **kwds):
        self._grid = grid
        if routing is not 'D8':
            assert routing is 'D4'
        self._routing = routing
        if ((type(self._grid) is landlab.grid.raster.RasterModelGrid) and
                (routing is 'D8')):
            self._D8 = True
            self.num_nbrs = 8
        else:
            self._D8 = False  # useful shorthand for thia test we do a lot
            if type(self._grid) is landlab.grid.raster.RasterModelGrid:
                self.num_nbrs = 4
        self._fill_slope = fill_slope
        self._apply_slope = apply_slope
        self.initialize()

    def initialize(self, input_stream=None):
        """
        The BMI-style initialize method takes an optional input_stream
        parameter, which may be either a ModelParameterDictionary object or
        an input stream from which a ModelParameterDictionary can read values.
        """
        # Create a ModelParameterDictionary for the inputs
        if input_stream is None:
            inputs = None
        elif type(input_stream) == ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)

        # Make sure the grid includes elevation data. This means either:
        #  1. The grid has a node field called 'topographic__elevation', or
        #  2. The input file has an item called 'ELEVATION_FIELD_NAME' *and*
        #     a field by this name exists in the grid.
        try:
            self._elev = self._grid.at_node['topographic__elevation']
        except FieldError:
            try:
                self.topo_field_name = inputs.read_string('ELEVATION_' +
                                                          'FIELD_NAME')
            except AttributeError:
                print('Error: Because your grid does not have a node field')
                print('called "topographic__elevation", you need to pass the')
                print('name of a text input file or ModelParameterDictionary,')
                print('and this file or dictionary needs to include the name')
                print('of another field in your grid that contains your')
                print('elevation data.')
                raise AttributeError
            except MissingKeyError:
                print('Error: Because your grid does not have a node field')
                print('called "topographic__elevation", your input file (or')
                print('ModelParameterDictionary) must include an entry with')
                print('the key "ELEVATION_FIELD_NAME", which gives the name')
                print('of a field in your grid that contains your elevation')
                print('data.')
                raise MissingKeyError('ELEVATION_FIELD_NAME')
            try:
                self._elev = self._grid.at_node[self.topo_field_name]
            except AttributeError:
                print('Your grid does not seem to have a node field called',
                      self.topo_field_name)
        else:
            self.topo_field_name = 'topographic__elevation'
        # create the only new output field:
        self.sed_fill_depth = self._grid.add_zeros('node',
                                                   'sediment_fill__depth',
                                                   noclobber=False)

        self._lf = DepressionFinderAndRouter(self._grid, routing=self._routing)
        self._fr = FlowRouter(self._grid, method=self._routing)

    def fill_pits(self, **kwds):
        """
        This is a synonym for the main method :func:`run_one_step`.
        """
        self.run_one_step(**kwds)

    def run_one_step(self, **kwds):
        """
        This is the main method. Call it to fill depressions in a starting
        topography.
        """
        # added for back-compatibility with old formats
        try:
            self._apply_slope = kwds['apply_slope']
        except KeyError:
            pass
        self.original_elev = self._elev.copy()
        # We need this, as we'll have to do ALL this again if we manage
        # to jack the elevs too high in one of the "subsidiary" lakes.
        # We're going to implement the lake_mapper component to do the heavy
        # lifting here, then delete its fields. This means we first need to
        # test if these fields already exist, in which case, we should *not*
        # delete them!
        existing_fields = {}
        spurious_fields = set()
        set_of_outputs = (set(self._lf.output_var_names) |
                          set(self._fr.output_var_names))
        try:
            set_of_outputs.remove(self.topo_field_name)
        except KeyError:
            pass
        for field in set_of_outputs:
            try:
                existing_fields[field] = self._grid.at_node[field].copy()
            except FieldError:  # not there; good!
                spurious_fields.add(field)

        self._fr.route_flow()
        self._lf.map_depressions(pits=self._grid.at_node['flow__sink_flag'],
                                 reroute_flow=True)
        # add the depression depths to get up to flat:
        self._elev += self._grid.at_node['depression__depth']
        # if apply_slope is none, we're now done! But if not...
        if self._apply_slope:
            # new way of doing this - use the upstream structure! Should be
            # both more general and more efficient
            for (outlet_node, lake_code) in zip(self._lf.lake_outlets,
                                                self._lf.lake_codes):
                lake_nodes = np.where(self._lf.lake_map == lake_code)[0]
                lake_perim = self._get_lake_ext_margin(lake_nodes)
                perim_elevs = self._elev[lake_perim]
                out_elev = self._elev[outlet_node]
                lowest_elev_perim = perim_elevs[perim_elevs != out_elev].min()
                # note we exclude the outlet node
                elev_increment = ((lowest_elev_perim-self._elev[outlet_node]) /
                                  (lake_nodes.size + 2.))
                assert elev_increment > 0.
                all_ordering = self._grid.at_node['flow__upstream_node_order']
                upstream_order_bool = np.in1d(all_ordering, lake_nodes,
                                              assume_unique=True)
                lake_upstream_order = all_ordering[upstream_order_bool]
                argsort_lake = np.argsort(lake_upstream_order)
                elevs_to_add = (np.arange(lake_nodes.size, dtype=float) +
                                1.) * elev_increment
                sorted_elevs_to_add = elevs_to_add[argsort_lake]
                self._elev[lake_nodes] += sorted_elevs_to_add
        # now put back any fields that were present initially, and wipe the
        # rest:
        for delete_me in spurious_fields:
            self._grid.delete_field('node', delete_me)
        for update_me in existing_fields.keys():
            self.grid.at_node[update_me][:] = existing_fields[update_me]
        # fill the output field
        self.sed_fill_depth[:] = self._elev - self.original_elev

    @deprecated(use='fill_pits', version=1.0)
    def _fill_pits_old(self, apply_slope=None):
        """

        .. deprecated:: 0.1.38
            Use :func:`fill_pits` instead.

        This is the main method. Call it to fill depressions in a starting
        topography.

        **Output fields**

        *  `topographic__elevation` : the updated elevations
        *  `sediment_fill__depth` : the depth of sediment added at each node

        Parameters
        ----------
        apply_slope : None, bool, or float
            If a float is provided this is the slope of the surface down
            towards the lake outlet. Supply a small positive number, e.g.,
            1.e-5 (or True, to use this default value).
            A test is performed to ensure applying this slope will not alter
            the drainage structure at the edge of the filled region (i.e.,
            that we are not accidentally reversing the flow direction far
            from the outlet.) The component will automatically decrease the
            (supplied or default) gradient a number of times to try to
            accommodate this, but will eventually raise an OverflowError
            if it can't deal with it. If you pass True, the method will use
            the default value of 1.e-5.
        """
        self.original_elev = self._elev.copy()
        # We need this, as we'll have to do ALL this again if we manage
        # to jack the elevs too high in one of the "subsidiary" lakes.
        # We're going to implement the lake_mapper component to do the heavy
        # lifting here, then delete its fields. This means we first need to
        # test if these fields already exist, in which case, we should *not*
        # delete them!
        existing_fields = {}
        spurious_fields = set()
        set_of_outputs = self._lf.output_var_names | self._fr.output_var_names
        try:
            set_of_outputs.remove(self.topo_field_name)
        except KeyError:
            pass
        for field in set_of_outputs:
            try:
                existing_fields[field] = self._grid.at_node[field].copy()
            except FieldError:  # not there; good!
                spurious_fields.add(field)

        self._fr.route_flow()
        self._lf.map_depressions(pits=self._grid.at_node['flow__sink_flag'],
                                 reroute_flow=False)
        # add the depression depths to get up to flat:
        self._elev += self._grid.at_node['depression__depth']
        # if apply_slope is none, we're now done! But if not...
        if apply_slope is True:
            apply_slope = self._fill_slope
        elif type(apply_slope) in (float, int):
            assert apply_slope >= 0.
        if apply_slope:
            # this isn't very efficient, but OK as we're only running this
            # code ONCE in almost all use cases
            sublake = False
            unstable = True
            stability_increment = 0
            self.lake_nodes_treated = np.array([], dtype=int)
            while unstable:
                while 1:
                    for (outlet_node, lake_code) in zip(self._lf.lake_outlets,
                                                        self._lf.lake_codes):
                        self._apply_slope_current_lake(apply_slope,
                                                       outlet_node,
                                                       lake_code, sublake)
                    # Call the mapper again here. Bail out if no core pits are
                    # found.
                    # This is necessary as there are some configs where adding
                    # the slope could create subsidiary pits in the topo
                    self._lf.map_depressions(pits=None, reroute_flow=False)
                    if len(self._lf.lake_outlets) == 0.:
                        break
                    self._elev += self._grid.at_node['depression__depth']
                    sublake = True
                    self.lake_nodes_treated = np.array([], dtype=int)
                # final test that all lakes are not reversing flow dirs
                all_lakes = np.where(self._lf.flood_status <
                                     BAD_INDEX_VALUE)[0]
                unstable = self.drainage_directions_change(all_lakes,
                                                           self.original_elev,
                                                           self._elev)
                if unstable:
                    apply_slope *= 0.1
                    sublake = False
                    self.lake_nodes_treated = np.array([], dtype=int)
                    self._elev[:] = original_elev  # put back init conds
                    stability_increment += 1
                    if stability_increment == 10:
                        raise OverflowError('Filler could not find a stable ' +
                                            'condition with a sloping ' +
                                            'surface!')
        # now put back any fields that were present initially, and wipe the
        # rest:
        for delete_me in spurious_fields:
            self._grid.delete_field('node', delete_me)
        for update_me in existing_fields.keys():
            self.grid.at_node[update_me] = existing_fields[update_me]
        # fill the output field
        self.sed_fill_depth[:] = self._elev - self.original_elev

    def _add_slopes(self, slope, outlet_node, lake_code):
        """
        Assuming you have already run the lake_mapper, adds an incline towards
        the outlet to the nodes in the lake.
        """
        new_elevs = self._elev.copy()
        outlet_coord = (self._grid.node_x[outlet_node],
                        self._grid.node_y[outlet_node])
        lake_nodes = np.where(self._lf.lake_map == lake_code)[0]
        lake_nodes = np.setdiff1d(lake_nodes, self.lake_nodes_treated)
        lake_ext_margin = self._get_lake_ext_margin(lake_nodes)
        d = self._grid.calc_distances_of_nodes_to_point(outlet_coord,
                                                       node_subset=lake_nodes)
        add_vals = slope*d
        new_elevs[lake_nodes] += add_vals
        self.lake_nodes_treated = np.union1d(self.lake_nodes_treated,
                                             lake_nodes)
        return new_elevs, lake_nodes

    def _get_lake_ext_margin(self, lake_nodes):
        """
        Returns the nodes forming the external margin of the lake, honoring
        the *routing* method (D4/D8) if applicable.
        """
        if self._D8 is True:
            all_poss = np.union1d(self._grid.active_neighbors_at_node(lake_nodes),
                                  self._grid._get_diagonal_list(lake_nodes))
        else:
            all_poss = np.unique(self._grid.active_neighbors_at_node(
                lake_nodes))
        lake_ext_edge = np.setdiff1d(all_poss, lake_nodes)
        return lake_ext_edge[lake_ext_edge != BAD_INDEX_VALUE]

    def _get_lake_int_margin(self, lake_nodes, lake_ext_edge):
        """
        Returns the nodes forming the internal margin of the lake, honoring
        the *routing* method (D4/D8) if applicable.
        """
        lee = lake_ext_edge
        if self._D8 is True:
            all_poss_int = np.union1d(self._grid.active_neighbors_at_node(lee),
                                      self._grid._get_diagonal_list(lee))
        else:
            all_poss_int = np.unique(self._grid.active_neighbors_at_node(lee))
        lake_int_edge = np.intersect1d(all_poss_int, lake_nodes)
        return lake_int_edge[lake_int_edge != BAD_INDEX_VALUE]

    def _apply_slope_current_lake(self, apply_slope, outlet_node, lake_code,
                                 sublake):
        """
        Wraps the _add_slopes method to allow handling of conditions where the
        drainage structure would be changed or we're dealing with a sublake.
        """
        while 1:
            starting_elevs = self._elev.copy()
            self._elev[:], lake_nodes = self._add_slopes(apply_slope,
                                                         outlet_node,
                                                         lake_code)
            ext_edge = self._get_lake_ext_margin(lake_nodes)
            if sublake:
                break
            else:
                if not self.drainage_directions_change(lake_nodes,
                                                       starting_elevs,
                                                       self._elev):
                    break
                else:
                    # put the elevs back...
                    self._elev[lake_nodes] = starting_elevs[lake_nodes]
                    # the slope was too big. Reduce it.
                    apply_slope *= 0.1
        # if we get here, either sublake, or drainage dirs are stable

    def drainage_directions_change(self, lake_nodes, old_elevs, new_elevs):
        """
        True if the drainage structure at lake margin changes, False otherwise.
        """
        ext_edge = self._get_lake_ext_margin(lake_nodes)
        if self._D8:
            edge_neighbors = np.hstack((self._grid.active_neighbors_at_node(ext_edge),
                                        self._grid._get_diagonal_list(
                                            ext_edge)))
        else:
            edge_neighbors = self._grid.active_neighbors_at_node(ext_edge).copy()
        edge_neighbors[edge_neighbors == BAD_INDEX_VALUE] = -1
        # ^value irrelevant
        old_neighbor_elevs = old_elevs[edge_neighbors]
        new_neighbor_elevs = new_elevs[edge_neighbors]
        # enforce the "don't change drainage direction" condition:
        edge_elevs = old_elevs[ext_edge].reshape((ext_edge.size, 1))
        cond = np.allclose((edge_elevs >= old_neighbor_elevs),
                           (edge_elevs >= new_neighbor_elevs))
        # if True, we're good, the tilting didn't mess with the fr
        return not cond
