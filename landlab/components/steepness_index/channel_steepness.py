# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19.

@author: dejh
"""
from __future__ import print_function
from six.moves import range  # this is Python 3's generator, not P2's list

import landlab
from landlab import ModelParameterDictionary, Component, FieldError, \
                    FIXED_VALUE_BOUNDARY
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.components.flow_routing.lake_mapper import \
    DepressionFinderAndRouter
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.grid.base import BAD_INDEX_VALUE
from landlab.utils.decorators import use_file_name_or_kwds
import numpy as np


class SteepnessFinder(Component):
    """
    This component calculates steepness indices, sensu Wobus et al. 2006,
    for a Landlab landscape.
    Follows broadly the approach used in GeomorphTools, geomorphtools.org.

    Construction::

        SteepnessFinder(grid, reference_concavity=0.5, min_drainage_area=1.e6,
                        elev_step=0., discretization_length=0.)

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    reference_concavity : float
        The reference concavity to use in the calculation.
    min_drainage_area : float (m**2; default 1.e6)
        The minimum drainage area above which steepness indices are 
        calculated.
        Defaults to 1.e6 m**2, per Wobus et al. 2006.
    elev_step : float (m; default 0.)
        If >0., becomes a vertical elevation change step to use to
        discretize the data (per Wobus). If 0., all nodes are used and
        no discretization happens.
    discretization_length : float (m; default 0.)
        If >0., becomes the lengthscale over which to segment the profiles -
        i.e., one different steepness index value is calculated every
        discretization_length. If only one (or no) points are present in a
        segment, it will be lumped together with the next segment.
        If zero, one value is assigned to each channel node.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
    >>> from landlab.components import FlowRouter, FastscapeEroder
    >>> mg = RasterModelGrid((3, 10), (100., 100.))
    >>> for nodes in (mg.nodes_at_right_edge, mg.nodes_at_bottom_edge,
    ...               mg.nodes_at_top_edge):
    ...     mg.status_at_node[nodes] = CLOSED_BOUNDARY
    >>> _ = mg.add_zeros('node', 'topographic__elevation')
    >>> mg.at_node['topographic__elevation'][mg.core_nodes] = mg.node_x[
    ...     mg.core_nodes]/1000.
    >>> fr = FlowRouter(mg)
    >>> sp = FastscapeEroder(mg, K_sp=0.01)
    >>> sf = SteepnessFinder(mg, min_drainage_area=10000.)
    >>> for i in range(10):
    ...     mg.at_node['topographic__elevation'][mg.core_nodes] += 10.
    ...     _ = fr.route_flow()
    ...     sp.run_one_step(1000.)
    >>> sf.calculate_steepnesses()
    >>> mg.at_node['channel__steepness_index'].reshape((3, 10))[1, :]
    array([  0.        ,  29.28427125,   1.        ,   1.        ,
             1.        ,   1.        ,   1.        ,   1.        ,
             0.99999997,   0.        ])
    >>> sf.hillslope_mask
    array([ True,  True,  True,  True,  True,  True,  True,  True,  True,
            True, False, False, False, False, False, False, False, False,
           False,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True], dtype=bool)

    >>> sf.calculate_steepnesses(discretization_length=350.)
    >>> mg.at_node['channel__steepness_index'].reshape((3, 10))[1, :]
    array([ 0.        ,  3.08232295,  3.08232295,  3.08232295,  1.        ,
            1.        ,  1.        ,  1.        ,  0.        ,  0.        ])

    >>> sf.calculate_steepnesses(elev_step=1.5)
    >>> mg.at_node['channel__steepness_index'].reshape((3, 10))[1, :]
    array([ 0.        ,  1.22673541,  1.2593727 ,  1.27781936,  1.25659369,
            1.12393156,  0.97335328,  0.79473963,  0.56196578,  0.        ])

    """
    _name = 'SteepnessFinder'

    _input_var_names = (
        'topographic__elevation',
        'drainage_area',
        'topographic__steepest_slope',
        'flow__receiver_node',
        'flow__upstream_node_order',
        'flow__link_to_receiver_node',
    )

    _output_var_names = (
        'channel__steepness_index',
    )

    _var_units = {'topographic__elevation': 'm',
                  'drainage_area': 'm**2',
                  'topographic__steepest_slope': '-',
                  'flow__receiver_node': '-',
                  'flow__upstream_node_order': '-',
                  'flow__link_to_receiver_node': '-',
                  'channel__steepness_index': 'variable',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'drainage_area': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow__receiver_node': 'node',
                    'flow__upstream_node_order': 'node',
                    'flow__link_to_receiver_node': 'node',
                    'channel__steepness_index': 'node',
                    }

    _var_doc = {'topographic__elevation': 'Surface topographic elevation',
                'drainage_area': 'upstream drainage area',
                'topographic__steepest_slope': ('the steepest downslope ' +
                                                'rise/run leaving the node'),
                'flow__receiver_node': ('the downstream node at the end of the ' +
                                  'steepest link'),
                'flow__upstream_node_order': ('node order such that nodes must ' +
                                        'appear in the list after all nodes ' +
                                        'downstream of them'),
                'flow__link_to_receiver_node':
                    ('ID of link downstream of each node, which carries the ' +
                     'discharge'),
                'channel__steepness_index': 'the local steepness index',
                }

    def __init__(self, grid, reference_concavity=0.5, min_drainage_area=1.e6,
                 elev_step=0., discretization_length=0., **kwds):
        """
        Constructor for the component.
        """
        self._grid = grid
        self._reftheta = reference_concavity
        self.min_drainage = min_drainage_area
        assert elev_step >= 0., "elev_step must be >= 0!"
        self._elev_step = elev_step
        self._discretization = discretization_length
        self.ksn = self._grid.add_zeros('node', 'channel__steepness_index')
        self._mask = self.grid.ones('node', dtype=bool)
        # this one needs modifying if smooth_elev
        self._elev = self.grid.at_node['topographic__elevation']

    def calculate_steepnesses(self, **kwds):
        """
        This is the main method. Call it to calculate local steepness indices
        at all points with drainage areas greater than *min_drainage_area*.

        This "run" method can optionally take the same parameter set as
        provided at instantiation. If they are provided, they will override
        the existing values from instantiation.

        Normalized steepness of any node without a defined value is reported
        as 0. These nodes are also identified in the mask retrieved with
        :func:`hillslope_mask`.
        """
        self._mask.fill(True)
        self.ksn.fill(0.)
        # test for new kwds:
        reftheta = kwds.get('reference_concavity', self._reftheta)
        min_drainage = kwds.get('min_drainage_area', self.min_drainage)
        elev_step = kwds.get('elev_step', self._elev_step)
        discretization_length = kwds.get('discretization_length',
                                         self._discretization)

        upstr_order = self.grid.at_node['flow__upstream_node_order']
        # get an array of only nodes with A above threshold:
        valid_dstr_order = (upstr_order[self.grid.at_node['drainage_area'][
            upstr_order] >= min_drainage])[::-1]
        valid_dstr_elevs = self.grid.at_node['topographic__elevation'][
            valid_dstr_order]
        valid_dstr_areas = self.grid.at_node['drainage_area'][valid_dstr_order]
        # note elevs are guaranteed to be in order, UNLESS a fill
        # algorithm has been used.
        nodes_incorporated = self.grid.zeros('node', dtype=bool)
        # now do each poss channel in turn
        # get the head of the first (longest!) channel:
        for dstr_order_index in range(valid_dstr_order.size):
            this_ch_top_node = valid_dstr_order[dstr_order_index]  # top node
            if not nodes_incorporated[this_ch_top_node]:
                nodes_incorporated[this_ch_top_node] = True
                nodes_in_channel = [this_ch_top_node, ]
                penultimate_node = this_ch_top_node
                current_node_incorporated = False
                while not current_node_incorporated:
                    next_node = self.grid.at_node['flow__receiver_node'][
                        penultimate_node]
                    if next_node == penultimate_node:  # end of flow path
                        break
                    nodes_in_channel.append(next_node)
                    current_node_incorporated = nodes_incorporated[next_node]
                    # ^ this is a COPY op, so we're free to update the array
                    nodes_incorporated[next_node] = True
                    penultimate_node = next_node
                # by here, we have a full, unique reach in nodes_in_channel
                # it incorporates a single, duplicate node at the lower end
                # Now, if this segment long enough?
                if elev_step:
                    top_elev = self._elev[nodes_in_channel[0]]
                    base_elev = self._elev[nodes_in_channel[-1]]
                    # work up the channel from the base to make new interp pts
                    interp_pt_elevs = np.arange(base_elev, top_elev,
                                                elev_step)
                    if interp_pt_elevs.size <= 1:
                        # <1 step; bail on this whole segment
                        break
                    # now we can fairly closely follow the Geomorphtools
                    # algorithm:
                    ch_nodes = np.array(nodes_in_channel)
                    # ^ this is top-to-bottom
                    ch_A = self.grid.at_node['drainage_area'][ch_nodes]
                    ch_dists = self.channel_distances_downstream(ch_nodes)
                    ch_S = self.interpolate_slopes_with_step(ch_nodes,
                                                             ch_dists,
                                                             interp_pt_elevs)
                else:
                    # all the nodes; much easier as links work
                    ch_nodes = np.array(nodes_in_channel)
                    ch_dists = self.channel_distances_downstream(ch_nodes)
                    ch_A = self.grid.at_node['drainage_area'][ch_nodes]
                    ch_S = self.grid.at_node['topographic__steepest_slope'][
                        ch_nodes]
                    assert np.all(ch_S >= 0.)
                # if we're doing spatial discretization, do it here:
                if discretization_length:
                    ch_ksn = self.calc_ksn_discretized(ch_dists, ch_A, ch_S,
                                                       reftheta,
                                                       discretization_length)
                else:  # not discretized
                    # also chopping off the final node, as above
                    log_A = np.log10(ch_A[:-1])
                    log_S = np.log10(ch_S[:-1])
                    # we're potentially propagating nans here if S<=0
                    log_ksn = log_S + reftheta * log_A
                    ch_ksn = 10.**log_ksn
                # save the answers into the main arrays:
                assert np.all(self._mask[ch_nodes[:-1]])
                # Final node gets trimmed off...
                self.ksn[ch_nodes[:-1]] = ch_ksn
                self._mask[ch_nodes] = False
        # now a final sweep to remove any undefined ksn values:
        self._mask[self.ksn == -1.] = True
        self.ksn[self.ksn == -1.] = 0.

    def channel_distances_downstream(self, ch_nodes):
        """
        Calculates distances downstream from top node of a defined flowpath.

        Parameters
        ----------
        ch_nodes : array of ints
            The nodes along a single defined flow path, starting upstream.

        Returns
        -------
        ch_dists : array of floats
            Distances downstream from top node of ch_nodes.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import FlowRouter
        >>> mg = RasterModelGrid((4,5), (5., 10.))
        >>> for nodes in (mg.nodes_at_right_edge, mg.nodes_at_bottom_edge,
        ...               mg.nodes_at_top_edge):
        ...     mg.status_at_node[nodes] = CLOSED_BOUNDARY
        >>> mg.status_at_node[[6, 12, 13, 14]] = CLOSED_BOUNDARY
        >>> _ = mg.add_field('node', 'topographic__elevation', mg.node_x)
        >>> fr = FlowRouter(mg)
        >>> sf = SteepnessFinder(mg)
        >>> _ = fr.route_flow()
        >>> ch_nodes = np.array([8, 7, 11, 10])
        >>> sf.channel_distances_downstream(ch_nodes)
        array([  0.        ,  10.        ,  21.18033989,  31.18033989])
        """
        ch_links = self.grid.at_node['flow__link_to_receiver_node'][
            ch_nodes]
        ch_dists = np.empty_like(ch_nodes, dtype=float)
        # dists from ch head, NOT drainage divide
        ch_dists[0] = 0.
        np.cumsum(self.grid._length_of_link_with_diagonals[ch_links[:-1]],
                  out=ch_dists[1:])
        return ch_dists

    def interpolate_slopes_with_step(self, ch_nodes, ch_dists,
                                     interp_pt_elevs):
        """
        Maps slopes to nodes, interpolating withing defined vertical intervals.

        This follows Geomorphtools' discretization methods. It is essentially a
        downwind map of the slopes.

        Parameters
        ----------
        ch_nodes : array of ints
            The nodes along a single defined flow path, starting upstream.
        ch_dists : array of floats
            Distances downstream from top node of ch_nodes.
        interp_pt_elevs : array of floats
            Elevations at the discretizing points along the profile, in order
            of increasing elevation.

        Returns
        -------
        ch_S : array of floats
            Interpolated slopes at each node in the flowpath (always positive).

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import FlowRouter
        >>> mg = RasterModelGrid((3,10), (5., 10.))
        >>> for nodes in (mg.nodes_at_right_edge, mg.nodes_at_bottom_edge,
        ...               mg.nodes_at_top_edge):
        ...     mg.status_at_node[nodes] = CLOSED_BOUNDARY
        >>> _ = mg.add_field('node', 'topographic__elevation', mg.node_x**1.1)
        >>> fr = FlowRouter(mg)
        >>> sf = SteepnessFinder(mg)
        >>> _ = fr.route_flow()
        >>> ch_nodes = np.arange(18, 9, -1)
        >>> ch_dists = sf.channel_distances_downstream(ch_nodes)
        >>> interp_pt_elevs = np.array([0., 30., 60., 90., 120.])
        >>> sf.interpolate_slopes_with_step(ch_nodes, ch_dists,
        ...                                 interp_pt_elevs)
        array([ 1.67970205,  1.67970205,  1.67970205,  1.65129294,  1.62115336,
        1.5811951 ,  1.53157521,  1.44240187,  1.36442227])
        >>> mg.at_node['topographic__steepest_slope'][ch_nodes]
        array([ 1.69383001,  1.66972677,  1.64200694,  1.60928598,  1.56915472,
        1.51678178,  1.43964028,  1.25892541,  0.        ])
        >>> mg.at_node['topographic__elevation'][:] = mg.node_x
        >>> interp_pt_elevs = np.array([0., 25., 50., 75., 80.])
        >>> sf.interpolate_slopes_with_step(ch_nodes, ch_dists,
        ...                                 interp_pt_elevs)
        array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
        """
        ch_z = self.grid.at_node['topographic__elevation'][ch_nodes]
        assert ch_z[0] >= interp_pt_elevs[-1], \
            'Highest interp_pt_elev must be below top channel node'
        interp_pt_x = np.interp(interp_pt_elevs, ch_z[::-1], ch_dists[::-1])
        interp_pt_S = np.empty_like(interp_pt_elevs)
        # now a downwind map of the slopes onto the nodes
        # slopes are defined positive
        z_diff = interp_pt_elevs[:-1] - interp_pt_elevs[1:]
        x_diff = interp_pt_x[1:] - interp_pt_x[:-1]
        np.divide(z_diff, x_diff, out=interp_pt_S[:-1])
        interp_pt_S[-1] = interp_pt_S[-2]
        # Map S back onto nodes
        ch_S = np.interp(ch_z, interp_pt_elevs, interp_pt_S)

        return ch_S

    def calc_ksn_discretized(self, ch_dists, ch_A, ch_S,
                             ref_theta, discretization_length):
        """
        Calculate normalized steepness index on defined channel segments.

        Every segment must have at least 2 nodes along it. If not, segments
        will be automatically merged to achieve this. The channel will be
        segmented starting at the *downstream* end.

        NB: The final node in the channel does not receive an index, as it
        either belongs to a longer, existing flow path, or it is a boundary
        node with S = 0. Neither works.

        Parameters
        ----------
        ch_dists : array of floats
            Distances downstream from top node of a single stream path.
        ch_A : array of floats
            Drainage areas at each node in the flowpath.
        ch_S : array of floats
            Slope at each node in the flowpath (defined as positive).
        ref_theta : float
            The reference concavity; must be positive.
        discretization_length : float (m)
            The streamwise length of each segment.

        Returns
        -------
        ch_ksn : array of floats
            The normalized steepness index at each node in the flowpath,
            EXCEPT THE LAST. (i.e., length is (ch_dists.size - 1)). Values
            will be the same within each defined segment.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import FlowRouter
        >>> mg = RasterModelGrid((3,10), (5., 10.))
        >>> for nodes in (mg.nodes_at_right_edge, mg.nodes_at_bottom_edge,
        ...               mg.nodes_at_top_edge):
        ...     mg.status_at_node[nodes] = CLOSED_BOUNDARY
        >>> _ = mg.add_field('node', 'topographic__elevation', mg.node_x)
        >>> fr = FlowRouter(mg)
        >>> sf = SteepnessFinder(mg)
        >>> _ = fr.route_flow()
        >>> ch_nodes = np.arange(18, 9, -1)
        >>> ch_dists = sf.channel_distances_downstream(ch_nodes)
        >>> ch_A = mg.at_node['drainage_area'][ch_nodes]
        >>> ch_S = mg.at_node['topographic__steepest_slope'][ch_nodes]

        >>> ksn_25 = sf.calc_ksn_discretized(ch_dists, ch_A, ch_S, 0.5, 25.)
        >>> ksn_25.size == ch_dists.size - 1
        True
        >>> ksn_25
        array([ -1.        ,  11.0668192 ,  11.0668192 ,  15.70417802,
                15.70417802,  15.70417802,  19.3433642 ,  19.3433642 ])

        >>> ksn_10 = sf.calc_ksn_discretized(ch_dists, ch_A, ch_S, 0.5, 10.)
        >>> ksn_10
        array([  8.40896415,   8.40896415,  13.16074013,  13.16074013,
                16.5487546 ,  16.5487546 ,  19.3433642 ,  19.3433642 ])

        >>> ch_ksn_overdiscretized = sf.calc_ksn_discretized(
        ...     ch_dists, ch_A, ch_S, 0.5, 10.)
        >>> np.allclose(ch_ksn_overdiscretized, ksn_10)
        True
        """
        ch_ksn = np.empty_like(ch_A)
        # need to remove the influence of the final node in the seg,
        # as it reflects either the edge of the grid (S=0) or a point
        # after a confluence - hence the 0.000001
        seg_ends = np.arange(ch_dists[-1]-0.000001, 0.,
                             -discretization_length)[::-1]
        # ^ counts up from 0, but terminates at the far end cleanly
        pts_in_each_seg = np.searchsorted(seg_ends, ch_dists)
        num_segs = pts_in_each_seg[-1]
        i = num_segs - 1  # the final pt is no longer included
        while i >= 0:
            old_i = i
            pts_in_seg = pts_in_each_seg == i
            num_pts_in_seg = int(pts_in_seg.sum())
            # if i == num_segs:
            #     true_pts_in_seg = pts_in_each_seg.copy()
            #     pts_in_each_seg[-1] = False
            # else:
            #     true_pts_in_seg = pts_in_each_seg
            # make sure there's always 2 pts in the seg...
            while num_pts_in_seg < 2:
                i -= 1
                pts_in_seg = np.logical_and(pts_in_each_seg <= old_i,
                                            pts_in_each_seg >= i)
                num_pts_in_seg = int(pts_in_seg.sum())
                if i < 0:
                    break
            if num_pts_in_seg < 2:
                # must be at the end of the seg...
                # nodes in invalid segs at the end get ksn = -1.
                ch_ksn[pts_in_seg] = -1.
                break
            seg_A = ch_A[pts_in_seg]
            seg_S = ch_S[pts_in_seg]
            logseg_A = np.log10(seg_A)
            logseg_S = np.log10(seg_S)
            meanlogseg_A = np.mean(logseg_A)
            meanlogseg_S = np.mean(logseg_S)
            logseg_ksn = meanlogseg_S + ref_theta*meanlogseg_A
            ch_ksn[pts_in_seg] = 10.**logseg_ksn
            i -= 1

        return ch_ksn[:-1]

    @property
    def steepness_indices(self):
        """
        Return the array of channel steepness indices.
        Nodes not in the channel receive zeros.
        """
        return self.ksn

    @property
    def hillslope_mask(self):
        """
        Return a boolean array, False where steepness indices exist.
        """
        return self._mask

    @property
    def masked_steepness_indices(self):
        """
        Returns a masked array version of the 'channel__steepness_index' field.
        This enables easier plotting of the values with
        :func:`landlab.imshow_grid_at_node` or similar.

        Examples
        --------
        Make a topographic map with an overlay of steepness values:

        >>> from landlab import imshow_grid_at_node
        >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
        >>> from landlab.components import FlowRouter, FastscapeEroder
        >>> mg = RasterModelGrid((5, 5), 100.)
        >>> for nodes in (mg.nodes_at_right_edge, mg.nodes_at_bottom_edge,
        ...               mg.nodes_at_top_edge):
        ...     mg.status_at_node[nodes] = CLOSED_BOUNDARY
        >>> _ = mg.add_zeros('node', 'topographic__elevation')
        >>> mg.at_node['topographic__elevation'][mg.core_nodes] = mg.node_x[
        ...     mg.core_nodes]/1000.
        >>> np.random.seed(0)
        >>> mg.at_node['topographic__elevation'][
        ...     mg.core_nodes] += np.random.rand(mg.number_of_core_nodes)
        >>> fr = FlowRouter(mg)
        >>> sp = FastscapeEroder(mg, K_sp=0.01)
        >>> cf = SteepnessFinder(mg, min_drainage_area=20000.)
        >>> for i in range(10):
        ...     mg.at_node['topographic__elevation'][mg.core_nodes] += 10.
        ...     _ = fr.route_flow()
        ...     sp.run_one_step(1000.)
        >>> _ = fr.route_flow()
        >>> cf.calculate_steepnesses()

        >>> imshow_grid_at_node(mg, 'topographic__elevation',
        ...                     allow_colorbar=False)
        >>> imshow_grid_at_node(mg, cf.masked_steepness_indices,
        ...                     color_for_closed=None, cmap='winter')
        """
        return np.ma.array(self.steepness_indices, mask=self.hillslope_mask)
