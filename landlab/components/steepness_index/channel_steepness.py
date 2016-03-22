# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19.

@author: dejh
"""
from __future__ import print_function

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
                        elev_step=0., smooth_elev=False, smooth_slope=False)

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    input_stream : str, file_like, or ModelParameterDictionary, optional
        ModelParameterDictionary that holds the input parameters.
    reference_concavity : float
        The reference concavity to use in the calculation.
    min_drainage_area : float (m**2; default 1.e6)
        The drrainage area down to which to calculate steepness indices.
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
    smooth_elev : bool (default False)
        Controls whether to smooth elevations. geomorphtools do this.
    smooth_slope : bool (default False)
        Controls whether to use a mean local slope (all surrounding nodes),
        or whether to use only the steepest downstream slope (default).
        geomorphtools do the former.
    """
    _name = 'SteepnessFinder'

    _input_var_names = (
        'topographic__elevation',
        'drainage_area',
        'topographic__steepest_slope',
        'flow_receiver',
        'upstream_node_order',
        'links_to_flow_receiver',
    )

    _output_var_names = (
        'channel__steepness_index',
    )

    _var_units = {'topographic__elevation': 'm',
                  'drainage_area': 'm**2',
                  'topographic__steepest_slope': '-',
                  'flow_receiver': '-',
                  'upstream_node_order': '-',
                  'links_to_flow_receiver': '-',
                  'channel__steepness_index': 'variable',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'drainage_area': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow_receiver': 'node',
                    'upstream_node_order': 'node',
                    'links_to_flow_receiver': 'node',
                    'channel__steepness_index': 'node',
                    }

    _var_doc = {'topographic__elevation': 'Surface topographic elevation',
                'drainage_area': 'upstream drainage area',
                'topographic__steepest_slope': ('the steepest downslope ' +
                                                'rise/run leaving the node'),
                'flow_receiver': ('the downstream node at the end of the ' +
                                  'steepest link'),
                'upstream_node_order': ('node order such that nodes must ' +
                                        'appear in the list after all nodes ' +
                                        'downstream of them'),
                'links_to_flow_receiver':
                    ('ID of link downstream of each node, which carries the ' +
                     'discharge'),
                'channel__steepness_index': 'the local steepness index',
                }

    def __init__(self, grid, reference_concavity=0.5, min_drainage_area=1.e6,
                 elev_step=0., discretization_length=0.,
                 smooth_elev=False, smooth_slope=False, **kwds):
        """
        Constructor for the component.
        """
        self._grid = grid
        self._reftheta = reference_concavity
        self.min_drainage = min_drainage_area
        assert elev_step >= 0., "elev_step must be >= 0!"
        self._elev_step = elev_step
        self._discritization = discretization_length
        self._smooth_elev = smooth_elev
        self._smooth_slope = smooth_slope
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
        """
        self._mask.fill(True)
        self.ksn.fill(0.)
        # test for new kwds:
        try:
            reftheta = kwds['reference_concavity']
        except KeyError:
            reftheta = self._reftheta
        try:
            min_drainage = kwds['min_drainage_area']
        except KeyError:
            min_drainage = self.min_drainage
        try:
            elev_step = kwds['elev_step']
        except KeyError:
            elev_step = self._elev_step
        try:
            smooth_elev = kwds['smooth_elev']
        except KeyError:
            smooth_elev = self._smooth_elev
        try:
            smooth_slope = kwds['smooth_slope']
        except KeyError:
            smooth_slope = self._smooth_slope
        upstr_order = self.grid.at_node['upstream_node_order']
        # get an array of only nodes with A above threshold:
        valid_dstr_order = (upstr_order[self.grid.at_node['drainage_area'][
            upstr_order] > min_drainage_area])[::-1]
        valid_dstr_elevs = self.grid.at_node['topographic__elevation'][
            valid_dstr_order]
        valid_dstr_areas = self.grid.at_node['drainage_area'][valid_dstr_order]
        # note elevs are guaranteed to be in order, UNLESS a fill
        # algorithm has been used.
        nodes_incorporated = self.grid.zeros('node', dtype=bool)
        # setup the array for the steepnesses
        self._steepnesses = self.grid.ones
        # now do each poss channel in turn
        # get the head of the first (longest!) channel:
        for dstr_order_index in xrange(valid_dstr_order.size):
            this_ch_top_node = valid_dstr_order[dstr_order_index]  # top node
            if not nodes_incorporated[this_ch_top_node]:
                nodes_incorporated[this_ch_top_node] = True
                nodes_in_channel = [this_ch_top_node, ]
                penultimate_node = this_ch_top_node
                current_node_incorporated = False
                while not current_node_incorporated:
                    next_node = self.grid.at_node['flow_receiver'][
                        penultimate_node]
                    if next_node == penultimate_node:  # end of flow path
                        break
                    nodes_in_channel.append(next_node)
                    current_node_incorporated = nodes_incorporated[next_node]
                    # ^ this is a COPY op, so we're free to update the array
                    nodes_incorporated[next_node] = True
                    penultimate_node = current_node
                # by here, we have a full, unique reach in nodes_in_channel
                # it incorporates a single, duplicate node at the lower end
                # Now, if this segment long enough?
                if self._elev_step:
                    top_elev = self._elev[nodes_in_channel[0]]
                    base_elev = self._elev[nodes_in_channel[-1]]
                    # work up the channel from the base to make new interp pts
                    interp_pt_elevs = np.arange(base_elev, top_elev,
                                                self._elev_step)
                    if interp_pt_elevs.size <= 1:
                        # <1 step; bail on this whole segment
                        break
                    # now we can fairly closely follow the Topotools algorithm:
                    ch_nodes = np.array(nodes_in_channel)
                    # ^ this is top-to-bottom
                    ch_z = self.grid.at_node['topographic__elevation'][
                        ch_nodes]
                    ch_A = self.grid.at_node['drainage_area'][ch_nodes]
                    ch_links = self.grid.at_node['links_to_flow_receiver'][
                        ch_nodes]
                    ch_dists = np.empty_like(ch_A)
                    # dists from ch head, NOT drainage divide
                    ch_dists[0] = 0.
                    ch_dists = np.cumsum(self.grid.link_length[ch_links[:-1]],
                                         out=ch_dists[1:])
                    interp_pt_x = np.interp(interp_pt_elevs, ch_z, ch_dists)
                    interp_pt_S = np.empty_like(interp_pt_A)
                    # now a downwind map of the slopes onto the nodes
                    # slopes are defined positive
                    z_diff = interp_pt_elevs[:-1] - interp_pt_elevs[1:]
                    x_diff = interp_pt_x[1:] - interp_pt_x[:-1]
                    np.divide(z_diff, x_diff, out=interp_pt_S[:-1])
                    interp_pt_S[-1] = interp_pt_S[-2]
                    # Map S back onto nodes
                    ch_S = np.interp(ch_z, interp_pt_elevs, interp_pt_S)
                else:
                    # all the nodes; much easier as links work
                    ch_nodes = np.array(nodes_in_channel)
                    ch_A = self.grid.at_node['drainage_area'][ch_nodes]
                    ch_S = self.grid.at_node['topographic__steepest_slope'][
                        ch_nodes]
                    assert np.all(ch_S >= 0.)
                # if we're doing spatial discritization, do it here:
                if self._discritization:
                    # num_segs = ch_dists[-1]//self._discritization
                    ch_ksn = np.empty_like(ch_A)
                    seg_ends = np.arange(ch_dists[-1], 0.,
                                         -self._discritization)[::-1]
                    # ^ counts up from 0, but terminates at the far end cleanly
                    pts_in_each_seg = np.searchsorted(seg_ends, ch_dists)
                    num_segs = np.searchsorted[-1] + 1
                    i = num_segs - 1
                    while i >= 0:
                        pts_in_seg = pts_in_each_seg == i
                        num_pts_in_seg = int(pts_in_seg.sum())
                        # make sure there's always 2 pts in the seg...
                        while num_pts_in_seg < 2:
                            old_i = i
                            i -= 1
                            pts_in_seg = np.intersect1d(
                                pts_in_each_seg <= old_i,
                                pts_in_each_seg >= i)
                            num_pts_in_seg = int(pts_in_seg.sum())
                            if i < 0:
                                break
                        if num_pts_in_seg < 2:
                            # must be at the end of the seg...
                            break
                        seg_A = ch_A[pts_in_seg]
                        seg_S = ch_S[pts_in_seg]
                        logseg_A = np.log10(seg_A)
                        logseg_S = np.log10(seg_S)
                        meanlogseg_A = np.mean(logseg_A)
                        meanlogseg_S = np.mean(logseg_S)
                        logseg_ksn = meanlogseg_S - self._reftheta*meanlogseg_A
                        ch_ksn[pts_in_seg] = 10.**logseg_ksn
                        i -= 1
                else:  # not discritized
                    log_A = np.log10(ch_A)
                    log_S = np.log10(ch_S)
                    meanlog_A = np.mean(log_A)
                    meanlog_S = np.mean(log_S)
                    # we're potentially propagating nans here if S<=0
                    log_ksn = meanlog_S - self._reftheta*meanlog_A
                    ch_ksn = 10.**log_ksn
                # save the answers into the main arrays:
                assert np.all(self._mask[ch_nodes])
                # ...perhaps we need to trim off the final node?
                self.ksn[ch_nodes] = ch_ksn
                self._mask[ch_nodes] = False

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
