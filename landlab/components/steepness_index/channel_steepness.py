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
    min_drainage_area : float (default 1.e6)
        The drrainage area down to which to calculate steepness indices.
        Defaults to 1.e6 m**2, per Wobus et al. 2006.
    elev_step : float (default 0.)
        If >0., becomes a vertical elevation change step to use to
        discretize the data (per Wobus). If 0., all nodes are used and
        no discretization happens.
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
    )

    _output_var_names = (
        'channel__steepness_index',
    )

    _var_units = {'topographic__elevation': 'm',
                  'drainage_area': 'm**2',
                  'topographic__steepest_slope': '-',
                  'flow_receiver': '-',
                  'upstream_node_order': '-',
                  'channel__steepness_index': 'variable',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'drainage_area': 'node',
                    'topographic__steepest_slope': 'node',
                    'flow_receiver': 'node',
                    'upstream_node_order': 'node',
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
                'channel__steepness_index': 'the local steepness index',
                }

    def __init__(self, grid, reference_concavity=0.5, min_drainage_area=1.e6,
                 elev_step=0., smooth_elev=False, smooth_slope=False, **kwds):
        """
        Constructor for the component.
        """
        self._grid = grid
        self._reftheta = reference_concavity
        self.ksn = self._grid.add_zeros('node', 'channel__steepness_index')

    def calculate_steepnesses(self, **kwds):
        """
        This is the main method. Call it to calculate local steepness indices
        at all points with drainage areas greater than *min_drainage_area*.

        This "run" method can optionally take the same parameter set as
        provided at instantiation. If they are provided, they will override
        the existing values from instantiation.
        """
        upstr_order = mg.at_node['upstream_node_order']
        valid_dstr_order = (upstr_order[mg.at_node['drainage_area'][
            upstr_order] > min_drainage_area])[::-1]
        if elev_step:
            # use the full Wobus style method
            nodes_in_channel = np.array([])
            top_node = np.setdiff1d(valid_dstr_order, nodes_in_channel)[0]
            penultimate_node = top_node
            top_elev = self._elev[top_node]
            current_node = mg.at_node['flow_receiver'][penultimate_node]
