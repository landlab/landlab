# -*- coding: utf-8 -*-
"""
Unit tests for landlab.clast_tracker
Tests run_on_step output

Last updated 02/11/2019

"""

import pytest
import numpy as np
from landlab.components import FlowAccumulator, LinearDiffuser







def test_run_one_step(cc_south, grid_south):
    fa = FlowAccumulator(grid_south,
                         'topographic__elevation',
                         flow_director='FlowDirectorSteepest')
    kappa = 0.001
    ld = LinearDiffuser(grid_south,
                        linear_diffusivity=kappa)
    cc_south.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    # attributes
    assert cc_south.attrs['kappa'] == 0.001
    # coordinates
    assert cc_south.time_coordinates == [0.0, 10.0]
    # slope
    assert np.isclose(cc_south.slope__steepest_azimuth.values[0], 4.71238898)
    assert np.isclose(cc_south.slope__steepest_dip.values[0], 0.2914567944)
    # travel distance
    assert (cc_south.total_travelled_dist.values[0] == \
            cc_south.hop_length.values[0,1])