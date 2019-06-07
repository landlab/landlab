# -*- coding: utf-8 -*-
"""
Unit tests for landlab.clast_tracker
Tests output of building a clast collection (ClastCollection)

Last updated 02/11/2019

"""

import numpy as np


def test_grid_shape(cc_south):
    assert cc_south._grid.number_of_node_rows == 5
    assert cc_south._grid.number_of_node_columns == 5


def test_cc_name(cc_south):
    assert cc_south._name == 'ClastCollection'


def test_cc_dim(cc_south):
    assert len(cc_south.dims) == 2
    assert cc_south.dims['item_id'] == 2
    assert cc_south.dims['time'] == 1
    assert np.allclose(cc_south.time.values, [0.])
    assert np.allclose(cc_south.item_id.values, [0, 1])


def test_cc_inherited_from_DR(cc_south):
    assert cc_south.earliest_time == 0.


def test_var_names(cc_south):
    assert cc_south.variable_names == ['grid_element',
                                       'element_id',
                                       'clast__x',
                                       'clast__y',
                                       'clast__elev',
                                       'clast__node',
                                       'clast__initial_radius',
                                       'clast__radius',
                                       'lambda_0',
                                       'lambda_mean',
                                       'slope__WE',
                                       'slope__SN',
                                       'slope__steepest_azimuth',
                                       'slope__steepest_dip',
                                       'distance__to_exit',
                                       'target_node',
                                       'target_node_flag',
                                       'distance__to_travel',
                                       'change_x',
                                       'change_y',
                                       'hop_length',
                                       'total_travelled_dist',
                                       'close2boundary']


def test_var_val(cc_south):
    assert (cc_south.clast__elev.values == ([[0.6]])).all()
    assert max(cc_south.total_travelled_dist.values) == 0.
    assert np.isnan(max(cc_south.slope__SN.values))


def test_var_dim(cc_south):
    assert cc_south.clast__x.dims == ('item_id', 'time')
    assert cc_south.slope__SN.dims == ('item_id',)


def test_phantom(cc_south):
    assert np.logical_not(cc_south.phantom(0))
