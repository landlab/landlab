# -*- coding: utf-8 -*-
"""
Unit tests for landlab.data_record.data_record.DataRecord
Dimension = time

Last updated 10/18/2018

"""

from landlab import RasterModelGrid

grid = RasterModelGrid((3, 3))


def test_dr_nodim_name(dr_nodim):
    assert dr_nodim._name == "DataRecord"
