# -*- coding: utf-8 -*-
"""
Unit tests for landlab.components.overland_flow.KinwaveOverlandFlowModel

last updated: 3/14/16
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.overland_flow import KinwaveOverlandFlowModel

(_SHAPE, _SPACING, _ORIGIN) = ((10, 10), (25, 25), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((10, 10), spacing=0.5)
    grid.add_zeros('node', 'topographic__elevation', dtype=float)
    grid.add_zeros('node', 'topographic__gradient')

    globals().update({
        'KinWaveOF': KinwaveOverlandFlowModel(grid)})


@with_setup(setup_grid)
def test_KinWaveOF_name():
    assert_equal(KinWaveOF.name, 'KinwaveOverlandFlowModel')


@with_setup(setup_grid)
def test_KinWaveOF_input_var_names():
    assert_equal(KinWaveOF.input_var_names,  ('topographic__elevation',
                                          'topographic__gradient',))


@with_setup(setup_grid)
def test_KinWaveOF_output_var_names():
    assert_equal(KinWaveOF.output_var_names, ('surface_water__depth',
                                          'water__velocity',
                                          'water__specific_discharge'))

@with_setup(setup_grid)
def test_KinWaveOF_var_units():
    assert_equal(set(KinWaveOF.input_var_names) |
                 set(KinWaveOF.output_var_names),
                 set(dict(KinWaveOF.units).keys()))

    assert_equal(KinWaveOF.var_units('topographic__elevation'), 'm')
    assert_equal(KinWaveOF.var_units('topographic__gradient'), 'm/m')
    assert_equal(KinWaveOF.var_units('surface_water__depth'), 'm')
    assert_equal(KinWaveOF.var_units('water__velocity'), 'm/s')
    assert_equal(KinWaveOF.var_units('water__specific_discharge'), 'm2/s')

@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(KinWaveOF.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(KinWaveOF.grid.number_of_node_columns, _SHAPE[1])

def test_run_one_step():
    from landlab import RasterModelGrid
    import numpy as np
    from landlab.components.overland_flow import KinwaveOverlandFlowModel

    grid = RasterModelGrid((10, 10), spacing=0.5)
    grid.add_zeros('node', 'topographic__elevation', dtype=float)
    grid.add_zeros('node', 'topographic__gradient')

    topo_arr = np.ones(100).reshape(10, 10)
    i=0
    while i <= 9:
        topo_arr[:, i]  = 5 + (0.002*i)
        i+=1
    topo_arr = topo_arr.flatten()
    grid['node']['topographic__elevation'] = topo_arr
    KinWaveOF = KinwaveOverlandFlowModel(grid, precip_rate=100.,
        precip_duration=1.0, roughness=0.02)

    KinWaveOF.run_one_step(60)

    # I'll admit this is very non-robust. Solution roughly based on plot #9
    # from Heng et. al, (2009): "Modeling overland flow and soil eroion on
    # non uniform hillslopes: A finite volume scheme." They do not provide the
    # numerical solution but the plots match...
    max_h_mm = max(grid['node']['surface_water__depth']) * 1000.
    np.testing.assert_almost_equal(max_h_mm, 1.66666666667)
