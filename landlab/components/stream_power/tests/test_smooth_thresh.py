from nose.tools import assert_raises

from landlab import RasterModelGrid
from landlab.components.stream_power import StreamPowerSmoothThresholdEroder as Spst

def test_bad_nsp():

    mg = RasterModelGrid(4, 4, 1)
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

    assert_raises(ValueError, Spst, mg, K_sp = 1.0, n_sp = 1.01)
