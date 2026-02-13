import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.sort.sort import sort_spokes_at_hub


@pytest.mark.skipif(
    np.lib.NumpyVersion(np.__version__) < "2.4.0",
    reason="warning only present in numpy >= 2.4",
)
def test_no_numpy_warning(recwarn):
    spokes_at_hub = np.asarray([[0, 1, 2, -1]])
    xy_of_hub = np.asarray([[0.0, 0.0]])
    xy_of_spokes = np.asarray(
        [
            [1.0, 0.0],
            [-1.0, 0.0],
            [0.0, 1.0],
        ]
    )
    expected = [[0, 2, 1, -1]]
    actual = sort_spokes_at_hub(spokes_at_hub, xy_of_hub, xy_of_spokes)
    assert not any("'where' used without 'out'" in str(w.message) for w in recwarn)
    assert_array_equal(actual, expected)
