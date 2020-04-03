import numpy as np
from numpy.testing import assert_array_equal

from landlab.components.flow_accum import find_drainage_area_and_discharge
from landlab.components.flow_accum.flow_accum_to_n import (
    find_drainage_area_and_discharge_to_n,
)


def test_boundary_to_n():
    r = np.array(
        [
            [1, 2],
            [4, 5],
            [1, 5],
            [6, 2],
            [4, -1],
            [4, -1],
            [5, 7],
            [4, 5],
            [6, 7],
            [7, 8],
        ]
    )
    p = np.array(
        [
            [0.6, 0.4],
            [0.85, 0.15],
            [0.65, 0.35],
            [0.9, 0.1],
            [1.0, 0.0],
            [1.0, 0.0],
            [0.75, 0.25],
            [0.55, 0.45],
            [0.8, 0.2],
            [0.95, 0.05],
        ]
    )
    s = np.array([4, 5, 1, 7, 2, 6, 0, 8, 3, 9])

    a, q = find_drainage_area_and_discharge_to_n(s, r, p, boundary_nodes=[0])
    true_a = np.array([0.0, 1.715, 1.1, 1.0, 9.0, 4.9775, 2.74, 2.845, 1.05, 1.0])
    assert np.allclose(a, true_a)


def test_boundary_bw():
    r = np.array([2, 5, 2, 7, 5, 5, 6, 5, 7, 8]) - 1
    s = np.array([4, 1, 0, 2, 5, 6, 3, 8, 7, 9])
    a, q = find_drainage_area_and_discharge(s, r, boundary_nodes=[0])
    true_a = np.array([0.0, 2.0, 1.0, 1.0, 9.0, 4.0, 3.0, 2.0, 1.0, 1.0])
    assert_array_equal(a, true_a)
