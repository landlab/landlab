"""Test non-linear diffuser.

This tester turns over the nonlinear diffuser a couple of times to ensure basic
functionality is working.

***
Note that the result of this test is not correct! In that, it does not honor
boundary conditions correctly. Nevertheless, this test stands for the moment
to prevent grid changes damaging core component functionality.
***
"""
from __future__ import print_function

import os

import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components.nonlinear_diffusion import PerronNLDiffuse

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))

INPUTS = os.path.join(_THIS_DIR, "drive_perron_params.txt")

nrows = 10
ncols = 20
dx = 1.0
dt = 0.1
time_to_run = 0.3
uplift = 0.01

t_z = np.array(
    [
        0.00113696,
        0.00113696,
        0.00121257,
        0.00126092,
        0.00129241,
        0.001313,
        0.00132635,
        0.00133475,
        0.00133964,
        0.00134185,
        0.00134175,
        0.00133932,
        0.00133416,
        0.00132538,
        0.0013115,
        0.00129011,
        0.00125741,
        0.00120719,
        0.00112856,
        0.00112856,
        0.00113696,
        0.00014565,
        0.00022846,
        0.00028308,
        0.00031976,
        0.00034447,
        0.00036094,
        0.00037156,
        0.00037785,
        0.00038073,
        0.00038059,
        0.00037743,
        0.00037079,
        0.00035971,
        0.00034261,
        0.00031699,
        0.00027899,
        0.00022238,
        0.00013644,
        0.0,
        0.00121048,
        0.00022585,
        0.00036385,
        0.00045798,
        0.00052224,
        0.00056588,
        0.00059508,
        0.00061395,
        0.00062515,
        0.00063027,
        0.00063003,
        0.0006244,
        0.00061258,
        0.0005929,
        0.00056258,
        0.00051736,
        0.00045082,
        0.00035338,
        0.00021051,
        0.0,
        0.00125268,
        0.00027286,
        0.00044594,
        0.00056669,
        0.00065019,
        0.00070732,
        0.0007457,
        0.00077055,
        0.00078532,
        0.00079207,
        0.00079175,
        0.00078433,
        0.00076875,
        0.00074283,
        0.00070298,
        0.0006438,
        0.00055738,
        0.0004325,
        0.00025361,
        0.0,
        0.00127173,
        0.00029436,
        0.00048404,
        0.00061779,
        0.00071091,
        0.00077486,
        0.00081791,
        0.00084583,
        0.00086243,
        0.00087002,
        0.00086967,
        0.00086132,
        0.0008438,
        0.00081469,
        0.00076999,
        0.00070375,
        0.00060741,
        0.00046916,
        0.00027327,
        0.0,
        0.00127086,
        0.00029335,
        0.00048223,
        0.00061533,
        0.00070795,
        0.00077154,
        0.00081434,
        0.00084209,
        0.00085859,
        0.00086613,
        0.00086578,
        0.00085748,
        0.00084007,
        0.00081113,
        0.0007667,
        0.00070083,
        0.00060501,
        0.00046742,
        0.00027235,
        0.0,
        0.00124992,
        0.00026971,
        0.00044033,
        0.00055911,
        0.00064114,
        0.0006972,
        0.00073484,
        0.0007592,
        0.00077367,
        0.00078029,
        0.00077998,
        0.0007727,
        0.00075743,
        0.00073202,
        0.00069295,
        0.00063487,
        0.00054997,
        0.00042711,
        0.00025074,
        0.0,
        0.00120541,
        0.00022016,
        0.00035387,
        0.00044471,
        0.00050656,
        0.0005485,
        0.00057653,
        0.00059462,
        0.00060536,
        0.00061027,
        0.00061004,
        0.00060464,
        0.00059331,
        0.00057443,
        0.00054533,
        0.00050187,
        0.00043782,
        0.00034377,
        0.0002053,
        0.0,
        0.00112861,
        0.00013651,
        0.000213,
        0.00026306,
        0.00029655,
        0.00031906,
        0.00033403,
        0.00034368,
        0.0003494,
        0.00035201,
        0.00035188,
        0.00034901,
        0.00034298,
        0.00033292,
        0.00031737,
        0.00029403,
        0.00025933,
        0.00020743,
        0.00012801,
        0.0,
        0.003,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
)


def test_sniff_Perron():
    mg = RasterModelGrid((nrows, ncols), xy_spacing=(dx, dx))
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)
    mg.add_zeros("topographic__elevation", at="node")
    diffusion_component = PerronNLDiffuse(mg, INPUTS)

    elapsed_time = 0.0
    while elapsed_time < time_to_run:
        diffusion_component.input_timestep(dt)
        mg.at_node["topographic__elevation"][mg.core_nodes] += uplift * dt
        mg.at_node["topographic__elevation"][mg.nodes_at_left_edge] += uplift * dt
        mg.at_node["topographic__elevation"][mg.nodes_at_bottom_edge] += uplift * dt
        mg = diffusion_component.diffuse(mg, elapsed_time)
        elapsed_time += dt

    assert_array_almost_equal(mg.at_node["topographic__elevation"], t_z)
