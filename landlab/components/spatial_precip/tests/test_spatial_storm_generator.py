import os

import pytest
from six.moves import range
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from matplotlib.pyplot import plot, show, figure
from landlab import imshow_grid_at_node

import landlab
from landlab import RasterModelGrid, VoronoiDelaunayGrid, RadialModelGrid, FieldError
from landlab.components import SpatialPrecipitationDistribution

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_MS_params():
    XYZ = np.loadtxt(_THIS_DIR + '/XYZ.txt')
    X = XYZ[:, 0]
    Y = XYZ[:, 1]
    z = XYZ[:, 2]
    vdg = VoronoiDelaunayGrid(X+np.random.rand(len(X)), Y+np.random.rand(len(Y)))
    z = vdg.add_field('node', 'topographic__elevation', z)
    rain = SpatialPrecipitationDistribution(vdg, number_of_years=10,
                                            orographic_scenario='Singer')
    [storm for (storm, istorm) in rain.yield_storms()]



    [istorm for (storm, istorm) in rain.yield_storms()]