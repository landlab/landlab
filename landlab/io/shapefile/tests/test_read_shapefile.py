
import os
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from landlab.io.shapefile import read_shapefile

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

def test_read_methow():
    file = os.path.join(_TEST_DATA_DIR, 'methow', 'Methow_Network.shp')
    grid = read_shapefile.read_shapefile(file)
    # TODO add assertions about the resulting grid.

    # # for plotting and testing
    # x_of_polylines = grid['link']['x_of_polyline']
    # y_of_polylines = grid['link']['y_of_polyline']
    #
    # segments = []
    #
    # for i in range(len(x_of_polylines)):
    #     x = np.array(x_of_polylines[i])
    #     y = np.array(y_of_polylines[i])
    #     segment = np.array((x, y)).T
    #     segments.append(segment)
    #
    #
    # from landlab.plot import graph
    # fig, ax = plt.subplots(figsize=(8,8), dpi=300)
    # graph.plot_links(grid, color='c', linestyle='solid', with_id=False,
    #                as_arrow=False, linewidth=1)
    # graph.plot_nodes(grid, color='r', with_id=False, markersize=1)
    #
    # line_segments = LineCollection(segments, color='b', linewidth=0.5)
    # ax.add_collection(line_segments)
    # plt.savefig('test.png')

def test_read_elwah_dhsvm():
    file = os.path.join(_TEST_DATA_DIR, 'elwah_dhsvm', 'elwha_example.shp')
    grid = read_shapefile.read_shapefile(file)
    # TODO add assertions about the resulting grid.
