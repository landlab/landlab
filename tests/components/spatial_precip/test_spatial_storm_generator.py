import os

import numpy as np

from landlab import RasterModelGrid
from landlab.components import SpatialPrecipitationDistribution

# from matplotlib.pyplot import plot, show, figure
# from landlab import imshow_grid_at_node
# import shapefile as shp


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
# _THIS_DIR = '.'
# datapath = os.path.join(_THIS_DIR, 'model_input')


def test_MS_params():

    # This commented code block ingests MS's grid setup to come close to a
    # version of this topo (DEJH believes there may be a bug in the STORM
    # code's allocation of elevations, which precludes a precise replication)
    # a = shp.Reader(os.path.join(datapath, 'boundary', 'boundary.shp'))
    # # This is shapefile of the watershed boundary
    # b = shp.Reader(os.path.join(
    #     datapath, 'point_elevations', 'point_elevations',
    #     'point_elevations.shp'))
    # # This is shapefile of the elevations at each 'gauging' point on the
    # # grid within the watershed boundary
    # coords = a.shapeRecords()[0].shape.__geo_interface__['coordinates'][0]
    # Yy = [Y for (X, Y) in coords]
    # Xx = [X for (X, Y) in coords]
    # polypts_xy = [[X, Y] for (X, Y) in coords]
    # (X1, Y1) = np.meshgrid(
    #     np.linspace(min(Xx), max(Xx),
    #                 int(round((max(Xx)-min(Xx))/1000.))),
    #     np.linspace(min(Yy), max(Yy),
    #                 int(round((max(Yy)-min(Yy))/1000.))))
    # # creates a mesh with 1000m spacings
    # isin = Path(polypts_xy).contains_points(
    #     [(zx, zy) for (zx, zy) in zip(X1.flat, Y1.flat)]).reshape(X1.shape)
    # # isin=inpolygon(X1(:),Y1(:),Xx,Yy)
    # Yin = Y1[isin]
    # Xin = X1[isin]
    #
    # Easting = np.loadtxt(os.path.join(_THIS_DIR, 'model_input', 'Easting.csv'))
    # # This is the Longitudinal data for each gauge.
    # Northing = np.loadtxt(os.path.join(_THIS_DIR, 'model_input', 'Northing.csv'))
    # # This is the Latitudinal data for each gauge. It will be sampled
    # # below.
    # vdg = VoronoiDelaunayGrid(Easting, Northing)
    # gauges = np.loadtxt(os.path.join(_THIS_DIR, 'model_input', 'gauges.csv'))
    # # This is the list of gauge numbers. It will be sampled below.
    # gauge_elev = np.loadtxt(os.path.join(_THIS_DIR, 'model_input', 'gauge_elev.csv'))
    # # This is the list of gauge numbers. It will be sampled below.
    # # put the elevs on the grid. Mind the ordering
    # vdg_z = vdg.add_field('node', 'topographic__elevation',
    #                       gauge_elev[np.argsort(Northing)])
    # numgauges = len(gauges)
    #
    # mg = RasterModelGrid((12, 26), (1042.3713, 1102.0973))
    # mg.status_at_node[:] = 4
    # mg.status_at_node[isin.flatten()] = 0
    # z = mg.add_zeros('node', 'topographic__elevation')
    #
    # closest_core_node_in_vdg = []
    # for E, N in zip(Xin, Yin):
    #     closest_core_node_in_vdg.append(
    #        np.argmin(vdg.calc_distances_of_nodes_to_point((E, N))))
    # z[mg.status_at_node == 0] = vdg_z[np.array(closest_core_node_in_vdg)]

    mg = RasterModelGrid((12, 26), xy_spacing=(1102.0973, 1042.3713))
    mg.status_at_node = np.loadtxt(os.path.join(_THIS_DIR, "BCs_Singer.txt"))
    mg.add_field(
        "node",
        "topographic__elevation",
        np.loadtxt(os.path.join(_THIS_DIR, "elevs_Singer.txt")),
    )

    np.random.seed(10)
    rain = SpatialPrecipitationDistribution(
        mg, number_of_years=2, orographic_scenario="Singer"
    )

    max_intensity = []
    storm_dur = []
    istorm_dur = []
    rec = []
    depth = []
    count = 0
    for (storm, istorm) in rain.yield_storms(
        # style='monsoonal', limit='total_rainfall'):
        style="whole_year",
        limit="total_rainfall",
    ):
        # print('storm dur:', storm, rain.storm_duration_last_storm)
        # print('istorm dur:', istorm)
        # print('intensity:', rain.storm_intensity_last_storm)
        # print('recession:', rain.storm_recession_value_last_storm)
        # print('depth:', rain.storm_depth_last_storm)
        # # print('accum depth:', rain.total_rainfall_this_season)
        # print('target depth:', rain.target_median_total_rainfall_this_season)
        # print('***')
        max_intensity.append(rain.storm_intensity_last_storm)
        storm_dur.append(storm)
        istorm_dur.append(istorm)
        rec.append(rain.storm_recession_value_last_storm)
        depth.append(rain.storm_depth_last_storm)
        count += 1
    # print('Total number of storms:', count)
    # print('Target_depth:', rain.target_median_total_rainfall_this_season)
    assert np.isclose(np.mean(max_intensity), 35.84806631859928)  # mm/hr
    assert np.isclose(np.mean(storm), 0.40865380457460571)  # hrs
    assert np.isclose(np.mean(istorm), 85.258871894485694)  # hrs

    # XYZ = np.loadtxt(_THIS_DIR + '/XYZ.txt')
    # X = XYZ[:, 0]
    # Y = XYZ[:, 1]
    # z = XYZ[:, 2]
    # vdg = VoronoiDelaunayGrid(X+np.random.rand(len(X)), Y+np.random.rand(len(Y)))
    # vdg = RasterModelGrid((12, 25), (1042.3713, 1102.0973))
    # z = vdg.add_field('node', 'topographic__elevation', z)
    # rain = SpatialPrecipitationDistribution(vdg, number_of_years=1,
    #                                         orographic_scenario='Singer')
    #
    # storms = [storm for (storm, istorm) in rain.yield_storms(
    #     style='monsoonal', monsoon_storm_interarrival_GEV={
    #                          'shape': -0.807971, 'sigma': 9.4957,
    #                          'mu': 10.6108, 'trunc_interval': (0., 720.)})]
    #
    # istorms = [istorm for (storm, istorm) in rain.yield_storms(
    #     style='monsoonal', monsoon_storm_interarrival_GEV={
    #                          'shape': -0.807971, 'sigma': 9.4957,
    #                          'mu': 10.6108, 'trunc_interval': (0., 720.)})]
