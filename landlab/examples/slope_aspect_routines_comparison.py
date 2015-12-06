

# 23 Sep 14 - SN & EI
# For SI2 meeting - comparisons of slope and aspect calculated
# by Burroughs, Horn and best fit plane methods available in Landlab

from landlab import RasterModelGrid
import matplotlib.pyplot as plt
import numpy as np
from landlab.plot.imshow import imshow_grid

# Elevation_NS.npy is predetermined elevation profile that has
# a North and a South facing slope with a flat valley in the middle
# Slopes were defined to range from 0 through 55 in steps of 5 degrees
# Aspect is supposed to be 0 deg for North facing and 180 deg for South facing

grid = RasterModelGrid(53, 67, 10.)
elev = np.load('elevation_NS.npy')
grid['node']['Elevation'] = elev

ids = grid.node_at_cell
# Burroughs
slope_burrough, aspect_burrough = \
    grid.calculate_slope_aspect_at_nodes_burrough(ids, vals='Elevation')
# Horn
slope_horn, aspect_horn = \
    grid.calculate_slope_aspect_at_nodes_horn(ids, vals='Elevation')
# BestFitPlane
slope_bfp, aspect_bfp = grid.calculate_slope_aspect_at_nodes_best_fit_plane(
    ids, elev)


# node_slopes_using_patches
slope_NSP = grid.node_slopes_using_patches(elevs='Elevation', unit='degrees')

pic = 0
plt.figure(pic)
imshow_grid(grid, 'Elevation', values_at='node', grid_units=('m', 'm'))
plt.title('Elevation in m')
# plt.savefig('Elevation_NS')

pic += 1
plt.figure(pic)
imshow_grid(grid, np.degrees(slope_burrough),
            values_at='cell', grid_units=('m', 'm'))
plt.title('Slope in degrees - Burrough 1998')
# plt.savefig('Slope_burrough_NS')

pic += 1
plt.figure(pic)
imshow_grid(grid, np.degrees(aspect_burrough),
            values_at='cell', grid_units=('m', 'm'))
plt.title('Aspect in degrees - Burrough 1998')
# plt.savefig('aspect_burrough_NS')

pic += 1
plt.figure(pic)
imshow_grid(grid, np.degrees(slope_horn),
            values_at='cell', grid_units=('m', 'm'))
plt.title('Slope in degrees - Horn')
# plt.savefig('Slope_Horn_NS')

pic += 1
plt.figure(pic)
imshow_grid(grid, np.degrees(aspect_horn),
            values_at='cell', grid_units=('m', 'm'))
plt.title('Aspect in degrees - Horn')
# plt.savefig('aspect_Horn_NS')

pic += 1
plt.figure(pic)
imshow_grid(grid, np.degrees(slope_bfp),
            values_at='cell', grid_units=('m', 'm'))
plt.title('Slope in degrees - bestFitPlane')
# plt.savefig('Slope_bfp_NS')

pic += 1
plt.figure(pic)
imshow_grid(grid, np.degrees(aspect_bfp),
            values_at='cell', grid_units=('m', 'm'))
plt.title('Aspect in degrees - bestFitPlane')
# plt.savefig('aspect_bfp_NS')

pic += 1
plt.figure(pic)
imshow_grid(grid, slope_NSP, values_at='node', grid_units=('m', 'm'))
plt.title('Slope in degrees - Node Slopes using Patches')
# plt.savefig('aspect_bfp_NS')

plt.show()
