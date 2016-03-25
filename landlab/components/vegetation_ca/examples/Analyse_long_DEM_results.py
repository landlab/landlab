
GRASS = 0
SHRUB = 1
TREE = 2
BARE = 3
SHRUBSEEDLING = 4
TREESEEDLING = 5
import os
from landlab import RasterModelGrid as rmg
from landlab.io   import read_esri_ascii
from landlab.plot import imshow_grid
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

## Point to the input DEM
_DEFAULT_INPUT_FILE_1 = os.path.join(os.path.dirname(__file__),
                                 'DEM_10m.asc')
## Point to the input elevation file
_DEFAULT_INPUT_FILE_2 = os.path.join(os.path.dirname(__file__),
                                 'elevation_NS.npy')
USE_DEM = 1    # Make this 0 to use a custom grid
if USE_DEM == 1:
    ## Importing Grid and Elevations from DEM
    (grid,elevation) = read_esri_ascii(_DEFAULT_INPUT_FILE_1)
    grid['node']['Elevation'] = elevation
else:
    ## Use an open book like grid - custom grid
    grid = rmg(53,67,10.)
    elev = np.load(_DEFAULT_INPUT_FILE_2)
    grid['node']['Elevation'] = elev

sim = 'long_DEM_'
CumWaterStress = np.load(sim+'CumWaterStress.npy')
P = np.load(sim+'P.npy')
Tb = np.load(sim+'Tb.npy')
Tr = np.load(sim+'Tr.npy')
yrs = np.load(sim+'Years.npy')
VegType = np.load(sim+'VegType.npy')

n = P.shape[0]    # Number of iterations
Time = np.empty(n)
Time[0] = 0
for x in range(1,n):
    Time[x] = Time[x-1]+(Tb[x]+Tr[x])/(24.*365)

grass_cov = np.empty(yrs)
shrub_cov = np.empty(yrs)
tree_cov = np.empty(yrs)
grid_size = float(VegType.shape[1])
for x in range(0,yrs):
    grass_cov[x] = (VegType[x][VegType[x] == GRASS].size/grid_size) * 100
    shrub_cov[x] = (VegType[x][VegType[x] == SHRUB].size/grid_size) * 100 + \
                 (VegType[x][VegType[x] == SHRUBSEEDLING].size/grid_size) * 100
    tree_cov[x] = (VegType[x][VegType[x] == TREE].size/grid_size) * 100 + \
                 (VegType[x][VegType[x] == TREESEEDLING].size/grid_size) * 100


years = range(0,yrs)
pic = 0
plt.figure(pic)
plt.plot(years, grass_cov, '-g', label = 'Grass')
plt.hold(True)
plt.plot(years, shrub_cov, '-r', label = 'Shrub')
plt.hold(True)
plt.plot(years, tree_cov, '-k', label = 'Tree')
plt.ylabel(' % Coverage ')
plt.xlabel('Time in years' )
plt.legend(loc = 0)

## Plotting
cmap = mpl.colors.ListedColormap(   \
                    [ 'green', 'red', 'black', 'white', 'red', 'black' ] )
bounds = [-0.5,0.5,1.5,2.5,3.5,4.5,5.5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

## Plot images to make gif.
for year in range(0,yrs,10):
    filename = 'Year = ' +  "%05d" % year
    pic += 1
    plt.figure(pic)
    imshow_grid(grid,VegType[year],values_at = 'cell', cmap=cmap, norm=norm, limits = [0,5])
    plt.title(filename)
    plt.savefig(filename)

plt.show()
