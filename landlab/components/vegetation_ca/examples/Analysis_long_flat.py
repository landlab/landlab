# Jai Sri Sainath!

# GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
# SHRUBSEEDLING = 4; TREESEEDLING = 5
GRASS = 0
SHRUB = 1
TREE = 2
BARE = 3
SHRUBSEEDLING = 4
TREESEEDLING = 5
import numpy as np
import matplotlib.pyplot as plt

sim = 'sim_4_'
CumWaterStress = np.load(sim+'CumWaterStress.npy')
P = np.load(sim+'P.npy')
Tb = np.load(sim+'Tb.npy')
Tr = np.load(sim+'Tr.npy')
yrs = np.load(sim+'Years.npy')
VegType = np.load(sim+'VegType.npy')
#MAP = np.load(sim+'MAP.npy')

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
#plt.figure(pic)
#plt.plot(years, MAP[0:yrs])
#plt.xlim([0, yrs])
#plt.xlabel('Time in years')
#plt.ylabel('Mean Annual Precipitation in mm')
#plt.title('MAP')

pic += 1
plt.figure(pic)
plt.plot(years, CumWaterStress[0:yrs,10])
plt.xlim([0, yrs])
plt.xlabel('Years')
plt.ylabel('Cumulative Water Stress')

pic += 1
plt.figure(pic)
plt.plot(years, grass_cov, '-g', label = 'Grass')
plt.hold(True)
plt.plot(years, shrub_cov, '-r', label = 'Shrub')
plt.hold(True)
plt.plot(years, tree_cov, '-k', label = 'Tree')
plt.ylabel(' % Coverage ')
plt.xlabel('Time in years' )
plt.legend(loc = 0)

plt.show()
