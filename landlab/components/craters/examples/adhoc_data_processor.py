import numpy as np
import os
from landlab import RasterModelGrid
from pylab import *
#list_of_dirs = ['0to50k','50to100k','100to150k','150to200k','200to250k','250to300k','300to350k','350to400k','400to450k','450to500k']
#let's focus only on the big chance crater to start:
list_of_dirs = ['50to100k','100to150k','150to200k','200to250k','250to300k','300to350k','350to400k','400to450k','450to500k']
home_dir = '/Users/danhobley/Xcodesvn/PyLL/trunk/landlab/components/craters/examples'
case = 'no_angle_small_mimic'
mg = RasterModelGrid(1000,1000,0.002)
elev_list = []
elev_list_raster = []
xsecs = []
for i in list_of_dirs:
    os.chdir(home_dir)
    os.chdir(i)
    os.chdir(case)
    elev_list.append(np.load('craterssave10000.npy'))

os.chdir(home_dir)
os.chdir('data_for_'+case)

for i in xrange(len(elev_list)):
    this_raster = mg.node_vector_to_raster(elev_list[i])
    elev_list_raster.append(this_raster)
    xsecs.append(elev_list_raster[i][:,550])

min_elevs = np.amin(elev_list_raster)
max_elevs = np.amax(elev_list_raster)
count = 0
for i in elev_list_raster:
    imshow(i,vmin=min_elevs,vmax=max_elevs)
    colorbar()
    mystring = 'map_'+str(count)+'.png'
    savefig(mystring)
    clf()
    count+=1

count = 0
num_lines = len(elev_list)
hoz_range = list(np.arange(1000)*0.002)
fig = figure(1)
for i in xsecs:
    if count == 0:
        colormap = 'r'
    else:
        colormap = str(count/(1.5*num_lines))
    plot(hoz_range,i,color=colormap,label=str(count))
    count+=1

xlabel('Horizontal distance (km)')
ylabel('Elevation (km)')
#legend()
savefig('xsec_550_v.png')
os.chdir(home_dir)
show()