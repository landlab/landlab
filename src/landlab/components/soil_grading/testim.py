results=[]
slide_depth=5
soil_depth=3
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components.soil_grading import SoilGrading
from landlab.components import ExtendedGravelBedrockEroder
import numpy as np

xy_spacing=100
dz_erosion = slide_depth

grid = RasterModelGrid((4, 3), xy_spacing=xy_spacing)
elev = grid.add_zeros("topographic__elevation", at="node")
sed_depth = soil_depth
porosity=0.5
sed_weight = sed_depth * 2650 * (1-porosity) # mass per grid node area per size
grains_weight =[sed_weight, 0]
grain_sizes = [0.001, 0.1]
sg = SoilGrading(grid,
                 meansizes=grain_sizes,
                 grains_weight=grains_weight,
                 phi=porosity)
fa = FlowAccumulator(grid)
fa.run_one_step()
eroder = ExtendedGravelBedrockEroder(grid)

grid.at_node['bedrock__elevation']+=100
grid.at_node['bed_grains__proportions'][:,0]=0
grid.at_node['bed_grains__proportions'][:,1]=1

erosion = np.zeros_like(grid.nodes.flatten()).astype(float)
deposition = np.zeros_like(grid.nodes.flatten()).astype(float)

grid.add_field('landslide__erosion',erosion,at='node')
grid.add_field('landslide__deposition',deposition,at='node')

grid.at_node['landslide__erosion'][4] = dz_erosion
grid.at_node['landslide__deposition'][7] = dz_erosion

grid.at_node['soil__depth'][4]
np.sum(grid.at_node['grains__weight'],1)[4]

g0 = np.copy(grid.at_node['grains__weight'])
grid.at_node['grains__weight']


sg.update_mass_based_on_outsource_dz()
results.append(grid.at_node['grains__weight'][7,1] / np.sum(grid.at_node['grains__weight'][7,:]))