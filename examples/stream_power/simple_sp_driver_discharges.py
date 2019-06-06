"""
simple_sp_driver.py

A simple driver implementing Braun-Willett flow routing and then a
(non-fastscape) stream power component.
DEJH, 09/15/14
"""


import numpy as np
import pylab


from landlab import ModelParameterDictionary, RasterModelGrid
from landlab.components import FlowAccumulator, StreamPowerEroder
from landlab.components.stream_power import FastscapeEroder as Fsc
from landlab.plot.imshow import imshow_grid

inputs = ModelParameterDictionary("./drive_sp_params_discharge.txt")
nrows = 5
ncols = 5
dx = inputs.read_float("dx")
dt = inputs.read_float("dt")
time_to_run = inputs.read_float("run_time")
# nt needs defining
uplift = inputs.read_float("uplift_rate")
init_elev = inputs.read_float("init_elev")

mg = RasterModelGrid((nrows, ncols), xy_spacing=dx)

# create the fields in the grid
mg.add_zeros("topographic__elevation", at="node")
z = np.array(
    [
        5.,
        5.,
        0.,
        5.,
        5.,
        5.,
        2.,
        1.,
        2.,
        5.,
        5.,
        3.,
        2.,
        3.,
        5.,
        5.,
        4.,
        4.,
        4.,
        5.,
        5.,
        5.,
        5.,
        5.,
        5.,
    ]
)
mg["node"]["topographic__elevation"] = z

print("Running ...")

# instantiate the components:
fr = FlowAccumulator(mg, flow_director="D8")
sp = StreamPowerEroder(mg, "./drive_sp_params_discharge.txt")
# load the Fastscape module too, to allow direct comparison
fsp = Fsc(mg, "./drive_sp_params_discharge.txt")

# perform the loop (once!)
for i in range(1):
    fr.route_flow(method="D8")
    my_Q = mg.at_node["surface_water__discharge"] * 1.
    sp.erode(
        mg,
        dt,
        node_drainage_areas="drainage_area",
        slopes_at_nodes="topographic__steepest_slope",
        Q_if_used=my_Q,
    )
    # no uplift

# print the stream power that was calculated:
print("stream power values:")
print(mg.at_node["stream_power_erosion"])

# Finalize and plot
elev = mg["node"]["topographic__elevation"]
elev_r = mg.node_vector_to_raster(elev)

# Clear previous plots
pylab.figure(1)
pylab.close()

# Plot topography
pylab.figure(1)
im = imshow_grid(mg, "topographic__elevation")  # display a colored image

pylab.figure(2)
im = pylab.plot(dx * np.arange(nrows), elev_r[:, int(ncols // 2)])
pylab.title("Vertical cross section")

pylab.show()

print("Done.")
