import matplotlib.pyplot as plt
import numpy as np
import pytest

# from landlab.components import NetworkSedimentTransporter
from landlab import BAD_INDEX_VALUE
from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.plot import graph

_OUT_OF_NETWORK = BAD_INDEX_VALUE - 1


# Example Grid

y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

example_nmg = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

# add variables to grid
example_nmg.at_node["topographic__elevation"] = [0.0, 0.1, 0.3, 0.2, 0.35, 0.45, 0.5, 0.6]
example_nmg.at_node["bedrock__elevation"] = [0.0, 0.1, 0.3, 0.2, 0.35, 0.45, 0.5, 0.6]
area = example_nmg.add_ones("cell_area_at_node", at="node")
example_nmg.at_link["drainage_area"] = [100e6, 10e6, 70e6, 20e6, 70e6, 30e6, 40e6]  # m2
example_nmg.at_link["channel_slope"] = [0.01, 0.02, 0.01, 0.02, 0.02, 0.03, 0.03]
example_nmg.at_link["link_length"] = [10000, 10000, 10000, 10000, 10000, 10000, 10000]  # m

example_nmg.at_link["channel_width"] = 15 * np.ones(np.size(example_nmg.at_link["drainage_area"]))

# Example director
example_flow_director = FlowDirectorSteepest(example_nmg)
example_flow_director.run_one_step()

# Example flow depth
timesteps = 30

Qgage = 8000.0  # (m3/s)
dt = 60 * 60 * 24  # (seconds) daily timestep

Hgage = 1.703 * Qgage ** 0.3447
# (m)
Agage = 4.5895e9
# (m2)

example_flow_depth = (
    np.tile(Hgage, (example_nmg.number_of_links)) / (Agage ** 0.4)
) * np.tile(example_nmg.at_link["drainage_area"], (timesteps + 1, 1)) ** 0.4



# Test
time = [0.0]  # probably not the sensible way to do this...

items = {"grid_element": "link",
       "element_id": np.array([[6]])}

initial_volume = np.array([[1]])
abrasion_rate = np.array([0.0001])

variables = {
  "starting_link": (["item_id"], np.array([6])),
  "abrasion_rate": (["item_id"], abrasion_rate),
  "density": (["item_id"], np.array([2650])),
  "time_arrival_in_link": (["item_id", "time"], np.array([[0.71518937]])),
  "active_layer": (["item_id", "time"], np.array([[1]])),
  "location_in_link": (["item_id", "time"], np.array([[0]])),
  "D": (["item_id", "time"], np.array([[0.05]])),
  "volume": (["item_id", "time"], initial_volume),
}

one_parcel = DataRecord(
  example_nmg,
  items=items,
  time=time,
  data_vars=variables,
  dummy_elements={"link": [_OUT_OF_NETWORK]},
)

timesteps = 4

example_flow_depth = example_flow_depth*10# outrageously high transport rate

nst = NetworkSedimentTransporter(
      example_nmg,
      one_parcel,
      example_flow_director,
      example_flow_depth,
      bed_porosity=0.03,
      g=9.81,
      fluid_density=1000,
      channel_width="channel_width",
      transport_method="WilcockCrowe",
  )

dt = 60 * 60 * 24  # (seconds) daily timestep

# Parcel should
with pytest.raises(RuntimeError):
  for t in range(0, (timesteps * dt), dt):
      nst.run_one_step(dt)
