

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import pytest

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.grid.network import NetworkModelGrid
from landlab import BAD_INDEX_VALUE
from landlab.data_record import DataRecord
_OUT_OF_NETWORK = BAD_INDEX_VALUE - 1

def test_basic_init(
    example_nmg, example_parcels, example_flow_director, example_flow_depth
):

    nst = NetworkSedimentTransporter(
        example_nmg,
        example_parcels,
        example_flow_director,
        example_flow_depth,
        bed_porosity=0.3,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
    )


def test_bad_flow_director(
    example_nmg, example_parcels, example_flow_depth
):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            example_parcels,
            "bad_fd",
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="WilcockCrowe",
        )


def test_bad_network_model_grid(
    example_parcels, example_flow_director, example_flow_depth
):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            "bad_nmg",
            example_parcels,
            example_flow_director,
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="WilcockCrowe",
        )

def test_bad_parcels(
    example_nmg, example_flow_director, example_flow_depth
):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            "bad_parcels",
            example_flow_director,
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="WilcockCrowe",
        )


def test_bad_porosity(
    example_nmg, example_parcels, example_flow_director, example_flow_depth
):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            example_parcels,
            example_flow_director,
            example_flow_depth,
            bed_porosity= -0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="WilcockCrowe",
        )

def test_bad_transport_method(
    example_nmg, example_parcels, example_flow_director, example_flow_depth
):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            example_parcels,
            example_flow_director,
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="bad_transport_method",
        )

def test_no_flow_no_transport(
    example_nmg, example_parcels, example_flow_director
):

    timesteps = 3
    no_flow = np.zeros([(timesteps+1),(example_nmg.number_of_links)])

    nst = NetworkSedimentTransporter(
        example_nmg,
        example_parcels,
        example_flow_director,
        no_flow,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
    )

    dt = 60 * 60 * 24  # (seconds) daily timestep



    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)

        # Need to define original_node_elev after a timestep has passed.
        if t/(60*60*24) == 1:
            original_node_elev = example_nmg.at_node["topographic__elevation"]

    # No flow? Parcels should stay in same locations in links
    assert_array_almost_equal(
            example_parcels.dataset.location_in_link[:,0],
            example_parcels.dataset.location_in_link[:,1])

    # No flow? Total parcel volume should stay the same
    assert_array_equal(
            np.sum(example_parcels.dataset["volume"].values, axis=0)[0],
            np.sum(example_parcels.dataset["volume"].values, axis=0)[1])

    # No flow? Elevations should stay the same
    assert_array_equal(
            original_node_elev,
            example_nmg.at_node["topographic__elevation"])

def test_parcel_leaves(
    example_nmg, example_flow_depth, example_flow_director
):

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

    timesteps = 3

    example_flow_depth = example_flow_depth*20# outrageously high transport rate

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




def test_abrasion(
    example_nmg, example_parcels, example_flow_depth,example_flow_director
):
    time = [0.0]  # probably not the sensible way to do this...
    
    items = {"grid_element": "link",
             "element_id": np.array([[6],[6]])}
    
    initial_volume = np.array([[1],[1]])
    abrasion_rate = np.array([0.0001,0])
    
    variables = {
        "starting_link": (["item_id"], np.array([6,6])),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], np.array([2650,90000])),
        "time_arrival_in_link": (["item_id", "time"], np.array([[0.5],[0]])),
        "active_layer": (["item_id", "time"], np.array([[1],[1]])),
        "location_in_link": (["item_id", "time"], np.array([[0],[0]])),
        "D": (["item_id", "time"], np.array([[0.05],[0.05]])),
        "volume": (["item_id", "time"], initial_volume),
    }
    
    two_parcels = DataRecord(
        example_nmg,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [_OUT_OF_NETWORK]},
    )
        
    timesteps = 8
    
    example_flow_depth = example_flow_depth*5 # outrageously high transport rate
    
    nst = NetworkSedimentTransporter(
            example_nmg,
            two_parcels,
            example_flow_director,
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="WilcockCrowe",
        )
    
    dt = 60 * 60 * 24  # (seconds) daily timestep
    
    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)
        print("Successfully completed a timestep")
        # Need to define original_node_elev after a timestep has passed.
        if t/(60*60*24) == 1:
            original_parcel_vol = example_parcels.dataset.volume[:,0][0]
            
    
    # Parcel volume should decrease according to abrasion rate
    volume_after_transport = (np.squeeze(np.transpose(initial_volume)) 
                              * np.exp(nst._distance_traveled_cumulative
                                       *-abrasion_rate)
                              )
    
    print("volume_after_transport", volume_after_transport)
                              
    assert_array_almost_equal(
            volume_after_transport,
            two_parcels.dataset.volume[0:2,-1]
            )

def test_defined_parcel_transport():
    
    y_of_node = (0, 0, 0, 0)
    x_of_node = (0, 100, 200, 300)
    nodes_at_link = ((0,1), (1,2), (2,3))
    
    nmg_constant_slope = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    
    # add variables to nmg
    nmg_constant_slope.at_node["topographic__elevation"] = [3., 2., 1., 0.]
    nmg_constant_slope.at_node["bedrock__elevation"] = [3., 2., 1., 0.]
    area = nmg_constant_slope.add_ones("cell_area_at_node", at="node")
    nmg_constant_slope.at_link["drainage_area"] = [10e6, 10e6, 10e6]  # m2
    nmg_constant_slope.at_link["channel_slope"] = [0.001, 0.001, 0.001]
    nmg_constant_slope.at_link["link_length"] = [100, 100, 100]  # m
    nmg_constant_slope.at_link["channel_width"] = 15 * np.ones(np.size(nmg_constant_slope.at_link["drainage_area"]))
    
    flow_director = FlowDirectorSteepest(nmg_constant_slope)
    flow_director.run_one_step()
    
    timesteps = 11
    
    example_flow_depth = (
        np.tile(2, (nmg_constant_slope.number_of_links))
    ) * np.tile(1, (timesteps + 1, 1)) # 2 meter flow depth
    
    time = [0.0]
    
    items = {"grid_element": "link",
             "element_id": np.array([[0]])}
    
    initial_volume = np.array([[1]])
    abrasion_rate = np.array([0])
    
    variables = {
        "starting_link": (["item_id"], np.array([0])),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], np.array([2650])),
        "time_arrival_in_link": (["item_id", "time"], np.array([[0.71518937]])),
        "active_layer": (["item_id", "time"], np.array([[1]])),
        "location_in_link": (["item_id", "time"], np.array([[0]])),
        "D": (["item_id", "time"], np.array([[0.05]])),
        "volume": (["item_id", "time"], initial_volume),
    }
    
    one_parcel = DataRecord(
        nmg_constant_slope,
        items=items,
        time=time,
        data_vars=variables,
        dummy_elements={"link": [_OUT_OF_NETWORK]},
    )

    nst = NetworkSedimentTransporter(
            nmg_constant_slope,
            one_parcel,
            flow_director,
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="WilcockCrowe",
        )
    
    dt = 60  # (seconds) 1 min timestep
    
    distance_traveled = np.arange(0.,timesteps)
    active_layer_thickness_array = np.arange(0.,timesteps)
    #distance_traveled = np.arange(0.,timesteps)
    
    for t in range(0, (timesteps * dt), dt):
            nst.run_one_step(dt)
            distance_traveled[np.int(t/dt)]= nst._distance_traveled_cumulative
            active_layer_thickness_array[np.int(t/dt)]=nst.active_layer_thickness_array
    # NEED TO CALCULATE THINGS HERE. 
    # Transport distance should match?
    Distance_Traveled_Should_Be = ([  21.61998527,   43.23997054,   64.85995581,
         86.47994109,  108.09992636,  129.71991163,  151.3398969 ,
        172.95988217,  194.57986744,  216.19985272,  237.81983799])
        
    assert_array_almost_equal(Distance_Traveled_Should_Be,distance_traveled, decimal = -1)
