import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import pytest

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter


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
