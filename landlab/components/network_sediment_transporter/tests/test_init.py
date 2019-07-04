import numpy as np
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
        active_layer_thickness=0.5,
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
            active_layer_thickness=0.5,
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
            active_layer_thickness=0.5,
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
            active_layer_thickness=0.5,
            bed_porosity=-0.03,
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
            active_layer_thickness=0.5,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            channel_width="channel_width",
            transport_method="bad_transport_method",
        )
