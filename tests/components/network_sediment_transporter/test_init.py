import pytest

from landlab.components import NetworkSedimentTransporter


def test_basic_init(example_nmg, example_parcels, example_flow_director):

    _ = NetworkSedimentTransporter(
        example_nmg,
        example_parcels,
        example_flow_director,
        bed_porosity=0.3,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )


def test_bad_flow_director(example_nmg, example_parcels):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            example_parcels,
            "bad_fd",
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            transport_method="WilcockCrowe",
        )


def test_bad_network_model_grid(example_parcels, example_flow_director):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            "bad_nmg",
            example_parcels,
            example_flow_director,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            transport_method="WilcockCrowe",
        )


def test_bad_parcels(example_nmg, example_flow_director):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            "bad_parcels",
            example_flow_director,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            transport_method="WilcockCrowe",
        )


def test_bad_porosity(example_nmg, example_parcels, example_flow_director):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            example_parcels,
            example_flow_director,
            bed_porosity=-0.03,
            g=9.81,
            fluid_density=1000,
            transport_method="WilcockCrowe",
        )


def test_bad_transport_method(example_nmg, example_parcels, example_flow_director):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            example_parcels,
            example_flow_director,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
            transport_method="bad_transport_method",
        )
