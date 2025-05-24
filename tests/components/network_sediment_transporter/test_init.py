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


@pytest.mark.parametrize(
    "bad_args",
    ({"flow_director": "bad_fd"}, {"parcels": "bad_parcels"}, {"grid": "bad_nmg"}),
)
def test_bad_arg_raises_type_error(
    example_nmg, example_parcels, example_flow_director, bad_args
):
    good_args = {
        "grid": example_nmg,
        "parcels": example_parcels,
        "flow_director": example_flow_director,
    }
    with pytest.raises(TypeError):
        NetworkSedimentTransporter(
            **{**good_args, **bad_args},
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
