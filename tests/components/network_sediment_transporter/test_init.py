import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid

_OUT_OF_NETWORK = NetworkModelGrid.BAD_INDEX - 1


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
        transport_method="WilcockCrowe",
    )


def test_bad_flow_director(example_nmg, example_parcels, example_flow_depth):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            example_parcels,
            "bad_fd",
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
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
            transport_method="WilcockCrowe",
        )


def test_bad_parcels(example_nmg, example_flow_director, example_flow_depth):

    with pytest.raises(ValueError):
        NetworkSedimentTransporter(
            example_nmg,
            "bad_parcels",
            example_flow_director,
            example_flow_depth,
            bed_porosity=0.03,
            g=9.81,
            fluid_density=1000,
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
            bed_porosity=-0.03,
            g=9.81,
            fluid_density=1000,
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
            transport_method="bad_transport_method",
        )
