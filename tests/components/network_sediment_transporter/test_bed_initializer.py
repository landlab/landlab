import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab.components import NetworkSedimentTransporter
from landlab.data_record import DataRecord


def test_bed_initializerSOMETHING(example_nmg,  example_flow_director):

    BPI_created_parcels = BedParcelInitializer(
        example_nmg,
        flow_depth_at_link = example_nmg.at_link["flow_depth"]
        )

    nst = NetworkSedimentTransporter(
        example_nmg,
        BPI_created_parcels,
        example_flow_director,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        transport_method="WilcockCrowe",
    )

    dt = 60 * 60 * 24  # (seconds) daily timestep

    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)

    # Some criteria we want to test
    thing = np.squeeze(np.transpose(initial_volume)) * np.exp(
        nst._distance_traveled_cumulative * -abrasion_rate
    )

    assert_array_almost_equal(
        thing, two_parcels.dataset.volume[0:2, -1]
    )
