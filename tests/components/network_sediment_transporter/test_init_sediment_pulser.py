import pytest
import numpy as np


from landlab.components.network_sediment_transporter.sediment_pulser_at_links import SedimentPulserAtLinks
from landlab.components.network_sediment_transporter.sediment_pulser_each_parcel import SedimentPulserEachParcel

from landlab import RasterModelGrid

def test_basic_init_each_parcel(example_nmg, example_parcels):
    """test SedimentPulserEachParcel initialization"""
    _ = SedimentPulserEachParcel(example_nmg, parcels = example_parcels)


def test_basic_init_at_link(example_nmg, example_parcels):
    """test SedimentPulserAtLinks initialization, no time to pulse given"""
    _ = SedimentPulserAtLinks(example_nmg, parcels = example_parcels)

def test_grid_is_nmg():
    """test ValueError exception is raised when nmg is a mg"""
    nmg = RasterModelGrid((5, 5))
    num_starting_parcels = 2
    with pytest.raises(ValueError):
        initialize_parcels = SedimentPulserAtLinks(nmg)                                                
    with pytest.raises(ValueError):
        initialize_parcels = SedimentPulserEachParcel(nmg)

