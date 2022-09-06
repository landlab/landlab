import pytest


from landlab.components import SedimentPulserAtLinks, SedimentPulserEachParcel
from landlab import RasterModelGrid


def test_basic_init_each_parcel(example_nmg, example_parcels):
    """test SedimentPulserEachParcel initialization"""
    _ = SedimentPulserEachParcel(example_nmg, parcels=example_parcels)


def test_basic_init_at_link(example_nmg, example_parcels):
    """test SedimentPulserAtLinks initialization, no time to pulse given"""
    _ = SedimentPulserAtLinks(example_nmg, parcels=example_parcels)


def test_grid_is_nmg():
    """test ValueError exception is raised when nmg is a mg"""
    nmg = RasterModelGrid((5, 5))
    with pytest.raises(ValueError):
        _ = SedimentPulserAtLinks(nmg)
    with pytest.raises(ValueError):
        _ = SedimentPulserEachParcel(nmg)
