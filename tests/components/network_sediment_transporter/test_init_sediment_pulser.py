import pytest


from landlab.components import SedimentPulserAtLinks, SedimentPulserEachParcel
from landlab import RasterModelGrid


variable_list = [
    "_grid",
    "_parcels",
    "_D50",
    "_D84_D50",
    "_rho_sediment",
    "_parcel_volume",
    "_abrasion_rate",
]


def test_basic_init_each_parcel(example_nmg, example_parcels):
    """test all class variables are present after SedimentPulserEachParcel initializes"""
    pulser = SedimentPulserEachParcel(example_nmg, parcels=example_parcels)
    TF = [True if val in list(vars(pulser).keys()) else False for val in variable_list]
    assert all(TF)


def test_basic_init_at_link(example_nmg, example_parcels):
    """test all class variables are present after SedimentPulserAtLinks initializes"""
    pulser = SedimentPulserAtLinks(example_nmg, parcels=example_parcels)
    TF = [True if val in list(vars(pulser).keys()) else False for val in variable_list]
    assert all(TF)


def test_grid_is_nmg():
    """test ValueError exception is raised when nmg is a mg"""
    nmg = RasterModelGrid((5, 5))
    with pytest.raises(ValueError):
        _ = SedimentPulserAtLinks(nmg)
    with pytest.raises(ValueError):
        _ = SedimentPulserEachParcel(nmg)
