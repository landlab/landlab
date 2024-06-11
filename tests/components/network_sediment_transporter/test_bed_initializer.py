"""
Tests written by Jeff Keck and Allison Pfeiffer
"""

import numpy as np
import pytest

from landlab.components import BedParcelInitializerArea
from landlab.components import BedParcelInitializerDepth
from landlab.components import BedParcelInitializerDischarge
from landlab.components import BedParcelInitializerUserD50
from landlab.components.network_sediment_transporter.bed_parcel_initializers import (
    _determine_approx_parcel_volume,
)
from landlab.components.network_sediment_transporter.bed_parcel_initializers import (
    calc_d50_dArea_scaling,
)
from landlab.components.network_sediment_transporter.bed_parcel_initializers import (
    calc_d50_depth,
)
from landlab.components.network_sediment_transporter.bed_parcel_initializers import (
    calc_d50_discharge,
)
from landlab.components.network_sediment_transporter.bed_parcel_initializers import (
    calc_total_parcel_volume,
)
from landlab.data_record import DataRecord


# Basic test: that you can call the four initializers, create DataRecord
def test_call_area_BPI(example_nmg2):
    initialize_parcels = BedParcelInitializerArea(
        example_nmg2, drainage_area_coefficient=0.18, drainage_area_exponent=-0.12
    )
    parcels = initialize_parcels()
    assert isinstance(parcels, DataRecord)


def test_call_discharge_BPI(example_nmg2):
    initialize_parcels = BedParcelInitializerDischarge(
        example_nmg2, discharge_at_link=np.full(example_nmg2.number_of_links, 10.0)
    )
    parcels = initialize_parcels()
    assert isinstance(parcels, DataRecord)


def test_call_depth_BPI(example_nmg2):
    depth_at_link = np.full(example_nmg2.number_of_links, 10.0)  # m
    initialize_parcels = BedParcelInitializerDepth(
        example_nmg2,
        flow_depth_at_link=depth_at_link,
    )
    parcels = initialize_parcels()
    assert isinstance(parcels, DataRecord)


def test_call_userD50_BPI(example_nmg2):
    d50_at_link = np.full(example_nmg2.number_of_links, 0.08)  # m
    initialize_parcels = BedParcelInitializerUserD50(
        example_nmg2,
        user_d50=d50_at_link,
    )
    parcels = initialize_parcels()
    assert isinstance(parcels, DataRecord)


# %% Tests of BPI base class helper functions


def test_calc_total_parcel_volume():
    vole = 12.5
    vol = calc_total_parcel_volume(5, 5, 0.5)
    np.testing.assert_allclose(vol, vole, rtol=1e-4)


def test_det_approx_parcel_volume():
    total_parcel_volume_at_link = 10
    median_number_of_starting_parcels = 50
    pvole = 0.2
    pvol = _determine_approx_parcel_volume(
        total_parcel_volume_at_link, median_number_of_starting_parcels
    )
    np.testing.assert_allclose(pvol, pvole, rtol=1e-4)


# %% Test to make sure we're getting correct correct d50 correct_values


def test_calc_D50_discharge():
    """test calc D50 grain size give correct values"""
    width = 20
    slope = 0.01
    discharge = 100
    mannings_n = 0.05
    gravity = 9.81
    rho_water = 1000
    rho_sediment = 3000
    tau_c_50 = 0.05

    D50e = 0.173286

    D50 = calc_d50_discharge(
        width,
        slope,
        discharge,
        mannings_n,
        gravity,
        rho_water,
        rho_sediment,
        tau_c_50,
    )

    np.testing.assert_allclose(D50, D50e, rtol=1e-4)


class Test_calc_d50_dArea_scaling:
    def test_normal_d50_dArea_scaling(self):
        """test normal values 1"""
        D50e = 0.13182567
        a = 0.1
        n = 0.12
        CA = 10  # km^2
        D50 = calc_d50_dArea_scaling(drainage_area=CA, a=a, n=n)
        np.testing.assert_allclose(D50, D50e, rtol=1e-4)

    def test_normal_d50_dArea_scaling_2(self):
        """test normal values 2"""
        D50e = 0.1
        a = 0.1
        n = 0
        CA = 10  # km^2
        D50 = calc_d50_dArea_scaling(drainage_area=CA, a=a, n=n)
        np.testing.assert_allclose(D50, D50e, rtol=1e-4)


def test_calc_d50_depth():
    """test normal values for calc_d50_depth"""
    D50e = 0.1
    D50 = calc_d50_depth(
        slope=0.01,
        flow_depth=1,
        tau_c_multiplier=1,
        rho_water=1000.0,
        rho_sediment=3000.0,
        tau_c_50=0.05,
    )
    np.testing.assert_allclose(D50, D50e, rtol=1e-4)


# %% Test for errors raised
def test_D50_not_specified(example_nmg2):
    """test ValueError exception raised if D50 input is None"""
    with pytest.raises(ValueError):
        BedParcelInitializerUserD50(example_nmg2, user_d50=None)


# %% Test for expected correct values


class Test_BedParcelInitializer:
    def test_normal_BPI(self, example_nmg2):
        """
        Minimum attributes specified, most attributes should use
        defaults specified at instantiation. uses BedParcelInitializerDischarge.
        """
        discharge = np.full(example_nmg2.number_of_links, 1.0)
        initialize_parcels = BedParcelInitializerDischarge(
            example_nmg2,
            discharge_at_link=discharge,
            median_number_of_starting_parcels=1,
            rng=5,
        )

        parcels = initialize_parcels()

        # create expected values and get values from datarecord
        GEe = np.expand_dims(["link"] * 8, axis=1).astype(object)
        GE = parcels.dataset["grid_element"]
        EIe = np.expand_dims(np.array([0, 1, 2, 3, 4, 5, 5, 6]), axis=1)
        EI = parcels.dataset["element_id"]
        SLe = np.array([0, 1, 2, 3, 4, 5, 5, 6])
        SL = parcels.dataset["starting_link"]
        ARe = np.ones(8) * 0
        AR = parcels.dataset["abrasion_rate"]
        De = np.ones(8) * 2650
        D = parcels.dataset["density"]
        TAe = np.expand_dims(np.zeros(8), axis=1)
        #   TAe = np.array([[ 0.08074127], [ 0.7384403 ], [ 0.44130922], [ 0.15830987],
        #                 [ 0.87993703], [ 0.27408646], [ 0.41423502], [ 0.29607993]])
        TA = parcels.dataset["time_arrival_in_link"]
        ALe = np.expand_dims(np.ones(8), axis=1)
        AL = parcels.dataset["active_layer"]
        LL = parcels.dataset["location_in_link"]
        De = [
            [0.02704727],
            [0.02982003],
            [0.05161594],
            [0.10882225],
            [0.15829283],
            [0.09817898],
            [0.06006139],
            [0.0445011],
        ]
        D = parcels.dataset["D"]
        V = parcels.dataset["volume"]

        assert list(GE.values) == list(GEe)
        np.testing.assert_allclose(EI, EIe, rtol=1e-4)
        np.testing.assert_allclose(SL, SLe, rtol=1e-4)
        np.testing.assert_allclose(AR, ARe, rtol=1e-4)
        np.testing.assert_allclose(D, De, rtol=1e-4)
        np.testing.assert_array_less(TAe, TA)
        np.testing.assert_allclose(AL, ALe, rtol=1e-4)
        np.testing.assert_array_less(LL, 1)
        np.testing.assert_array_less(0, D)
        np.testing.assert_array_less(0, V)

    def test_normal_BPI_abrasion(self, example_nmg2):
        """
        Test normal outputs of bed parcel initializer with abrasion.
        uses BedParcelInitializerDischarge.
        """
        discharge = np.full(example_nmg2.number_of_links, 1.0)
        initialize_parcels = BedParcelInitializerDischarge(
            example_nmg2,
            discharge_at_link=discharge,
            median_number_of_starting_parcels=1,
            abrasion_rate=0.1,
            rng=5,
        )

        parcels = initialize_parcels()

        # create expected values and get values from datarecord
        GEe = np.expand_dims(["link"] * 8, axis=1).astype(object)
        GE = parcels.dataset["grid_element"]
        EIe = np.expand_dims(np.array([0, 1, 2, 3, 4, 5, 5, 6]), axis=1)
        EI = parcels.dataset["element_id"]
        SLe = np.array([0, 1, 2, 3, 4, 5, 5, 6])
        SL = parcels.dataset["starting_link"]
        ARe = np.ones(8) * 0.1
        AR = parcels.dataset["abrasion_rate"]
        De = np.ones(8) * 2650
        D = parcels.dataset["density"]
        # TAe = np.array([[ 0.08074127], [ 0.7384403 ], [ 0.44130922], [ 0.15830987],
        #                 [ 0.87993703], [ 0.27408646], [ 0.41423502], [ 0.29607993]])
        TAe = np.expand_dims(np.zeros(8), axis=1)
        TA = parcels.dataset["time_arrival_in_link"]
        ALe = np.expand_dims(np.ones(8), axis=1)
        AL = parcels.dataset["active_layer"]
        De = [
            [0.02704727],
            [0.02982003],
            [0.05161594],
            [0.10882225],
            [0.15829283],
            [0.09817898],
            [0.06006139],
            [0.0445011],
        ]
        D = parcels.dataset["D"]
        V = parcels.dataset["volume"]

        assert list(GE.values) == list(GEe)
        np.testing.assert_allclose(EI, EIe, rtol=1e-4)
        np.testing.assert_allclose(SL, SLe, rtol=1e-4)
        np.testing.assert_allclose(AR, ARe, rtol=1e-4)
        np.testing.assert_allclose(D, De, rtol=1e-4)
        # np.testing.assert_allclose(TA, TAe, rtol = 1e-4)
        np.testing.assert_array_less(TAe, TA)
        np.testing.assert_allclose(AL, ALe, rtol=1e-4)
        # np.testing.assert_allclose(LL, LLe, rtol = 1e-4)
        np.testing.assert_array_less(0, D)
        np.testing.assert_array_less(0, V)
