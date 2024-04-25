import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from pytest import approx

from landlab.components.network_sediment_transporter.sediment_pulser_at_links import (
    SedimentPulserAtLinks,
)
from landlab.components.network_sediment_transporter.sediment_pulser_base import (
    SedimentPulserBase,
)
from landlab.components.network_sediment_transporter.sediment_pulser_each_parcel import (
    SedimentPulserEachParcel,
)


def always_time_to_pulse(time):
    return True


def time_to_pulse_list(time):
    Ptime = [19, 20, 22, 23, 24, 75, 76]
    return time in Ptime


def test_call_SedimentPulserBase(example_nmg2):
    """test exception raised if SedimentPulserBase is called"""
    grid = example_nmg2
    make_pulse = SedimentPulserBase(grid)

    with pytest.raises(NotImplementedError) as exc_info:
        make_pulse()
    assert exc_info.match("the base component has no call method")


class Test_SedimentPulserAtLinks:
    def test_normal_1(self, example_nmg2):
        """only time specified, links and number parcels specified,
        should use defaults in base class"""
        grid = example_nmg2
        time_to_pulse = always_time_to_pulse

        make_pulse = SedimentPulserAtLinks(grid, time_to_pulse=time_to_pulse, rng=1947)

        time = 11.0
        links = [0]
        n_parcels_at_link = [10]
        parcels = make_pulse(
            time=time, links=links, n_parcels_at_link=n_parcels_at_link
        )

        D50 = 0.05  # default grain size parameters
        D84_D50 = 2.1
        D_sd = D50 * D84_D50 - D50

        assert np.all(parcels.dataset["grid_element"] == "link")
        assert np.all(parcels.dataset["element_id"] == 0)
        assert np.all(parcels.dataset["starting_link"] == 0)
        assert np.all(parcels.dataset["abrasion_rate"] == 0)
        assert np.all(parcels.dataset["density"] == approx(2650.0))
        assert np.all(parcels.dataset["time_arrival_in_link"] == approx(time))
        assert np.all(parcels.dataset["active_layer"] == approx(1.0))
        assert np.all(parcels.dataset["volume"] == approx(0.5))

        assert np.all(
            (parcels.dataset["location_in_link"] >= 0.0)
            & (parcels.dataset["location_in_link"] <= 1.0)
        )
        # check mean randomly selected grain size with within 3 standard deviations
        # of specified mean
        assert np.abs(parcels.dataset["D"].mean() - D50) < 3 * D_sd

    def test_normal_2(self, example_nmg2):
        """only D50 specified, all other parcel attributes should
        use defaults in base class"""

        grid = example_nmg2
        time_to_pulse = always_time_to_pulse

        make_pulse = SedimentPulserAtLinks(grid, time_to_pulse=time_to_pulse, rng=2001)
        time = 11
        links = [2, 6]
        n_parcels_at_link = [2, 3]
        D50 = [0.3, 0.12]
        parcels = make_pulse(
            time=time, links=links, n_parcels_at_link=n_parcels_at_link, D50=D50
        )

        D50_1 = D50[0]  # grain size
        D84_D50_1 = 2.1
        D_sd_1 = D50_1 * D84_D50_1 - D50_1

        D50_2 = D50[1]  # grain size
        D84_D50_2 = 2.1
        D_sd_2 = D50_2 * D84_D50_2 - D50_2

        D = parcels.dataset["D"]

        assert np.all(parcels.dataset["grid_element"] == "link")
        assert np.all(parcels.dataset["element_id"] == [[2], [2], [6], [6], [6]])
        assert np.all(parcels.dataset["starting_link"] == [2, 2, 6, 6, 6])
        assert np.all(parcels.dataset["abrasion_rate"] == 0)
        assert np.all(parcels.dataset["density"] == approx(2650.0))
        assert np.all(parcels.dataset["time_arrival_in_link"] == approx(time))
        assert np.all(parcels.dataset["active_layer"] == approx(1.0))
        assert np.all(parcels.dataset["volume"] == approx(0.5))

        assert np.all(
            (parcels.dataset["location_in_link"] >= 0.0)
            & (parcels.dataset["location_in_link"] <= 1.0)
        )
        # check mean randomly selected grain size with within 3 standard deviations
        # of specified mean
        assert np.abs(D[0:2].mean() - D50_1) < 3 * D_sd_1
        assert np.abs(D[2:].mean() - D50_2) < 3 * D_sd_2

    def test_normal_3(self, example_nmg2):
        """two pulses. First, only time, links and number of parcels specified,
        uses defaults in base class for all other parcel attributes
        second, two parcels in link two and three parcels in link six
        are added and all attributes specified"""

        grid = example_nmg2
        time_to_pulse = always_time_to_pulse

        make_pulse = SedimentPulserAtLinks(grid, time_to_pulse=time_to_pulse, rng=1776)

        time = 11
        links = [0]
        n_parcels_at_link = [2]
        parcels = make_pulse(
            time=time, links=links, n_parcels_at_link=n_parcels_at_link
        )

        time = 12
        links = [2, 6]
        n_parcels_at_link = [2, 3]
        D50 = [0.3, 0.12]
        D84_D50 = [2.1, 1.5]
        parcel_volume = [1, 0.5]
        rho_sediment = [2650, 2500]
        abrasion_rate = [0.1, 0.3]
        parcels = make_pulse(
            time=time,
            links=links,
            n_parcels_at_link=n_parcels_at_link,
            D50=D50,
            D84_D50=D84_D50,
            parcel_volume=parcel_volume,
            rho_sediment=rho_sediment,
            abrasion_rate=abrasion_rate,
        )

        assert_array_equal(
            parcels.dataset["grid_element"].values.tolist(),
            [
                ["link", np.nan],
                ["link", np.nan],
                [np.nan, "link"],
                [np.nan, "link"],
                [np.nan, "link"],
                [np.nan, "link"],
                [np.nan, "link"],
            ],
        )
        assert_allclose(
            parcels.dataset["element_id"],
            [
                [0.0, np.nan],
                [0.0, np.nan],
                [np.nan, 2.0],
                [np.nan, 2.0],
                [np.nan, 6.0],
                [np.nan, 6.0],
                [np.nan, 6.0],
            ],
        )
        assert_allclose(
            parcels.dataset["starting_link"], [0.0, 0.0, 2.0, 2.0, 6.0, 6.0, 6.0]
        )
        assert_allclose(
            parcels.dataset["abrasion_rate"], [0.0, 0.0, 0.1, 0.1, 0.3, 0.3, 0.3]
        )
        assert_allclose(
            parcels.dataset["density"],
            [2650.0, 2650.0, 2650.0, 2650.0, 2500.0, 2500.0, 2500.0],
        )
        assert_allclose(
            parcels.dataset["time_arrival_in_link"],
            [
                [11.0, np.nan],
                [11.0, np.nan],
                [np.nan, 12.0],
                [np.nan, 12.0],
                [np.nan, 12.0],
                [np.nan, 12.0],
                [np.nan, 12.0],
            ],
        )
        assert_allclose(
            parcels.dataset["active_layer"],
            [
                [1.0, np.nan],
                [1.0, np.nan],
                [np.nan, 1.0],
                [np.nan, 1.0],
                [np.nan, 1.0],
                [np.nan, 1.0],
                [np.nan, 1.0],
            ],
        )
        assert_allclose(
            parcels.dataset["volume"],
            [
                [0.5, np.nan],
                [0.5, np.nan],
                [np.nan, 1.0],
                [np.nan, 1.0],
                [np.nan, 0.5],
                [np.nan, 0.5],
                [np.nan, 0.5],
            ],
        )

        assert np.all(
            (
                (parcels.dataset["location_in_link"] >= 0.0)
                & (parcels.dataset["location_in_link"] <= 1.0)
            )
            | np.isnan(parcels.dataset["location_in_link"])
        )
        assert np.all((parcels.dataset["D"] > 0.0) | np.isnan(parcels.dataset["D"]))

    def test_special_1(self, example_nmg2):
        """user entered time is not a pulse time, calling instance returns
        the original parcels datarecord, which is None if there is no
        original datarecord"""

        grid = example_nmg2
        time_to_pulse = time_to_pulse_list

        make_pulse = SedimentPulserAtLinks(grid, time_to_pulse=time_to_pulse, rng=1492)

        time = 11
        links = [0]
        n_parcels_at_link = [2]
        parcels = make_pulse(
            time=time, links=links, n_parcels_at_link=n_parcels_at_link
        )

        assert parcels is None


# @pytest.mark.xfail(reason = "TDD, test class is not yet implemented")
class Test_SedimentPulserEachParcel:
    def test_normal_1(self, example_nmg2):
        """minimum attributes specified in Pulse, attributes should use
        defaults specified at instantiation"""

        grid = example_nmg2

        make_pulse = SedimentPulserEachParcel(grid, rng=1975)

        PulseDF = pd.DataFrame(
            {
                "pulse_volume": [0.2, 1, 1.1, 0.5],
                "link_#": [1, 3, 5, 2],
                "normalized_downstream_distance": [0.8, 0.7, 0.5, 0.2],
            }
        )
        time = 7
        parcels = make_pulse(time, PulseDF)

        D50 = 0.05  # default grain size parameters
        D84_D50 = 2.1
        D_sd = D50 * D84_D50 - D50
        D = parcels.dataset["D"]

        assert np.all(parcels.dataset["grid_element"] == "link")
        assert np.all(
            parcels.dataset["element_id"] == [[1], [3], [3], [5], [5], [5], [2]]
        )
        assert np.all(parcels.dataset["starting_link"] == [1, 3, 3, 5, 5, 5, 2])
        assert np.all(parcels.dataset["abrasion_rate"] == approx(0.0))
        assert np.all(parcels.dataset["density"] == approx(2650.0))
        assert np.all(parcels.dataset["time_arrival_in_link"] == approx(time))
        assert np.all(parcels.dataset["active_layer"] == approx(1.0))
        assert_allclose(
            parcels.dataset["location_in_link"],
            [[0.8], [0.7], [0.7], [0.5], [0.5], [0.5], [0.2]],
        )
        assert_allclose(
            parcels.dataset["volume"], [[0.2], [0.5], [0.5], [0.5], [0.5], [0.1], [0.5]]
        )

        assert np.abs(D.mean() - D50) < 3 * D_sd

    def test_normal_2(self, example_nmg2):
        """all attributes specified in Pulse"""

        grid = example_nmg2

        make_pulse = SedimentPulserEachParcel(grid, rng=1066)

        PulseDF = pd.DataFrame(
            {
                "pulse_volume": [0.2, 1, 1.1, 0.5],
                "link_#": [1, 3, 5, 2],
                "normalized_downstream_distance": [0.8, 0.7, 0.5, 0.2],
                "D50": [0.15, 0.2, 0.22, 0.1],
                "D84_D50": [1, 1, 1, 1],
                "abrasion_rate": [0.01, 0.02, 0.005, 0.03],
                "rho_sediment": [2650, 2300, 2750, 2100],
                "parcel_volume": [0.1, 1, 1, 0.2],
            }
        )
        time = 7.0
        parcels = make_pulse(time, PulseDF)

        assert np.all(parcels.dataset["grid_element"] == "link")
        assert np.all(
            parcels.dataset["element_id"] == [[1], [1], [3], [5], [5], [2], [2], [2]]
        )
        assert np.all(parcels.dataset["starting_link"] == [1, 1, 3, 5, 5, 2, 2, 2])
        assert_allclose(
            parcels.dataset["abrasion_rate"],
            [0.01, 0.01, 0.02, 0.005, 0.005, 0.03, 0.03, 0.03],
        )
        assert_allclose(
            parcels.dataset["density"],
            [2650.0, 2650.0, 2300.0, 2750.0, 2750.0, 2100.0, 2100.0, 2100.0],
        )
        assert np.all(parcels.dataset["time_arrival_in_link"] == approx(time))
        assert np.all(parcels.dataset["active_layer"] == approx(1.0))
        assert_allclose(
            parcels.dataset["location_in_link"],
            [[0.8], [0.8], [0.7], [0.5], [0.5], [0.2], [0.2], [0.2]],
        )
        assert_allclose(
            parcels.dataset["D"],
            [[0.15], [0.15], [0.2], [0.22], [0.22], [0.1], [0.1], [0.1]],
        )
        assert_allclose(
            parcels.dataset["volume"],
            [[0.1], [0.1], [1], [1], [0.1], [0.2], [0.2], [0.1]],
        )

    def test_normal_3(self, example_nmg2):
        """two pulses. First, only minimum attributes specified,
        second, two parcels in link two and three parcels in link six
        are added and all attributes specified"""

        grid = example_nmg2
        np.random.seed(seed=5)

        make_pulse = SedimentPulserEachParcel(grid, rng=1945)

        PulseDF = pd.DataFrame(
            {
                "pulse_volume": [0.2, 1],
                "link_#": [1, 3],
                "normalized_downstream_distance": [0.8, 0.7],
            }
        )
        time = 7
        parcels = make_pulse(time, PulseDF)

        PulseDF = pd.DataFrame(
            {
                "pulse_volume": [1.1, 0.5],
                "link_#": [5, 2],
                "normalized_downstream_distance": [0.5, 0.2],
                "D50": [0.22, 0.1],
                "D84_D50": [2.1, 1.5],
                "abrasion_rate": [0.005, 0.03],
                "density": [2750, 2100],
                "parcel_volume": [1, 0.2],
            }
        )
        time = 8
        parcels = make_pulse(time, PulseDF)

        assert_array_equal(
            parcels.dataset["grid_element"].values.tolist(),
            [
                ["link", np.nan],
                ["link", np.nan],
                ["link", np.nan],
                [np.nan, "link"],
                [np.nan, "link"],
                [np.nan, "link"],
                [np.nan, "link"],
                [np.nan, "link"],
            ],
        )
        assert_allclose(
            parcels.dataset["element_id"],
            [
                [1.0, np.nan],
                [3.0, np.nan],
                [3.0, np.nan],
                [np.nan, 5.0],
                [np.nan, 5.0],
                [np.nan, 2.0],
                [np.nan, 2.0],
                [np.nan, 2.0],
            ],
        )
        assert_array_equal(parcels.dataset["starting_link"], [1, 3, 3, 5, 5, 2, 2, 2])
        assert_allclose(
            parcels.dataset["abrasion_rate"],
            [0.0, 0.0, 0.0, 0.005, 0.005, 0.03, 0.03, 0.03],
        )
        assert_allclose(parcels.dataset["density"], 2650.0)
        assert_allclose(
            parcels.dataset["time_arrival_in_link"],
            [
                [7.0, np.nan],
                [7.0, np.nan],
                [7.0, np.nan],
                [np.nan, 8.0],
                [np.nan, 8.0],
                [np.nan, 8.0],
                [np.nan, 8.0],
                [np.nan, 8.0],
            ],
        )
        assert_allclose(
            parcels.dataset["active_layer"],
            [
                [1.0, np.nan],
                [1.0, np.nan],
                [1.0, np.nan],
                [np.nan, 1.0],
                [np.nan, 1.0],
                [np.nan, 1.0],
                [np.nan, 1.0],
                [np.nan, 1.0],
            ],
        )
        assert_allclose(
            parcels.dataset["location_in_link"],
            [
                [0.8, np.nan],
                [0.7, np.nan],
                [0.7, np.nan],
                [np.nan, 0.5],
                [np.nan, 0.5],
                [np.nan, 0.2],
                [np.nan, 0.2],
                [np.nan, 0.2],
            ],
        )
        assert_allclose(
            parcels.dataset["volume"],
            [
                [0.2, np.nan],
                [0.5, np.nan],
                [0.5, np.nan],
                [np.nan, 1.0],
                [np.nan, 0.1],
                [np.nan, 0.2],
                [np.nan, 0.2],
                [np.nan, 0.1],
            ],
        )

        assert np.all((parcels.dataset["D"] > 0.0) | np.isnan(parcels.dataset["D"]))

    def test_normal_4(self, example_nmg2):
        """Series of pulses using both the SedimentPulserAtlinks and
        SedimentPulserEachParcel."""

        grid = example_nmg2

        # define the initial parcels datarecord
        make_pulse_links = SedimentPulserAtLinks(
            grid, time_to_pulse=always_time_to_pulse, rng=5
        )

        # pulse 1
        time = 7
        links = [0]
        n_parcels_at_link = [2]
        parcels = make_pulse_links(
            time=time, links=links, n_parcels_at_link=n_parcels_at_link
        )

        # datarecord is input to SedimentPulserEachParcel, now both pulser
        # instances will point to the same datarecord when called
        make_pulse_pulseDF = SedimentPulserEachParcel(
            parcels=parcels, grid=grid, rng=1863
        )

        # pulse 2
        PulseDF = pd.DataFrame(
            {
                "pulse_volume": [1.1],
                "link_#": [5],
                "normalized_downstream_distance": [0.5],
                "D50": [0.22],
                "D84_D50": [1],
                "abrasion_rate": [0.005],
                "density": [2750],
                "parcel_volume": [1],
            }
        )
        time = 8
        parcels = make_pulse_pulseDF(time, PulseDF)

        # pulse 3
        time = 9
        links = [5]
        n_parcels_at_link = [1]
        parcels = make_pulse_links(
            time=time, links=links, n_parcels_at_link=n_parcels_at_link
        )

        # pulse 4
        PulseDF = pd.DataFrame(
            {
                "pulse_volume": [0.5],
                "link_#": [6],
                "normalized_downstream_distance": [0.2],
            }
        )
        time = 10
        parcels = make_pulse_pulseDF(time, PulseDF)

        assert_array_equal(
            parcels.dataset["grid_element"].values.tolist(),
            [
                ["link", np.nan, np.nan, np.nan],
                ["link", np.nan, np.nan, np.nan],
                [np.nan, "link", np.nan, np.nan],
                [np.nan, "link", np.nan, np.nan],
                [np.nan, np.nan, "link", np.nan],
                [np.nan, np.nan, np.nan, "link"],
            ],
        )
        assert_allclose(
            parcels.dataset["element_id"],
            [
                [0.0, np.nan, np.nan, np.nan],
                [0.0, np.nan, np.nan, np.nan],
                [np.nan, 5.0, np.nan, np.nan],
                [np.nan, 5.0, np.nan, np.nan],
                [np.nan, np.nan, 5.0, np.nan],
                [np.nan, np.nan, np.nan, 6.0],
            ],
        )
        assert_array_equal(parcels.dataset["starting_link"], [0, 0, 5, 5, 5, 6])
        assert_allclose(
            parcels.dataset["abrasion_rate"], [0.0, 0.0, 0.005, 0.005, 0.0, 0.0]
        )
        assert_allclose(parcels.dataset["density"], 2650.0)
        assert_allclose(
            parcels.dataset["time_arrival_in_link"],
            [
                [7.0, np.nan, np.nan, np.nan],
                [7.0, np.nan, np.nan, np.nan],
                [np.nan, 8.0, np.nan, np.nan],
                [np.nan, 8.0, np.nan, np.nan],
                [np.nan, np.nan, 9.0, np.nan],
                [np.nan, np.nan, np.nan, 10.0],
            ],
        )
        assert_allclose(
            parcels.dataset["active_layer"],
            [
                [1.0, np.nan, np.nan, np.nan],
                [1.0, np.nan, np.nan, np.nan],
                [np.nan, 1.0, np.nan, np.nan],
                [np.nan, 1.0, np.nan, np.nan],
                [np.nan, np.nan, 1.0, np.nan],
                [np.nan, np.nan, np.nan, 1.0],
            ],
        )
        assert_allclose(
            parcels.dataset["location_in_link"],
            [
                [0.51532556, np.nan, np.nan, np.nan],
                [0.28580138, np.nan, np.nan, np.nan],
                [np.nan, 0.5, np.nan, np.nan],
                [np.nan, 0.5, np.nan, np.nan],
                [np.nan, np.nan, 0.38336888, np.nan],
                [np.nan, np.nan, np.nan, 0.2],
            ],
        )
        assert_allclose(
            parcels.dataset["volume"],
            [
                [0.5, np.nan, np.nan, np.nan],
                [0.5, np.nan, np.nan, np.nan],
                [np.nan, 1.0, np.nan, np.nan],
                [np.nan, 0.1, np.nan, np.nan],
                [np.nan, np.nan, 0.5, np.nan],
                [np.nan, np.nan, np.nan, 0.5],
            ],
        )

        # grain size must be greater than 0
        assert np.all((parcels.dataset["D"] > 0.0) | np.isnan(parcels.dataset["D"]))

    def test_bad_1(self, example_nmg2):
        """test exception raised if instance is called without specifying the
        pulserDF"""

        grid = example_nmg2
        make_pulse = SedimentPulserEachParcel(grid)
        PulseDF = None
        time = 7

        with pytest.raises(ValueError) as exc_info:
            make_pulse(time, PulseDF)
        assert exc_info.match("PulseDF was not specified")

    def test_special_1(self, example_nmg2):
        """test that calling with an empty PulseDF returns the existing
        datarecord"""

        grid = example_nmg2

        make_pulse = SedimentPulserEachParcel(grid, rng=5)

        PulseDF = pd.DataFrame(
            {
                "pulse_volume": [0.2, 1, 1.1, 0.5],
                "link_#": [1, 3, 5, 2],
                "normalized_downstream_distance": [0.8, 0.7, 0.5, 0.2],
            }
        )
        # call SedimentPulserEachParcel using setup in normal_test_1
        time1 = 7.0
        parcels = make_pulse(time1, PulseDF)

        # call again, using an empty PulserDF
        time2 = 8.0
        parcels = make_pulse(time2, pd.DataFrame([]))

        assert np.all(parcels.dataset["grid_element"] == "link")
        assert np.all(
            parcels.dataset["element_id"] == [[1], [3], [3], [5], [5], [5], [2]]
        )
        assert np.all(parcels.dataset["starting_link"] == [1, 3, 3, 5, 5, 5, 2])
        assert np.all(parcels.dataset["abrasion_rate"] == approx(0.0))
        assert np.all(parcels.dataset["time_arrival_in_link"] == approx(time1))
        assert np.all(parcels.dataset["active_layer"] == approx(1.0))
        assert_allclose(
            parcels.dataset["location_in_link"],
            [[0.8], [0.7], [0.7], [0.5], [0.5], [0.5], [0.2]],
        )
        assert_allclose(
            parcels.dataset["D"],
            [
                [0.0275786],
                [0.01871699],
                [0.04158561],
                [0.06830391],
                [0.11615185],
                [0.05423998],
                [0.03318153],
            ],
        )
        assert_allclose(
            parcels.dataset["volume"],
            np.array([[0.2], [0.5], [0.5], [0.5], [0.5], [0.1], [0.5]]),
        )
