# ! /usr/env/python
"""
Created on Tue Feb 27 16:25:11 2018

@author: barnhark
"""
import matplotlib
import numpy as np
import pytest

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import ChannelProfiler
from landlab.components import DepressionFinderAndRouter
from landlab.components import FastscapeEroder
from landlab.components import FlowAccumulator
from landlab.components import LinearDiffuser

matplotlib.use("agg")


def test_assertion_error():
    """Test that the correct assertion error will be raised."""
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += 200 + mg.x_of_node + mg.y_of_node + np.random.randn(mg.size("node"))

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(0, z, -9999)
    fa = FlowAccumulator(
        mg, flow_director="D8", depression_finder=DepressionFinderAndRouter
    )
    sp = FastscapeEroder(mg, K_sp=0.0001, m_sp=0.5, n_sp=1, erode_flooded_nodes=True)
    ld = LinearDiffuser(mg, linear_diffusivity=0.0001)

    dt = 100
    for _ in range(200):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        ld.run_one_step(dt=dt)
        mg.at_node["topographic__elevation"][0] -= 0.001  # Uplift

    with pytest.raises(ValueError):
        ChannelProfiler(mg, outlet_nodes=[0], number_of_watersheds=2)


def test_asking_for_too_many_watersheds():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += 200 + mg.x_of_node + mg.y_of_node
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(0, z, -9999)
    fa = FlowAccumulator(mg, flow_director="D8")
    sp = FastscapeEroder(mg, K_sp=0.0001, m_sp=0.5, n_sp=1)

    dt = 100
    for _ in range(200):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        mg.at_node["topographic__elevation"][0] -= 0.001

    with pytest.raises(ValueError):
        ChannelProfiler(mg, number_of_watersheds=3)

    with pytest.raises(ValueError):
        ChannelProfiler(mg, number_of_watersheds=None, minimum_outlet_threshold=200)


def test_no_minimum_channel_threshold():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += 200 + mg.x_of_node + mg.y_of_node + np.random.randn(mg.size("node"))

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(0, z, -9999)
    fa = FlowAccumulator(
        mg, flow_director="D8", depression_finder=DepressionFinderAndRouter
    )
    fa.run_one_step()

    profiler = ChannelProfiler(mg)

    assert profiler._minimum_channel_threshold == 0.0


def test_no_flow__link_to_receiver_node():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.add_zeros("flow__receiver_node", at="node")
    with pytest.raises(FieldError):
        ChannelProfiler(mg)


def test_no_flow__receiver_node():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.add_zeros("flow__link_to_receiver_node", at="node")
    with pytest.raises(FieldError):
        ChannelProfiler(mg)


@pytest.fixture()
def profile_example_grid():
    mg = RasterModelGrid((40, 60))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += 200 + mg.x_of_node + mg.y_of_node
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(0, z, -9999)
    fa = FlowAccumulator(mg, flow_director="D8")
    sp = FastscapeEroder(mg, K_sp=0.0001, m_sp=0.5, n_sp=1)

    dt = 100
    for _ in range(200):
        fa.run_one_step()
        sp.run_one_step(dt=dt)
        mg.at_node["topographic__elevation"][0] -= 0.001
    return mg


def test_start_away_from_boundary(profile_example_grid):
    mg = profile_example_grid
    profiler = ChannelProfiler(mg, outlet_nodes=[64], cmap="viridis_r")
    profiler.run_one_step()
    assert profiler.distance_along_profile[0][0] == 0.0


def test_plotting_and_structure(profile_example_grid):
    mg = profile_example_grid
    profiler = ChannelProfiler(
        mg,
        number_of_watersheds=1,
        main_channel_only=False,
        minimum_channel_threshold=50,
    )
    profiler.run_one_step()

    profiler.plot_profiles()
    profiler.plot_profiles_in_map_view()

    # hard to test plotting... but in April 2019 KRB visually verified that the
    # plots were correct and has hard coded in what the profile structure was.
    # correct_structure = np.array(
    correct_structure = [
        np.array([0, 61]),
        np.array([61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]),
        np.array([61, 121, 181, 241, 301, 361, 421, 481, 541, 601, 661, 721, 781]),
        np.array([71, 72, 73, 74, 75]),
        np.array([71, 131, 191, 251, 311, 371, 431, 491, 551]),
        np.array([781, 841, 901, 961, 1021]),
        np.array([781, 842, 843, 844, 845, 846, 847]),
        np.array([75, 76, 77, 78, 79]),
        np.array(
            [
                75,
                135,
                195,
                255,
                315,
                375,
                435,
                495,
                555,
                615,
                675,
                735,
                795,
                855,
                915,
                975,
                1035,
            ]
        ),
        np.array([1021, 1081, 1141, 1201, 1261, 1321]),
        np.array([1021, 1082, 1083, 1084, 1085, 1086, 1087, 1088]),
        np.array([79, 80, 81, 82, 83]),
        np.array(
            [
                79,
                139,
                199,
                259,
                319,
                379,
                439,
                499,
                559,
                619,
                679,
                739,
                799,
                859,
                919,
                979,
                1039,
                1099,
            ]
        ),
        np.array([1321, 1322, 1323, 1324]),
        np.array([1321, 1381, 1441, 1501, 1561, 1621, 1681, 1741, 1801]),
        np.array([83, 84, 85, 86]),
        np.array(
            [
                83,
                143,
                203,
                263,
                323,
                383,
                443,
                503,
                563,
                623,
                683,
                743,
                803,
                863,
                923,
                983,
            ]
        ),
        np.array([86, 87, 88, 89, 90]),
        np.array(
            [
                86,
                147,
                207,
                267,
                327,
                387,
                447,
                507,
                567,
                627,
                687,
                747,
                807,
                867,
                927,
                987,
                1047,
            ]
        ),
        np.array([90, 91, 92, 93, 94, 95]),
        np.array([90, 151, 211, 271, 331, 391, 451, 511, 571, 631, 691, 751]),
        np.array([95, 96, 97, 98]),
        np.array([95, 155, 215, 275, 335, 395, 455, 515, 575, 635]),
        np.array([98, 99, 100, 101, 102, 103]),
        np.array([98, 159, 219, 279, 339, 399, 459]),
        np.array([103, 104, 105, 106, 107, 108, 109]),
        np.array([103, 163]),
    ]
    for actual, expected in zip(profiler.nodes, correct_structure):
        np.testing.assert_array_equal(actual, expected)


def test_end_nodes_only(profile_example_grid):
    mg = profile_example_grid
    # with the same grid, test some other profiler options.
    profiler2 = ChannelProfiler(
        mg,
        number_of_watersheds=None,
        main_channel_only=True,
        minimum_outlet_threshold=3,
        minimum_channel_threshold=50,
    )
    profiler2.run_one_step()

    profiler2.plot_profiles()
    profiler2.plot_profiles_in_map_view(endpoints_only=True)


def test_different_kwargs(profile_example_grid):
    mg = profile_example_grid
    # with the same grid, test some other profiler options.
    profiler2 = ChannelProfiler(
        mg,
        number_of_watersheds=None,
        main_channel_only=True,
        minimum_outlet_threshold=3,
        minimum_channel_threshold=50,
    )
    profiler2.run_one_step()

    profiler2.plot_profiles()
    profiler2.plot_profiles_in_map_view()

    correct_structure = np.array(
        [
            0,
            61,
            62,
            63,
            64,
            65,
            66,
            67,
            68,
            69,
            70,
            71,
            72,
            73,
            74,
            75,
            76,
            77,
            78,
            79,
            80,
            81,
            82,
            83,
            84,
            85,
            86,
            87,
            88,
            89,
            90,
            91,
            92,
            93,
            94,
            95,
            96,
            97,
            98,
            99,
            100,
            101,
            102,
            103,
            104,
            105,
            106,
            107,
            108,
            109,
        ]
    )
    np.testing.assert_array_equal(profiler2.nodes[0], correct_structure)


def test_re_calculating_nodes_and_distance():
    mg = RasterModelGrid((20, 20), xy_spacing=100)
    z = mg.add_zeros("topographic__elevation", at="node")
    z += np.random.rand(z.size)
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=False,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fa = FlowAccumulator(mg, flow_director="D8")
    sp = FastscapeEroder(mg, K_sp=0.0001, m_sp=0.5, n_sp=1)

    dt = 1000
    uplift_per_step = 0.001 * dt

    for _ in range(10):
        z[mg.core_nodes] += uplift_per_step
        fa.run_one_step()
        sp.run_one_step(dt=dt)

    profiler = ChannelProfiler(mg)
    profiler.run_one_step()
    assert len(profiler.distance_along_profile) == 1  # result: 1
    profiler.run_one_step()
    # here nathan originally found result: 2, a bug!
    assert len(profiler.distance_along_profile) == 1

    # make the most complicated profile structure
    profiler = ChannelProfiler(mg, main_channel_only=False, number_of_watersheds=2)
    profiler.run_one_step()
    p1 = list(profiler.nodes)
    d1 = list(profiler.distance_along_profile)

    profiler.run_one_step()
    p2 = list(profiler.nodes)
    d2 = list(profiler.distance_along_profile)

    # assert that these are copies, not pointers to same thing
    assert p1 is not p2
    assert d1 is not d2

    # test that structures are the same.
    for idx_watershed in range(len(p1)):
        p1_w = p1[idx_watershed]
        p2_w = p2[idx_watershed]

        d1_w = d1[idx_watershed]
        d2_w = d2[idx_watershed]

        for idx_segment in range(len(p1_w)):
            np.testing.assert_array_equal(p1_w[idx_segment], p2_w[idx_segment])
            np.testing.assert_array_equal(d1_w[idx_segment], d2_w[idx_segment])


@pytest.mark.parametrize("main", [True, False])
@pytest.mark.parametrize("nshed", [1, None, 3])
def test_getting_all_the_way_to_the_divide(main, nshed):
    np.random.seed(42)
    mg = RasterModelGrid((10, 12))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += np.random.rand(z.size)

    fa = FlowAccumulator(mg, flow_director="D8")
    sp = FastscapeEroder(mg, K_sp=0.0001, m_sp=0.5, n_sp=1)

    dt = 1000
    uplift_per_step = 0.001 * dt

    for _ in range(100):
        z[mg.core_nodes] += uplift_per_step
        fa.run_one_step()
        sp.run_one_step(dt=dt)

    profiler = ChannelProfiler(
        mg,
        number_of_watersheds=nshed,
        minimum_outlet_threshold=0,
        main_channel_only=main,
        minimum_channel_threshold=0,
    )
    profiler.run_one_step()

    # assert that with minimum_channel_threshold set to zero, we get all the way
    # to the top of the divide.
    for outlet_id in profiler._data_struct:
        seg_tuples = profiler._data_struct[outlet_id].keys()

        wshd_ids = [profiler._data_struct[outlet_id][seg]["ids"] for seg in seg_tuples]

        nodes = np.concatenate(wshd_ids).ravel()
        da = mg.at_node["drainage_area"][nodes]

        # if "profile" is just bits of the edge, then da is 0.
        assert (mg.area_of_cell.min() in da) or (0.0 in da)
