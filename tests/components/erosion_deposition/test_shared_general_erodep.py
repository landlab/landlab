import numpy as np
import pytest
from numpy import testing

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import SharedStreamPower


@pytest.mark.parametrize("frac_fines", (2.0, -0.5))
def test_fraction_fines_out_of_range(frac_fines):
    grid = RasterModelGrid((5, 5), xy_spacing=10.0)
    with pytest.raises(ValueError):
        SharedStreamPower(grid, F_f=frac_fines)


def test_q_as_field():
    """
    Test that passing in water discharge as a grid field results in self.q
    holding correct values
    """
    mg = RasterModelGrid((5, 5), xy_spacing=10.0)
    mg.add_ones("topographic__elevation", at="node")

    mg.at_node["user_imposed_discharge"] = np.random.rand(mg.number_of_nodes)
    expected = mg.at_node["user_imposed_discharge"].copy()

    FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    ed = SharedStreamPower(mg, discharge_field="user_imposed_discharge")

    testing.assert_array_equal(
        ed._q, expected, err_msg="E/D discharge field test failed", verbose=True
    )


def test_q_as_array():
    """
    Test that passing in water discharge as an array results in self.q
    holding correct values
    """
    mg = RasterModelGrid((5, 5), xy_spacing=10.0)
    mg.add_ones("topographic__elevation", at="node")

    discharge = np.random.rand(mg.number_of_nodes)
    expected = discharge.copy()

    FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    ed = SharedStreamPower(mg, discharge_field=discharge)

    testing.assert_array_equal(
        ed._q,
        expected,
        err_msg="E/D discharge array test failed",
        verbose=True,
    )


def test_sediment__outflux_already_created():
    """
    Test that an existing sediment flux grid field is not changed by
    instantiating SharedStreamPower.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=10.0)

    mg.add_ones("topographic__elevation", at="node")
    mg.at_node["sediment__outflux"] = np.random.rand(mg.number_of_nodes)
    expected = mg.at_node["sediment__outflux"].copy()

    FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    ed = SharedStreamPower(mg)

    testing.assert_array_equal(
        ed._qs,
        expected,
        err_msg="E/D sediment flux field test failed",
        verbose=True,
    )
