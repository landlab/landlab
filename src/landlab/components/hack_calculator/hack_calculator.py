"""Calculate Hack parameters."""

import collections
from itertools import chain

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from landlab import Component
from landlab.components.profiler.channel_profiler import ChannelProfiler
from landlab.utils.distance_to_divide import calculate_distance_to_divide


def _hacks_law(A, C, h):
    """Given A, C, and h calculate L.

    Where L = C * A**h

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.hack_calculator.hack_calculator import _hacks_law
    >>> _hacks_law(1, 1, 1)
    1
    >>> np.testing.assert_array_equal(_hacks_law([1, 2], 1, 1), np.array([1, 2]))
    >>> np.testing.assert_array_equal(_hacks_law([1, 2], 3, 2), np.array([3, 12]))
    """
    assert isinstance(C, (np.number, int, float))
    assert isinstance(h, (np.number, int, float))
    L = C * np.power(np.asarray(A), h)
    return L


def _estimate_hack_coeff(A, L):
    """Estimate Hack parameters.

    Given A and L, estimate C and h Where

    L = C * A**h

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> from landlab.components.hack_calculator.hack_calculator import (
    ...     _estimate_hack_coeff,
    ...     _hacks_law,
    ... )
    >>> C = 0.5
    >>> h = 0.75
    >>> A = np.arange(1, 1000)
    >>> L = _hacks_law(A, C, h) + np.random.randn(A.size)
    >>> C_hat, h_hat = _estimate_hack_coeff(A, L)
    >>> np.round(C_hat, decimals=3)
    0.497
    >>> np.round(h_hat, decimals=3)
    0.751
    """
    popt, pcov = curve_fit(_hacks_law, A, L, (0.5, 0.7))
    return popt


def _flatten(list_):
    """
    Examples
    --------
    >>> from landlab.components.hack_calculator.hack_calculator import _flatten
    >>> struct = [[1, 2, 3, 4], [[5, 6, 7, 9], [9, 10, 11, 12]], [13, 14, 15, 16]]
    >>> out = _flatten(struct)
    >>> np.testing.assert_array_equal(
    ...     out, np.array([1, 2, 3, 4, 5, 6, 7, 9, 9, 10, 11, 12, 13, 14, 15, 16])
    ... )
    >>> assert _flatten(None) is None
    """
    if list_ is None:
        return None
    if not hasattr(list_, "__iter__"):
        return [list_]
    else:
        return list(chain(*map(_flatten, list_)))


class HackCalculator(Component):
    """This component calculates Hack's law parameters for drainage basins.

    Hacks law is given as

    ..:math:
        L = C * A**h

    Where :math:`L` is the distance to the drainage divide along the channel,
    :math:`A` is the drainage area, and :math:`C`and :math:`h` are
    parameters.

    The HackCalculator uses a ChannelProfiler to determine the nodes on which
    to calculate the parameter fit.

    Examples
    --------
    >>> import pandas as pd
    >>> pd.set_option("display.max_columns", None)
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, FastscapeEroder, HackCalculator
    >>> np.random.seed(42)
    >>> mg = RasterModelGrid((50, 100), xy_spacing=100)
    >>> z = mg.add_zeros("node", "topographic__elevation")
    >>> z[mg.core_nodes] += np.random.randn(mg.core_nodes.size)
    >>> fa = FlowAccumulator(mg)
    >>> fs = FastscapeEroder(mg, K_sp=0.001)
    >>> for i in range(100):
    ...     fa.run_one_step()
    ...     fs.run_one_step(1000)
    ...     z[mg.core_nodes] += 0.01 * 1000
    ...
    >>> hc = HackCalculator(mg)
    >>> hc.calculate_hack_parameters()
    >>> largest_outlet = mg.boundary_nodes[
    ...     np.argsort(mg.at_node["drainage_area"][mg.boundary_nodes])[-1:]
    ... ][0]
    >>> largest_outlet
    4978
    >>> hc.hack_coefficient_dataframe.loc[largest_outlet, "A_max"]
    2830000.0
    >>> hc.hack_coefficient_dataframe.round(2)
    A_max     C          h     basin_outlet_id
    4978      2830000.0  0.31  0.62

    >>> hc = HackCalculator(
    ...     mg, number_of_watersheds=3, main_channel_only=False, save_full_df=True
    ... )
    >>> hc.calculate_hack_parameters()
    >>> hc.hack_coefficient_dataframe.round(2)
    A_max     C          h     basin_outlet_id
    39        2170000.0  0.13  0.69
    4929      2350000.0  0.13  0.68
    4978      2830000.0  0.23  0.64
    >>> hc.full_hack_dataframe.head().round(2)
    basin_outlet_id    A     L_obs      L_est   node_id
    39                 39.0  2170000.0  3200.0  2903.43
    139                39.0  2170000.0  3100.0  2903.43
    238                39.0    10000.0     0.0    71.61
    239                39.0  2160000.0  3000.0  2894.22
    240                39.0    10000.0     0.0    71.61

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Hack, J. T. Studies of longitudinal stream profiles in Virginia and
    Maryland (Vol. 294). U.S. Geological Survey Professional Paper 294-B (1957).
    https://doi.org/10.3133/pp294B

    """

    _name = "HackCalculator"

    _unit_agnostic = True

    _info = {
        "distance_to_divide": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Distance from drainage divide.",
        },
        "drainage_area": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(self, grid, save_full_df=False, **kwds):
        """
        Parameters
        ----------
        grid : Landlab Model Grid instance, required
        save_full_df: bool
            Flag indicating whether to create the ``full_hack_dataframe``.
        **kwds :
            Values to pass to the ChannelProfiler.
        """
        super().__init__(grid)
        super().initialize_output_fields()
        self._dist = grid.at_node["distance_to_divide"]

        self.profiler = ChannelProfiler(grid, **kwds)
        self._save_full_df = save_full_df

    @property
    def hack_coefficient_dataframe(self):
        """Hack coefficient dataframe.

        This dataframe is created and stored on the component.

        It is a pandas dataframe with one row for each basin for which Hack
        parameters are calculated. Thus, there are as many rows as the
        number of watersheds identified by the ChannelProfiler.

        The dataframe has the following index and columns.

            * Index
                * **basin_outlet_id**: The node ID of the watershed outlet
                  where each set of Hack parameters was estimated.

            * Columns
                * **A_max**: The drainage area of the watershed outlet.
                * **C**: The Hack coefficient as defined in the equations above.
                * **h**: The Hack exponent as defined in the equations above.
        """
        if hasattr(self, "_df"):
            return self._df
        else:
            raise RuntimeError(
                "The hack_coefficient_dataframe does not yet exist. "
                "Try running calculate_hack_parameters"
            )

    @property
    def full_hack_dataframe(self):
        """Full Hack calculation dataframe.

        This dataframe is optionally created and stored on the component when
        the keyword argument ``full_hack_dataframe=True`` is passed to the
        component init.

        It is pandas dataframe with a row for every model grid cell used to
        estimate the Hack parameters. It has the following index and columns.

            * Index
                * *node_id**: The node ID of the model grid cell.

            * Columns
                * **basin_outlet_id**: The node IDs of watershed outlet
                * **A**: The drainage are of the model grid cell.
                * **L_obs**: The observed distance to the divide.
                * **L_est**: The predicted distance to divide based on the
                  Hack coefficient fit.
        """
        if not self._save_full_df:
            raise NotImplementedError(
                "This instance of a HackCalculator was not set up to save "
                "the full_hack_dataframe. Try recreating it with "
                "save_full_df=True."
            )
        else:
            if hasattr(self, "_full_df"):
                return self._full_df
            else:
                raise RuntimeError(
                    "The full_hack_dataframe does not yet exist. "
                    "Try running calculate_hack_parameters"
                )

    def calculate_hack_parameters(self):
        """Calculate Hack parameters for desired watersheds."""
        out = collections.OrderedDict()
        self.profiler.run_one_step()

        self._dist[:] = calculate_distance_to_divide(self._grid, longest_path=True)

        if self._save_full_df:
            internal_df = []

        # for watershed in watersheds (in profile structure)
        for outlet_node in self.profiler._data_struct:
            seg_tuples = self.profiler._data_struct[outlet_node].keys()

            watershed = [
                self.profiler._data_struct[outlet_node][seg]["ids"]
                for seg in seg_tuples
            ]

            A_max = self._grid.at_node["drainage_area"][outlet_node]

            nodes = np.unique(_flatten(watershed))

            A = self._grid.at_node["drainage_area"][nodes]
            L = self._dist[nodes]
            C, h = _estimate_hack_coeff(A, L)

            out[outlet_node] = {"A_max": A_max, "C": C, "h": h}

            if self._save_full_df:
                internal_df.append(
                    pd.DataFrame.from_dict(
                        {
                            "basin_outlet_id": outlet_node * np.ones(A.shape),
                            "A": A,
                            "L_obs": L,
                            "L_est": C * A**h,
                            "node_id": nodes,
                        }
                    )
                )

        self._df = pd.DataFrame.from_dict(
            out, orient="index", columns=["A_max", "C", "h"]
        )
        self._df.index.name = "basin_outlet_id"

        if self._save_full_df:
            hdf = (
                pd.concat(internal_df, ignore_index=True)
                .set_index("node_id")
                .sort_index()
            )
            self._full_df = hdf
