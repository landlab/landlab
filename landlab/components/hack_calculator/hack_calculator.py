# -*- coding: utf-8 -*-
"""

@author: krb
"""
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
    >>> np.testing.assert_array_equal(
    ...     _hacks_law([1, 2], 1, 1),
    ...     np.array([1, 2]))
    >>> np.testing.assert_array_equal(
    ...     _hacks_law([1, 2], 3, 2),
    ...     np.array([ 3, 12]))
    """
    assert isinstance(C, (np.number, int, float))
    assert isinstance(h, (np.number, int, float))
    L = C * np.power(np.asarray(A), h)
    return L


def _estimate_hack_coeff(A, L):
    """Estimate Hack coefficients.

    Given A and L, estimate C and h Where

    L = C * A**h

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> from landlab.components.hack_calculator.hack_calculator import (
    ...     _estimate_hack_coeff, _hacks_law)
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


def _flatten(l):
    """
    Examples
    --------
    >>> from landlab.components.hack_calculator.hack_calculator import _flatten
    >>> struct = [[1, 2, 3, 4],
    ...           [[5, 6, 7, 9],
    ...            [9, 10, 11, 12]],
    ...           [13, 14, 15, 16]]
    >>> out = _flatten(struct)
    >>> np.testing.assert_array_equal(
    ...     out,
    ...     np.array([1, 2, 3, 4, 5, 6, 7, 9, 9, 10, 11, 12, 13, 14, 15, 16]))
    >>> assert _flatten(None) is None
    """
    if l is None:
        return None
    if not hasattr(l, "__iter__"):
        return [l]
    else:
        return list(chain(*map(_flatten, l)))


class HackCalculator(Component):
    """
    This component calculates Hack's law coefficients for drainage basins.

    Hacks law is given as

    ..:math:
        L = C * A**h

    Where :math:`L` is the distance to the drainage divide along the channel,
    :math:`A` is the drainage area, and :math:`C`and :math:`h` are
    coefficients.

    Examples
    --------
    >>> import pandas as pd
    >>> pd.set_option('display.max_columns', None)
    >>> import numpy as np
    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
    >>> from landlab.components import (
    ...     FlowAccumulator,
    ...     FastscapeEroder,
    ...     HackCalculator)
    >>> np.random.seed(42)
    >>> mg = RasterModelGrid((50, 100), xy_spacing=100)
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> z[mg.core_nodes] += np.random.randn(mg.core_nodes.size)
    >>> fa = FlowAccumulator(mg)
    >>> fs = FastscapeEroder(mg, K_sp=0.001)
    >>> for i in range(100):
    ...     fa.run_one_step()
    ...     fs.run_one_step(1000)
    ...     z[mg.core_nodes] += 0.01 * 1000
    >>> hc = HackCalculator(mg)
    >>> df = hc.calculate_hack_coefficients()
    >>> largest_outlet = mg.boundary_nodes[
    ...     np.argsort(mg.at_node['drainage_area'][mg.boundary_nodes])[-1:]][0]
    >>> largest_outlet
    4978
    >>> df.loc[largest_outlet, "A_max"]
    2830000.0
    >>> df.round(2)  # doctest: +NORMALIZE_WHITESPACE
              A_max     C     h
    4978  2830000.0  0.33  0.62

    >>> hc = HackCalculator(mg, number_of_watersheds=3, main_channel_only=False)
    >>> df = hc.calculate_hack_coefficients()
    >>> df.round(2)  # doctest: +NORMALIZE_WHITESPACE
              A_max     C     h
    39    2170000.0  0.24  0.64
    4929  2350000.0  0.25  0.63
    4978  2830000.0  0.46  0.60
    """

    _name = "HackCalculator"

    _input_var_names = (
        "topographic__elevation",
        "drainage_area",
        "flow__receiver_node",
        "flow__link_to_receiver_node",
        "flow__upstream_node_order",
    )

    _output_var_names = "distance_to_divide"

    _var_units = {
        "topographic__elevation": "m",
        "flow__receiver_node": "-",
        "drainage_area": "m**2",
        "flow__link_to_receiver_node": "-",
        "distance_to_divide": "m",
        "flow__upstream_node_order": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "flow__receiver_node": "node",
        "drainage_area": "node",
        "flow__link_to_receiver_node": "node",
        "distance_to_divide": "node",
        "flow__upstream_node_order": "node",
    }
    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current node)",
        "drainage_area": "Upstream accumulated surface area contributing to the node's discharge",
        "flow__link_to_receiver_node": "Node array containing ID of link that leads from each node to its receiver, or BAD_INDEX_VALUE if no link",
        "distance_to_divide": "Distance from drainage divide.",
        "flow__upstream_node_order": "node order such that nodes must appear in the list after all nodes downstream of them",
    }

    def __init__(
        self,
        grid,
        number_of_watersheds=1,
        main_channel_only=True,
        starting_nodes=None,
        threshold=None,
        **kwds
    ):
        """
        Parameters
        ----------
        grid : Landlab Model Grid instance, required
        number_of_watersheds : int, optional
            Total number of watersheds to calculate the Hack coefficients for.
            Default value is 1. If value is greater than 1 and starting_nodes
            is not specified, then the number_of_watersheds largest watersheds
            based on the drainage area at the model grid boundary.
        main_channel_only : Boolean, optional
            Use only the longest channel to calculate the Hack coefficients (if
            True, or use all the pixels in each watershed with drainage area
            above the threshold value).
        starting_nodes : length number_of_watersheds iterable, optional
            Length number_of_watersheds itterable containing the node IDs of
            nodes to start the channel profiles from. If not provided, the
            default is the number_of_watersheds node IDs on the model grid
            boundary with the largest terminal drainage area.
        threshold : float, optional
            Value to use for the minimum drainage area associated with a
            plotted channel segment. Default values is 2.0 x minimum grid cell
            area.
        """
        super(HackCalculator, self).__init__(grid)
        self._grid = grid
        self._profiler = ChannelProfiler(
            grid,
            number_of_watersheds=number_of_watersheds,
            main_channel_only=main_channel_only,
            starting_nodes=starting_nodes,
            threshold=threshold,
        )

    def calculate_hack_coefficients(self):
        """Calculate Hack coefficients for desired watersheds.

        Returns
        -------
        hack_structure : pandas DataFrame
            Index are the node IDs of watershed outlets where the Hack
            coefficients were estimated. They coorespond to the
            number_of_watersheds largest drainages on the model grid boundaries
            or the nodes indicated with starting_nodes. Column values are
            "A_max" for the drainage area of the watershed outlet, "C" for the
            coefficient, and "h" for the exponent.

        """
        out = collections.OrderedDict()
        self._profiler.run_one_step()

        dist = calculate_distance_to_divide(self._grid, longest_path=True)

        internal_df = []
        # for watershed in watersheds (in profile structure)
        for watershed in self._profiler._profile_structure:
            outlet_node = watershed[0][0]
            A_max = self._grid.at_node["drainage_area"][outlet_node]

            nodes = np.unique(_flatten(watershed))

            A = self._grid.at_node["drainage_area"][nodes]
            L = dist[nodes]
            C, h = _estimate_hack_coeff(A, L)

            out[outlet_node] = {"A_max": A_max, "C": C, "h": h}

            internal_df.append(
                pd.DataFrame.from_dict(
                    {
                        "basin_id": outlet_node * np.ones(A.shape),
                        "A": A,
                        "L_obs": L,
                        "L_est": C * A ** h,
                    }
                )
            )

        df = pd.DataFrame.from_dict(out, orient="index", columns=["A_max", "C", "h"])

        hdf = pd.concat(internal_df, ignore_index=True)
        amax_df = (
            hdf.drop(["L_obs", "L_est"], axis=1)
            .groupby("basin_id")
            .max()
            .sort_values(by="A")
        )
        s = pd.Series(hdf["basin_id"], dtype="category", index=hdf.index)
        s = s.cat.set_categories(amax_df.index.values, ordered=True)
        hdf["basin_id"] = s

        self._df = hdf

        return df
