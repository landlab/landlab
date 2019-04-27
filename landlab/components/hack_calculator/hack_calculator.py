# -*- coding: utf-8 -*-
"""

@author: krb
"""
import collections
from itertools import chain

import numpy as np
from scipy.optimize import curve_fit

from landlab import Component
from landlab.components.profiler.channel_profiler import ChannelProfiler
from landlab.utils.distance_to_divide import calculate_distance_to_divide


def _hacks_law(A, C, h):
    """Given A, C, and h calculate L.

    Where L = C * A**h

    Examples
    --------
    >>> from landlab.components.hack_calculator.hack_calculator import _hacks_law
    >>> _hacks_law(1, 1, 1)
    1
    >>> _hacks_law([1, 2], 1, 1)
    array([1, 2])
    >>> _hacks_law([1, 2], 3, 2)
    array([ 3, 12])
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
    >>> out
    [1, 2, 3, 4, 5, 6, 7, 9, 9, 10, 11, 12, 13, 14, 15, 16]
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
    >>> struct = hc.calculate_hack_coefficients()
    >>> largest_outlet = mg.boundary_nodes[
    ...     np.argsort(mg.at_node['drainage_area'][mg.boundary_nodes])[-1:]][0]
    >>> largest_outlet
    4978
    >>> struct[largest_outlet]["A_max"]
    2830000.0
    >>> np.testing.assert_almost_equal(struct[largest_outlet]["C"], 0.32742, decimal=4)
    >>> np.testing.assert_almost_equal(struct[largest_outlet]["h"], 0.61849, decimal=4)
    >>> hc = HackCalculator(mg, number_of_watersheds=3, main_channel_only=False)
    >>> struct = hc.calculate_hack_coefficients()
    >>> struct
    OrderedDict([(39, {'A_max': 2170000.0,
                       'C': 0.23102443230088521,
                       'h': 0.64388645022463153}),
                 (4929, {'A_max': 2350000.0,
                         'C': 0.16200093663042728,
                         'h': 0.66097209310822003}),
                 (4978, {'A_max': 2830000.0,
                         'C': 0.50944655173089382,
                         'h': 0.58851304479736888})])
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
        starting_nodes : length number_of_watersheds itterable, optional
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
        hack_structure : dict
            Dictionary keys are the node IDs of watersheds where Hack
            coefficients were estimated. Values are a dictionary of form:
            {"A_max": A_max, "C": C, "h": h} where A_max is the drainage area
            of the watershed outlet, C is the estimate of the Hack coefficient
            and h is the estimate of the Hack outlet (for that watershed).

        """
        out = collections.OrderedDict()
        self._profiler.run_one_step()

        dist = calculate_distance_to_divide(self._grid, longest_path=True)

        # for watershed in watersheds (in profile structure)
        for watershed in self._profiler._profile_structure:
            outlet_node = watershed[0][0]
            A_max = self._grid.at_node["drainage_area"][outlet_node]

            nodes = _flatten(watershed)

            C, h = _estimate_hack_coeff(
                self._grid.at_node["drainage_area"][nodes], dist[nodes]
            )
            out[outlet_node] = {"A_max": A_max, "C": C, "h": h}

        return out
