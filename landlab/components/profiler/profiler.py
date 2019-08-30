# coding: utf8
# ! /usr/env/python
"""profiler.py component to create user-defined profiles."""
from collections import OrderedDict

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from landlab import FieldError
from landlab.components.profiler.base_profiler import _BaseProfiler


class Profiler(_BaseProfiler):
    """Extract and plot profiles defined by points within a grid."""

    _name = "Profiler"

    def __init__(self, grid, field, trace_nodes, sample_spacing=None,
                 cmap="viridis"):
        """Instantiate Profiler.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        field : string
            The `grid` field at nodes in which to extract the profile.
        trace_nodes : tuple list or integer list
            Each element of the list is a node along the profile trace. The
            first and last nodes in this list are the trace endpoints. The
            second to penultimate nodes in this list are trace vertices.
            elements can be tuples of coordinates (x, y) and/or node ids.
        sample_spacing : float, integer, or none
            The distance along the profile to sample `field`. The default value
            of `none` will set the sample spacing to the `dx` of `grid`.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy import array
        >>> mg = RasterModelGrid((7, 7), 10)
        >>> z = array([
        ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
        ...     -9999.,    26.,     0.,    30.,    32.,    34., -9999.,
        ...     -9999.,    28.,     1.,    25.,    28.,    32., -9999.,
        ...     -9999.,    30.,     3.,     3.,    11.,    34., -9999.,
        ...     -9999.,    32.,    11.,    25.,    18.,    38., -9999.,
        ...     -9999.,    34.,    32.,    34.,    36.,    40., -9999.,
        ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.])
        >>> mg.at_node['topographic__elevation'] = z
        >>> mg.set_watershed_boundary_condition_outlet_id(2, z,
        ...     nodata_value=-9999.)

        Create profile with tuples of endpoint (ep) coordinates.
        >>> ep0 = (9.5, 28.1)
        >>> ep1 = (48.0, 29.5)
        >>> field = 'topographic__elevation'
        >>> fp = FieldProfiler(mg, field, [ep0, ep1])

        Arrays of profile and grid field values will match at the nodes where
        the profile crosses.
        >>> row_mask = mg.y_of_node[mg.core_nodes] == 3 * mg.dy
        >>> fp.field_value == mg.at_node[field][mg.core_nodes][row_mask]
        array([ True,  True,  True,  True,  True], dtype=bool)
        """
        super(_BaseProfiler, self).__init__(grid)

        if field not in grid.at_node.keys():
            raise FieldError('the field, {} must be a field of the input grid '
                             'at nodes to create a profile'.format(field))

        if not isinstance(trace_nodes, list) or len(trace_nodes) < 2:
            msg = '`trace_nodes` must be a list of at least 2 nodes'
            raise ValueError(msg)

        # Store inputs.

        self._grid = grid
        self._field = field
        self._trace_nodes = trace_nodes
        self._cmap = plt.get_cmap(cmap)

        if sample_spacing == None:
            self._spacing = self._grid.dx
        else:
            self._spacing = sample_spacing

    def _create_profile_structure(self):
        """ """
        # Track profile cumulative distance by iteratively adding segment
        # lengths.

        cum_distance = 0

        # The sample length remainder is tracked to apply it to the beginning
        # of the next segment next for the purpose of sample spacing.

        seg_sample_len_remainder = 0

        # Process the profile by segment. Segments are bound by successive
        # elements in `points`.

        self._net_struct = OrderedDict()
        n_points = len(self._trace_nodes)

        for i, node in enumerate(self._trace_nodes):
            x, y = self._grid.xy_of_node[node]

            if i == 0:
                continue

            # `h` is segment start index and `i` is segment end index.
            h = i - 1
            xh, yh = self._grid.xy_of_node[self._trace_nodes[h]]

            # Calculate total segment length.

            seg_len = np.sqrt((x - xh)**2 + (y - yh)**2)

            # Calculate `d_samples`: the distance values along the segment.

            if seg_sample_len_remainder == 0:
                d_min = 0
                n_samples = int(np.floor((seg_len - d_min) / self._spacing)) + 1
                d_max = (n_samples - 1) * self._spacing
                print(1, d_min, n_samples, d_max)

            else:
                d_min = self._spacing - seg_sample_len_remainder
                n_samples = int(np.floor((seg_len - d_min) / self._spacing)) + 1
                d_max = n_samples * self._spacing - seg_sample_len_remainder
                print(2, d_min, n_samples, d_max)

            d_samples = np.linspace(d_min, d_max, n_samples)

            # Ensure segment endpoint is included even when its distance is not
            # a multiple of sample spacing.

            if seg_len not in d_samples:
                d_samples = np.append(d_samples, seg_len)

            # Get segment sample coordinates.

            d_norm = d_samples / seg_len
            x_sample = xh + d_norm * (x - xh)
            y_sample = yh + d_norm * (y - yh)

            # Adjust first segment sample given the remainder sample length
            # from prior segment. Then update segment remainder.

            d_samples += seg_sample_len_remainder

            seg_sample_len_remainder = seg_len - d_max

            # Set the network structure for the segment.

            seg_nodes = self._grid.find_nearest_node((x_sample, y_sample))

            self._net_struct[i] = {'ids': seg_nodes,
                                   'distances': d_samples + cum_distance}

            cum_distance += d_max

        self._assign_colors()
        self._create_flat_structures()

    def _assign_colors(self, color_mapping=None):
        """Assign a unique color for each watershed.

        Parameters
        ----------
        color_mapping : str
            Color map name.
        """

        if color_mapping is None:
            segment_count = len(self._net_struct)
            norm = mpl.colors.Normalize(vmin=0, vmax=segment_count)
            mappable = cm.ScalarMappable(norm=norm, cmap=self._cmap)
            color_mapping = {
                segment_id: mappable.to_rgba(idx)
                for idx, segment_id in enumerate(self._net_struct)
            }

        for segment_id in self._net_struct:
            self._net_struct[segment_id]["color"] = color_mapping[
                segment_id
            ]

    def _create_flat_structures(self):
        """Create expected flattened structures for ids, distances, and colors.
        """
        self._net_ids = []
        self._distance_along_profile = []
        self._colors = []

        for segment_id in self._net_struct:
            self._net_ids.append(self._net_struct[segment_id]["ids"])
            self._distance_along_profile.append(
                    self._net_struct[segment_id]["distances"])
            self._colors.append(self._net_struct[segment_id]["color"])
