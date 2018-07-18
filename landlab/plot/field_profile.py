from landlab import FieldError
from landlab.plot.imshow import imshow_grid
import matplotlib.pyplot as plt
import numpy as np


class FieldProfiler:
    """Create a profile of field values at model grid nodes.

    The profile is laid out by `points` that act as the vertices, including
    endpoints, of the profile trace. Sections of the profile between the points
    are referred to as `segments`. The profile trace will be one segment
    when `points` has exactly two elements (making a straight line) . The
    profile will have multiple segments when `points` has greater than two
    elements and the trace can be more complex than a straight line.

    Boundary nodes will be included in the profile if these nodes are samples
    along the trace. Avoid the boundary when selecting inputs points if these
    nodes should not be included in the profile.

    The `field` of `grid` is sampled either at the resolution of `grid` (the
    default), or if set, at the length interval of `sample_spacing`. For
    example, a sample spacing of 4 will have this structure:

    start     end   , where `-` is a length equal to sample spacing and `+` is
      +---+---+       the location of a sample along the profile.

    The spacing between the last two points may be less than the spacing
    between all other points when the profile length is not long enough for a
    complete final segment:

    start       end
      +---+---+--+

    The sample length extends across segments. If the sample length does not
    coincide with the end of a segment, it will continue into the next segment:

    start
      +---+--
            |
            + end
    """

    def __init__(self, grid, field, points, sample_spacing=None):
        """Instantiate FieldProfiler.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        field : string
            The `grid` field at node in which to extract the profile.
        points : tuple list or integer list
            Each element of the list is a vertex of the profile trace. The
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
        if field not in grid.at_node.keys():
            raise FieldError('the field, {} must be a field of the input grid '
                             'to create a profile'.format(field))

        if not isinstance(points, list) or len(points) < 2:
            raise ValueError('`points` must be a list of at least 2 elements')

        # Store inputs.

        self._grid = grid
        self._field = field

        # Set sample delta distance, dd.

        if sample_spacing == None:
            dd = grid.dx
        else:
            dd = sample_spacing

        # Create vectors for x and y coordinates of input points.

        n_points = len(points)
        x = np.empty(n_points)
        y = np.empty(n_points)

        # Track profile cumulative distance by iteratively adding segment
        # lengths.

        cum_distance = 0

        # The sample length remainder is tracked in order to apply it to the
        # subsequent segment.

        seg_sample_len_remainder = 0

        # Create lists for profile data.

        self._distance = []
        self._field_value = []
        self._sample_x = []
        self._sample_y = []
        self._nodes = []

        # Process the profile by segment. Segments are bound by successive
        # elements in `points`.

        for i, p in enumerate(points):
            # Parse segment point x and y coordinates by input type.

            if isinstance(p, (float, int)):
                x[i] = grid.x_of_node[p]
                y[i] = grid.y_of_node[p]
            elif isinstance(p, tuple) and len(p) == 2:
                x[i] = p[0]
                y[i] = p[1]
            else:
                raise TypeError('each element of `points` must be a number '
                                'representing a node id or a tuple of node x, '
                                'y coordinates')

            if i > 0:
                # `h` is segment start index and `i` is segment end index.
                h = i - 1

                # Calculate total segment length.

                seg_len = np.sqrt((x[i] - x[h])**2 + (y[i] - y[h])**2)

                # Calculate profile-long distance for segment samples,
                # `d_samples`.

                if seg_sample_len_remainder == 0:
                    d_min = 0
                    n_samples = int(np.floor((seg_len - d_min) / dd)) + 1
                    d_max = (n_samples - 1) * dd

                else:
                    d_min = dd - seg_sample_len_remainder
                    n_samples = int(np.floor((seg_len - d_min) / dd)) + 1
                    d_max = n_samples * dd - seg_sample_len_remainder

                d_samples = np.linspace(d_min, d_max, n_samples)

                # Remove duplicate sample. This may be necessary when the prior
                # segment had a sample at its endpoint.

                if d_samples[0] + cum_distance in self._distance:
                    d_samples = np.delete(d_samples, 0)

                # Include profile endpoint.

                processing_final_segment = i == n_points - 1
                point_not_yet_included = seg_len not in d_samples

                if processing_final_segment and point_not_yet_included:
                    d_samples = np.append(d_samples, seg_len)

                # Get segment sample coordinates.

                d_norm = d_samples / seg_len
                x_seg_samples = x[h] + d_norm * (x[i] - x[h])
                y_seg_samples = y[h] + d_norm * (y[i] - y[h])

                # Adjust first segment sample given the remainder sample
                # length from prior segment. Also, update remainder variable.

                d_samples += seg_sample_len_remainder

                seg_sample_len_remainder = seg_len - d_max

                # Get the nodes nearest to the sample coordinates in order to
                # get the field values at these nodes.

                seg_nodes = grid.find_nearest_node((x_seg_samples,
                                                    y_seg_samples))

                # Update profile lists.

                self._distance.extend(d_samples + cum_distance)
                self._field_value.extend(grid.at_node[field][seg_nodes])
                self._sample_x.extend(x_seg_samples)
                self._sample_y.extend(y_seg_samples)
                self._nodes.extend(seg_nodes)

                # Add the segment length to the profile cumulative distance.

                cum_distance += d_max

        self._trace_x = x
        self._trace_y = y

    def plot(self, distance_unit_label, field_y_axis_label, **kwds):
        """Plot the grid field profile.

        The figure has two subplots:
        1.  The trace of the profile plotted on the grid shaded by the field
            values.
        2.  The profile illustrating the profile samples (distance vs the
            field value).

        Parameters
        ----------
        distance_unit_label : string
            The distance unit label for the profile figure x-axis.
        field_y_axis_label : string
            The field value label for the profile figure y-axis.
        **kwds : tuple or integer
            The same parameters in :func:`imshow_grid_at_node` can be used
            in this method.
        """
        fig, axes = plt.subplots(2, 1)

        # Plot field values as grid.

        z = self._grid.at_node[self._field]
        z_min = min(z[z != -9999])
        z_max = max(z[z != -9999])

        # Set colorbar label to the profile y axis label unless specified in
        # kwds.

        if 'colorbar_label' not in kwds:
            kwds['colorbar_label'] = field_y_axis_label

        plt.sca(axes[0])
        imshow_grid(self._grid, self._field, limits=(z_min, z_max), **kwds)

        # Plot profile trace on grid.

        axes[0].plot(self.trace_x, self.trace_y, 'k')

        axes[0].plot(self.sample_x[0], self.sample_y[0], 'mo', label='start')
        axes[0].plot(self.sample_x[-1], self.sample_y[-1], 'mx', label='end')

        axes[0].legend(bbox_to_anchor=(0.5, 1.3), loc=9, borderaxespad=0,
            ncol=2)

        # Plot profile.

        d = self.distance
        v = self.field_value

        axes[1].plot(d, v, 'k')
        axes[1].plot(d[0], v[0], 'mo', d[-1], v[-1], 'mx')

        axes[1].set_xlabel('Distance ({})'.format(distance_unit_label))
        axes[1].set_ylabel(field_y_axis_label)

        plt.tight_layout()

    @property
    def distance(self):
        """Get distance along the profile.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy.random import random_sample
        >>> shape = (3, 3)
        >>> mg = RasterModelGrid(shape)
        >>> mg.at_node['random_value'] = random_sample(shape)
        >>> fp = FieldProfiler(mg, 'random_value', [3, 5])

        The profile distance values and grid column x-coordinates will be equal
        at the nodes where the profile crosses given that the profile sample
        distance grid resolution and are 1.
        >>> fp.distance == mg.x_of_node[mg.y_of_node == 1]
        array([ True,  True,  True], dtype=bool)
        """
        return np.array(self._distance)

    @property
    def field_value(self):
        """Get the field value along the profile.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy.random import random_sample
        >>> shape = (3, 3)
        >>> mg = RasterModelGrid(shape)
        >>> mg.at_node['random_value'] = random_sample(shape)
        >>> fp = FieldProfiler(mg, 'random_value', [3, 5])

        Arrays of profile and grid field values will match at the nodes where
        the profile crosses.
        >>> fp.field_value == mg.at_node['random_value'][mg.y_of_node == 1]
        array([ True,  True,  True], dtype=bool)
        """
        return np.array(self._field_value)

    @property
    def sample_x(self):
        """Get the x-coordinate along the profile.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy.random import random_sample
        >>> shape = (3, 3)
        >>> mg = RasterModelGrid(shape)
        >>> mg.at_node['random_value'] = random_sample(shape)
        >>> fp = FieldProfiler(mg, 'random_value', [3, 5])

        Arrays of profile sample and grid column x-coordinates will match at
        the nodes where the profile crosses.
        >>> fp.sample_x == mg.x_of_node[mg.y_of_node == 1]
        array([ True,  True,  True], dtype=bool)
        """
        return np.array(self._sample_x)

    @property
    def sample_y(self):
        """Get the y-coordinate along the profile.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy.random import random_sample
        >>> shape = (3, 3)
        >>> mg = RasterModelGrid(shape)
        >>> mg.at_node['random_value'] = random_sample(shape)
        >>> fp = FieldProfiler(mg, 'random_value', [3, 5])

        Arrays of profile sample and grid column y-coordinates will match at
        the nodes where the profile crosses.
        >>> fp.sample_y == mg.y_of_node[mg.y_of_node == 1]
        array([ True,  True,  True], dtype=bool)
        """
        return np.array(self._sample_y)

    @property
    def nodes(self):
        """Get the nodes along the profile.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy.random import random_sample
        >>> shape = (3, 3)
        >>> mg = RasterModelGrid(shape)
        >>> mg.at_node['random_value'] = random_sample(shape)
        >>> fp = FieldProfiler(mg, 'random_value', [3, 5])

        Arrays of profile and grid node ids will match at the nodes where the
        profile crosses.
        >>> grid_nodes = mg.nodes.flatten()
        >>> fp.nodes == grid_nodes[mg.y_of_node == 1]
        array([ True,  True,  True], dtype=bool)
        """
        return np.array(self._nodes)

    @property
    def trace_x(self):
        """Get the x-coordinate along the profile trace.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy.random import random_sample
        >>> shape = (3, 3)
        >>> mg = RasterModelGrid(shape)
        >>> mg.at_node['random_value'] = random_sample(shape)
        >>> profile_nodes = [3, 5]
        >>> fp = FieldProfiler(mg, 'random_value', profile_nodes)

        The x-coordinate of the profile trace and grid coordinates of the input
        profile nodes will be equal.
        >>> fp.trace_x == mg.x_of_node[profile_nodes]
        array([ True,  True], dtype=bool)
        """
        return np.array(self._trace_x)

    @property
    def trace_y(self):
        """Get the y-coordinate along the profile trace.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.plot import FieldProfiler
        >>> from numpy.random import random_sample
        >>> shape = (3, 3)
        >>> mg = RasterModelGrid(shape)
        >>> mg.at_node['random_value'] = random_sample(shape)
        >>> profile_nodes = [3, 5]
        >>> fp = FieldProfiler(mg, 'random_value', profile_nodes)

        The y-coordinate of the profile trace and grid coordinates of the input
        profile nodes will be equal.
        >>> fp.trace_y == mg.y_of_node[profile_nodes]
        array([ True,  True], dtype=bool)
        """
        return np.array(self._trace_y)
