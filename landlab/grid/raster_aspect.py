#! /usr/bin/env python
"""Calculate slope aspects on a :any:`RasterModelGrid`."""
import numpy as np

from landlab.utils.decorators import deprecated


def _one_line_slopes(input_array, grid, vals):
    node = input_array[0]
    diagonals = input_array[5:]
    neighbors = input_array[1:5]

    if not grid.status_at_node[node] == 0:
        raise IndexError("One or more of the provided nodes was closed!")

    try:
        slope_we = (
            (vals[diagonals[1]] + 2.0 * vals[neighbors[2]] + vals[diagonals[2]])
            - (vals[diagonals[0]] + 2.0 * vals[neighbors[0]] + vals[diagonals[3]])
        ) / (8.0 * grid.dx)
        slope_sn = (
            (vals[diagonals[2]] + 2.0 * vals[neighbors[3]] + vals[diagonals[3]])
            - (vals[diagonals[1]] + 2.0 * vals[neighbors[:, 1]] + vals[diagonals[0]])
        ) / (8.0 * grid.dy)
        return slope_we, slope_sn
    except IndexError:
        C = vals[node]
        weighting_verticals = 4.0
        weighting_horizontals = 4.0
        try:
            vertical_grad = (vals[neighbors[3]] - vals[neighbors[1]]) / (2.0 * grid.dy)
        except IndexError:
            try:
                vertical_grad = (C - vals[neighbors[1]]) / grid.dy
            except IndexError:
                try:
                    vertical_grad = (vals[neighbors[3]] - C) / grid.dy
                except IndexError:
                    vertical_grad = 0.0
                    weighting_verticals -= 2.0
        try:
            horizontal_grad = (vals[neighbors[2]] - vals[neighbors[0]]) / (
                2.0 * grid.dx
            )
        except IndexError:
            try:
                horizontal_grad = (C - vals[neighbors[0]]) / grid.dx
            except IndexError:
                try:
                    horizontal_grad = (vals[neighbors[2]] - C) / grid.dx
                except IndexError:
                    horizontal_grad = 0.0
                    weighting_horizontals -= 2.0
        try:
            left_grad = (vals[diagonals[2]] - vals[diagonals[1]]) / (2.0 * grid.dx)
        except IndexError:
            try:
                C = vals[neighbors[2]]
            except IndexError:
                left_grad = 0.0
                weighting_verticals -= 1.0
            else:
                try:
                    left_grad = (C - vals[diagonals[1]]) / grid.dx
                except IndexError:
                    left_grad = (vals[diagonals[2]] - C) / grid.dx
        try:
            right_grad = (vals[diagonals[3]] - vals[diagonals[0]]) / (2.0 * grid.dx)
        except IndexError:
            try:
                C = vals[neighbors[0]]
            except IndexError:
                right_grad = 0.0
                weighting_verticals -= 1.0
            else:
                try:
                    right_grad = (C - vals[diagonals[0]]) / grid.dx
                except IndexError:
                    right_grad = (vals[diagonals[3]] - C) / grid.dx
        try:
            top_grad = (vals[diagonals[1]] - vals[diagonals[0]]) / (2.0 * grid.dy)
        except IndexError:
            try:
                C = vals[neighbors[1]]
            except IndexError:
                top_grad = 0.0
                weighting_horizontals -= 1.0
            else:
                try:
                    top_grad = (C - vals[diagonals[0]]) / grid.dy
                except IndexError:
                    top_grad = (vals[diagonals[1]] - C) / grid.dy
        try:
            bottom_grad = (vals[diagonals[2]] - vals[diagonals[3]]) / (2.0 * grid.dy)
        except IndexError:
            try:
                C = vals[neighbors[3]]
            except IndexError:
                bottom_grad = 0.0
                weighting_horizontals -= 1.0
            else:
                try:
                    bottom_grad = (C - vals[diagonals[3]]) / grid.dy
                except IndexError:
                    bottom_grad = (vals[diagonals[2]] - C) / grid.dy

        slope_we = (
            top_grad + 2.0 * horizontal_grad + bottom_grad
        ) / weighting_horizontals
        slope_sn = (left_grad + 2.0 * vertical_grad + right_grad) / weighting_verticals

        return slope_we, slope_sn


@deprecated(use="grid.calc_slope_at_node", version=1.0)
def calc_slope_aspect_of_nodes_horn(grid, ids=None, vals="topographic__elevation"):
    r"""Calculate slope and aspect.

    .. note::

        THIS CODE HAS ISSUES (SN 25-Sept-14): This code didn't perform well
        on a NS facing elevation profile. Please check
        slope_aspect_routines_comparison.py under landlab\examples before
        using this.  Suggested alternative:
        calculate_slope_aspect_at_nodes_burrough

    Calculates the local topographic slope (i.e., the down-dip slope, and
    presented as positive), and the aspect (dip direction in radians
    clockwise from north), at the given nodes, *ids*. All *ids* must be of
    core nodes.

    This method uses the Horn 1981 algorithm, the one employed by many GIS
    packages. It should be significantly faster than alternative slope
    methods.

    If *ids* is not provided, the slope will be returned for all core
    nodes.

    *vals* is either the name of an existing grid field from which to draw
    topographic data, or an array of values to use. If an array of values
    is passed, it must be nnodes long. If *vals* is not provided, this
    method will default to trying to use the field
    "topographic__elevation".

    Parameters
    ----------
    grid : RasterModelGrid
        Grid on which to calculate slopes and aspects.
    ids : array_like of int, optional
        Nodes on which to calculate slope and aspect.
    vals : str or ndarray, optional
        Node values used to measure slope and aspect.

    Returns
    -------
    (slope, aspect) : tuple of float
        *slope*: a len(ids) array of slopes at each node provided.
        *aspect*: a len(ids) array of aspects at each node provided.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> grid = RasterModelGrid((4, 5))

    Create a south-facing slope.

    >>> elevation = np.array([
    ...     0., 0., 0., 0., 0,
    ...     1., 1., 1., 1., 1,
    ...     2., 2., 2., 2., 2,
    ...     3., 3., 3., 3., 3])
    >>> (slope, aspect) = calc_slope_aspect_of_nodes_horn(grid,
    ...     vals=elevation)
    >>> len(slope) == grid.number_of_core_nodes
    True
    >>> len(aspect) == grid.number_of_core_nodes
    True
    >>> slope
    array([ 1.,  1.,  1.,  1.,  1.,  1.])
    >>> aspect * 180. / np.pi
    array([ 180.,  180.,  180.,  180.,  180.,  180.])

    Make the slope north-facing by multiplying elevations by -1. The slopes
    will still be positive but the aspects will change.

    >>> elevation *= -1
    >>> (slope, aspect) = calc_slope_aspect_of_nodes_horn(grid,
    ...     vals=elevation)
    >>> slope
    array([ 1.,  1.,  1.,  1.,  1.,  1.])
    >>> aspect * 180. / np.pi
    array([ 0.,  0.,  0.,  0.,  0.,  0.])

    Double the slope and make it west-facing.

    >>> elevation = np.array([
    ...     0., 1., 2., 3., 4.,
    ...     0., 1., 2., 3., 4.,
    ...     0., 1., 2., 3., 4.,
    ...     0., 1., 2., 3., 4.])
    >>> elevation *= 2.
    >>> (slope, aspect) = calc_slope_aspect_of_nodes_horn(grid,
    ...     vals=elevation)
    >>> slope
    array([ 2.,  2.,  2.,  2.,  2.,  2.])
    >>> aspect * 180. / np.pi
    array([ 270.,  270.,  270.,  270.,  270.,  270.])

    Make the slope east-facing by multiplying elevations by -1. The slopes
    will still be positive but the aspects will change.

    >>> elevation *= -1.
    >>> (slope, aspect) = calc_slope_aspect_of_nodes_horn(grid,
    ...     vals=elevation)
    >>> slope
    array([ 2.,  2.,  2.,  2.,  2.,  2.])
    >>> aspect * 180. / np.pi
    array([ 90.,  90.,  90.,  90.,  90.,  90.])

    LLCATS: DEPR NINF SURF
    """
    if ids is None:
        ids = grid.core_nodes
    if isinstance(vals, str):
        vals = grid.at_node[vals]
    else:
        if len(vals) != grid.number_of_nodes:
            raise IndexError("*vals* was not of a compatible length!")

    # [right, top, left, bottom]
    neighbors = grid.active_adjacent_nodes_at_node[ids]
    # [topright, topleft, bottomleft, bottomright]
    diagonals = grid.diagonal_adjacent_nodes_at_node[ids]

    input_array = np.empty((len(ids), 9), dtype=int)
    input_array[:, 0] = ids
    input_array[:, 1:5] = neighbors
    input_array[:, 5:] = diagonals

    slopes_array = np.apply_along_axis(_one_line_slopes, 1, input_array, grid, vals)
    slope_we = slopes_array[:, 0]
    slope_sn = slopes_array[:, 1]

    slope = np.sqrt(slope_we * slope_we + slope_sn * slope_sn)
    # aspect = np.empty_like(slope)
    aspect = np.ones(slope.size, dtype=float)
    simple_cases = slope_we != 0.0
    complex_cases = np.logical_not(simple_cases)

    complex_aspects = aspect[complex_cases]
    complex_aspects[slope_sn[complex_cases] < 0.0] = np.pi
    complex_aspects[slope_sn[complex_cases] >= 0.0] = 0.0
    aspect[complex_cases] = complex_aspects

    # +ve is CCW rotation from x axis
    angle_to_xaxis = np.arctan(slope_sn[simple_cases] / slope_we[simple_cases])
    aspect[simple_cases] = ((1.0 - np.sign(slope_we[simple_cases])) * 0.5) * np.pi + (
        0.5 * np.pi - angle_to_xaxis
    )

    return slope.ravel(), aspect.ravel()
