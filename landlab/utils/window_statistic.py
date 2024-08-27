"""Function to calculate node statistics in a moving window."""

import numpy as np

from landlab import FieldError


def calculate_window_statistic(
    grid, field, func, search_radius, calc_on_closed_nodes=True, **kwargs
):
    """Calculate a statistic using a function within a search window.

    .. note::

        This only works on grid **nodes** (not other grid elements e.g. links) for
        any :class:`~.ModelGrid` type.

    This utility outputs an array of length equal to the grid's number of
    nodes. Each element of the output array represents the node location in
    the grid. The value of each element is a function of the nodes within the
    search window surrounding that node location (see the model grid diagram
    below).

    The grid below contains six columns and five rows with cell spacing set
    to 10 distance units. This utility iteratively evaluates all nodes in the
    grid. The diagram shows evaluation of node ID 15 (marked ``x``). If the
    search radius is set to 20, twice the cell spacing, each node marked with
    a ``*`` is within the search window.
    ::

        · · * · · ·
        · * * * · ·
        * * x * * ·
        · * * * · ·
        · · * · · ·

    Increasing the search radius to 25 results in the following search window.
    ::

        · * * * · ·
        * * * * * ·
        * * x * * ·
        * * * * * ·
        · * * * · ·

    Decreasing the search radius to 15 results in the following search window.
    ::

        · · · · · ·
        · * * * · ·
        · * x * · ·
        · * * * · ·
        · · · · · ·

    The input field can be any field assigned to grid nodes (e.g.
    "topographic__elevation") and the input function can be any function that
    acts on the input field (e.g. "np.min" to find the minimum). The input
    function may be user defined and may contain any number of inputs, which
    are input as ``kwargs``.

    For example, if the input field is "topographic__elevation" and the input
    function is ``np.ptp`` (peak-to-peak, meaning max minus min value), then the
    output at node 15 will be the maximum elevation within the search window
    minus the minimum elevation within the search window (also known as relief).
    The ``np.percentile`` function, however, requires not only the input field,
    but also an input value to define the "q-th percentile" to be calculated.
    This second input would be added as a ``kwarg`` (e.g. ``q=90``) at the end of
    the inputs for :func:`~calculate_window_statistic`. Both of these scenarios are
    shown in the examples below.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab ModelGrid.
    field : string
        An existing grid field on which to calculate the statistic of interest.
        Must exist in grid.
    func : function
        The function that calculates the window statistic of *field*.
        The first parameter of the function must be the values at nodes within
        the window, which are used used to calculate the statistic for the
        node under evaluation. Additional parameters of the function can be
        passed with ``kwargs``.
    search_radius : float
        Radius of window within which the statistic is calculated.
    calc_on_closed_nodes : boolean, optional
        Toggle calculation over all nodes including closed nodes (``True``) or all
        nodes except closed nodes (``False``).
    kwargs : optional
        Keyword arguments passed to *func* that are additional to the array of
        node values within the search window.

    Returns
    -------
    output : ndarray
        Output array containing the calculated values of the statistic. Same
        length as input field.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.utils import window_statistic

    >>> grid = RasterModelGrid((5, 6), xy_spacing=10.0)
    >>> grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    >>> z = grid.add_zeros("topographic__elevation", at="node")
    >>> z += np.arange(len(z))

    Calculate relief using ``np.ptp`` function.

    >>> relief = calculate_window_statistic(
    ...     grid, "topographic__elevation", np.ptp, search_radius=15
    ... )
    >>> grid.at_node["topographic__elevation"]
    array([ 0.,   1.,   2.,   3.,   4.,   5.,
            6.,   7.,   8.,   9.,  10.,  11.,
           12.,  13.,  14.,  15.,  16.,  17.,
           18.,  19.,  20.,  21.,  22.,  23.,
           24.,  25.,  26.,  27.,  28.,  29.])
    >>> relief
    array([ 7.,   8.,   8.,   8.,   8.,   7.,
           13.,  14.,  14.,  14.,  14.,  13.,
           13.,  14.,  14.,  14.,  14.,  13.,
           13.,  14.,  14.,  14.,  14.,  13.,
            7.,   8.,   8.,   8.,   8.,   7.])

    Calculate relief using ``np.ptp`` function excluding closed nodes.

    >>> relief = calculate_window_statistic(
    ...     grid,
    ...     "topographic__elevation",
    ...     np.ptp,
    ...     search_radius=15,
    ...     calc_on_closed_nodes=False,
    ... )
    >>> grid.at_node["topographic__elevation"]
    array([ 0.,   1.,   2.,   3.,   4.,   5.,
            6.,   7.,   8.,   9.,  10.,  11.,
           12.,  13.,  14.,  15.,  16.,  17.,
           18.,  19.,  20.,  21.,  22.,  23.,
           24.,  25.,  26.,  27.,  28.,  29.])
    >>> relief
    array([nan,  nan,  nan,  nan,  nan,  nan,
            7.,   8.,   8.,   8.,   8.,   7.,
           13.,  14.,  14.,  14.,  14.,  13.,
            7.,   8.,   8.,   8.,   8.,   7.,
           nan,  nan,  nan,  nan,  nan,  nan])

    Calculate 90th percentile using ``np.percentile`` function and ``kwargs``.

    >>> perc_90 = calculate_window_statistic(
    ...     grid,
    ...     "topographic__elevation",
    ...     np.percentile,
    ...     search_radius=15,
    ...     calc_on_closed_nodes=False,
    ...     q=90,
    ... )
    >>> grid.at_node["topographic__elevation"]
    array([ 0.,   1.,   2.,   3.,   4.,   5.,
            6.,   7.,   8.,   9.,  10.,  11.,
           12.,  13.,  14.,  15.,  16.,  17.,
           18.,  19.,  20.,  21.,  22.,  23.,
           24.,  25.,  26.,  27.,  28.,  29.])
    >>> perc_90
    array([ nan,  nan,  nan,  nan,  nan,  nan,
           12.7, 13.5, 14.5, 15.5, 16.5, 16.7,
           18.5, 19.2, 20.2, 21.2, 22.2, 22.5,
           18.7, 19.5, 20.5, 21.5, 22.5, 22.7,
            nan,  nan,  nan,  nan,  nan,  nan])

    Calculate relief above 90th percentile elevation using a user-defined
    function and ``kwargs``.

    >>> def max_minus_percentile(elev, q):
    ...     output = np.max(elev) - np.percentile(elev, q)
    ...     return output
    ...
    >>> rel_above_90th_perc = calculate_window_statistic(
    ...     grid,
    ...     "topographic__elevation",
    ...     max_minus_percentile,
    ...     search_radius=15,
    ...     calc_on_closed_nodes=False,
    ...     q=90,
    ... )
    >>> grid.at_node["topographic__elevation"]
    array([ 0.,   1.,   2.,   3.,   4.,   5.,
            6.,   7.,   8.,   9.,  10.,  11.,
           12.,  13.,  14.,  15.,  16.,  17.,
           18.,  19.,  20.,  21.,  22.,  23.,
           24.,  25.,  26.,  27.,  28.,  29.])
    >>> rel_above_90th_perc
    array([nan,  nan,  nan,  nan,  nan,  nan,
           0.3,  0.5,  0.5,  0.5,  0.5,  0.3,
           0.5,  0.8,  0.8,  0.8,  0.8,  0.5,
           0.3,  0.5,  0.5,  0.5,  0.5,  0.3,
           nan,  nan,  nan,  nan,  nan,  nan])
    """
    if field not in grid.at_node:
        raise FieldError(f"A {field} field is required at the nodes of the input grid.")

    # Create output array
    output = np.zeros(grid.number_of_nodes)

    # Create arrays of x and y coords for input to "distance to point' calc
    x_coord = grid.x_of_node
    y_coord = grid.y_of_node

    nodes_in_loop = grid.nodes.flatten()
    nodes_to_include = np.ones(grid.number_of_nodes, dtype=bool)

    if calc_on_closed_nodes is False:
        closed_nodes = grid.status_at_node == grid.BC_NODE_IS_CLOSED
        nodes_in_loop = nodes_in_loop[~closed_nodes]
        nodes_to_include[closed_nodes] = False
        output[closed_nodes] = np.NaN

    # Calculate "dist to point" then local value at nodes within window.
    for node in nodes_in_loop:
        node_dist_to_point = grid.calc_distances_of_nodes_to_point(
            (x_coord[node], y_coord[node])
        )
        nodes_in_window = np.all(
            [node_dist_to_point <= search_radius, nodes_to_include], 0
        )
        values_in_window = grid.at_node[field][nodes_in_window]
        output[node] = func(values_in_window, **kwargs)

    return output
