def set_status_at_node_on_edges(grid, status_at_edge=None):
    """Set node status on grid edges.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    status_at_edge : dict_like
        Status of grid edges. This can be either a dict (or an OrderedDict)
        or a list of tuples.

    Examples
    --------
    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY

    >>> grid = RasterModelGrid((3, 4))
    >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([1, 1, 1, 1,
           1, 0, 0, 1,
           1, 1, 1, 1], dtype=int8)

    >>> grid.set_status_at_node_on_edges(
    ...     status_at_edge=dict(right=CLOSED_BOUNDARY))
    >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([1, 1, 1, 4,
           1, 0, 0, 4,
           1, 1, 1, 4], dtype=int8)

    >>> from landlab import FIXED_GRADIENT_BOUNDARY
    >>> grid = RasterModelGrid((3, 4))

    The statuses are set in the order they are provided in the *status_at_edge*
    keyword. This is useful for status values on corners. If two edges
    share a corner, it is the value of the last edge set that gives the value
    of the corner node.

    >>> status_at_edge = [('left', CLOSED_BOUNDARY),
    ...                   ('bottom', FIXED_GRADIENT_BOUNDARY)]
    >>> grid.set_status_at_node_on_edges(status_at_edge=status_at_edge)
    >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([2, 2, 2, 2,
           4, 0, 0, 1,
           4, 1, 1, 1], dtype=int8)

    In this case we set the left edge last, so the bottom left corner is
    assigned the value for the left edge.

    >>> status_at_edge = [('bottom', FIXED_GRADIENT_BOUNDARY),
    ...                   ('left', CLOSED_BOUNDARY)]
    >>> grid.set_status_at_node_on_edges(status_at_edge=status_at_edge)
    >>> grid.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([4, 2, 2, 2,
           4, 0, 0, 1,
           4, 1, 1, 1], dtype=int8)
    """
    status_at_edge = status_at_edge or []
    try:
        status_at_edge = status_at_edge.items()
    except AttributeError:
        pass

    for edge, val in status_at_edge:
        nodes = grid.nodes_at_edge(edge)
        grid.status_at_node[nodes] = val
