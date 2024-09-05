def set_status_at_node_on_edges(grid, right=None, top=None, left=None, bottom=None):
    """Set node status on grid edges.

    Parameters
    ----------
    grid : RasterModelGrid
        A grid.
    right : int, optional
        Node status along right edge.
    top : int, optional
        Node status along top edge.
    left : int, optional
        Node status along left edge.
    bottom : int, optional
        Node status along bottom edge.

    Examples
    --------
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 4))
    >>> grid.status_at_node.reshape(grid.shape)
    array([[1, 1, 1, 1],
           [1, 0, 0, 1],
           [1, 1, 1, 1]], dtype=uint8)

    >>> grid.set_status_at_node_on_edges(right=grid.BC_NODE_IS_CLOSED)
    >>> grid.status_at_node.reshape(grid.shape)
    array([[1, 1, 1, 4],
           [1, 0, 0, 4],
           [1, 1, 1, 4]], dtype=uint8)

    >>> grid = RasterModelGrid((3, 4))

    The status of a corner is set along with its clockwise edge. That is,
    if setting the status for the top and right edges, the upper-right corner
    has the status of the right edge.

    >>> grid.set_status_at_node_on_edges(
    ...     top=grid.BC_NODE_IS_CLOSED, right=grid.BC_NODE_IS_FIXED_GRADIENT
    ... )
    >>> grid.status_at_node.reshape(grid.shape)
    array([[1, 1, 1, 2],
           [1, 0, 0, 2],
           [4, 4, 4, 2]], dtype=uint8)

    In the above example, if you wanted the corner to have the status of the
    top edge, you need to make two calls to `set_status_at_node_on_edges`,

    >>> grid = RasterModelGrid((3, 4))
    >>> grid.set_status_at_node_on_edges(right=grid.BC_NODE_IS_FIXED_GRADIENT)
    >>> grid.set_status_at_node_on_edges(top=grid.BC_NODE_IS_CLOSED)
    >>> grid.status_at_node.reshape(grid.shape)
    array([[1, 1, 1, 2],
           [1, 0, 0, 2],
           [4, 4, 4, 4]], dtype=uint8)

    An example that sets all of the edges shows how corners are set.

    >>> grid.set_status_at_node_on_edges(right=1, top=2, left=3, bottom=4)
    >>> grid.status_at_node.reshape(grid.shape)
    array([[3, 4, 4, 4],
           [3, 0, 0, 1],
           [2, 2, 2, 1]], dtype=uint8)

    :meta landlab: boundary-condition
    """
    status_at_edge = (
        ("bottom", bottom),
        ("left", left),
        ("top", top),
        ("right", right),
    )

    for edge, val in status_at_edge:
        if val is not None:
            nodes = grid.nodes_at_edge(edge)
            grid.status_at_node[nodes] = val

    if right is not None and bottom is not None:
        lr = grid.nodes_at_right_edge[0]
        grid.status_at_node[lr] = bottom
