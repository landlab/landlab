import numpy as np

from ...core.utils import as_id_array
from ...graph.structured_quad.structured_quad import StructuredQuadGraphTopology
from . import _neighbors_at_link


def neighbors_at_link(shape, links):
    """Get neighbor links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.overland_flow._links import neighbors_at_link

    >>> neighbors_at_link((3, 2), np.arange(7))
    array([[-1,  3, -1, -1],
           [ 2,  4, -1, -1], [-1,  5,  1, -1],
           [-1,  6, -1,  0],
           [ 5,  7, -1,  1], [-1, -1,  4,  2],
           [-1, -1, -1,  3]])
    """
    links = np.asarray(links, dtype=int)
    out = np.full((links.size, 4), -1, dtype=int)
    _neighbors_at_link.neighbors_at_link(links, shape, out)
    return out


def vertical_link_ids(shape):
    """Vertical links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    (M, N) ndarray :
        Array of link IDs.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import vertical_link_ids
    >>> vertical_link_ids((3, 4))
    array([[ 3,  4,  5,  6],
           [10, 11, 12, 13]])
    """
    layout = StructuredQuadGraphTopology(shape)
    return layout.vertical_links.reshape((shape[0] - 1, shape[1]))


def horizontal_link_ids(shape):
    """Horizontal links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    (M, N) ndarray :
        Array of link IDs.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import horizontal_link_ids
    >>> horizontal_link_ids((3, 4))
    array([[ 0,  1,  2],
           [ 7,  8,  9],
           [14, 15, 16]])
    """
    layout = StructuredQuadGraphTopology(shape)
    return layout.horizontal_links.reshape((shape[0], shape[1] - 1))


def vertical_south_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of south, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *south* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     vertical_link_ids,
    ...     vertical_south_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> vertical_links = vertical_link_ids(rmg.shape)
    >>> vertical_south_link_neighbor(rmg.shape, vertical_links).reshape((3, 5))
    array([[-1, -1, -1, -1, -1],
           [ 4,  5,  6,  7,  8],
           [13, 14, 15, 16, 17]])
    """
    vertical_links = StructuredQuadGraphTopology(shape).vertical_links
    vertical_links[shape[1] :] = vertical_links[: -shape[1]]
    vertical_links[: shape[1]] = bad_index_value

    return vertical_links


def vertical_west_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of west, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids- MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *west* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     vertical_link_ids,
    ...     vertical_west_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> vertical_links = vertical_link_ids(rmg.shape)
    >>> vertical_west_link_neighbor(rmg.shape, vertical_links).reshape((3, 5))
    array([[-1,  4,  5,  6,  7],
           [-1, 13, 14, 15, 16],
           [-1, 22, 23, 24, 25]])
    """
    vertical_links = StructuredQuadGraphTopology(shape).vertical_links.reshape(
        (shape[0] - 1, shape[1])
    )
    vertical_links[:, 1:] = vertical_links[:, :-1]
    vertical_links[:, 0] = bad_index_value

    return vertical_links.reshape(-1)


def vertical_north_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of north, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids- MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *north* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     vertical_link_ids,
    ...     vertical_north_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> vertical_ids = vertical_link_ids(rmg.shape)
    >>> vertical_north_link_neighbor(rmg.shape, vertical_ids).reshape((3, 5))
    array([[13, 14, 15, 16, 17],
           [22, 23, 24, 25, 26],
           [-1, -1, -1, -1, -1]])
    """
    vertical_links = StructuredQuadGraphTopology(shape).vertical_links
    vertical_links[: -shape[1]] = vertical_links[shape[1] :]
    vertical_links[-shape[1] :] = bad_index_value

    return vertical_links


def vertical_east_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of east, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *east* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     vertical_link_ids,
    ...     vertical_east_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> vertical_links = vertical_link_ids(rmg.shape)
    >>> vertical_east_link_neighbor(rmg.shape, vertical_links).reshape((3, 5))
    array([[ 5,  6,  7,  8, -1],
           [14, 15, 16, 17, -1],
           [23, 24, 25, 26, -1]])
    """
    vertical_links = StructuredQuadGraphTopology(shape).vertical_links.reshape(
        (shape[0] - 1, shape[1])
    )
    vertical_links[:, :-1] = vertical_links[:, 1:]
    vertical_links[:, -1] = bad_index_value

    return vertical_links.base


def active_link_ids(shape, node_status):
    """Get active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    node_status : array_link
        Status of nodes in grid.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------
    >>> from landlab.grid import RasterModelGrid
    >>> from landlab.components.overland_flow._links import active_link_ids

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    >>> status = rmg.status_at_node
    >>> status.reshape(rmg.shape)
    array([[4, 4, 4, 4],
           [4, 0, 0, 4],
           [4, 4, 4, 4]], dtype=uint8)

    >>> active_link_ids((3, 4), status)
    array([8])
    """
    return as_id_array(np.where(is_active_link(shape, node_status))[0])


def is_active_link(shape, node_status):
    """Link IDs of active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    node_status : array_link
        Status of nodes in grid.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import is_active_link
    >>> from landlab.grid.nodestatus import NodeStatus

    >>> status = [
    ...     [NodeStatus.CLOSED, NodeStatus.CLOSED, NodeStatus.CLOSED],
    ...     [NodeStatus.CLOSED, NodeStatus.CORE, NodeStatus.CLOSED],
    ...     [NodeStatus.CLOSED, NodeStatus.CORE, NodeStatus.CLOSED],
    ...     [NodeStatus.CLOSED, NodeStatus.CLOSED, NodeStatus.CLOSED],
    ... ]
    >>> is_active_link((4, 3), status)
    array([False, False,
           False, False, False,
           False, False,
           False, True, False,
           False, False,
           False, False, False,
           False, False])
    """
    from ...grid.linkstatus import is_active_link

    node_status = np.asarray(node_status).reshape(-1)

    if np.prod(shape) != node_status.size:
        raise ValueError(
            "node status array does not match size of grid "
            "(%d != %d)" % (np.prod(shape), len(node_status))
        )

    # status_at_link_start = node_status.flat[node_id_at_link_start(shape)]
    # status_at_link_end = node_status.flat[node_id_at_link_end(shape)]

    # status_at_link = node_status[StructuredQuadGraphTopology(shape).nodes_at_link]

    return is_active_link(node_status[StructuredQuadGraphTopology(shape).nodes_at_link])


def vertical_active_link_ids(shape, active_ids, bad_index_value=-1):
    """ID of vertical active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    active_ids : array of int
        Array of all active link ids
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the VERTICAL active links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

          *---I-->*---I-->*---I-->*---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->o---H-->o---H-->o---I-->*
          ^       ^       ^       ^       ^
          I       6       7       8       I
          |       |       |       |       |
          *---I-->o---H-->o---H-->o---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->*---I-->*---I-->*---I-->*

    .. note::

        ``*`` indicates the nodes that are set to `NodeStatus.CLOSED`

        ``o`` indicates the nodes that are set to `NodeStatus.CORE`

        ``I`` indicates the links that are set to `LinkStatus.INACTIVE`

        ``H`` indicates horizontal active ids, which are ignored by this
        function

        Numeric values correspond to the vertical `LinkStatus.ACTIVE` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     active_link_ids,
    ...     vertical_active_link_ids,
    ... )

    >>> rmg = RasterModelGrid((4, 5))
    >>> active_ids = active_link_ids((4, 5), rmg.status_at_node)
    >>> active_ids
    array([ 5,  6,  7,
            9, 10, 11, 12,
           14, 15, 16,
           18, 19, 20, 21,
           23, 24, 25])

    >>> vertical_active_link_ids((4, 5), active_ids).reshape((3, 5))
    array([[-1,  5,  6,  7, -1],
           [-1, 14, 15, 16, -1],
           [-1, 23, 24, 25, -1]])

    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.status_at_node
    >>> active_ids = active_link_ids((4, 5), status)
    >>> vertical_active_link_ids((4, 5), active_ids).reshape((3, 5))
    array([[-1, -1, -1, -1, -1],
           [-1, 14, 15, 16, -1],
           [-1, -1, -1, -1, -1]])
    """
    number_of_vertical_links = (shape[0] - 1) * shape[1]
    out = np.full(number_of_vertical_links, bad_index_value, dtype=int)
    vertical_ids = active_ids[np.where(is_vertical_link(shape, active_ids))]

    out[nth_vertical_link(shape, vertical_ids)] = vertical_ids
    return out


def _number_of_links(shape):
    """Number of links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of links in grid.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import _number_of_links
    >>> _number_of_links((3, 4))
    17
    """
    return (shape[0] - 1) * shape[1] + shape[0] * (shape[1] - 1)
    # return number_of_vertical_links(shape) + number_of_horizontal_links(shape)


def number_of_vertical_links(shape):
    """Number of vertical links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of vertical links in grid.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import number_of_vertical_links
    >>> number_of_vertical_links((3, 4))
    8
    """
    return (shape[0] - 1) * shape[1]


def number_of_horizontal_links(shape):
    """Number of horizontal links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of horizontal links in grid.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import number_of_horizontal_links
    >>> number_of_horizontal_links((3, 4))
    9
    """
    return shape[0] * (shape[1] - 1)


def is_vertical_link(shape, links):
    """Test if links are vertical.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of bool
        `True` for links that are vertical.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import (
    ...     is_vertical_link,
    ...     _number_of_links,
    ... )
    >>> import numpy as np
    >>> shape = (3, 4)
    >>> links = np.arange(_number_of_links(shape))
    >>> is_vertical_link(shape, links)
    array([False, False, False,  True,  True,  True,  True,
           False, False, False,  True,  True,  True,  True,
           False, False, False])
    """
    return ((links % (2 * shape[1] - 1)) >= shape[1] - 1) & (
        links < _number_of_links(shape)
    )


def nth_vertical_link(shape, links):
    """Convert link ID to vertical link ID.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of int
        The link ID as the nth vertical links.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import nth_vertical_link
    >>> shape = (3, 4)
    >>> nth_vertical_link(shape, 4)
    1
    >>> nth_vertical_link(shape, (3, 4, 11))
    array([0, 1, 5])
    """
    links = np.asarray(links, dtype=int)
    return as_id_array(
        (links // (2 * shape[1] - 1)) * shape[1]
        + links % (2 * shape[1] - 1)
        - (shape[1] - 1)
    )


def horizontal_active_link_ids(shape, active_ids, bad_index_value=-1):
    """ID of horizontal active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    active_ids : array of int
        Array of all active link ids
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the HORIZONTAL active links. Length of
        number_of_horizontal_links.

    Examples
    --------

    The following example uses this grid::

          *---I-->*---I-->*---I-->*---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->o--24-->o--25-->o---I-->*
          ^       ^       ^       ^       ^
          I       V       V       V       I
          |       |       |       |       |
          *---I-->o--20-->o--21-->o---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->*---I-->*---I-->*---I-->*

    .. note::

        ``*`` indicates the nodes that are set to `NodeStatus.CLOSED`

        ``o`` indicates the nodes that are set to `NodeStatus.CORE`

        ``I`` indicates the links that are set to `LinkStatus.INACTIVE`

        ``V`` indicates vertical active ids, which are ignored by this
        function.

        Numeric values correspond to the horizontal `LinkStatus.ACTIVE`  ID.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     active_link_ids,
    ...     horizontal_active_link_ids,
    ... )

    >>> rmg = RasterModelGrid((4, 5))
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    >>> status = rmg.status_at_node
    >>> status.reshape(rmg.shape)
    array([[4, 4, 4, 4, 4],
           [4, 0, 0, 0, 4],
           [4, 0, 0, 0, 4],
           [4, 4, 4, 4, 4]], dtype=uint8)
    >>> active_ids = active_link_ids((4, 5), status)

    >>> horizontal_active_link_ids((4, 5), active_ids).reshape((4, 4))
    array([[-1, -1, -1, -1],
           [-1, 10, 11, -1],
           [-1, 19, 20, -1],
           [-1, -1, -1, -1]])
    """
    number_of_horizontal_links = shape[0] * (shape[1] - 1)
    out = np.full(number_of_horizontal_links, bad_index_value, dtype=int)
    horizontal_ids = active_ids[np.where(~is_vertical_link(shape, active_ids))]

    out[nth_horizontal_link(shape, horizontal_ids)] = horizontal_ids
    return out


def nth_horizontal_link(shape, links):
    """Convert link ID to horizontal link ID.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of int
        The link ID as the nth horizontal links.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import nth_horizontal_link
    >>> shape = (3, 4)
    >>> nth_horizontal_link(shape, 16)
    8
    >>> nth_horizontal_link(shape, (1, 7, 8))
    array([1, 3, 4])
    """
    links = np.asarray(links, dtype=int)
    return as_id_array(
        (links // (2 * shape[1] - 1)) * (shape[1] - 1) + links % (2 * shape[1] - 1)
    )


def is_horizontal_link(shape, links):
    """Test if a link is horizontal.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of bool
        `True` for links that are horizontal.

    Examples
    --------
    >>> from landlab.components.overland_flow._links import (
    ...     is_horizontal_link,
    ...     _number_of_links,
    ... )
    >>> import numpy as np
    >>> shape = (3, 4)
    >>> links = np.arange(_number_of_links(shape))
    >>> is_horizontal_link(shape, links)
    array([ True,  True,  True, False, False, False, False,
            True,  True,  True, False, False, False, False,
            True,  True,  True])
    """
    return (~is_vertical_link(shape, links)) & (links < _number_of_links(shape))


def horizontal_west_link_neighbor(shape, horizontal_ids, bad_index_value=-1):
    """ID of west, horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of west horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     horizontal_link_ids,
    ...     horizontal_west_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_west_link_neighbor(rmg.shape, horizontal_links)
    array([-1,  0,  1,  2, -1,  9, 10, 11, -1, 18, 19, 20, -1, 27, 28, 29])
    """
    links = np.roll(horizontal_ids.reshape((shape[0], shape[1] - 1)), 1, axis=1)
    links[:, 0] = bad_index_value

    return links.reshape(-1)


def horizontal_east_link_neighbor(shape, horizontal_ids, bad_index_value=-1):
    """IDs of east, horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of east horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal `LinkStatus.ACTIVE` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     horizontal_link_ids,
    ...     horizontal_east_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_east_link_neighbor(rmg.shape, horizontal_links)
    array([ 1,  2,  3, -1, 10, 11, 12, -1, 19, 20, 21, -1, 28, 29, 30, -1])
    """
    links = np.roll(horizontal_ids.reshape((shape[0], shape[1] - 1)), -1, axis=1)
    links[:, -1] = -1

    return links.reshape(-1)


def horizontal_north_link_neighbor(shape, horizontal_ids, bad_index_value=-1):
    """ID of north, horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of north horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal `LinkStatus.ACTIVE` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     horizontal_link_ids,
    ...     horizontal_north_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_north_link_neighbor(rmg.shape, horizontal_links)
    array([ 9, 10, 11, 12, 18, 19, 20, 21, 27, 28, 29, 30, -1, -1, -1, -1])
    """
    links = np.roll(horizontal_ids.reshape((shape[0], shape[1] - 1)), -1, axis=0)
    links[-1, :] = bad_index_value

    return links.reshape(-1)


def horizontal_south_link_neighbor(shape, horizontal_ids, bad_index_value=-1):
    """ID of south horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids *must be of len(horizontal_links)*.
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of south horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.overland_flow._links import (
    ...     horizontal_link_ids,
    ...     horizontal_north_link_neighbor,
    ... )
    >>> rmg = RasterModelGrid((4, 5))
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_south_link_neighbor(rmg.shape, horizontal_links)
    array([-1, -1, -1, -1,  0,  1,  2,  3,  9, 10, 11, 12, 18, 19, 20, 21])
    """
    links = np.roll(horizontal_ids.reshape((shape[0], shape[1] - 1)), 1, axis=0)
    links[0, :] = bad_index_value

    return links.reshape(-1)
