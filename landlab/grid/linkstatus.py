#! /usr/bin/env python
import numpy as np

from .nodestatus import (
    CLOSED_BOUNDARY,
    CORE_NODE,
    FIXED_GRADIENT_BOUNDARY,
    FIXED_VALUE_BOUNDARY,
)

# Define the link types

#: Indicates a link is *active*, and can carry flux
ACTIVE_LINK = 0

#: Indicates a link has a fixed (gradient) value, & behaves as a boundary
FIXED_LINK = 2

#: Indicates a link is *inactive*, and cannot carry flux
INACTIVE_LINK = 4

LINK_STATUS_FLAGS_LIST = [ACTIVE_LINK, FIXED_LINK, INACTIVE_LINK]
LINK_STATUS_FLAGS = set(LINK_STATUS_FLAGS_LIST)


def is_fixed_link(node_status_at_link):
    """Find links that are fixed.

    A link is fixed if it connects a core node with a fixed value
    boundary node.

    Parameters
    ----------
    node_status_at_link : ndarray of int, shape `(n_links, 2)`
        Node status a link tail and head.

    Returns
    -------
    ndarray of bool, shape `(n_links, )`
        True if link is fixed.

    Examples
    --------
    >>> from landlab.grid.linkstatus import is_fixed_link
    >>> from landlab import CORE_NODE, FIXED_GRADIENT_BOUNDARY
    >>> is_fixed_link([CORE_NODE, FIXED_GRADIENT_BOUNDARY])
    array([ True], dtype=bool)

    >>> from landlab import FIXED_VALUE_BOUNDARY
    >>> is_fixed_link([CORE_NODE, FIXED_VALUE_BOUNDARY])
    array([False], dtype=bool)

    >>> is_fixed_link([[FIXED_GRADIENT_BOUNDARY, CORE_NODE],
    ...                [CORE_NODE, CORE_NODE]])
    array([ True, False], dtype=bool)
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core_node = node_status_at_link == CORE_NODE
    is_fixed_gradient_node = node_status_at_link == FIXED_GRADIENT_BOUNDARY

    return (is_core_node[:, 0] & is_fixed_gradient_node[:, 1]) | (
        is_fixed_gradient_node[:, 0] & is_core_node[:, 1]
    )


def is_inactive_link(node_status_at_link):
    """Find links that are inactive.

    A link is inactive if it connects two boundary nodes or one of
    its nodes is closed.

    Parameters
    ----------
    node_status_at_link : ndarray of int, shape `(n_links, 2)`
        Node status a link tail and head.

    Returns
    -------
    ndarray of bool, shape `(n_links, )`
        True if link is isactive.

    Examples
    --------
    >>> from landlab.grid.linkstatus import is_inactive_link
    >>> from landlab import CORE_NODE, FIXED_GRADIENT_BOUNDARY
    >>> is_inactive_link([CORE_NODE, CLOSED_BOUNDARY])
    array([ True], dtype=bool)

    >>> from landlab import FIXED_VALUE_BOUNDARY
    >>> is_inactive_link([FIXED_GRADIENT_BOUNDARY, FIXED_VALUE_BOUNDARY])
    array([ True], dtype=bool)

    >>> is_inactive_link([[FIXED_GRADIENT_BOUNDARY, CLOSED_BOUNDARY],
    ...                   [CORE_NODE, CORE_NODE]])
    array([ True, False], dtype=bool)
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core = node_status_at_link == CORE_NODE
    is_fixed_value = node_status_at_link == FIXED_VALUE_BOUNDARY
    is_fixed_gradient = node_status_at_link == FIXED_GRADIENT_BOUNDARY
    is_closed = node_status_at_link == CLOSED_BOUNDARY
    is_boundary_node = is_fixed_value | is_fixed_gradient | is_closed

    return (
        (is_boundary_node[:, 0] & is_boundary_node[:, 1])
        | (is_closed[:, 0] & is_core[:, 1])
        | (is_core[:, 0] & is_closed[:, 1])
    )


def is_active_link(node_status_at_link):
    """Find links that are active.

    A link is active if it connects a core node with another core
    node or a fixed value boundary.

    Parameters
    ----------
    node_status_at_link : ndarray of int, shape `(n_links, 2)`
        Node status a link tail and head.

    Returns
    -------
    ndarray of bool, shape `(n_links, )`
        True if link is isactive.

    Examples
    --------
    >>> from landlab.grid.linkstatus import is_active_link
    >>> from landlab import CORE_NODE, FIXED_GRADIENT_BOUNDARY
    >>> is_active_link([CORE_NODE, FIXED_GRADIENT_BOUNDARY])
    array([False], dtype=bool)

    >>> from landlab import FIXED_VALUE_BOUNDARY
    >>> is_active_link([CORE_NODE, FIXED_VALUE_BOUNDARY])
    array([ True], dtype=bool)

    >>> is_active_link([[FIXED_GRADIENT_BOUNDARY, CORE_NODE],
    ...                 [CORE_NODE, CORE_NODE]])
    array([False, True], dtype=bool)
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core_node = node_status_at_link == CORE_NODE
    is_fixed_value_node = node_status_at_link == FIXED_VALUE_BOUNDARY
    return (
        (is_core_node[:, 0] & is_core_node[:, 1])
        | (is_core_node[:, 0] & is_fixed_value_node[:, 1])
        | (is_fixed_value_node[:, 0] & is_core_node[:, 1])
    )


def set_status_at_link(node_status_at_link, out=None):
    n_links = len(node_status_at_link)

    if out is None:
        out = np.full(n_links, 255, dtype=np.uint8)

    _is_fixed_link = is_fixed_link(node_status_at_link)
    _is_active_link = is_active_link(node_status_at_link)
    _is_inactive_link = is_inactive_link(node_status_at_link)

    assert np.all(
        np.sum(np.vstack((_is_active_link, _is_inactive_link, _is_fixed_link)), axis=0)
        == 1
    )

    out[_is_inactive_link] = INACTIVE_LINK
    out[_is_active_link] = ACTIVE_LINK
    out[_is_fixed_link] = FIXED_LINK

    return out
