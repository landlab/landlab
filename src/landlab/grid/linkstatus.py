#! /usr/bin/env python
from enum import IntEnum
from enum import unique

import numpy as np

from .nodestatus import NodeStatus


@unique
class LinkStatus(IntEnum):
    """Define the link types"""

    #: Indicate a link is *active*, and can carry flux
    ACTIVE = 0
    #: Indicate a link has a fixed (gradient) value, & behaves as a boundary
    FIXED = 2
    #: Indicate a link is *inactive*, and cannot carry flux
    INACTIVE = 4


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
    >>> from landlab import NodeStatus
    >>> is_fixed_link([NodeStatus.CORE, NodeStatus.FIXED_GRADIENT])
    array([ True])

    >>> is_fixed_link([NodeStatus.CORE, NodeStatus.FIXED_VALUE])
    array([False])

    >>> is_fixed_link(
    ...     [
    ...         [NodeStatus.FIXED_GRADIENT, NodeStatus.CORE],
    ...         [NodeStatus.CORE, NodeStatus.CORE],
    ...     ]
    ... )
    array([ True, False])
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core_node = node_status_at_link == NodeStatus.CORE
    is_fixed_gradient_node = node_status_at_link == NodeStatus.FIXED_GRADIENT

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
    >>> from landlab import NodeStatus
    >>> is_inactive_link([NodeStatus.CORE, NodeStatus.CLOSED])
    array([ True])

    >>> is_inactive_link([NodeStatus.FIXED_GRADIENT, NodeStatus.FIXED_VALUE])
    array([ True])

    >>> is_inactive_link(
    ...     [
    ...         [NodeStatus.FIXED_GRADIENT, NodeStatus.CLOSED],
    ...         [NodeStatus.CORE, NodeStatus.CORE],
    ...     ]
    ... )
    array([ True, False])
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core = node_status_at_link == NodeStatus.CORE
    is_fixed_value = node_status_at_link == NodeStatus.FIXED_VALUE
    is_fixed_gradient = node_status_at_link == NodeStatus.FIXED_GRADIENT
    is_closed = node_status_at_link == NodeStatus.CLOSED
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
    >>> from landlab import NodeStatus
    >>> is_active_link([NodeStatus.CORE, NodeStatus.FIXED_GRADIENT])
    array([False])

    >>> is_active_link([NodeStatus.CORE, NodeStatus.FIXED_VALUE])
    array([ True])

    >>> is_active_link(
    ...     [
    ...         [NodeStatus.FIXED_GRADIENT, NodeStatus.CORE],
    ...         [NodeStatus.CORE, NodeStatus.CORE],
    ...     ]
    ... )
    array([False, True])
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core_node = node_status_at_link == NodeStatus.CORE
    is_fixed_value_node = node_status_at_link == NodeStatus.FIXED_VALUE
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

    out[_is_inactive_link] = LinkStatus.INACTIVE
    out[_is_active_link] = LinkStatus.ACTIVE
    out[_is_fixed_link] = LinkStatus.FIXED

    return out
