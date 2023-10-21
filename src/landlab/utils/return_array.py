#! /usr/bin/env python
"""Return array with same shape as grid elements."""
from landlab.utils.decorators import use_field_name_array_or_value


@use_field_name_array_or_value("node")
def return_array_at_node(grid, value):
    """Function to return an array stored at node or of shape `(n_nodes,)`.

    This function exists to take advantage of the use_field_name_array_or_value
    decorator which permits providing the surface as a field name or array.

    Parameters
    ----------
    grid : ModelGrid
    value : field name, ndarray of shape `(n_nodes, )`, or single value.

    Returns
    -------
    array : ndarray of shape `(n_nodes, )`
    """
    return value


@use_field_name_array_or_value("link")
def return_array_at_link(grid, value):
    """Function to return an array stored at node or of shape `(n_nodes,)`.

    This function exists to take advantage of the use_field_name_array_or_value
    decorator which permits providing the surface as a field name or array.

    Parameters
    ----------
    grid : ModelGrid
    value : field name, ndarray of shape `(n_nodes, )`, or single value.

    Returns
    -------
    array : ndarray of shape `(n_nodes, )`
    """
    return value
