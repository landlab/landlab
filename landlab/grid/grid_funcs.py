"""
Utility functions that operate on landlab grids.
------------------------------------------------

"""
import numpy as np


def resolve_values_on_active_links(grid, active_link_values):
    """Resolve active-link values into x and y directions.

    Takes a set of values defined on active links, and returns those values
    resolved into the x and y directions.  Two link arrays are returned:
    x, then y.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    active_link_values : ndarray
        Values on active links.

    Returns
    -------
    tuple of ndarray
        Values resolved into x-component and y-component.
    """
    link_lengths = grid.length_of_link[grid.active_links]
    return (
        np.multiply(
            (
                (
                    grid.node_x[grid._activelink_tonode]
                    - grid.node_x[grid._activelink_fromnode]
                )
                / link_lengths
            ),
            active_link_values,
        ),
        np.multiply(
            (
                (
                    grid.node_y[grid._activelink_tonode]
                    - grid.node_y[grid._activelink_fromnode]
                )
                / link_lengths
            ),
            active_link_values,
        ),
    )


def resolve_values_on_links(grid, link_values):
    """Resolve link values into x and y directions.

    Takes a set of values defined on active links, and returns those values
    resolved into the x and y directions.  Two link arrays are returned:
    x, then y.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    link_values : ndarray
        Values on links.

    Returns
    -------
    tuple of ndarray
        Values resolved into x-component and y-component.
    """
    return (
        np.multiply(
            (
                (
                    grid.node_x[grid.node_at_link_head]
                    - grid.node_x[grid.node_at_link_tail]
                )
                / grid.length_of_link
            ),
            link_values,
        ),
        np.multiply(
            (
                (
                    grid.node_y[grid.node_at_link_head]
                    - grid.node_y[grid.node_at_link_tail]
                )
                / grid.length_of_link
            ),
            link_values,
        ),
    )
