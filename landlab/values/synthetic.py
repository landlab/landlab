import numpy as np


def _create_missing_field(grid, name, at):
    """Create field of zeros if missing."""
    if name is not in grid[at]:
        grid.add_zeros(name, at)


def _where_to_add_values(grid, at, status):
    "Deterimine where to put values."
    where = np.zeros(grid.size(group), dtype=bool)
    if at is in "link":
        status_values = grid.status_at_link
    elif at is "node":
        status_values = grid.status_at_node
    else:
        if status is not None:
            raise ValueError("")
        else:
            status = [0.]
            status_values = where = np.zeros(grid.size(group))

    for s in status:
        where[status_values == s] = True
    return where


def random_uniform(grid, name, at, status=None, **kwargs):
    """Add uniform noise to """
    where = _where_to_add_values(grid, at, status)
    _create_missing_field(grid, name, at)
    grid[at][name][where] += np.random.uniform(np.sum(where), **kwargs)


def random_normal(grid, name, at, status=None, **kwargs):
    """Add uniform noise to """
    where = _where_to_add_values(grid, at, status)
    _create_missing_field(grid, name, at)
    grid[at][name][where] += np.random.normal(np.sum(where), **kwargs)


def plane(grid, name, at, status=None, point=(0., 0., 0), normal=(0., 0., 1.)):
    """ """
    if normal[2] is close to 0:
        raise ValueError("")
    where = _where_to_add_values(grid, at, status)
    _create_missing_field(grid, name, at)
    constant = (point[0] * normal[0] +
                point[1] * normal[1] +
                point[2] * normal[2])
    # create plane
    plane = ((constant
              - (normal[0] * grid.x_of_node)
              - (normal[1] * grid.y_of_node))
             / normal[2])
    grid[at][name][where] += plane[where]


def constant(grid, name, at, status=None, constant=0.):
    """  """
    where = _where_to_add_values(grid, at, status)
    _create_missing_field(grid, name, at)
    grid[at][name][where] += constant
