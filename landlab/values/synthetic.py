"""
synthetic.py provides functions to add synthetic values to a model grid.

Values can be added to any valid grid element (e.g. link or node). If no field
exists, a field of float zeros will be initialized.

All functions add values to the field---this means that multiple functions can
be chained together.

All functions support adding values to only portions of the grid, based on the
``status_at_link`` and ``status_at_node`` attributes.

For example, if one wanted to construct an initial topographic elevation
represented by a tetrahedron and add normally distributed noise only to core
nodes, this could be acomplished as follows:

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab import CORE_NODE
>>> from landlab.values import random, plane
>>> mg = RasterModelGrid((11, 11))

"""
import numpy as np


def _create_missing_field(grid, name, at):
    """Create field of zeros if missing."""
    if name not in grid[at]:
        grid.add_zeros(at, name)


def _where_to_add_values(grid, at, where):
    "Determine where to put values."
    where_to_place = np.zeros(grid.size(at), dtype=bool)

    if at is "link":
        status_values = grid.status_at_link
    elif at is "node":
        status_values = grid.status_at_node
    else:
        if where_to_place is not None:
            raise ValueError(("No status information exists for grid elements "
                              "that are not nodes or links."))

    # based on status, set where to true. support value or iterable.
    if where_to_place is None:
        where_to_place = np.ones(grid.size(at), dtype=bool)
    else:
        try:
            where.size == grid.size(at)
            where_to_place = where
        except AttributeError:
            try:
                for w in where:
                    where_to_place[status_values == w] = True
            except TypeError:
                where_to_place[status_values == where] = True

    return where_to_place


def random(grid,
           name,
           at='node',
           where=None,
           distribution="uniform",
           **kwargs):
    """Add random values to a grid.

    Parameters
    ----------
    grid : ModelGrid
    name : str
        Name of the field.
    at : str, optional
        Grid location to store values. If not given, values are
        assumed to be on `node`.
    where : optional
    distribution : str, optional
        Name of the distribution provided by the np.random
        package.
    kwargs : dict
        Keyword arguments to pass to the ``np.random``.

    Returns
    -------
    values : array
        Array of the values added to the field.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab import CORE_NODE
    >>> from landlab.values import random
    >>> np.random.seed(42)
    >>> mg = RasterModelGrid((4, 4))
    >>> values = random(mg,
    ...                 'soil__depth',
    ...                 'node',
    ...                 status=CORE_NODE,
    ...                 distribution='uniform',
    ...                 high=3.,
    ...                 low=2.)
    >>> mg.at_node['soil__depth']
    array([ 0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  2.37454012,  2.95071431,  0.        ,
            0.        ,  2.73199394,  2.59865848,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ])
    """
    where = _where_to_add_values(grid, at, where)
    _create_missing_field(grid, name, at)
    values = np.zeros(grid.size(at))

    function = np.random.__dict__[distribution]
    values[where] += function(size=np.sum(where), **kwargs)
    grid[at][name][:] += values
    return values


def plane(grid, name, at, where=None, point=(0., 0., 0), normal=(0., 0., 1.)):
    """Add a single plane defined by a point and a normal to a grid.

    Parameters
    ----------
    grid : ModelGrid
    name : str
        Name of the field.
    at : str, optional
        Grid location to store values. If not given, values are
        assumed to be on `node`.
    where : status-at-grid-element or list, optional
        A value or list of the grid element status at which values
        are added. By default, values are added to all elements.
    point : tuple, optional
        A tuple defining a point the plane goes through in the
        format (x, y, z). Default is (0., 0., 0.)
    normal : tuple, optional
        A tuple defining the normal to the plane in the format
        (dx, dy, dz). Must not be verticaly oriented. Default
        is a horizontal plane (0., 0., 1.).

    Returns
    -------
    values : array
        Array of the values added to the field.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.values import plane
    >>> mg = RasterModelGrid((4, 4))
    >>> values = plane(mg,
    ...                'soil__depth',
    ...                'node',
    ...                point=(0., 0., 0.),
    ...                normal=(-1., -1., 1.))
    >>> mg.at_node['soil__depth']
    array([ 0.,  1.,  2.,  3.,
            1.,  2.,  3.,  4.,
            2.,  3.,  4.,  5.,
            3.,  4.,  5.,  6.])

    """
    x, y = _get_x_and_y(grid, at)

    where = _where_to_add_values(grid, at, where)
    _create_missing_field(grid, name, at)
    values = _plane_function(x, y, point, normal)[where]
    grid[at][name][where] += values[where]

    return values


def _plane_function(x, y, point, normal):
    """calculate the plane function"""
    if np.isclose(normal[2], 0):
        raise ValueError("")

    constant = (point[0] * normal[0] +
                point[1] * normal[1] +
                point[2] * normal[2])
    values = ((constant
              - (normal[0] * (x - point[0]))
              - (normal[1] * (y - point[1])))
              / normal[2])

    return values


def _get_x_and_y(grid, at):
    if at is "node":
        x = grid.x_of_node
        y = grid.y_of_node
    elif at is "link":
        x = grid.x_of_link
        y = grid.y_of_link
    elif at is "cell":
        x = grid.x_of_cell
        y = grid.y_of_cell
    elif at is "face":
        x = grid.x_of_face
        y = grid.y_of_face
    else:
        raise ValueError("")
    return x, y


def constant(grid, name, at, where=None, constant=0.):
    """Add a constant to a grid.

    Parameters
    ----------
    grid : ModelGrid
    name : str
        Name of the field.
    at : str, optional
        Grid location to store values. If not given, values are
        assumed to be on `node`.
    status : status-at-grid-element or list, optional
        A value or list of the grid element status at which values
        are added. By default, values are added to all elements.
    constant : float, optional
        Constant value to add to the grid. Default is 0.

    Returns
    -------
    values : array
        Array of the values added to the field.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab import ACTIVE_LINK
    >>> from landlab.values import constant
    >>> mg = RasterModelGrid((4, 4))
    >>> values = constant(mg,
    ...                  'some_flux',
    ...                  'link',
    ...                  status=ACTIVE_LINK,
    ...                  constant=10)
    >>> mg.at_link['some_flux']
    array([  0.,   0.,   0.,   0.,  10.,  10.,   0.,  10.,  10.,  10.,   0.,
            10.,  10.,   0.,  10.,  10.,  10.,   0.,  10.,  10.,   0.,   0.,
             0.,   0.])

    """
    where = _where_to_add_values(grid, at, where)
    _create_missing_field(grid, name, at)
    values = np.zeros(grid.size(at))
    values[where] += constant
    grid[at][name][:] += values
    return values


def sine(grid, name, at, where):
    """Add a two dimentional sin wave to a grid.

    .. math::
        z = a sin ((x-x0)/px) + b sin((y-y0)/py)


    Parameters
    ----------
    grid : ModelGrid
    name : str
        Name of the field.
    at : str, optional
        Grid location to store values. If not given, values are
        assumed to be on `node`.
    status : status-at-grid-element or list, optional
        A value or list of the grid element status at which values
        are added. By default, values are added to all elements.
    point : tuple, optional
        (x0, y0)
    amplitudes :
        (a, b)
    periods :
        (px, py)

    Returns
    -------
    values : array
        Array of the values added to the field.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab import ACTIVE_LINK
    >>> from landlab.values import constant
    >>> mg = RasterModelGrid((4, 4))
    >>> values = constant(mg,
    ...                  'some_flux',
    ...                  'link',
    ...                  status=ACTIVE_LINK,
    ...                  constant=10)
    >>> mg.at_link['some_flux']
    array([  0.,   0.,   0.,   0.,  10.,  10.,   0.,  10.,  10.,  10.,   0.,
            10.,  10.,   0.,  10.,  10.,  10.,   0.,  10.,  10.,   0.,   0.,
             0.,   0.])

    """
    pass
