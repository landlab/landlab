import numpy as np

from ..voronoi.voronoi import VoronoiGraph


def number_of_nodes(shape):
    """Get the number of nodes in a hex graph.

    Parameters
    ----------
    shape : tuple of int
        Number of rows and columns of the hex grid. The first value
        is the number of nodes in the first column and the second the
        number of nodes in the first column.

    Examples
    --------
    >>> from landlab.graph.hex.hex import number_of_nodes
    >>> number_of_nodes((3, 2))
    7
    >>> number_of_nodes((2, 4))
    9
    >>> number_of_nodes((4, 2))
    10
    """
    if shape[0] % 2 == 0:
        return (2 * shape[1] + 1) * (shape[0] // 2)
    else:
        return (2 * shape[1] + 1) * (shape[0] // 2) + shape[1]


def create_xy_of_node(shape, spacing=1., origin=(0., 0.),
                      orientation='horizontal'):
    """Create arrays of coordinates of a node on a hex grid.

    Parameters
    ----------
    shape : tuple of int
        Number of rows and columns of the hex grid. The first value
        is the number of nodes in the first column and the second the
        number of nodes in the first column.
    spacing : float, optional
        Length of links.
    origin : tuple of float, optional
        Coordinates of lower-left corner of the grid.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph.hex.hex import create_xy_of_node

    >>> x, y = create_xy_of_node((3, 2))
    >>> x
    array([ 0.5,  1.5,  0. ,  1. ,  2. ,  0.5,  1.5])
    >>> y / (np.sqrt(3) / 2.)
    array([ 0.,  0.,  1.,  1.,  1.,  2.,  2.])
    >>> x, y = create_xy_of_node((2, 2))
    >>> x
    array([ 0.5,  1.5,  0. ,  1. ,  2. ])
    >>> y / (np.sqrt(3) / 2.)
    array([ 0.,  0.,  1.,  1.,  1.])

    >>> x, y = create_xy_of_node((2, 2), spacing=2, origin=(1, 2))
    >>> x
    array([ 3.,  5.,  2.,  4.,  6.])
    >>> (y - 1) / (np.sqrt(3) / 2.)
    array([ 0.,  0.,  2.,  2.,  2.])
    """
    from .ext.hex import get_xy_of_node

    if orientation == 'vertical':
        return create_xy_of_node((shape[1], shape[0]), spacing=spacing,
                                 origin=(origin[1], origin[0]),
                                 orientation='horizontal')[::-1]

    n_nodes = number_of_nodes(shape)

    x_of_node = np.empty((n_nodes, ), dtype=float)
    y_of_node = np.empty((n_nodes, ), dtype=float)

    get_xy_of_node(shape, x_of_node, y_of_node)

    x_of_node *= spacing
    x_of_node += origin[1]
    y_of_node *= (spacing * np.sin(np.pi / 3.))
    y_of_node += origin[0]

    return (x_of_node, y_of_node)


class HexGraph(VoronoiGraph):

    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import HexGraph

    >>> graph = HexGraph((3, 2))
    >>> graph.number_of_nodes
    7
    >>> np.round(graph.y_of_node * 2. / np.sqrt(3))
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  0.,  1.,  1.,  1.,  2.,  2.])
    >>> graph.x_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5,  1.5,  0. ,  1. ,  2. ,  0.5,  1.5])
    """

    def __init__(self, shape, spacing=1., origin=(0., 0.),
                 orientation='horizontal'):
        """Create a structured grid of triangles.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of the hex grid. The first value
            is the number of nodes in the first column and the second the
            number of nodes in the first column.
        spacing : float, optional
            Length of links.
        origin : tuple of float, optional
            Coordinates of lower-left corner of the grid.
        """
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError('spacing must be a float')

        if orientation not in ('horizontal', 'vertical'):
            raise ValueError('orientation not understood')

        x_of_node, y_of_node = create_xy_of_node(shape, spacing=spacing,
                                                 origin=origin,
                                                 orientation=orientation)

        super(HexGraph, self).__init__(
            (y_of_node, x_of_node), xy_sort=True, rot_sort=True,
            max_node_spacing=shape[1] + 1)
