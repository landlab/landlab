#! /usr/bin/env python
"""Calculate gradients of quantities over links.

Gradient calculation functions
+++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.grid.gradients.calc_grad_at_active_link
    ~landlab.grid.gradients.calc_grad_at_link
    ~landlab.grid.gradients.calculate_gradients_at_faces
    ~landlab.grid.gradients.calculate_diff_at_links
    ~landlab.grid.gradients.calculate_diff_at_active_links

"""

import numpy as np
from landlab.utils.decorators import use_field_name_or_array, deprecated
from landlab.core.utils import radians_to_degrees


@use_field_name_or_array('node')
def calc_grad_at_link(grid, node_values, out=None):
    """Calculate gradients of node values at links.

    Calculates the gradient in `node_values` at each link in the grid,
    returning an array of length `number_of_links`.

    Construction::

        calc_grad_at_link(grid, node_values, out=None)

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name (x number of nodes)
        Values at grid nodes.
    out : ndarray, optional (x number of links)
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Gradients across active links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid(3, 4, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> z[5] = 50.0
    >>> z[6] = 36.0
    >>> calc_grad_at_link(rg, z)  # there are 17 links
    array([ 0. ,  0. ,  0. ,  0. ,  5. ,  3.6,  0. ,  5. , -1.4, -3.6,  0. ,
           -5. , -3.6,  0. ,  0. ,  0. ,  0. ])

    >>> from landlab import HexModelGrid
    >>> hg = HexModelGrid(3, 3, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation', noclobber=False)
    >>> z[4] = 50.0
    >>> z[5] = 36.0
    >>> calc_grad_at_link(hg, z)  # there are 11 faces
    array([ 0. ,  0. ,  0. ,  5. ,  5. ,  3.6,  3.6,  0. ,  5. , -1.4, -3.6,
            0. , -5. , -5. , -3.6, -3.6,  0. ,  0. ,  0. ])
    """
    if out is None:
        out = grid.empty(at='link')
    return np.divide(node_values[grid.node_at_link_head] -
                     node_values[grid.node_at_link_tail],
                     grid.length_of_link, out=out)


@deprecated(use='calc_grad_at_link', version='1.0beta')
def calc_grad_of_active_link(grid, node_values, out=None):
    """Calculate gradients at active links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 4))
    >>> z = np.array([0., 0., 0., 0.,
    ...               1., 1., 1., 1.,
    ...               3., 3., 3., 3.])
    >>> grid.calc_grad_of_active_link(z)
    array([ 1.,  1.,  0.,  0.,  0.,  2.,  2.])

    This method is *deprecated*. Instead, use ``calc_grad_at_link``.

    >>> vals = grid.calc_grad_at_link(z)
    >>> vals[grid.active_links]
    array([ 1.,  1.,  0.,  0.,  0.,  2.,  2.])
    """
    return calc_grad_at_active_link(grid, node_values, out)


@deprecated(use='calc_grad_at_link', version='1.0beta')
@use_field_name_or_array('node')
def calc_grad_at_active_link(grid, node_values, out=None):
    """Calculate gradients of node values over active links.

    Calculates the gradient in *quantity* node values at each active link in
    the grid.

    Construction::

        calc_grad_at_active_link(grid, node_values, out=None)

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Gradients across active links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 4))
    >>> z = np.array([0., 0., 0., 0.,
    ...               1., 1., 1., 1.,
    ...               3., 3., 3., 3.])
    >>> grid.calc_grad_at_active_link(z)
    array([ 1.,  1.,  0.,  0.,  0.,  2.,  2.])

    This method is *deprecated*. Instead, use ``calc_grad_at_link``.

    >>> vals = grid.calc_grad_at_link(z)
    >>> vals[grid.active_links]
    array([ 1.,  1.,  0.,  0.,  0.,  2.,  2.])
    """
    if out is None:
        out = grid.empty(at='active_link')
    return np.divide(node_values[grid._activelink_tonode] -
                     node_values[grid._activelink_fromnode],
                     grid.length_of_link[grid.active_links], out=out)


@deprecated(use='calc_grad_at_link', version='1.0beta')
@use_field_name_or_array('node')
def calculate_gradients_at_faces(grid, node_values, out=None):
    """Calculate gradients of node values over faces.

    Calculate and return gradient in *node_values* at each face in the grid.
    Gradients are calculated from the nodes at either end of the link that
    crosses each face.

    Construction::

        calculate_gradients_at_faces(grid, node_values, out=None)

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray (x number of faces)
        Gradients across faces.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid(3, 4, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> z[5] = 50.0
    >>> z[6] = 36.0
    >>> calculate_gradients_at_faces(rg, z)  # there are 7 faces
    array([ 5. ,  3.6,  5. , -1.4, -3.6, -5. , -3.6])

    >>> from landlab import HexModelGrid
    >>> hg = HexModelGrid(3, 3, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation', noclobber=False)
    >>> z[4] = 50.0
    >>> z[5] = 36.0
    >>> calculate_gradients_at_faces(hg, z)  # there are 11 faces
    array([ 5. ,  5. ,  3.6,  3.6,  5. , -1.4, -3.6, -5. , -5. , -3.6, -3.6])
    """
    if out is None:
        out = grid.empty(at='face')
    laf = grid.link_at_face
    return np.divide(node_values[grid.node_at_link_head[laf]] -
                     node_values[grid.node_at_link_tail[laf]],
                     grid.length_of_link[laf], out=out)


@use_field_name_or_array('node')
def calc_diff_at_link(grid, node_values, out=None):
    """Calculate differences of node values over links.

    Calculates the difference in quantity *node_values* at each link in the
    grid.

    Construction::

        calc_diff_at_link(grid, node_values, out=None)

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Differences across links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> rmg = RasterModelGrid((3, 3))
    >>> z = np.zeros(9)
    >>> z[4] = 1.
    >>> rmg.calc_diff_at_link(z)
    array([ 0.,  0.,  0.,  1.,  0.,  1., -1.,  0., -1.,  0.,  0.,  0.])
    """
    if out is None:
        out = grid.empty(at='link')
    node_values = np.asarray(node_values)
    return np.subtract(node_values[grid.node_at_link_head],
                       node_values[grid.node_at_link_tail], out=out)


@deprecated(use='calc_diff_at_link', version='1.0beta')
@use_field_name_or_array('node')
def calculate_diff_at_links(grid, node_values, out=None):
    """Calculate differences of node values over links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 3))
    >>> z = np.zeros(9)
    >>> z[4] = 1.

    >>> grid.calculate_diff_at_links(z)
    array([ 0.,  0.,  0.,  1.,  0.,  1., -1.,  0., -1.,  0.,  0.,  0.])

    >>> grid.calc_diff_at_link(z)
    array([ 0.,  0.,  0.,  1.,  0.,  1., -1.,  0., -1.,  0.,  0.,  0.])
    """
    return calc_diff_at_link(grid, node_values, out)


@deprecated(use='calc_diff_at_link', version='1.0beta')
@use_field_name_or_array('node')
def calculate_diff_at_active_links(grid, node_values, out=None):
    """Calculate differences of node values over active links.

    Calculates the difference in quantity *node_values* at each active link
    in the grid.

    Construction::

        calculate_diff_at_active_links(grid, node_values, out=None)

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    node_values : ndarray or field name
        Values at grid nodes.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Differences across active links.
    """
    if out is None:
        out = grid.empty(at='active_link')
    node_values = np.asarray(node_values)
    return np.subtract(node_values[grid._activelink_tonode],
                       node_values[grid._activelink_fromnode], out=out)


def calc_unit_normal_at_patch(grid, elevs='topographic__elevation'):
    """Calculate and return the unit normal vector <a, b, c> to a patch.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    elevs : str or ndarray, optional
        Field name or array of node values.

    Returns
    -------
    nhat : num-patches x length-3 array
        The unit normal vector <a, b, c> to each patch.

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> mg = HexModelGrid(3, 3)
    >>> z = mg.node_x * 3. / 4.
    >>> mg.calc_unit_normal_at_patch(z)
    array([[-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8],
           [-0.6,  0. ,  0.8]])
    """
    try:
        z = grid.at_node[elevs]
    except TypeError:
        z = elevs
    # conceptualize patches as sets of 3 nodes, PQR
    diff_xyz_PQ = np.empty((grid.number_of_patches, 3))
    # ^this is the vector (xQ-xP, yQ-yP, zQ-yP)
    diff_xyz_PR = np.empty((grid.number_of_patches, 3))
    P = grid.nodes_at_patch[:, 0]
    Q = grid.nodes_at_patch[:, 1]
    R = grid.nodes_at_patch[:, 2]
    x_P = grid.node_x[P]
    y_P = grid.node_y[P]
    z_P = z[P]
    diff_xyz_PQ[:, 0] = grid.node_x[Q] - x_P
    diff_xyz_PQ[:, 1] = grid.node_y[Q] - y_P
    diff_xyz_PQ[:, 2] = z[Q] - z_P
    diff_xyz_PR[:, 0] = grid.node_x[R] - x_P
    diff_xyz_PR[:, 1] = grid.node_y[R] - y_P
    diff_xyz_PR[:, 2] = z[R] - z_P
    # cross product is orthogonal to both vectors, and is the normal
    # n = <a, b, c>, where plane is ax + by + cz = d
    nhat = np.cross(diff_xyz_PQ, diff_xyz_PR)  # <a, b, c>
    nmag = np.sqrt(np.square(nhat).sum(axis=1))

    return nhat / nmag.reshape(grid.number_of_patches, 1)


def calc_slope_at_patch(grid, elevs='topographic__elevation',
                        unit_normal=None):
    """Calculate the slope (positive magnitude of gradient) at patches.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    elevs : str or ndarray, optional
        Field name or array of node values.
    unit_normal : array with shape (num_patches, 3) (optional)
        The unit normal vector to each patch, if already known.

    Returns
    -------
    slopes_at_patch : n_patches-long array
        The slope (positive gradient magnitude) of each patch.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((4, 5))
    >>> z = mg.node_x
    >>> S = mg.calc_slope_at_patch(elevs=z)
    >>> S.size == mg.number_of_patches
    True
    >>> np.allclose(S, np.pi / 4.)
    True
    """
    if unit_normal is not None:
        assert unit_normal.shape[1] == 3
        nhat = unit_normal
    else:
        nhat = grid.calc_unit_normal_at_patch(elevs)
    dotprod = nhat[:, 2]  # by definition
    cos_slopes_at_patch = dotprod  # ...because it's now a unit vector
    slopes_at_patch = np.arccos(cos_slopes_at_patch)

    return slopes_at_patch


def calc_grad_at_patch(grid, elevs='topographic__elevation',
                       unit_normal=None, slope_magnitude=None):
    """Calculate the components of the gradient at each patch.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    elevs : str or ndarray, optional
        Field name or array of node values.
    unit_normal : array with shape (num_patches, 3) (optional)
        The unit normal vector to each patch, if already known.
    slope_magnitude : array with size num_patches (optional)
        The slope of each patch, if already known.

    Returns
    -------
    gradient_tuple : (x_component_at_patch, y_component_at_patch)
        Len-2 tuple of arrays giving components of gradient in the x and y
        directions, in the units of *units*.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((4, 5))
    >>> z = mg.node_y
    >>> (x_grad, y_grad) = mg.calc_grad_at_patch(elevs=z)
    >>> np.allclose(y_grad, np.pi / 4.)
    True
    >>> np.allclose(x_grad, 0.)
    True
    """
    if unit_normal is not None:
        assert unit_normal.shape[1] == 3
        nhat = unit_normal
    else:
        nhat = grid.calc_unit_normal_at_patch(elevs)
    if slope_magnitude is not None:
        assert slope_magnitude.size == grid.number_of_patches
        slopes_at_patch = slope_magnitude
    else:
        slopes_at_patch = grid.calc_slope_at_patch(elevs=elevs,
                                                   unit_normal=nhat)
    theta = np.arctan2(- nhat[:, 1], - nhat[:, 0])
    x_slope_patches = np.cos(theta) * slopes_at_patch
    y_slope_patches = np.sin(theta) * slopes_at_patch

    return (x_slope_patches, y_slope_patches)


def calc_slope_at_node(grid, elevs='topographic__elevation',
                       method='patch_mean', return_components=False, **kwds):
    """Array of slopes at nodes, averaged over neighboring patches.

    Produces a value for node slope (i.e., mean gradient magnitude)
    at each node in a manner analogous to a GIS-style slope map.
    It averages the gradient on each of the
    patches surrounding the node, creating a value for node slope that
    better incorporates nonlocal elevation information. Directional
    information can still be returned through use of the return_components
    keyword.

    Note that under these definitions, it is not always true that::

        mag, cmp = mg.calc_slope_at_node(z)
        mag ** 2 == cmp[0] ** 2 + cmp[1] ** 2  # not always true

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    elevs : str or ndarray, optional
        Field name or array of node values.
    method : {'patch_mean', 'Horn'}
        By equivalence to the raster version, `'patch_mean'` returns a scalar
        mean on the patches; `'Horn'` returns a vector mean on the patches.
    return_components : bool
        If True, return a tuple, (array_of_magnitude,
        (array_of_slope_x_radians, array_of_slope_y_radians)).
        If false, return an array of floats of the slope magnitude.

    Returns
    -------
    float array or length-2 tuple of float arrays
        If return_components, returns (array_of_magnitude,
        (array_of_slope_x_radians, array_of_slope_y_radians)).
        If not return_components, returns an array of slope magnitudes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RadialModelGrid, RasterModelGrid
    >>> mg = RasterModelGrid((4, 5), 1.)
    >>> z = mg.node_x
    >>> slopes = mg.calc_slope_at_node(elevs=z)
    >>> np.allclose(slopes, 45. / 180. * np.pi)
    True

    >>> mg = RasterModelGrid((4, 5), 1.)
    >>> z = - mg.node_y
    >>> slope_mag, cmp = mg.calc_slope_at_node(elevs=z,
    ...                                        return_components=True)
    >>> np.allclose(slope_mag, np.pi / 4.)
    True
    >>> np.allclose(cmp[0], 0.)
    True
    >>> np.allclose(cmp[1], - np.pi / 4.)
    True

    >>> mg = RadialModelGrid(num_shells=9)
    >>> z = mg.radius_at_node
    >>> slopes = mg.calc_slope_at_node(elevs=z)
    >>> mean_ring_slope = []
    >>> for i in range(10):
    ...     mean_ring_slope.append(
    ...         slopes[np.isclose(mg.radius_at_node, i)].mean())

    Notice the small amounts of numerical error here:

    >>> target_mean_ring_slope = [0.85707194785013108, 0.79363155567711452,
    ...                           0.77922185867135429, 0.78359813570962411,
    ...                           0.78433070957439543, 0.78452745144699965,
    ...                           0.78477643475446901, 0.78506472422668094,
    ...                           0.78505793680521629, 0.78661256633611021]
    >>> np.allclose(mean_ring_slope, target_mean_ring_slope)
    True
    """
    if method not in ('patch_mean', 'Horn'):
        raise ValueError('method name not understood')

    patches_at_node = np.ma.masked_where(
        grid.patches_at_node == -1, grid.patches_at_node, copy=False)

    nhat = grid.calc_unit_normal_at_patch(elevs=elevs)
    slopes_at_patch = grid.calc_slope_at_patch(elevs=elevs,
                                               unit_normal=nhat)

    # now CAREFUL - patches_at_node is MASKED
    slopes_at_node_unmasked = slopes_at_patch[patches_at_node]
    slopes_at_node_masked = np.ma.array(slopes_at_node_unmasked,
                                        mask=patches_at_node.mask)
    slope_mag = np.mean(slopes_at_node_masked, axis=1).data

    if return_components or method == 'Horn':
        (x_slope_patches, y_slope_patches) = grid.calc_grad_at_patch(
            elevs=elevs, unit_normal=nhat,
            slope_magnitude=slopes_at_patch)
        x_slope_unmasked = x_slope_patches[patches_at_node]
        x_slope_masked = np.ma.array(x_slope_unmasked,
                                     mask=patches_at_node.mask)
        x_slope = np.mean(x_slope_masked, axis=1).data
        y_slope_unmasked = y_slope_patches[patches_at_node]
        y_slope_masked = np.ma.array(y_slope_unmasked,
                                     mask=patches_at_node.mask)
        y_slope = np.mean(y_slope_masked, axis=1).data
        mean_grad_x = x_slope
        mean_grad_y = y_slope

        if method == 'Horn':
            slope_mag = np.arctan(np.sqrt(np.tan(y_slope_masked) ** 2 +
                                          np.tan(x_slope_masked) ** 2))
            return slope_mag
        else:
            return slope_mag, (mean_grad_x, mean_grad_y)

    else:
        return slope_mag


def calc_aspect_at_node(grid, slope_component_tuple=None,
                        elevs='topographic__elevation', unit='degrees'):
    """Get array of aspect of a surface.

    Calculates at returns the aspect of a surface. Aspect is returned as
    radians clockwise of north, unless input parameter units is set to
    'degrees'.

    If slope_component_tuple is provided, i.e., (slope_x, slope_y), the
    aspect will be calculated from these data.

    If it is not, it will be derived from elevation data at the nodes,
    which can either be a string referring to a grid field (default:
    'topographic__elevation'), or an nnodes-long numpy array of the
    values themselves.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    slope_component_tuple : (slope_x_array, slope_y_array) (optional)
        Tuple of components of slope in the x and y directions, defined
        on nodes, if already known. If not, provide *elevs*.
    elevs : str or array (optional)
        Node field name or node array of elevations.
        If *slope_component_tuple* is not provided, must be set, but unused
        otherwise.
    unit : {'degrees', 'radians'}
        Controls the unit that the aspect is returned as.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> mg = RasterModelGrid((4, 4))
    >>> z = mg.node_x ** 2 + mg.node_y ** 2
    >>> mg.calc_aspect_at_node(elevs=z)
    array([ 225.        ,  240.16585039,  255.2796318 ,  258.69006753,
            209.83414961,  225.        ,  243.54632481,  248.77808974,
            194.7203682 ,  206.45367519,  225.        ,  231.94498651,
            191.30993247,  201.22191026,  218.05501349,  225.        ])
    >>> z = z.max() - z
    >>> mg.calc_aspect_at_node(elevs=z)
    array([ 45.        ,  60.16585039,  75.2796318 ,  78.69006753,
            29.83414961,  45.        ,  63.54632481,  68.77808974,
            14.7203682 ,  26.45367519,  45.        ,  51.94498651,
            11.30993247,  21.22191026,  38.05501349,  45.        ])

    >>> mg = RasterModelGrid((4, 4), (2., 3.))
    >>> z = mg.node_x ** 2 + mg.node_y ** 2
    >>> mg.calc_aspect_at_node(elevs=z)
    array([ 236.30993247,  247.52001262,  259.97326008,  262.40535663,
            220.75264634,  234.41577266,  251.13402374,  255.29210302,
            201.54258265,  215.47930877,  235.73541937,  242.24162456,
            196.69924423,  209.43534223,  229.19345757,  236.30993247])

    Note that a small amount of asymmetry arises at the grid edges due
    to the "missing" nodes beyond the edge of the grid.
    """
    if slope_component_tuple:
        if not isinstance(slope_component_tuple, (tuple, list)):
            raise TypeError('slope_component_tuple must be tuple')
        if len(slope_component_tuple) != 2:
            raise ValueError('slope_component_tuple must be of length 2')
    else:
        try:
            elev_array = grid.at_node[elevs]
        except (KeyError, TypeError):
            assert elevs.size == grid.number_of_nodes
            elev_array = elevs

        _, slope_component_tuple = grid.calc_slope_at_node(
            elevs=elev_array, return_components=True)

    angle_from_x_ccw = np.arctan2(
        - slope_component_tuple[1], - slope_component_tuple[0])

    if unit == 'degrees':
        return radians_to_degrees(angle_from_x_ccw)
    elif unit == 'radians':
        angle_from_north_cw = (5. * np.pi / 2. -
                               angle_from_x_ccw) % (2. * np.pi)
        return angle_from_north_cw
    else:
        raise TypeError("unit must be 'degrees' or 'radians'")
