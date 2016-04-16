#! /usr/bin/env python
"""Calculate vector divergence and related quantities at nodes or cells."""
import numpy as np
from landlab.utils.decorators import use_field_name_or_array


@use_field_name_or_array('face')
def calc_net_face_flux_at_cells(grid, unit_flux_at_faces, out=None):
    """Calculate net face fluxes at cells.

    Given a flux per unit width across each face in the grid, calculate the net
    outflux (or influx, if negative) at each cell.

    Construction::

        calc_net_face_flux_at_cells(grid, unit_flux_at_faces, out=None)

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    unit_flux_at_faces : ndarray or field name
        Flux per unit width associated with faces.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray (x number of cells)
        Net flux at cells.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid(3, 4, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> z[5] = 50.0
    >>> z[6] = 36.0
    >>> fg = rg.calculate_gradients_at_faces(z)  # there are 7 faces
    >>> fg
    array([ 5. ,  3.6,  5. , -1.4, -3.6, -5. , -3.6])
    >>> calc_net_face_flux_at_cells(rg, -fg)
    array([-164.,  -94.])

    >>> from landlab import HexModelGrid
    >>> hg = HexModelGrid(3, 3, 10.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation', noclobber=False)
    >>> z[4] = 50.0
    >>> z[5] = 36.0
    >>> fg = hg.calculate_gradients_at_faces(z)  # there are 11 faces
    >>> fg
    array([ 5. ,  5. ,  3.6,  3.6,  5. , -1.4, -3.6, -5. , -5. , -3.6, -3.6])
    >>> nffc = calc_net_face_flux_at_cells(hg, -fg)
    >>> np.round(nffc)
    array([-152.,  -96.])

    Notes
    -----
    This is essentially a line integral for the fluxes along the boundaries of
    each cell. Hence, the resulting output has dimensions of total flux (so,
    if the unit flux happens to be mass per time per face width, the output
    will be in mass per unit time).
    """
    if out is None:
        out = grid.empty(centering='cell')
    total_flux = unit_flux_at_faces * grid.face_width
    out = np.zeros(grid.number_of_cells)
    fac = grid.faces_at_cell
    for c in range(grid.link_dirs_at_node.shape[1]):
        out += total_flux[fac[:,c]] \
        * grid.link_dirs_at_node[grid.node_at_cell,c]
    return out

