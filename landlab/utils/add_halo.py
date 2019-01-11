import numpy as np

def add_halo(data, halo, shape, nodata_value):
    """Add a halo of no data value to data.

    Parameters
    ----------
    data : array-like
    halo : int
    shape : tuple
        Shape of the data
    nodata_value : float

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils import add_halo
    >>> data = np.array([1, 2, 3, 4, 5, 6])
    >>> halo_data, new_shape = add_halo(data, 1, (2,3), 9)
    >>> halo_data
    array([9, 9, 9, 9, 9, 9, 1, 2, 3, 9, 9, 4, 5, 6, 9, 9, 9, 9, 9, 9])
    >>> new_shape
    (4, 5)


    """
    new_shape = (shape[0] + 2 * halo, shape[1] + 2 * halo)

    helper_row = np.ones(new_shape[1]) * nodata_value
    # for the first halo row(s), add num cols worth of nodata vals to data
    for i in range(0, halo):
        data = np.insert(data, 0, helper_row)
    # then for header['nrows'] add halo number nodata vals, header['ncols']
    # of data, then halo number of nodata vals
    helper_row_ends = np.ones(halo) * nodata_value
    for i in range(halo, shape[0] + halo):
        # this adds at the beginning of the row
        data = np.insert(data, i * new_shape[1], helper_row_ends)
        # this adds at the end of the row
        data = np.insert(data, (i + 1) * new_shape[1] - halo, helper_row_ends)
    # for the last halo row(s), add num cols worth of nodata vals to data
    for i in range(shape[0] + halo, new_shape[0]):
        data = np.insert(data, data.size, helper_row)
    return data, new_shape
