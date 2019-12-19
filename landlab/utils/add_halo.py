import numpy as np


def add_halo(data, halo=1, halo_value=None):
    """Add a halo of no data value to data.

    Parameters
    ----------
    data : array-like
        Array to add the halo to.
    halo : int, optional
        The size of the halo.
    halo_value : float, optional
        Value to fill the halo with. If not provided, the new data will is not
        initialized.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils import add_halo
    >>> data = np.array([[1, 2, 3], [4, 5, 6]])
    >>> add_halo(data, halo_value=9)
    array([[9, 9, 9, 9, 9], [9, 1, 2, 3, 9], [9, 4, 5, 6, 9], [9, 9, 9, 9, 9]])
    """
    if halo < 0:
        raise ValueError("halo must be greater than or equal to zero")
    elif halo == 0:
        return data.copy()

    data_with_halo = np.empty([dim + 2 * halo for dim in data.shape], dtype=data.dtype)
    data_with_halo[halo:-halo, halo:-halo] = data
    if halo_value is not None:
        data_with_halo[:halo, :] = halo_value
        data_with_halo[-halo:, :] = halo_value
        data_with_halo[:, :halo] = halo_value
        data_with_halo[:, -halo:] = halo_value

    return data_with_halo
