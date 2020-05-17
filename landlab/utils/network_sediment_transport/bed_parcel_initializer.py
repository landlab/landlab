from landlab.data_record import DataRecord


def synthetic_bed_parcel_initializer(grid):
    """Initialize bed sediment on a network for the NetworkSedimentTransporter


    More description here.

    Parameters
    ----------
    grid : NetworkModelGrid

    Returns
    -------
    parcels : DataRecord

    Examples
    --------
    >>> from landlab.utils.network_sediment_transport import (
    ...     synthetic_bed_parcel_initializer)

    """
    if not isinstance(grid, NetworkModelGrid):
        msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
        raise ValueError(msg)

    # given the input arguments, keyword arguments create bed sediment.
    # return the parcel datastructure.

    parcels = None
    return parcels
