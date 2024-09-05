from landlab import Component
from landlab.grid.network import NetworkModelGrid


class SedimentPulserBase(Component):
    """Base class of :class:`~.SedimentPulserAtLinks` and :class:`~.SedimentPulserEachParcel`.

    :class:`~.SedimentPulserAtLinks` and :class:`~.SedimentPulserEachParcel` run the
    landlab :class:`~.DataRecord` :meth:`~.DataRecord.add_item` method on a
    :class:`~.DataRecord` configured for :class:`~.NetworkSedimentTransporter`.


    .. codeauthor: Jeff Keck, Allison Pfeiffer, Shelby Ahrendt
                   (with help from Eric Hutton and Katy Barnhart)

    Parameters
    ----------
    grid : ModelGrid
        landlab *ModelGrid* to place sediment parcels on.
    parcels: landlab DataRecord
        Tracks parcel location and variables
    D50: float, optional
        median grain size [m]
    D84_D50: float, optional
        ratio of 84th percentile grain size to the median grain size
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    parcel_volume : float, optional
        parcel volume used for all parcels that do not have a specified volume
    abrasion_rate: float, optional
        volumetric abrasion exponent [1/m]


    Examples
    --------
    >>> import numpy as np
    >>> from landlab import NetworkModelGrid

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> grid.at_link["channel_width"] = np.full(grid.number_of_links, 1.0)  # m
    >>> grid.at_link["channel_slope"] = np.full(grid.number_of_links, 0.01)  # m / m
    >>> grid.at_link["reach_length"] = np.full(grid.number_of_links, 100.0)  # m
    >>> make_pulse_base = SedimentPulserBase(grid)
    >>> make_pulse_base._parcels

    SedimentPulserBase does not have any methods for adding a pulse

    >>> a_pulse = make_pulse_base()
    Traceback (most recent call last):
    ...
    NotImplementedError: the base component has no call method
    """

    _name = "SedimentPulserBase"

    _unit_agnostic = False

    _info = {}  # works with the DataRecord

    def __init__(
        self,
        grid,
        parcels=None,
        D50=0.05,
        D84_D50=2.1,
        rho_sediment=2650.0,
        parcel_volume=0.5,
        abrasion_rate=0.0,
    ):
        self._grid = grid
        self._parcels = parcels
        self._D50 = D50
        self._D84_D50 = D84_D50
        self._rho_sediment = rho_sediment
        self._parcel_volume = parcel_volume
        self._abrasion_rate = abrasion_rate

        if not isinstance(grid, NetworkModelGrid):
            raise ValueError(
                "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            )

    def __call__(self):
        """__call__ is not implemented for this component."""
        raise NotImplementedError("the base component has no call method")
