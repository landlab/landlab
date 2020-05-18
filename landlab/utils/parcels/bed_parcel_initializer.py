import numpy as np
import scipy.constants

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
    >>> from landlab.utils.parcels import synthetic_bed_parcel_initializer

    """
    if not isinstance(grid, NetworkModelGrid):
        msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
        raise ValueError(msg)

    # given the input arguments, keyword arguments create bed sediment.
    # return the parcel datastructure.

    parcels = None
    return parcels


_OUT_OF_NETWORK = -2


class BedParcelInitializer:

    """

    Parameters
    ----------
    grid : ModelGrid
        landlab *ModelGrid* to place sediment parcels on.
    mannings_n : float, optional
        Mannings's n.
    tau_critical : float, optional
        Critical shear stress.
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    rho_water : float, optional
        Density of water [kg / m^3].
    gravity : float, optional
        Accelertion due to gravity [m / s^2].
    std_dev : float, optional
        Standard deviation of lognormal distribution of grain size.

    Examples
    --------
    >>> from landlab import NetworkModelGrid
    >>> from landlab.utils.parcels import BedParcelInitializer

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> grid.at_link["channel_width"] = np.full(grid.number_of_links, 1.0)  # m
    >>> grid.at_link["channel_slope"] = np.full(grid.number_of_links, .01)  # m / m
    >>> grid.at_link["reach_length"] = np.full(grid.number_of_links, 100.0)  # m

    >>> initialize_parcels = BedParcelInitializer(grid)
    >>> discharge_at_link = np.full(grid.number_of_links, 10.0)  # m^3 / s
    >>> parcels = initialize_parcels(discharge_at_link)
    """
    def __init__(
        self,
        grid,
        mannings_n=0.035,
        tau_critical=0.04,
        rho_sediment=2650.0,
        rho_water=1000.0,
        gravity=scipy.constants.g,
        std_dev=2.1,
    ):
        self._grid = grid
        self._mannings_n = mannings_n
        self._tau_critical = tau_critical
        self._rho_sediment = rho_sediment
        self._rho_water = rho_water
        self._gravity = gravity
        self._std_dev = std_dev

    def __call__(self, discharge_at_link):
        d50 = calc_d50_grain_size(
            discharge_at_link,
            self._grid.at_link["channel_width"],
            self._grid.at_link["channel_slope"],
            mannings_n=self._mannings_n,
            gravity=self._gravity,
            rho_water=self._rho_water,
            rho_sediment=self._rho_sediment,
            tau_critical=self._tau_critical,
        )
        d84 = d50 * self._std_dev

        total_parcel_volume_at_link = calc_total_parcel_volume(
            self._grid.at_link["channel_width"],
            self._grid.at_link["reach_length"],
            d84 * 2.0 * 2.0,
        )
        max_parcel_volume = _determine_approx_parcel_volume(total_parcel_volume_at_link)

        variables, items = _parcel_characteristics(
            total_parcel_volume_at_link,
            max_parcel_volume,
            d50,
            self._std_dev,
            self._rho_sediment,
            0.0,
        )

        return DataRecord(
            self._grid,
            items=items,
            time=[0.0],
            data_vars=variables,
            dummy_elements={"link": [_OUT_OF_NETWORK]},
        )

def _parcel_characteristics(
    total_parcel_volume_at_link,
    max_parcel_volume,
    d50,
    std_dev,
    rho_sediment,
    abrasion_rate,
):
    n_parcels_at_link = np.ceil(total_parcel_volume_at_link / max_parcel_volume).astype(dtype=int)

    element_id = np.empty(np.sum(n_parcels_at_link), dtype=int)

    # volume = np.full(np.sum(n_parcels_at_link), max_parcel_volume, dtype=float)
    volume = np.full_like(element_id, max_parcel_volume, dtype=float)
    grain_size = np.empty_like(element_id, dtype=float)
    offset = 0
    for link, n_parcels in enumerate(n_parcels_at_link):
        element_id[offset:offset + n_parcels] = link
        grain_size[offset:offset + n_parcels] = np.random.lognormal(
            np.log(d50[link]), np.log(std_dev), n_parcels
        )
        volume[offset] = total_parcel_volume_at_link[link] % n_parcels
        offset += n_parcels
    starting_link = element_id.copy()
    abrasion_rate = np.full_like(element_id, abrasion_rate, dtype=float)
    density = np.full_like(element_id, rho_sediment, dtype=float)

    element_id = np.expand_dims(element_id, axis=1)
    volume = np.expand_dims(volume, axis=1)
    grain_size = np.expand_dims(grain_size, axis=1)

    time_arrival_in_link = np.expand_dims(np.random.rand(np.sum(n_parcels_at_link)), axis=1)
    location_in_link = np.expand_dims(np.random.rand(np.sum(n_parcels_at_link)), axis=1)
    # abrasion_rate = np.full(np.sum(n_parcels_at_link), abrasion_rate, dtype=float)
    # starting_link = element_id.copy()
    # active_layer = np.empty(np.sum(n_parcels_at_link), dtype=float)
    active_layer = np.empty_like(element_id, dtype=float)

    return {
        "starting_link": (["item_id"], starting_link),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], density),
        "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
        "active_layer": (["item_id", "time"], active_layer),
        "location_in_link": (["item_id", "time"], location_in_link),
        "D": (["item_id", "time"], grain_size),
        "volume": (["item_id", "time"], volume),
    }, {"grid_element": "link", "element_id": element_id}


def _determine_approx_parcel_volume(total_parcel_volume_at_link):
    min_link_volume = np.min(total_parcel_volume_at_link)
    min_number_of_starting_parcels = 100
    return min_link_volume / min_number_of_starting_parcels


def calc_total_parcel_volume(width, length, sediment_thickness):
    return width * length * sediment_thickness


def calc_dominant_discharge(discharge):
    """Calculate dominant discharge.

    Parameters
    ----------
    discharge : ndarray of float
        Time series of discharge.

    Returns
    -------
    ndarray of float
        Dominant discharge.
    """
    raise NotImplementedError("calc_dominant_discharge")


def calc_d50_grain_size(
    dominant_discharge,
    width,
    slope,
    mannings_n=0.45,
    gravity=scipy.constants.g,
    rho_water=1000.0,
    rho_sediment=2650.0,
    tau_critical=0.04,
):
    return (
        rho_water * gravity * mannings_n **  3 / 5 * dominant_discharge ** 3 / 5 * width ** (- 3 / 5) * slope ** (7 / 10)
    ) / (
        (rho_sediment - rho_water) * gravity * tau_critical
    )
