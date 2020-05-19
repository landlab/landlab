import numpy as np
import scipy.constants
import warnings

from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid

_OUT_OF_NETWORK = -2


class BedParcelInitializer:

    """
    This function creates a landlab DataRecord to represent parcels of sediment
    on a river network (represented by a NetworkModelGrid). The function takes
    discharge data for each link as input, as well as channel geometry
    (`channel_width`, `reach_length`, `channel_slope`) fields attached to the
    NetworkModelGrid.

    This function currently estimates median parcel grain size at a link
    according to Snyder et al. (2013), assuming a lognormal parcel grain size
    distribution.

    Wish list of future edits:
    - allow user to pass flow depth instead of discharge to determine d50
    -

    authors: Eric Hutton, Allison Pfeiffer, Muneer Ahammad

    last updated: May 2020

    Parameters
    ----------
    grid : ModelGrid
        landlab *ModelGrid* to place sediment parcels on.
    time : float, optional
        The initial time to add to the record.
    mannings_n : float, optional
        Manning's n value for all links, used to calculate median parcel grain
        size at a link.
    tau_50 : float, optional
        Shields stress for d50 at dominant discharge for all links, used to
        calculate median parcel grain size
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    rho_water : float, optional
        Density of water [kg / m^3].
    gravity : float, optional
        Accelertion due to gravity [m / s^2].
    std_dev : float, optional
        Standard deviation of lognormal distribution of grain size.
    sed_thickness : float, optional
        Sediment thickness in multiples of d84.
    abrasion_rate : float, optional
        Abrasion rate of parcels during transport in units of 1/m.
    median_number_of_starting_parcels : int, optional
        median number of parcels in a link.
    extra_parcel_attributes : str or list of str, optional
        name of user-defined parcel attribute to be added to parcel data record,
        which will be returned as an empty parcel attribute

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
        time = [0.0],
        mannings_n=0.035,
        tau_50=0.04,
        rho_sediment=2650.0,
        rho_water=1000.0,
        gravity=scipy.constants.g,
        std_dev=2.1,
        sed_thickness=4,
        abrasion_rate=0.0,
        median_number_of_starting_parcels=100,
        extra_parcel_attributes=None,

    ):
        self._time = time
        self._grid = grid
        self._mannings_n = mannings_n
        self._tau_50 = tau_50
        self._rho_sediment = rho_sediment
        self._rho_water = rho_water
        self._gravity = gravity
        self._std_dev = std_dev
        self._abrasion_rate = abrasion_rate
        self._extra_parcel_attributes = extra_parcel_attributes
        self._sed_thickness = sed_thickness
        self._median_number_of_starting_parcels = median_number_of_starting_parcels

        if not isinstance(grid, NetworkModelGrid):
            msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            raise ValueError(msg)

    def __call__(self, discharge_at_link):
        d50 = calc_d50_grain_size(
            discharge_at_link,
            self._grid.at_link["channel_width"],
            self._grid.at_link["channel_slope"],
            mannings_n=self._mannings_n,
            gravity=self._gravity,
            rho_water=self._rho_water,
            rho_sediment=self._rho_sediment,
            tau_50=self._tau_50,

        )
        self.D50 = d50
        d84 = d50 * self._std_dev

        total_parcel_volume_at_link = calc_total_parcel_volume(
            self._grid.at_link["channel_width"],
            self._grid.at_link["reach_length"],
            d84 * self._sed_thickness,
        )
        max_parcel_volume = _determine_approx_parcel_volume(
            total_parcel_volume_at_link,
            self._median_number_of_starting_parcels
            )

        variables, items = _parcel_characteristics(
            total_parcel_volume_at_link,
            max_parcel_volume,
            d50,
            self._std_dev,
            self._rho_sediment,
            self._abrasion_rate,
            self._extra_parcel_attributes
        )

        if np.min(self._sed_thickness) < 0.05:
            msg =(
            "BedParcelInitializer: Sediment thickness is unrealistically low. Minimum link sediment thickness = "
            + str(self._sed_thickness)
            +" m")
            warnings.warn(msg)

        if np.max(d50) > 0.5:
            msg =(
            "BedParcelInitializer: Median grain sizes are too large for physically realistic application of the NST. Maximum link D50 = "
            + str(d50)
            +" m")
            warnings.warn(msg)

        if np.min(d50) < 0.002:
            msg =(
            "BedParcelInitializer: The equations used in this initializer are intended for gravel bedded rivers."
            "Calculated median grain sizes are too small. Minimum link D50 = "
            + str(d50)
            +" m")
            warnings.warn(msg)

        if np.max(np.abs(self._mannings_n-0.035)) > 0.3:
            msg =(
            "BedParcelInitializer: Manning's n value is unrealistic. Value given = "
            + str(self._mannings_n))
            warnings.warn(msg)

        if np.max(np.abs(self._tau_50-0.055)) > 0.35:
            msg =(
            "BedParcelInitializer: Shields stress value is unrealistic. Value given = "
            + str(self._sed_thickness))
            warnings.warn(msg)

        if max_parcel_volume < 0.05:
            msg =(
            "BedParcelInitializer: Careful! Default parcel value is extremely small: "
            + str(max_parcel_volume)
            + ' m')
            warnings.warn(msg)

        return DataRecord(
            self._grid,
            items=items,
            time=self._time,
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
    extra_parcel_attributes
):
    n_parcels_at_link = np.ceil(total_parcel_volume_at_link / max_parcel_volume).astype(dtype=int)
    if np.min(n_parcels_at_link) <10:
        msg =(
        "BedParcelInitializer: At least one link has only "
        + str(n_parcels_at_link)
        + " parcels.")
        warnings.warn(msg)

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
        volume[offset] = (
                        total_parcel_volume_at_link[link]
                        - ((n_parcels-1)*max_parcel_volume)
                        ) # small remaining volume

        offset += n_parcels
    starting_link = element_id.copy()
    abrasion_rate = np.full_like(element_id, abrasion_rate, dtype=float)
    density = np.full_like(element_id, rho_sediment, dtype=float)

    element_id = np.expand_dims(element_id, axis=1)
    volume = np.expand_dims(volume, axis=1)
    grain_size = np.expand_dims(grain_size, axis=1)

    time_arrival_in_link = np.expand_dims(np.random.rand(np.sum(n_parcels_at_link)), axis=1)
    location_in_link = np.expand_dims(np.random.rand(np.sum(n_parcels_at_link)), axis=1)

    active_layer = np.empty_like(element_id, dtype=float)
    variables = {
        "starting_link": (["item_id"], starting_link),
        "abrasion_rate": (["item_id"], abrasion_rate),
        "density": (["item_id"], density),
        "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
        "active_layer": (["item_id", "time"], active_layer),
        "location_in_link": (["item_id", "time"], location_in_link),
        "D": (["item_id", "time"], grain_size),
        "volume": (["item_id", "time"], volume),
    }

    if extra_parcel_attributes is not None:

        for attrib in extra_parcel_attributes:
            variables[attrib] = (
                            ['item_id'],
                            np.nan*np.zeros_like(density)
                            )

    return variables, {"grid_element": "link", "element_id": element_id}

def _determine_approx_parcel_volume(
    total_parcel_volume_at_link,
    median_number_of_starting_parcels
    ):
    median_link_volume = np.median(total_parcel_volume_at_link)
    return median_link_volume / median_number_of_starting_parcels

def calc_total_parcel_volume(width, length, sediment_thickness):
    return width * length * sediment_thickness

def calc_d50_grain_size(
    dominant_discharge,
    width,
    slope,
    mannings_n,
    gravity,
    rho_water,
    rho_sediment,
    tau_50,
):
    """Calculate median grain size according to Snyder et al. (2013)

    Parameters
    ----------
    see above

    Returns
    -------
    ndarray of float
        d50.
    """

    return (
        rho_water * gravity * mannings_n ** (3 / 5) * dominant_discharge ** (3 / 5) * width ** (- 3 / 5) * slope ** (7 / 10)
    ) / (
        (rho_sediment - rho_water) * gravity * tau_50
    )
