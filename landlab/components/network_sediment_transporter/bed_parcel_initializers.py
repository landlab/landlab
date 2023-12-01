"""
Landlab components to initialize river bed sediment "parcels", represented as
items in a landlab :class:`~.DataRecord`, in each link in a river network (represented
by a landlab :class:`~.NetworkModelGrid`). The different *BedParcelInitializers* allow
the user to define the median grain size on a given link several different ways.

.. codeauthor:: Eric Hutton, Allison Pfeiffer, Muneer Ahammad, and Jon Czuba
"""

import warnings

import numpy as np
import scipy.constants

from landlab import Component
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid

_OUT_OF_NETWORK = -2


class BedParcelInitializerBase(Component):
    def __init__(
        self,
        grid,
        time=0.0,
        tau_c_50=0.04,
        rho_sediment=2650.0,
        rho_water=1000.0,
        gravity=scipy.constants.g,
        D84_D50=2.1,
        sed_thickness=2,
        abrasion_rate=0.0,
        median_number_of_starting_parcels=100,
        extra_parcel_attributes=None,
        rng=None,
    ):
        if rng is None:
            rng = np.random.default_rng()
        elif isinstance(rng, int):
            rng = np.random.default_rng(seed=rng)
        self._rng = rng

        if not isinstance(grid, NetworkModelGrid):
            raise TypeError("grid must be a NetworkModelGrid")

        if np.min(sed_thickness) < 0.05:
            warnings.warn(
                f"sed_thickness is unrealistically low ({sed_thickness} * d84)",
                stacklevel=2,
            )

        if np.max(np.abs(tau_c_50 - 0.055)) > 0.35:
            warnings.warn(f"tau_c_50 is unrealistic ({tau_c_50})", stacklevel=2)

        self._time = [time]
        self._grid = grid
        self._tau_c_50 = tau_c_50
        self._rho_sediment = rho_sediment
        self._rho_water = rho_water
        self._gravity = gravity
        self._D84_D50 = D84_D50
        self._abrasion_rate = abrasion_rate
        self._extra_parcel_attributes = extra_parcel_attributes
        self._sed_thickness = sed_thickness
        self._median_number_of_starting_parcels = median_number_of_starting_parcels

    def __call__(self):
        d50 = self.calc_d50()

        d84 = d50 * self._D84_D50

        total_parcel_volume_at_link = calc_total_parcel_volume(
            self._grid.at_link["channel_width"],
            self._grid.at_link["reach_length"],
            d84 * self._sed_thickness,
        )
        max_parcel_volume = _determine_approx_parcel_volume(
            total_parcel_volume_at_link, self._median_number_of_starting_parcels
        )

        variables, items = _parcel_characteristics(
            total_parcel_volume_at_link,
            max_parcel_volume,
            d50,
            self._D84_D50,
            self._rho_sediment,
            self._abrasion_rate,
            self._extra_parcel_attributes,
            rng=self._rng,
        )

        if np.max(d50) > 0.5:
            warnings.warn(
                f"calculated d50 is unrealistically large ({d50} m)", stacklevel=2
            )

        if np.min(d50) < 0.002:
            warnings.warn(
                f"calculated d50 is unrealistically low ({d50} m). The equations used "
                "in this initializer are intended for gravel bedded rivers.",
                stacklevel=2,
            )

        if max_parcel_volume < 0.05:
            warnings.warn(
                f"default parcel volume is extremely small ({max_parcel_volume} m)",
                stacklevel=2,
            )

        return DataRecord(
            self._grid,
            items=items,
            time=self._time,
            data_vars=variables,
            dummy_elements={"link": [_OUT_OF_NETWORK]},
        )

    def calc_d50(self):
        raise NotImplementedError("calc_d50")


class BedParcelInitializerDischarge(BedParcelInitializerBase):
    """Create a landlab :class:`~.DataRecord` to represent parcels of sediment on
    a river network.

    The function takes discharge data for each link as input, as well as channel
    geometry (``channel_width``, ``reach_length``, ``channel_slope``) fields attached
    to the :class:`~.NetworkModelGrid`.

    This function currently estimates median parcel grain size at a link
    according to Snyder et al. (2013), assuming a lognormal parcel grain size
    distribution.

    .. codeauthor:: Eric Hutton, Allison Pfeiffer, Muneer Ahammad

    Parameters
    ----------
    grid : NetworkModelGrid
        *landlab* :class:`~.NetworkModelGrid` to place sediment parcels on.
    time : float, optional
        The initial time to add to the record.
    discharge_at_link: float
        Dominant/formative discharge at each link in the network [m^3 / s].
    mannings_n : float, optional
        Manning's *n* value for all links, used to calculate median parcel grain
        size at a link.
    tau_c_50 : float, optional
        Critical Shields stress for *d50* at dominant discharge for all links, used to
        calculate median parcel grain size.
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    rho_water : float, optional
        Density of water [kg / m^3].
    gravity : float, optional
        Acceleration due to gravity [m / s^2].
    D84_D50 : float, optional
        Ratio of *D84:D50,* used to set lognormal distribution of grain size.
    sed_thickness : float, optional
        Sediment thickness in multiples of *d84*.
    abrasion_rate : float, optional
        Abrasion rate of parcels during transport [1/m].
    median_number_of_starting_parcels : int, optional
        Median number of parcels in a link.
    extra_parcel_attributes : str or list of str, optional
        Name of user-defined parcel attribute to be added to parcel data record,
        which will be returned as an empty parcel attribute.

    Examples
    --------
    >>> from landlab import NetworkModelGrid
    >>> from landlab.components.network_sediment_transporter import (
    ...     BedParcelInitializerDischarge,
    ... )

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    >>> _ = grid.add_full("channel_width", 1.0, at="link")  # m
    >>> _ = grid.add_full("channel_slope", 0.01, at="link")  # m / m
    >>> _ = grid.add_full("reach_length", 100.0, at="link")  # m

    >>> discharge = np.full(grid.number_of_links, 10.0)  # m^3 / s
    >>> initialize_parcels = BedParcelInitializerDischarge(
    ...     grid, discharge_at_link=discharge
    ... )
    >>> parcels = initialize_parcels()
    """

    _name = "BedParcelInitializerDischarge"

    _unit_agnostic = False

    __version__ = "1.0"

    _info = {
        "discharge_at_link": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m^3 / s",
            "mapping": "link",
            "doc": "Dominant/formative discharge at each link in the network",
        },
    }

    def __init__(
        self, grid, time=0.0, discharge_at_link=None, mannings_n=0.035, **kwds
    ):
        if np.max(np.abs(mannings_n - 0.035)) > 0.3:
            warnings.warn(
                f"Manning's n value is unrealistic ({mannings_n})", stacklevel=2
            )

        if discharge_at_link is None:
            raise ValueError("User must provide discharge_at_link")

        if np.size(discharge_at_link) != grid.number_of_links:
            raise ValueError(
                "discharge_at_link should be size number_of_links "
                f"({np.size(discharge_at_link)} != {grid.number_of_links})"
            )

        self._discharge = discharge_at_link
        self._mannings_n = mannings_n

        BedParcelInitializerBase.__init__(self, grid, time=time, **kwds)

    def calc_d50(self):
        return calc_d50_discharge(
            self._grid.at_link["channel_width"],
            self._grid.at_link["channel_slope"],
            discharge=self._discharge,
            mannings_n=self._mannings_n,
            gravity=self._gravity,
            rho_water=self._rho_water,
            rho_sediment=self._rho_sediment,
            tau_c_50=self._tau_c_50,
        )


class BedParcelInitializerDepth(BedParcelInitializerBase):
    """Create a *landlab* :class:`~.DataRecord` to represent parcels of sediment on
    a river network.

    The function takes dominant flow depth for each link as input, as well as channel
    geometry (*channel_width*, *reach_length*, *channel_slope*) fields attached to the
    :class:`~.NetworkModelGrid`.

    This function currently estimates median parcel grain size at a link
    using the formative Shields stress (as in Pfeiffer et al., 2017), assuming
    a lognormal parcel grain size distribution.

    .. codeauthor:: Eric Hutton, Allison Pfeiffer, Muneer Ahammad

    Parameters
    ----------
    grid : ModelGrid
        *landlab* :class:`~.ModelGrid` to place sediment parcels on.
    time : float, optional
        The initial time to add to the record.
    flow_depth_at_link: float, optional
        Dominant/formative flow depth at each link in the network.
    tau_c_multiplier: float, optional
        Coefficient to relate critical and dominant/bankfull/formative Shields
        stress. Dominant/formative/bankfull Shields stress is calculated as
        ``multiplier * critical``.
    tau_c_50 : float, optional
        Critical Shields stress for *d50* at dominant discharge for all links, used to
        calculate median parcel grain size
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    rho_water : float, optional
        Density of water [kg / m^3].
    gravity : float, optional
        Acceleration due to gravity [m / s^2].
    D84_D50 : float, optional
        Ratio of *D84:D50*, used to set lognormal distribution of grain size.
    sed_thickness : float, optional
        Sediment thickness in multiples of *d84*.
    abrasion_rate : float, optional
        Abrasion rate of parcels during transport [1 / m].
    median_number_of_starting_parcels : int, optional
        Median number of parcels in a link.
    extra_parcel_attributes : str or list of str, optional
        Name of user-defined parcel attribute to be added to parcel data record,
        which will be returned as an empty parcel attribute.

    Examples
    --------
    >>> from landlab import NetworkModelGrid
    >>> from landlab.components.network_sediment_transporter import (
    ...     BedParcelInitializerDepth,
    ... )

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> _ = grid.add_full("channel_width", 1.0, at="link")  # m
    >>> _ = grid.add_full("channel_slope", 0.01, at="link")  # m / m
    >>> _ = grid.add_full("reach_length", 100.0, at="link")  # m


    >>> depth = np.full(grid.number_of_links, 1.0)  # m
    >>> initialize_parcels = BedParcelInitializerDepth(grid, flow_depth_at_link=depth)
    >>> parcels = initialize_parcels()
    """

    _name = "BedParcelInitializerDepth"

    _unit_agnostic = False

    __version__ = "1.0"

    _info = {
        "flow_depth_at_link": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Dominant/formative flow depth at each link in the network",
        },
    }

    def __init__(
        self, grid, time=0.0, flow_depth_at_link=None, tau_c_multiplier=1.0, **kwds
    ):
        if flow_depth_at_link is None:
            raise ValueError("User must provide flow_depth_at_link")

        if np.size(flow_depth_at_link) != grid.number_of_links:
            raise ValueError("flow_depth_at_link should be size number_of_links")

        self._flow_depth = flow_depth_at_link
        self._tau_c_multiplier = tau_c_multiplier

        BedParcelInitializerBase.__init__(self, grid, time=time, **kwds)

    def calc_d50(self):
        return calc_d50_depth(
            self._grid.at_link["channel_slope"],
            flow_depth=self._flow_depth,
            tau_c_multiplier=self._tau_c_multiplier,
            rho_water=self._rho_water,
            rho_sediment=self._rho_sediment,
            tau_c_50=self._tau_c_50,
        )


class BedParcelInitializerArea(BedParcelInitializerBase):
    """Create a *landlab* :class:`~.DataRecord` to represent parcels of sediment on
    a river network.

    The function takes a coefficient and exponent in a grain size-drainage area power
    law scaling relationship, as well as channel attribute (`drainage_area`,
    *channel_width*, *reach_length*, *channel_slope*) fields attached to the
    :class:`~.NetworkModelGrid`.

    This function currently estimates median parcel grain size at a link
    using a power-law scaling relationship between drainage area and median
    grain size (:math:`d_{50} = c A^n`), assuming a lognormal parcel grain size
    distribution.

    .. codeauthor:: Eric Hutton, Allison Pfeiffer, Muneer Ahammad

    Parameters
    ----------
    grid : ModelGrid
        *landlab* :class:`~.ModelGrid` to place sediment parcels on.
    time : float, optional
        The initial time to add to the record.
    drainage_area_coefficient : float, optional
        Coefficient in a power law grain size-drainage area scaling relationship.
    drainage_area_exponent : float, optional
        Exponent in a power law grain size-drainage area scaling relationship.
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    rho_water : float, optional
        Density of water [kg / m^3].
    gravity : float, optional
        Acceleration due to gravity [m / s^2].
    D84_D50 : float, optional
        Ratio of *D84:D50*, used to set lognormal distribution of grain size.
    sed_thickness : float, optional
        Sediment thickness in multiples of *d84*.
    abrasion_rate : float, optional
        Abrasion rate of parcels during transport [1 / m].
    median_number_of_starting_parcels : int, optional
        Median number of parcels in a link.
    extra_parcel_attributes : str or list of str, optional
        Name of user-defined parcel attribute to be added to parcel data record,
        which will be returned as an empty parcel attribute.

    Examples
    --------
    >>> from landlab import NetworkModelGrid
    >>> from landlab.components.network_sediment_transporter import (
    ...     BedParcelInitializerArea,
    ... )

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    >>> _ = grid.add_full("channel_width", 1.0, at="link")  # m
    >>> _ = grid.add_full("channel_slope", 0.01, at="link")  # m / m
    >>> _ = grid.add_full("reach_length", 100.0, at="link")  # m
    >>> _ = grid.add_full("drainage_area", 100.0, at="link")


    >>> initialize_parcels = BedParcelInitializerArea(
    ...     grid, drainage_area_coefficient=0.1, drainage_area_exponent=0.3
    ... )
    >>> parcels = initialize_parcels()
    """

    _name = "BedParcelInitializerArea"

    _unit_agnostic = False

    __version__ = "1.0"

    _info = {
        "time": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "s",
            "mapping": "link",
            "doc": "The initial time to add to the record",
        },
        "drainage_area_coefficient": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "--",
            "mapping": "link",
            "doc": (
                "Coefficient in a power law grain size-drainage area scaling "
                "relationship"
            ),
        },
        "drainage_area_exponent": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "--",
            "mapping": "link",
            "doc": (
                "Exponent in a power law grain size-drainage area scaling "
                "relationship."
            ),
        },
    }

    def __init__(
        self,
        grid,
        time=0.0,
        drainage_area_coefficient=None,
        drainage_area_exponent=None,
        **kwds,
    ):
        if drainage_area_coefficient is None:
            raise ValueError("User must provide drainage_area_coefficient")

        if drainage_area_exponent is None:
            raise ValueError("User must provide drainage_area_exponent")

        self._drainage_area_coefficient = drainage_area_coefficient
        self._drainage_area_exponent = drainage_area_exponent

        BedParcelInitializerBase.__init__(self, grid, time=time, **kwds)

    def calc_d50(self):
        return calc_d50_dArea_scaling(
            self._grid.at_link["drainage_area"],
            a=self._drainage_area_coefficient,
            n=self._drainage_area_exponent,
        )


class BedParcelInitializerUserD50(BedParcelInitializerBase):
    """Create a *landlab* :class:`~.DataRecord` to represent parcels of sediment on
    a river network.

    The function takes either a scalar value or an array (of of length,
    *number_of_links*) to assign the median grain size for parcels on each link in
    the network grid.

    This function creates a lognormal grain size distribution for the parcels
    in the link.

    .. codeauthor:: Eric Hutton, Allison Pfeiffer, Muneer Ahammad

    Parameters
    ----------
    grid : ModelGrid
        landlab :class:`~.ModelGrid` to place sediment parcels on.
    time : float, optional
        The initial time to add to the record.
    user_d50 : float, optional
        Either an array of *d50* (of length *number_of_links*) or a scalar to be
        applied to all links in the network.
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    rho_water : float, optional
        Density of water [kg / m^3].
    gravity : float, optional
        Acceleration due to gravity [m / s^2].
    D84_D50 : float, optional
        Ratio of *D84:D50*, used to set lognormal distribution of grain size.
    sed_thickness : float, optional
        Sediment thickness in multiples of *d84*.
    abrasion_rate : float, optional
        Abrasion rate of parcels during transport [1 / m].
    median_number_of_starting_parcels : int, optional
        Median number of parcels in a link.
    extra_parcel_attributes : str or list of str, optional
        Name of user-defined parcel attribute to be added to parcel data record,
        which will be returned as an empty parcel attribute.

    Examples
    --------
    >>> from landlab import NetworkModelGrid
    >>> from landlab.components.network_sediment_transporter import (
    ...     BedParcelInitializerUserD50,
    ... )

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    >>> _ = grid.add_full("channel_width", 1.0, at="link")  # m
    >>> _ = grid.add_full("channel_slope", 0.01, at="link")  # m / m
    >>> _ = grid.add_full("reach_length", 100.0, at="link")  # m


    >>> initialize_parcels = BedParcelInitializerUserD50(grid, user_d50=0.05)
    >>> parcels = initialize_parcels()
    """

    _name = "BedParcelInitializerUserD50"

    _unit_agnostic = False

    __version__ = "1.0"

    _info = {
        "time": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "s",
            "mapping": "link",
            "doc": "The initial time to add to the record",
        },
        "user_d50": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Median grain size of the bed sediment in each link",
        },
    }

    def __init__(self, grid, time=0.0, user_d50=None, **kwds):
        if user_d50 is None:
            raise ValueError("User must provide user_d50")

        self._user_d50 = user_d50

        BedParcelInitializerBase.__init__(self, grid, time=time, **kwds)

    def calc_d50(self):
        if np.size(self._user_d50) == 1:  # one d50, all links
            d50 = np.full_like(self._grid.length_of_link, self._user_d50, dtype=float)

            return d50

        elif np.size(self._user_d50) == (
            self._grid.number_of_links
        ):  # different d50 each link
            d50 = self._user_d50

            return d50
        else:
            raise ValueError(
                "user defined D50 must be either a scalar or size(element_id)"
            )


# BedParcelInitializerBase helper functions
def _parcel_characteristics(
    total_parcel_volume_at_link,
    max_parcel_volume,
    d50,
    D84_D50,
    rho_sediment,
    abrasion_rate,
    extra_parcel_attributes,
    rng=None,
):
    if rng is None:
        rng = np.random.default_rng()

    n_parcels_at_link = np.ceil(total_parcel_volume_at_link / max_parcel_volume).astype(
        dtype=int
    )
    if np.min(n_parcels_at_link) < 10:
        warnings.warn(
            f"at least one link has only {n_parcels_at_link} parcels.", stacklevel=2
        )

    element_id = np.empty(np.sum(n_parcels_at_link), dtype=int)

    # volume = np.full(np.sum(n_parcels_at_link), max_parcel_volume, dtype=float)
    volume = np.full_like(element_id, max_parcel_volume, dtype=float)
    grain_size = np.empty_like(element_id, dtype=float)
    offset = 0
    for link, n_parcels in enumerate(n_parcels_at_link):
        element_id[offset : offset + n_parcels] = link
        grain_size[offset : offset + n_parcels] = rng.lognormal(
            np.log(d50[link]), np.log(D84_D50), n_parcels
        )
        volume[offset] = total_parcel_volume_at_link[link] - (
            (n_parcels - 1) * max_parcel_volume
        )  # small remaining volume

        offset += n_parcels
    starting_link = element_id.copy()
    abrasion_rate = np.full_like(element_id, abrasion_rate, dtype=float)
    density = np.full_like(element_id, rho_sediment, dtype=float)

    element_id = np.expand_dims(element_id, axis=1)
    volume = np.expand_dims(volume, axis=1)
    grain_size = np.expand_dims(grain_size, axis=1)

    time_arrival_in_link = np.expand_dims(
        rng.uniform(size=np.sum(n_parcels_at_link)), axis=1
    )
    location_in_link = np.expand_dims(
        rng.uniform(size=np.sum(n_parcels_at_link)), axis=1
    )

    active_layer = np.ones(np.shape(element_id))
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
            variables[attrib] = (["item_id"], np.nan * np.zeros_like(density))

    return variables, {"grid_element": "link", "element_id": element_id}


def _determine_approx_parcel_volume(
    total_parcel_volume_at_link, median_number_of_starting_parcels
):
    """What size parcels should we use?"""
    median_link_volume = np.median(total_parcel_volume_at_link)
    return median_link_volume / median_number_of_starting_parcels


def calc_total_parcel_volume(width, length, sediment_thickness):
    """Simple rectangular prism geometry. Total parcel vol in each link = L*W*H"""

    return width * length * sediment_thickness


# calc_d50 helper functions


def calc_d50_discharge(
    width,
    slope,
    discharge,
    mannings_n,
    gravity,
    rho_water,
    rho_sediment,
    tau_c_50,
):
    """Calculate median grain size via dominant discharge and channel width
    according to Snyder et al. (2013)

    Returns
    -------
    ndarray of float
        d50.

    Examples
    --------
    >>> from numpy.testing import assert_almost_equal

    >>> w = 20
    >>> S = 0.01
    >>> Q = 100
    >>> n = 0.05
    >>> g = 9.81
    >>> rho_w = 1000
    >>> rho_s = 3000
    >>> tau_c_50 = 0.05

    >>> expected_value = (
    ...     rho_w * g * n ** (3 / 5) * Q ** (3 / 5) * w ** (-3 / 5) * S ** (7 / 10)
    ... ) / ((rho_s - rho_w) * g * tau_c_50)
    >>> print(np.round(expected_value, decimals=3))
    0.173
    >>> assert_almost_equal(
    ...     calc_d50_discharge(20, 0.01, 100, 0.05, 9.81, 1000, 3000, 0.05),
    ...     expected_value,
    ... )
    """

    return (
        rho_water
        * gravity
        * mannings_n ** (3 / 5)
        * discharge ** (3 / 5)
        * width ** (-3 / 5)
        * slope ** (7 / 10)
    ) / ((rho_sediment - rho_water) * gravity * tau_c_50)


def calc_d50_depth(
    slope,
    flow_depth,
    tau_c_multiplier,
    rho_water,
    rho_sediment,
    tau_c_50,
):
    """Calculate median grain size via dominant flow depth according to
    Pfeiffer et al. (2017).

    Returns
    -------
    ndarray of float
        *d50*.

    Examples
    --------
    >>> from numpy.testing import assert_almost_equal

    >>> slope = 0.01
    >>> depth = 1
    >>> tau_c_multiplier = 1
    >>> rho_w = 1000
    >>> rho_s = 3000
    >>> tau_c_50 = 0.05

    >>> expected_value = (rho_w * depth * slope) / (
    ...     (rho_s - rho_w) * tau_c_50 * tau_c_multiplier
    ... )
    >>> print(np.round(expected_value, decimals=3))
    0.1
    >>> assert_almost_equal(
    ...     calc_d50_depth(0.01, 1, 1, 1000, 3000, 0.05), expected_value
    ... )
    """

    return (rho_water * flow_depth * slope) / (
        (rho_sediment - rho_water) * tau_c_multiplier * tau_c_50
    )


def calc_d50_dArea_scaling(drainage_area, a, n):
    """Calculate median grain size via power law scaling relationship with drainage
    area.

    Returns
    -------
    ndarray of float
        *d50*.

    Examples
    --------
    >>> from numpy.testing import assert_almost_equal

    >>> drainage_area = 10
    >>> drainage_area_coefficient = 1
    >>> drainage_area_exponent = -0.1

    >>> expected_value = drainage_area_coefficient * (
    ...     drainage_area**drainage_area_exponent
    ... )
    >>> print(np.round(expected_value, decimals=3))
    0.794
    >>> assert_almost_equal(calc_d50_dArea_scaling(10, 1, -0.1), expected_value)
    """
    d50 = a * drainage_area**n

    return d50
