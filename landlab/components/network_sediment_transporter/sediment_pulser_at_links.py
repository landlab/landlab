import numpy as np

from landlab.data_record import DataRecord

from .sediment_pulser_base import SedimentPulserBase

_OUT_OF_NETWORK = -2


class SedimentPulserAtLinks(SedimentPulserBase):
    """Send a pulse of parcels to specific links in a channel network

    :class:`~.SedimentPulserAtLinks` is instantiated by specifying the
    :class:`~.NetworkModelGrid` it will pulse the parcels into and the time(s) when
    a pulse is allowed to occur.  It inherits attributes and functions from the
    :class:`~.SedimentPulserBase`.

    :class:`~.SedimentPulserAtLinks` is run (adds parcels to ``DataRecord``) by
    calling the instance with a list of links and a list of the number of parcels
    added to each link.

    If parcel attributes are constant with time and uniform
    across the basin, these constant-uniform-attributes can be defined
    when :class:`~.SedimentPulserAtLinks` is instantiated. If parcel attributes vary
    with location and time, the user specifies the varying parcel attributes
    each time the instance is called with a list for each attribute of length
    equal to the number of links included in the pulse.


    .. codeauthor:: Jeff Keck, Allison Pfeiffer, Shelby Ahrendt
                    (with help from Eric Hutton and Katy Barnhart)


    Examples
    --------
    >>> import numpy as np
    >>> from landlab import NetworkModelGrid

    Create the network model grid the parcels will be added to.

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> grid.at_link["channel_width"] = np.full(grid.number_of_links, 1.0)  # m
    >>> grid.at_link["channel_slope"] = np.full(grid.number_of_links, 0.01)  # m / m
    >>> grid.at_link["reach_length"] = np.full(grid.number_of_links, 100.0)  # m

    Define a function that contains which times a pulse is allowed to occur.
    This function says a pulse can occur at any time

    >>> def time_to_pulse(time):
    ...     return True
    ...

    Instantiate :class:`~.SedimentPulserAtLinks`

    >>> make_pulse = SedimentPulserAtLinks(grid, time_to_pulse=time_to_pulse)

    Run the instance with inputs for the time, link location and number of
    parcels. Other attributes will use the default values in the base class

    >>> time = 11
    >>> links = [2, 6]
    >>> n_parcels_at_link = [2, 3]
    >>> parcels = make_pulse(
    ...     time=time, links=links, n_parcels_at_link=n_parcels_at_link
    ... )

    Check the element_id of each parcel

    >>> print(parcels.dataset["element_id"].values)
    [[2]
     [2]
     [6]
     [6]
     [6]]
    """

    _name = "SedimentPulserAtLinks"

    _unit_agnostic = False

    _info = {}  # works with the DataRecord

    def __init__(
        self,
        grid,
        time_to_pulse=None,
        parcels=None,
        D50=0.05,
        D84_D50=2.1,
        rho_sediment=2650.0,
        parcel_volume=0.5,
        abrasion_rate=0.0,
        rng=None,
    ):
        """Create :class:`~.SedimentPulserAtLinks`.

        Parameters
        ----------

        grid : ModelGrid
            landlab :class:`~.ModelGrid` to place sediment parcels on.
        time_to_pulse: function, optional
            The condition when a pulse occurs using the ``_pulse_characteristics``
            method. If not specified, a pulse occurs whenever the instance is run
        parcels: landlab DataRecord
            Tracks parcel location and variables.
        D50: float, optional
            Median grain size [m].
        D84_D50: float, optional
            Ratio of 84th percentile grain size to the median grain size.
        rho_sediment : float, optional
            Sediment grain density [kg / m^3].
        parcel_volume : float
            Parcel volume [m^3]
        abrasion_rate: float
            Volumetric abrasion exponent [1 / m]
        """
        if rng is None:
            rng = np.random.default_rng()
        elif isinstance(rng, int):
            rng = np.random.default_rng(seed=rng)
        self._rng = rng

        SedimentPulserBase.__init__(
            self,
            grid,
            parcels=parcels,
            D50=D50,
            D84_D50=D84_D50,
            rho_sediment=rho_sediment,
            parcel_volume=parcel_volume,
            abrasion_rate=abrasion_rate,
        )

        # set time_to_pulse to True if not specified
        if time_to_pulse is None:
            self._time_to_pulse = lambda time: True
        else:
            self._time_to_pulse = time_to_pulse

    def __call__(
        self,
        time,
        links=None,
        n_parcels_at_link=None,
        D50=None,
        D84_D50=None,
        rho_sediment=None,
        parcel_volume=None,
        abrasion_rate=None,
    ):
        """
        specify the time, link(s) and attributes of pulses added to a
        :class:`~.NetworkModelGrid` at stochastically determined locations within the
        link(s).

        Parameters
        ----------
        time : integer or datetime64
            Time that the pulse is occurs.
        links : list of int
            Link ID # that parcels are added to.
        n_parcels_at_link: list of int
            Number of parcels added to each link listed in links
        D50 : list of float, optional
            Median grain size of parcels added to each link listed in links [m].
        D84_D50 : list of float, optional
            Ratio of 84th percentile grain size to the median grain size.
        rho_sediment: list of float, optional
            Density of grains [kg / m^3].
        parcel_volume : list of float, optional
            Volume of each parcel added to link listed in links [m^3].
        abrasion_rate: list of float, optional
            Rate that grain size decreases with distance along channel [mm / km?].

        Returns
        -------
        parcels
            :class:`~.DataRecord` containing all information on each individual parcel.

        """

        # check user provided links and number of parcels sent to each link
        assert (
            links is not None and n_parcels_at_link is not None
        ), "must provide links and number of parcels entered into each link"

        links = np.array(links)
        n_parcels_at_link = np.array(n_parcels_at_link)

        # any parameters not specified with __Call__ method use default values
        # specified in the base class
        if D50 is None:
            D50 = np.full_like(links, self._D50, dtype=float)
        else:
            D50 = np.array(D50)

        if D84_D50 is None:
            D84_D50 = np.full_like(links, self._D84_D50, dtype=float)
        else:
            D84_D50 = np.array(D84_D50)

        if rho_sediment is None:
            rho_sediment = np.full_like(links, self._rho_sediment, dtype=float)
        else:
            rho_sediment = np.array(rho_sediment)

        if parcel_volume is None:
            parcel_volume = np.full_like(links, self._parcel_volume, dtype=float)
        else:
            parcel_volume = np.array(parcel_volume)

        if abrasion_rate is None:
            abrasion_rate = np.full_like(links, self._abrasion_rate, dtype=float)
        else:
            abrasion_rate = np.array(abrasion_rate)

        # before running, check that no inputs < 0
        # check for negative inputs
        if (
            np.array([D50, D84_D50, rho_sediment, parcel_volume, abrasion_rate]) < 0
        ).any():
            raise AssertionError("parcel attributes cannot be less than zero")
        # before running, check if time to pulse
        if not self._time_to_pulse(time):
            # if not time to pulse, return the existing parcels
            print("user provided time not a time-to-pulse, parcels have not changed")

            return self._parcels

        # create items and variables for DataRecord
        variables, items = self._sediment_pulse_stochastic(
            time,
            links,
            n_parcels_at_link,
            parcel_volume,
            D50,
            D84_D50,
            abrasion_rate,
            rho_sediment,
        )

        # if DataRecord does not exist, create one
        if self._parcels is None:
            self._parcels = DataRecord(
                self._grid,
                items=items,
                time=[time],
                data_vars=variables,
                dummy_elements={"link": [_OUT_OF_NETWORK]},
            )
        # else, add parcels to existing DataRecrod
        else:
            self._parcels.add_item(time=[time], new_item=items, new_item_spec=variables)

        return self._parcels

    def _sediment_pulse_stochastic(
        self,
        time,
        links,
        n_parcels_at_link,
        parcel_volume,
        D50,
        D84_D50,
        abrasion_rate,
        rho_sediment,
    ):
        """Convert lists of link ids and link parcel parameters to a dataset
        that describes the network location and attributes of each individual parcel

        Returns
        -------
        dict
            Dictionary with keys and data in format for :class:`~.DataRecord`.

        """

        # Create np array for each paracel attribute. Length of array is equal
        # to the number of parcels

        # link id, D50 and volume
        element_id = np.empty(np.sum(n_parcels_at_link), dtype=int)
        grain_size = np.empty(np.sum(n_parcels_at_link))
        volume = np.empty(np.sum(n_parcels_at_link))
        offset = 0
        for link, n_parcels in enumerate(n_parcels_at_link):
            element_id[offset : offset + n_parcels] = links[link]
            grain_size[offset : offset + n_parcels] = self._rng.lognormal(
                np.log(D50[link]), np.log(D84_D50[link]), n_parcels
            )
            volume[offset : offset + n_parcels] = parcel_volume[link] % n_parcels
            offset += n_parcels
        starting_link = element_id.copy()

        # abrasion rate and density
        abrasion_rate_L = []
        density_L = []
        for c, ei in enumerate(np.unique(element_id)):
            element_id_subset = element_id[element_id == ei]
            abrasion_rate_L = abrasion_rate_L + list(
                np.full_like(element_id_subset, abrasion_rate[c], dtype=float)
            )
            density_L = density_L + list(
                np.full_like(element_id_subset, rho_sediment[c], dtype=float)
            )
        abrasion_rate = np.array(abrasion_rate_L)
        density = np.array(density_L)

        element_id = np.expand_dims(element_id, axis=1)
        grain_size = np.expand_dims(grain_size, axis=1)
        volume = np.expand_dims(volume, axis=1)

        # time of arrivial (time instance called)
        time_arrival_in_link = np.full(np.shape(element_id), time, dtype=float)

        # link location (distance from link inlet / link length) is stochastically
        # determined
        location_in_link = np.expand_dims(
            self._rng.uniform(size=np.sum(n_parcels_at_link)), axis=1
        )

        # All parcels in pulse are in the active layer (1) rather than subsurface (0)
        active_layer = np.ones(np.shape(element_id))

        # specify that parcels are in the links of the network model grid
        grid_element = ["link"] * np.size(element_id)
        grid_element = np.expand_dims(grid_element, axis=1)

        return {
            "starting_link": (["item_id"], starting_link),
            "abrasion_rate": (["item_id"], abrasion_rate),
            "density": (["item_id"], density),
            "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
            "active_layer": (["item_id", "time"], active_layer),
            "location_in_link": (["item_id", "time"], location_in_link),
            "D": (["item_id", "time"], grain_size),
            "volume": (["item_id", "time"], volume),
        }, {"grid_element": grid_element, "element_id": element_id}
