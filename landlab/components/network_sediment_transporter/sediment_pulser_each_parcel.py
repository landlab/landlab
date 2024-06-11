import warnings

import numpy as np

from landlab.components.network_sediment_transporter.sediment_pulser_base import (
    SedimentPulserBase,
)
from landlab.data_record import DataRecord

_OUT_OF_NETWORK = -2


class SedimentPulserEachParcel(SedimentPulserBase):
    """Send pulses of sediment to specific point locations within the channel
    network and divide the pulses into parcels. Pulses may be any volume.
    Parcels must be less than or equal to a user specified maximum volume.

    SedimentPulserEachParcel is instantiated by specifying the network model grid
    it will pulse the parcels into

    SedimentPulserEachParcel is run (adds parcels to DataRecrod) by calling the
    SedimentPulserEachParcel instance with the time that pulses are added to
    the channel network and a sediment pulse table (PulseDF)

    PulseDF is a pandas dataframe. At a minimum, the dataframe must have columns 'Link#'
    'normalized_downstream_distance' and 'pulse_volume'. Optionally, the parcel
    volume that the pulse is divided into and grain characteristics of each pulse
    can also be specified in PulseDF.


    .. codeauthor: Jeff Keck, Allison Pfeiffer, Shelby Ahrendt
                   (with help from Eric Hutton and Katy Barnhart)


    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> from landlab import NetworkModelGrid

    Create the network model grid. Pulses are added to the links of the network
    model grid.

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> grid.at_link["channel_width"] = np.full(grid.number_of_links, 1.0)  # m
    >>> grid.at_link["channel_slope"] = np.full(grid.number_of_links, 0.01)  # m / m
    >>> grid.at_link["reach_length"] = np.full(grid.number_of_links, 100.0)  # m


    Instantiate 'SedimentPulserEachParcel'

    >>> make_pulse = SedimentPulserEachParcel(grid)

    Define the PulseDF and time of the pulse

    >>> PulseDF = pd.DataFrame(
    ...     {
    ...         "pulse_volume": [0.2, 1, 1.1, 0.5],
    ...         "link_#": [1, 3, 5, 2],
    ...         "normalized_downstream_distance": [0.8, 0.7, 0.5, 0.2],
    ...     }
    ... )
    >>> time = 7

    Run the instance

    >>> parcels = make_pulse(time, PulseDF)

    This should yield a UserWarning: Parcels not provided, created a new DataRecord

    check element_id of each parcel

    >>> print(parcels.dataset["element_id"].values)
    [[1]
    [3]
    [3]
    [5]
    [5]
    [5]
    [2]]

    """

    _name = "SedimentPulserEachParcel"

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
        rng=None,
    ):
        """
        instantiate SedimentPulserEachParcel

        Parameters
        ----------
        grid : ModelGrid
            landlab *ModelGrid* to place sediment parcels on.
        parcels: landlab DataRecord, optional
            Tracks parcel location and attributes
        D50: float, optional
            median grain size [m]
        D84_D50: float, optional
            ratio of 84th percentile grain size to the median grain size
        rho_sediment : float, optional
            Sediment grain density [kg/m^3].
        parcel_volume : float, optional
            parcel volume [m^3]
        abrasion_rate: float, optional
            volumetric abrasion exponent [1/m]
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

    def __call__(self, time, PulseDF=None):
        """specify the location and attributes of each pulse of material added to
        a Network Model Grid DataRecord

        Parameters
        ----------
        time : integer or datetime64 value equal to nst.time
            time that the pulse is triggered in the network sediment transporter
        PulseDF : pandas dataframe
            each row contains information on the deposition location and volume of
            a single pulse of sediment. The pulse is divided into 'n' number of
            parcels, where 'n' equals the np.ceil(pulse volume / parcel volume)
            For details on the format of the DataFrame, see the docstring for
            function _sediment_pulse_dataframe


        Returns
        -------
            self._parcels
                a DataRecord containing all information on each individual parcel
        """
        # If no PulseDF provided, raise error. Should at least provide an empty PulseDF
        if PulseDF is None:
            raise ValueError("PulseDF was not specified")

        if (
            PulseDF.empty is True
        ):  # if empty, pulser stops, returns the existing parcels, call stops
            warnings.warn("Pulse DataFrame is EMPTY", stacklevel=2)
            return self._parcels

        variables, items = self._sediment_pulse_dataframe(
            time,  # create variabels and and items needed to create the data record
            PulseDF,
        )

        if self._parcels is None:  # if no parcels, create parcels
            self._parcels = DataRecord(
                self._grid,
                items=items,
                time=[time],
                data_vars=variables,
                dummy_elements={"link": [_OUT_OF_NETWORK]},
            )

            warnings.warn(
                "Parcels not provided, created a new DataRecord", stacklevel=2
            )

        else:  # else use the add item method to add parcels
            self._parcels.add_item(time=[time], new_item=items, new_item_spec=variables)

        return self._parcels

    def _sediment_pulse_dataframe(self, time, PulseDF):
        """Convert PulseDF to a :class:`~.DataRecord` formatted for the
        :class:`~.NetworkSedimentTransporter`.

        Parameters
        ----------
        time : integer or datetime64

        PulseDF : pandas dataframe

            The PulseDF must include the following columns:
                'link_#', 'pulse_volume', 'normalized_downstream_distance'

            Optionally, PulseDF can include the following columns:
               'D50', 'D84_D50', 'abrasion_rate', 'rho_sediment', 'parcel_volume'

            Values in each columne are defined as follows:

            'link_#': int - link number pulse enters the channel network
            'pulse_volume: float - total volume of the pulse [m^3]
            'normalized_downstream_distance': float - distance from link inlet
                                                      divided by link length
            'D50': float - median grain size [m]
            'D84_D50': float - grain-size standard deviation [m]
            'abrasion_rate': float - rate that grain size decreases with
                                     distance along channel [mm/km?]
            'rho_sediment': float - density grains [kg/m^3]
            'parcel_volume': float - maximum volume of one parcel [m^3]


            if the optional columns are not included in PulseDF, those parameters
            are assumed uniform across the basin, constant with time and equal
            to the corrisponding class variables.

        Returns
        -------
        tuple: (variables, items)
            variables: dictionary, attribute values for all new parcels
            item_id: dictionary, model grid element and index of element of each parcel

        """
        # split pulse into parcels.
        p_np = []  # list of number of parcels in each pulse
        volume = np.array([])  # list of parcel volumes from all pulses
        for _index, row in PulseDF.iterrows():
            # set the maximum allowable parcel volume using either
            # the default value or value in PulseDF
            if "parcel_volume" in PulseDF:
                mpv = row["parcel_volume"]
            else:
                mpv = self._parcel_volume

            # split the pulse into parcels
            if row["pulse_volume"] < mpv:
                # only one partial parcel volume
                v_p = np.array([row["pulse_volume"]])
            else:
                # number of whole parcels
                n_wp = int(np.floor(row["pulse_volume"] / mpv))
                # array of volumes, whole parcels
                v_wp = np.ones(n_wp) * mpv
                # volume of last parcel, a partial parcel
                v_pp = np.array([row["pulse_volume"] % mpv])
                # array of all parcel volumes
                # partial parcel included if volume > 0
                if v_pp > 0:
                    v_p = np.concatenate((v_wp, v_pp))
                else:
                    v_p = v_wp
            volume = np.concatenate((volume, v_p))
            p_np.append(len(v_p))
        volume = np.expand_dims(volume, axis=1)

        # link location
        link_distance_ratio = np.array([])
        for i, val in enumerate(PulseDF["normalized_downstream_distance"].values):
            # parcels from the same pulse enter channel at the same point
            link_distance_ratio = np.concatenate(
                (link_distance_ratio, np.ones(p_np[i]) * val)
            )
        location_in_link = np.expand_dims(link_distance_ratio, axis=1)

        # element id and starting link
        element_id = np.array([])
        for i, row in PulseDF.iterrows():
            element_id = np.concatenate((element_id, np.ones(p_np[i]) * row["link_#"]))
        starting_link = element_id.copy()
        element_id = np.expand_dims(element_id.astype(int), axis=1)

        # specify that parcels are in the links of the network model grid
        grid_element = ["link"] * np.size(element_id)
        grid_element = np.expand_dims(grid_element, axis=1)

        # time of arrivial (time instance called)
        time_arrival_in_link = np.full(np.shape(element_id), time, dtype=float)

        # All parcels in pulse are in the active layer (1) rather than subsurface (0)
        active_layer = np.ones(np.shape(element_id))

        if "rho_sediment" in PulseDF.columns:
            density = np.array([])
            for i, row in PulseDF.iterrows():
                density = np.concatenate(
                    (density, np.ones(p_np[i]) * row["rho_sediment"])
                )
        else:
            density = self._rho_sediment * np.ones(np.shape(starting_link))

        if "abrasion_rate" in PulseDF.columns:
            abrasion_rate = np.array([])
            for i, row in PulseDF.iterrows():
                abrasion_rate = np.concatenate(
                    (abrasion_rate, np.ones(p_np[i]) * row["abrasion_rate"])
                )
        else:
            abrasion_rate = self._abrasion_rate * np.ones(np.shape(starting_link))

        if "D50" in PulseDF.columns and "D84_D50" in PulseDF.columns:
            grain_size = np.array([])
            for i, row in PulseDF.iterrows():
                # det D50 and D84_D50
                n_parcels = p_np[i]
                D50 = row["D50"]
                D84_D50 = row["D84_D50"]
                grain_size_pulse = self._rng.lognormal(
                    np.log(D50), np.log(D84_D50), n_parcels
                )
                grain_size = np.concatenate((grain_size, grain_size_pulse))
        else:
            n_parcels = sum(p_np)
            D50 = self._D50
            D84_D50 = self._D84_D50
            grain_size = self._rng.lognormal(np.log(D50), np.log(D84_D50), n_parcels)

        grain_size = np.expand_dims(grain_size, axis=1)

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
