#!/usr/env/python

"""
Landlab component that simulates the transport of bed material
sediment through a 1-D river network, while tracking the resulting changes
in bed material grain size and river bed elevation. Model framework
described in Czuba (2018). Additions include: particle abrasion, variable
active layer thickness (Wong et al., 2007).

.. codeauthor:: Allison Pfeiffer, Katy Barnhart, Jon Czuba, Eric Hutton

Created on Tu May 8, 2018
Last edit was sometime after February 2020
"""

import warnings

import numpy as np
import scipy.constants
import xarray as xr

from landlab import Component
from landlab.components import FlowDirectorSteepest
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid

_SUPPORTED_TRANSPORT_METHODS = ["WilcockCrowe"]
_SUPPORTED_ACTIVE_LAYER_METHODS = ["WongParker", "GrainSizeDependent", "Constant10cm"]

_REQUIRED_PARCEL_ATTRIBUTES = [
    "time_arrival_in_link",
    "abrasion_rate",
    "density",
    "active_layer",
    "location_in_link",
    "D",
    "volume",
]

_ACTIVE = 1
_INACTIVE = 0

_SAND_SIZE = 0.002
_INIT_ACTIVE_LAYER_THICKNESS = 0.03116362


class NetworkSedimentTransporter(Component):
    """Move sediment parcels on a river network.

    Landlab component that simulates the transport of bed material
    sediment through a 1-D river network, while tracking the resulting changes
    in bed material grain size and river bed elevation. Model framework
    described in Czuba (2018). Additions include: particle abrasion, variable
    active layer thickness (Wong et al., 2007).

    This component cares about units. Its time, length, and mass units are
    seconds, meters, and kilograms, by default. The horizontal unit of the
    grid, and the units of the parameters ``g`` and ``fluid_density`` are
    what specify the component units. In addition, the expressions used
    to calculate the transport have units (Wilcock and Crowe, 2003).

    There is a function that assists in plotting the output of this component.
    It is called
    :py:func:`~landlab.plot.network_sediment_transporter.plot_network_and_parcels`.
    Examples of its usage can be found in the NetworkSedimentTransporter
    notebooks (located in the "notebooks" folder).

    Attributes
    ----------
    OUT_OF_NETWORK : int
        Indicates a parcel is out of network.

    Examples
    ----------
    >>> import numpy as np
    >>> from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
    >>> from landlab import NetworkModelGrid
    >>> from landlab.data_record import DataRecord

    The NetworkSedimentTransporter moves "parcels" of sediment down a network
    based on a given flow and a given sediment transport formulation. The river
    network is represented by a landlab :py:class:`~landlab.grid.network.NetworkModelGrid`. Flow direction in the
    network is determined using a landlab flow director. Sediment parcels are
    represented as items within a landlab :py:class:`~landlab.data_record.data_record.DataRecord`. The landlab
    :py:class:`~landlab.data_record.data_record.DataRecord` is used to track the location, grain size, sediment density,
    and total volume of each parcel.

    Create a :py:class:`~landlab.grid.network.NetworkModelGrid` to represent the river channel network. In
    this case, the grid is a single line of 4 nodes connected by 3 links. Each
    link represents a reach of river.

    >>> y_of_node = (0, 0, 0, 0)
    >>> x_of_node = (0, 100, 200, 300)
    >>> nodes_at_link = ((0,1), (1,2), (2,3))
    >>> nmg = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

    Add required channel and topographic variables to the :py:class:`~landlab.grid.network.NetworkModelGrid`.

    >>> _ = nmg.add_field("bedrock__elevation", [3., 2., 1., 0.], at="node") # m
    >>> _ = nmg.add_field("reach_length", [100., 100., 100.], at="link")  # m
    >>> _ = nmg.add_field(
    ...     "channel_width",
    ...     (15 * np.ones(nmg.size("link"))),
    ...     at="link")
    >>> _ = nmg.add_field(
    ...     "flow_depth",
    ...     (2 * np.ones(nmg.size("link"))),
    ...     at="link") # m

    Add ``topographic__elevation`` to the grid because the
    :py:class:`~landlab.components.FlowDirectorSteepest` will look to it to determine the direction of
    sediment transport through the network. Each time we run the
    ``NetworkSedimentTransporter`` the topography will be updated based on the
    bedrock elevation and the distribution of alluvium.

    >>> _ = nmg.add_field("topographic__elevation", np.copy(nmg.at_node["bedrock__elevation"]), at="node")

    Run :py:class:`~landlab.components.FlowDirectorSteepest` to determine the direction of sediment transport through the network.

    >>> flow_director = FlowDirectorSteepest(nmg)
    >>> flow_director.run_one_step()

    Define the starting time and the number of timesteps for this model run.

    >>> timesteps = 10
    >>> time = [0.0]

    Define the sediment characteristics that will be used to create the parcels
    ``DataRecord``.

    >>> items = {"grid_element": "link",
    ...          "element_id": np.array([[0]])}

    >>> variables = {
    ...     "starting_link": (["item_id"], np.array([0])),
    ...     "abrasion_rate": (["item_id"], np.array([0])),
    ...     "density": (["item_id"], np.array([2650])),
    ...     "time_arrival_in_link": (["item_id", "time"], np.array([[0]])),
    ...     "active_layer": (["item_id", "time"], np.array([[1]])),
    ...     "location_in_link": (["item_id", "time"], np.array([[0]])),
    ...     "D": (["item_id", "time"], np.array([[0.05]])),
    ...     "volume": (["item_id", "time"], np.array([[1]])),
    ... }

    Create the sediment parcel :py:class:`~landlab.data_record.data_record.DataRecord`. In this case, we are creating a
    single sediment parcel with all of the required attributes.

    >>> one_parcel = DataRecord(
    ...     nmg,
    ...     items=items,
    ...     time=time,
    ...     data_vars=variables,
    ...     dummy_elements={
    ...         "link": [NetworkSedimentTransporter.OUT_OF_NETWORK]},
    ... )

    Instantiate the model run

    >>> nst = NetworkSedimentTransporter(
    ...         nmg,
    ...         one_parcel,
    ...         flow_director,
    ...         bed_porosity=0.03,
    ...         g=9.81,
    ...         fluid_density=1000,
    ...         transport_method="WilcockCrowe",
    ...         active_layer_method="WongParker"
    ...     )

    >>> dt = 60  # (seconds) 1 min timestep

    Run the model

    >>> for t in range(0, (timesteps * dt), dt):
    ...     nst.run_one_step(dt)

    We can the link location of the parcel at each timestep

    >>> print(one_parcel.dataset.element_id.values)
    [[ 0.  0.  0.  0.  0.  1.  1.  1.  1.  1.  2.]]

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    A JOSS submission has been prepared.

    **Additional References**

    Czuba, J. A. (2018). A Lagrangian framework for exploring complexities of mixed-size sediment transport in gravel-bedded river networks. Geomorphology, 321, 146-152.

    Wilcock, P. R., & Crowe, J. C. (2003). Surface-based transport model for mixed-size sediment. Journal of Hydraulic Engineering, 129(2), 120-128.

    Wong, M., Parker, G., DeVries, P., Brown, T. M., & Burges, S. J. (2007). Experiments on dispersion of tracer stones under lower‐regime plane‐bed equilibrium bed load transport. Water Resources Research, 43(3).
    """

    _name = "NetworkSedimentTransporter"

    _unit_agnostic = False

    __version__ = "1.0"

    _info = {
        "bedrock__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of the bedrock surface",
        },
        "channel_slope": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "Slope of the river channel through each reach",
        },
        "channel_width": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Flow width of the channel, assuming constant width",
        },
        "flow_depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Flow depth of the channel",
        },
        "reach_length": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Length of each reach",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    #: Indicates a parcel is out of network
    OUT_OF_NETWORK = NetworkModelGrid.BAD_INDEX - 1

    def __init__(
        self,
        grid,
        parcels,
        flow_director,
        bed_porosity=0.3,
        g=scipy.constants.g,
        fluid_density=1000.0,
        transport_method="WilcockCrowe",
        active_layer_method="WongParker",
        active_layer_d_multiplier=2,
    ):
        """
        Parameters
        ----------
        grid: NetworkModelGrid
            A :py:class:`~landlab.grid.network.NetworkModelGrid` in which links are stream channel
            segments.
        parcels: DataRecord
            A landlab :py:class:`~landlab.data_record.data_record.DataRecord` describing the characteristics and location of
            sediment "parcels".
            At any given timestep, each parcel is located at a specified point
            along (location_in_link) a particular link (element_id). Each
            parcel has a total sediment volume (volume), sediment grain size (D),
            sediment density (density), and bed material abrasion rate
            (abrasion_rate). During a given timestep, parcels may be in the
            "active layer" of most recently deposited sediment
            (active_layer = 1), or they may be buried and not subject to
            transport (active_layer = 0). Whether a sediment parcel is active
            or not is determined based on flow conditions and parcel attributes
            in 'run_one_step'
        flow_director: :py:class:`~landlab.components.FlowDirectorSteepest`
            A landlab flow director. Currently, must be :py:class:`~landlab.components.FlowDirectorSteepest`.
        bed_porosity: float, optional
            Proportion of void space between grains in the river channel bed.
            Default value is 0.3.
        g: float, optional
            Acceleration due to gravity. Default value is 9.81 (m/s^2)
        fluid_density: float, optional
            Density of the fluid (generally, water) in which sediment is
            moving. Default value is 1000 (kg/m^3)
        transport_method: string
            Sediment transport equation option. Default (and currently only)
            option is "WilcockCrowe".
        active_layer_method: string, optional
            Option for treating sediment active layer as a constant or variable
            (default, "WongParker")
        """
        if not isinstance(grid, NetworkModelGrid):
            msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            raise ValueError(msg)

        # run super. this will check for required inputs specified by _info
        super().__init__(grid)

        # check key information about the parcels, including that all required
        # attributes are present.
        if not isinstance(parcels, DataRecord):
            msg = (
                "NetworkSedimentTransporter: parcels must be an instance"
                "of DataRecord"
            )
            raise ValueError(msg)

        for rpa in _REQUIRED_PARCEL_ATTRIBUTES:
            if rpa not in parcels.dataset:
                msg = "NetworkSedimentTransporter: {rpa} must be assigned to the parcels".format(
                    rpa=rpa
                )
                raise ValueError(msg)

        # save key information about the parcels.
        self._parcels = parcels
        self._num_parcels = self._parcels.number_of_items
        self._time_variable_parcel_attributes = [
            "time_arrival_in_link",
            "active_layer",
            "location_in_link",
            "D",
            "volume",
        ]

        # assert that the flow director is a component and is of type
        # FlowDirectorSteepest
        if not isinstance(flow_director, FlowDirectorSteepest):
            msg = (
                "NetworkSedimentTransporter: flow_director must be "
                "FlowDirectorSteepest."
            )
            raise ValueError(msg)

        # save reference to flow director
        self._fd = flow_director

        # verify and save the bed porosity.
        if not 0 <= bed_porosity < 1:
            msg = "NetworkSedimentTransporter: bed_porosity must be" "between 0 and 1"
            raise ValueError(msg)
        self._bed_porosity = bed_porosity

        # save or create other key properties.
        self._g = g
        self._fluid_density = fluid_density
        self._time_idx = 0
        self._time = 0.0
        self._distance_traveled_cumulative = np.zeros(self._num_parcels)

        # check the transport method is valid.
        if transport_method in _SUPPORTED_TRANSPORT_METHODS:
            self._transport_method = transport_method
        else:
            msg = "NetworkSedimentTransporter: Valid transport method not supported."
            raise ValueError(msg)

        # update the update_transport_time function to be the correct function
        # for the transport method.
        if self._transport_method == "WilcockCrowe":
            self._update_transport_time = self._calc_transport_wilcock_crowe

        if active_layer_method in _SUPPORTED_ACTIVE_LAYER_METHODS:
            self._active_layer_method = active_layer_method
        else:
            msg = "NetworkSedimentTransporter: Active layer method not supported."
            raise ValueError(msg)

        if self._active_layer_method == "GrainSizeDependent":
            self._active_layer_d_multiplier = active_layer_d_multiplier

        # save reference to key fields
        self._width = self._grid.at_link["channel_width"]
        self._topographic__elevation = self._grid.at_node["topographic__elevation"]

        # create field for channel_slope and topographic__elevation if they
        # don't yet exist.
        self.initialize_output_fields()
        self._channel_slope = self._grid.at_link["channel_slope"]

        # Adjust topographic elevation based on the parcels present.
        # Note that at present FlowDirector is just used for network connectivity.
        # get alluvium depth and calculate topography from br+alluvium, then update slopes.

        self._create_new_parcel_time()
        self._calculate_mean_D_and_rho()

        self._partition_active_and_storage_layers()
        self._adjust_node_elevation()
        self._update_channel_slopes()

    @property
    def time(self):
        """Return current time."""
        return self._time

    @property
    def d_mean_active(self):
        """Mean parcel grain size of active parcels aggregated at link."""
        return self._d_mean_active

    @property
    def rhos_mean_active(self):
        """Mean parcel density of active parcels aggregated at link."""
        return self._rhos_mean_active

    def _create_new_parcel_time(self):
        """ If we are going to track parcels through time in :py:class:`~landlab.data_record.data_record.DataRecord`, we
        need to add a new time column to the parcels dataframe. This method simply
        copies over the attributes of the parcels from the former timestep.
        Attributes will be updated over the course of this step.
        """

        if self._time_idx != 0:

            self._parcels.add_record(time=[self._time])

            self._parcels.ffill_grid_element_and_id()

            # copy parcel attributes forward in time.
            for at in self._time_variable_parcel_attributes:
                self._parcels.dataset[at].values[
                    :, self._time_idx
                ] = self._parcels.dataset[at].values[:, self._time_idx - 1]

        self._this_timesteps_parcels = np.zeros_like(
            self._parcels.dataset.element_id, dtype=bool
        )
        self._this_timesteps_parcels[:, -1] = True

        parcels_off_grid = (
            self._parcels.dataset.element_id[:, -1] == self.OUT_OF_NETWORK
        )
        self._this_timesteps_parcels[parcels_off_grid, -1] = False

        self._num_parcels = self._parcels.number_of_items
        # ^ needs to run just in case we've added more parcels

    def _update_channel_slopes(self):
        """Re-calculate channel slopes during each timestep."""

        for i in range(self._grid.number_of_links):

            upstream_node_id = self._fd.upstream_node_at_link()[i]
            downstream_node_id = self._fd.downstream_node_at_link()[i]

            self._channel_slope[i] = _recalculate_channel_slope(
                self._grid.at_node["topographic__elevation"][upstream_node_id],
                self._grid.at_node["topographic__elevation"][downstream_node_id],
                self._grid.at_link["reach_length"][i],
            )

    def _calculate_mean_D_and_rho(self):
        """Calculate mean grain size and density on each link"""

        current_parcels = self._parcels.dataset.isel(time=self._time_idx)

        # In the first full timestep, we need to calc grain size & rho_sed.
        # Assume all parcels are in the active layer for the purposes of
        # grain size and mean sediment density calculations

        # FUTURE: make it possible to circumvent this if mean grain size
        # has already been calculated (e.g. during 'zeroing' runs)

        # Calculate mean values for density and grain size (weighted by volume).
        sel_parcels = current_parcels.where(
            current_parcels.element_id != self.OUT_OF_NETWORK
        )

        d_weighted = sel_parcels.D * sel_parcels.volume
        rho_weighted = sel_parcels.density * sel_parcels.volume
        d_weighted.name = "d_weighted"
        rho_weighted.name = "rho_weighted"

        grouped_by_element = xr.merge(
            (sel_parcels.element_id, sel_parcels.volume, d_weighted, rho_weighted)
        ).groupby("element_id")

        d_avg = grouped_by_element.sum().d_weighted / grouped_by_element.sum().volume
        rho_avg = (
            grouped_by_element.sum().rho_weighted / grouped_by_element.sum().volume
        )

        self._d_mean_active = np.zeros(self._grid.size("link"))
        self._d_mean_active[d_avg.element_id.values.astype(int)] = d_avg.values

        self._rhos_mean_active = np.zeros(self._grid.size("link"))
        self._rhos_mean_active[rho_avg.element_id.values.astype(int)] = rho_avg.values

    def _partition_active_and_storage_layers(self, **kwds):
        """For each parcel in the network, determines whether it is in the
        active or storage layer during this timestep, then updates node
        elevations.

        """
        self._vol_tot = self._parcels.calc_aggregate_value(
            np.sum,
            "volume",
            at="link",
            filter_array=self._this_timesteps_parcels,
            fill_value=0.0,
        )

        if self._active_layer_method == "WongParker":
            # Wong et al. (2007) approximation for active layer thickness.
            # NOTE: calculated using grain size and grain density calculated for
            # the active layer grains in each link at the **previous** timestep.
            # This circumvents the need for an iterative scheme to determine grain
            # size of the active layer before determining which grains are in the
            # active layer.

            # calculate tau
            tau = (
                self._fluid_density
                * self._g
                * self._grid.at_link["channel_slope"]
                * self._grid.at_link["flow_depth"]
            )

            # calcuate taustar
            taustar = tau / (
                (self._rhos_mean_active - self._fluid_density)
                * self._g
                * self._d_mean_active
            )

            # calculate active layer thickness
            self._active_layer_thickness = (
                0.515 * self._d_mean_active * (3.09 * (taustar - 0.0549) ** 0.56)
            )  # in units of m

        elif self._active_layer_method == "GrainSizeDependent":
            # Set all active layers to a multiple of the lnk mean grain size
            self._active_layer_thickness = (
                self._d_mean_active * self._active_layer_d_multiplier
            )

        elif self._active_layer_method == "Constant10cm":
            # Set all active layers to 10 cm thickness.
            self._active_layer_thickness = 0.1 * np.ones_like(self._d_mean_active)

        # If links have no parcels, we still need to assign them an active layer
        # thickness..
        links_with_no_active_layer = np.isnan(self._active_layer_thickness)
        self._active_layer_thickness[links_with_no_active_layer] = np.mean(
            self._active_layer_thickness[links_with_no_active_layer == 0]
        )  # assign links with no parcels an average value

        if np.sum(np.isfinite(self._active_layer_thickness)) == 0:
            self._active_layer_thickness.fill(_INIT_ACTIVE_LAYER_THICKNESS)
            # handles the case of the first timestep -- assigns a modest value

        capacity = (
            self._grid.at_link["channel_width"]
            * self._grid.at_link["reach_length"]
            * self._active_layer_thickness
        )  # in units of m^3

        active_inactive = _INACTIVE * np.ones(self._num_parcels)

        current_link = self._parcels.dataset.element_id.values[:, -1].astype(int)
        time_arrival = self._parcels.dataset.time_arrival_in_link.values[:, -1]
        volumes = self._parcels.dataset.volume.values[:, -1]

        for i in range(self._grid.number_of_links):

            if (
                self._vol_tot[i] > 0
            ):  # only do this check capacity if parcels are in link

                # First In Last Out.

                # Find parcels on this link.
                this_links_parcels = np.where(current_link == i)[0]

                # sort them by arrival time.
                time_arrival_sort = np.flip(
                    np.argsort(time_arrival[this_links_parcels], 0,)
                )
                parcel_id_time_sorted = this_links_parcels[time_arrival_sort]

                # calculate the cumulative volume (in sorted order).
                cumvol = np.cumsum(volumes[parcel_id_time_sorted])

                # determine which parcels are within capacity and set those to
                # active.
                make_active = parcel_id_time_sorted[cumvol <= capacity[i]]

                active_inactive[make_active] = _ACTIVE

        self._parcels.dataset.active_layer[:, -1] = active_inactive

        # set active here. reference it below in wilcock crowe
        self._active_parcel_records = (
            self._parcels.dataset.active_layer == _ACTIVE
        ) * (self._this_timesteps_parcels)

        self._vol_act = self._parcels.calc_aggregate_value(
            np.sum,
            "volume",
            at="link",
            filter_array=self._active_parcel_records,
            fill_value=0.0,
        )

        self._vol_stor = (self._vol_tot - self._vol_act) / (1 - self._bed_porosity)

    def _adjust_node_elevation(self):
        """Adjusts slope for each link based on parcel motions from last
        timestep and additions from this timestep.
        """

        number_of_contributors = np.sum(
            self._fd.flow_link_incoming_at_node() == 1, axis=1
        )
        downstream_link_id = self._fd.link_to_flow_receiving_node
        # USED TO BE      downstream_link_id = self._fd.link_to_flow_receiving_node[
        #            self._fd.downstream_node_at_link()
        #        ]
        upstream_contributing_links_at_node = np.where(
            self._fd.flow_link_incoming_at_node() == 1, self._grid.links_at_node, -1
        )

        # Update the node topographic elevations depending on the quantity of stored sediment
        for n in range(self._grid.number_of_nodes):

            if number_of_contributors[n] > 0:  # we don't update head node elevations

                upstream_links = upstream_contributing_links_at_node[n]
                real_upstream_links = upstream_links[
                    upstream_links != self._grid.BAD_INDEX
                ]
                width_of_upstream_links = self._grid.at_link["channel_width"][
                    real_upstream_links
                ]
                length_of_upstream_links = self._grid.at_link["reach_length"][
                    real_upstream_links
                ]

                #                ALERT: Moved this to the "else" statement below. AP 11/11/19
                #                length_of_downstream_link = self._grid.at_link["reach_length"][
                #                    downstream_link_id
                #                ][n]
                #                width_of_downstream_link = self._grid.at_link["channel_width"][
                #                    downstream_link_id
                #                ][n]

                if (
                    downstream_link_id[n] == self._grid.BAD_INDEX
                ):  # I'm sure there's a better way to do this, but...
                    length_of_downstream_link = 0
                    width_of_downstream_link = 0
                else:
                    length_of_downstream_link = self._grid.at_link["reach_length"][
                        downstream_link_id
                    ][n]
                    width_of_downstream_link = self._grid.at_link["channel_width"][
                        downstream_link_id
                    ][n]

                alluvium__depth = _calculate_alluvium_depth(
                    self._vol_stor[downstream_link_id][n],
                    width_of_upstream_links,
                    length_of_upstream_links,
                    width_of_downstream_link,
                    length_of_downstream_link,
                    self._bed_porosity,
                )

                self._grid.at_node["topographic__elevation"][n] = (
                    self._grid.at_node["bedrock__elevation"][n] + alluvium__depth
                )

    def _calc_transport_wilcock_crowe(self):
        """Method to determine the transport time for each parcel in the active
        layer using a sediment transport equation.

        Note: could have options here (e.g. Wilcock and Crowe, FLVB, MPM, etc)
        """
        # Initialize _pvelocity, the virtual velocity of each parcel (link length / link travel time)
        self._pvelocity = np.zeros(self._num_parcels)

        # parcel attribute arrays from DataRecord

        Darray = self._parcels.dataset.D[:, self._time_idx]
        Activearray = self._parcels.dataset.active_layer[:, self._time_idx].values
        Rhoarray = self._parcels.dataset.density.values
        Volarray = self._parcels.dataset.volume[:, self._time_idx].values
        Linkarray = self._parcels.dataset.element_id[
            :, self._time_idx
        ].values  # link that the parcel is currently in

        R = (Rhoarray - self._fluid_density) / self._fluid_density

        # parcel attribute arrays to populate below
        frac_sand_array = np.zeros(self._num_parcels)
        vol_act_array = np.zeros(self._num_parcels)
        Sarray = np.zeros(self._num_parcels)
        Harray = np.zeros(self._num_parcels)
        Larray = np.zeros(self._num_parcels)
        D_mean_activearray = np.zeros(self._num_parcels) * (np.nan)
        active_layer_thickness_array = np.zeros(self._num_parcels) * np.nan
        #        rhos_mean_active = np.zeros(self._num_parcels)
        #        rhos_mean_active.fill(np.nan)

        # find active sand
        findactivesand = (
            self._parcels.dataset.D < _SAND_SIZE
        ) * self._active_parcel_records  # since find active already sets all prior timesteps to False, we can use D for all timesteps here.

        vol_act_sand = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=findactivesand, fill_value=0.0
        )

        frac_sand = np.zeros_like(self._vol_act)
        frac_sand[self._vol_act != 0.0] = (
            vol_act_sand[self._vol_act != 0.0] / self._vol_act[self._vol_act != 0.0]
        )
        frac_sand[np.isnan(frac_sand)] = 0.0

        # Calc attributes for each link, map to parcel arrays
        for i in range(self._grid.number_of_links):

            active_here = np.nonzero(
                np.logical_and(Linkarray == i, Activearray == _ACTIVE)
            )[0]
            d_act_i = Darray[active_here]
            vol_act_i = Volarray[active_here]
            rhos_act_i = Rhoarray[active_here]
            vol_act_tot_i = np.sum(vol_act_i)

            self._d_mean_active[i] = np.sum(d_act_i * vol_act_i) / (vol_act_tot_i)
            if vol_act_tot_i > 0:
                self._rhos_mean_active[i] = np.sum(rhos_act_i * vol_act_i) / (
                    vol_act_tot_i
                )
            else:
                self._rhos_mean_active[i] = np.nan
            D_mean_activearray[Linkarray == i] = self._d_mean_active[i]
            frac_sand_array[Linkarray == i] = frac_sand[i]
            vol_act_array[Linkarray == i] = self._vol_act[i]
            Sarray[Linkarray == i] = self._grid.at_link["channel_slope"][i]
            Harray[Linkarray == i] = self._grid.at_link["flow_depth"][i]
            Larray[Linkarray == i] = self._grid.at_link["reach_length"][i]
            active_layer_thickness_array[Linkarray == i] = self._active_layer_thickness[
                i
            ]

        # Wilcock and Crowe calculate transport for all parcels (active and inactive)
        taursg = _calculate_reference_shear_stress(
            self._fluid_density, R, self._g, D_mean_activearray, frac_sand_array
        )

        frac_parcel = np.nan * np.zeros_like(Volarray)
        frac_parcel[vol_act_array != 0.0] = (
            Volarray[vol_act_array != 0.0] / vol_act_array[vol_act_array != 0.0]
        )

        b = 0.67 / (1.0 + np.exp(1.5 - Darray / D_mean_activearray))

        tau = self._fluid_density * self._g * Harray * Sarray
        tau = np.atleast_1d(tau)

        taur = taursg * (Darray / D_mean_activearray) ** b
        tautaur = tau / taur
        tautaur_cplx = tautaur.astype(np.complex128)
        # ^ work around needed b/c np fails with non-integer powers of negative numbers

        W = 0.002 * np.power(tautaur_cplx.real, 7.5)
        W[tautaur >= 1.35] = 14 * np.power(
            (1 - (0.894 / np.sqrt(tautaur_cplx.real[tautaur >= 1.35]))), 4.5
        )

        active_parcel_idx = Activearray == _ACTIVE

        # compute parcel virtual velocity, m/s
        self._pvelocity[active_parcel_idx] = (
            W.real[active_parcel_idx]
            * (tau[active_parcel_idx] ** (3.0 / 2.0))
            * frac_parcel[active_parcel_idx]
            / (self._fluid_density ** (3.0 / 2.0))
            / self._g
            / R[active_parcel_idx]
            / active_layer_thickness_array[active_parcel_idx]
        )

        self._pvelocity[np.isnan(self._pvelocity)] = 0.0

        if np.max(self._pvelocity) > 1:
            warnings.warn(
                "NetworkSedimentTransporter: Maximum parcel virtual velocity exceeds 1 m/s"
            )

        # Assign those things to the grid -- might be useful for plotting
        self._grid.at_link["sediment_total_volume"] = self._vol_tot
        self._grid.at_link["sediment__active__volume"] = self._vol_act
        self._grid.at_link["sediment__active__sand_fraction"] = frac_sand

    def _move_parcel_downstream(self, dt):
        """Method to update parcel location for each parcel in the active
        layer.
        """
        # determine where parcels are starting
        current_link = self._parcels.dataset.element_id.values[:, -1].astype(int)
        self.current_link = current_link

        # determine location within link where parcels are starting.
        location_in_link = self._parcels.dataset.location_in_link.values[:, -1]

        # determine how far each parcel needs to travel this timestep.
        distance_to_travel_this_timestep = self._pvelocity * dt
        # total distance traveled in dt at parcel virtual velocity
        # Note: movement in current and any DS links at this dt is at the same
        # velocity as in the current link perhaps modify in the future

        # Accumulate the total distance traveled by a parcel for abrasion rate
        # calculations.
        if np.size(self._distance_traveled_cumulative) != np.size(
            distance_to_travel_this_timestep
        ):
            dist_array = distance_to_travel_this_timestep
            dist_array[: self._num_parcels] += distance_to_travel_this_timestep
            self._distance_traveled_cumulative = dist_array
        else:
            self._distance_traveled_cumulative += distance_to_travel_this_timestep
            # ^ accumulates total distanced traveled for testing abrasion

        # active parcels on the network:
        in_network = (
            self._parcels.dataset.element_id.values[:, self._time_idx]
            != self.OUT_OF_NETWORK
        )
        active = distance_to_travel_this_timestep > 0.0
        active_parcel_ids = np.nonzero(in_network * active)[0]

        distance_left_to_travel = distance_to_travel_this_timestep.copy()
        while np.any(distance_left_to_travel > 0.0):

            # Step 1: Move parcels downstream.
            on_network = current_link != self.OUT_OF_NETWORK

            # Get current link lengths:
            current_link_lengths = self._grid.at_link["reach_length"][current_link]

            # Determine where they are in the current link.
            distance_to_exit_current_link = current_link_lengths * (
                1.0 - location_in_link
            )

            # Identify which ones will come to rest in the current link.
            rest_this_link = (
                (distance_left_to_travel < distance_to_exit_current_link)
                * on_network
                * (distance_left_to_travel > 0.0)
            )

            # Deal with those staying in the current link.
            if np.any(rest_this_link):
                # print('  {x} coming to rest'.format(x=np.sum(rest_this_link)))

                # for those staying in this link, calculate the location in link
                # (note that this is a proportional distance). AND change distance_left_to_travel to 0.0
                location_in_link[rest_this_link] = 1.0 - (
                    (
                        distance_to_exit_current_link[rest_this_link]
                        - distance_left_to_travel[rest_this_link]
                    )
                    / current_link_lengths[rest_this_link]
                )

                distance_left_to_travel[rest_this_link] = 0.0

            # Deal with those moving to a downstream link.
            moving_downstream = (
                (distance_left_to_travel >= distance_to_exit_current_link)
                * on_network
                * (distance_left_to_travel > 0.0)
            )
            if np.any(moving_downstream):
                # print('  {x} next link'.format(x=np.sum(moving_downstream)))
                # change location in link to 0
                location_in_link[moving_downstream] = 0.0

                # decrease distance to travel.
                distance_left_to_travel[
                    moving_downstream
                ] -= distance_to_exit_current_link[moving_downstream]

                # change current link to the downstream link.

                # get the downstream link at link:
                downstream_node = self._fd.downstream_node_at_link()[current_link]
                downstream_link = self._fd.link_to_flow_receiving_node[downstream_node]

                # assign new values to current link.
                current_link[moving_downstream] = downstream_link[moving_downstream]

                # find and address those links who have moved out of network.
                moved_oon = downstream_link == self._grid.BAD_INDEX

                if np.any(moved_oon):
                    # print('  {x} exiting network'.format(x=np.sum(moved_oon)))

                    current_link[moved_oon] = self.OUT_OF_NETWORK
                    # assign location in link of np.nan to those which moved oon
                    location_in_link[moved_oon] = np.nan
                    distance_left_to_travel[moved_oon] = 0.0

        # Step 2: Parcel is at rest... Now update its information.

        # reduce D and volume due to abrasion
        vol = _calculate_parcel_volume_post_abrasion(
            self._parcels.dataset.volume[active_parcel_ids, self._time_idx],
            distance_to_travel_this_timestep[active_parcel_ids],
            self._parcels.dataset.abrasion_rate[active_parcel_ids],
        )

        D = _calculate_parcel_grain_diameter_post_abrasion(
            self._parcels.dataset.D[active_parcel_ids, self._time_idx],
            self._parcels.dataset.volume[active_parcel_ids, self._time_idx],
            vol,
        )

        # update parcel attributes

        # arrival time in link
        self._parcels.dataset.time_arrival_in_link[
            active_parcel_ids, self._time_idx
        ] = self._time_idx

        # location in link
        self._parcels.dataset.location_in_link[
            active_parcel_ids, self._time_idx
        ] = location_in_link[active_parcel_ids]

        self._parcels.dataset.element_id[
            active_parcel_ids, self._time_idx
        ] = current_link[active_parcel_ids]
        #                self._parcels.dataset.active_layer[p, self._time_idx] = 1
        # ^ reset to 1 (active) to be recomputed/determined at next timestep
        self._parcels.dataset.D[active_parcel_ids, self._time_idx] = D
        self._parcels.dataset.volume[active_parcel_ids, self._time_idx] = vol

    def run_one_step(self, dt):
        """Run NetworkSedimentTransporter forward in time.

        When the NetworkSedimentTransporter runs forward in time the following
        steps occur:

            1. A new set of records is created in the Parcels that corresponds to the new time
            2. If parcels are on the network then:
                a. Active parcels are identifed based on entrainment critera.
                b. Effective bed slope is calculated based on inactive parcel volumes.
                c. Transport rate is calculated.
                d. Active parcels are moved based on the tranport rate.

        Parameters
        ----------
        dt : float
            Duration of time to run the NetworkSedimentTransporter forward.

        Returns
        -------
        RuntimeError if no parcels remain on the grid.

        """
        self._time += dt
        self._time_idx += 1
        self._create_new_parcel_time()

        if self._this_timesteps_parcels.any():
            self._partition_active_and_storage_layers()
            self._adjust_node_elevation()
            self._update_channel_slopes()
            self._update_transport_time()
            self._move_parcel_downstream(dt)

        else:
            msg = "No more parcels on grid"
            raise RuntimeError(msg)


# %% Methods referenced above, separated for purposes of testing


def _recalculate_channel_slope(z_up, z_down, dx, threshold=1e-4):
    """Recalculate channel slope based on elevation.

    Parameters
    ----------
    z_up : float
        Upstream elevation.
    z_down : float
        Downstream elevation.
    dz : float
        Distance.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _recalculate_channel_slope
    >>> import pytest
    >>> _recalculate_channel_slope(10., 0., 10.)
    1.0
    >>> _recalculate_channel_slope(0., 0., 10.)
    0.0001
    >>> with pytest.warns(UserWarning):
    ...     _recalculate_channel_slope(0., 10., 10.)
    0.0

    """
    chan_slope = (z_up - z_down) / dx

    if chan_slope < 0.0:
        chan_slope = 0.0
        warnings.warn(
            "NetworkSedimentTransporter: Negative channel slope encountered.",
            UserWarning,
        )

    elif chan_slope < threshold:
        chan_slope = threshold

    return chan_slope


def _calculate_alluvium_depth(
    stored_volume,
    width_of_upstream_links,
    length_of_upstream_links,
    width_of_downstream_link,
    length_of_downstream_link,
    porosity,
):
    """Calculate alluvium depth based on adjacent link inactive parcel volumes.

    Parameters
    ----------
    stored_volume : float
        Total volume of inactive parcels in this link.
    width_of_upstream_links : float
        Channel widths of upstream links.
    length_of_upstream_link : float
        Channel lengths of upstream links.
    width_of_downstream_link : float
        Channel widths of downstream links.
    length_of_downstream_link : float
        Channel lengths of downstream links.
    porosity: float
        Channel bed sediment porosity.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_alluvium_depth
    >>> import pytest
    >>> _calculate_alluvium_depth(100,np.array([0.5,1]),np.array([10,10]), 1, 10, 0.2)
    10.0
    >>> _calculate_alluvium_depth(24,np.array([0.1,3]),np.array([10,10]), 1, 1, 0.5)
    3.0
    >>> with pytest.raises(ValueError):
    ...     _calculate_alluvium_depth(24,np.array([0.1,3]),np.array([10,10]), 1, 1, 2)

    """

    alluvium__depth = (
        2
        * stored_volume
        / (
            np.sum(width_of_upstream_links * length_of_upstream_links)
            + width_of_downstream_link * length_of_downstream_link
        )
        / (1 - porosity)
    )

    if alluvium__depth < 0.0:
        raise ValueError("NST Alluvium Depth Negative")

    return alluvium__depth


def _calculate_reference_shear_stress(
    fluid_density, R, g, mean_active_grain_size, frac_sand
):
    """Calculate reference Shields stress (taursg) using the sand content of
    the bed surface, as per Wilcock and Crowe (2003).

    Parameters
    ----------
    fluid_density : float
        Density of fluid (generally, water).
    R: float
        Specific weight..?
    g: float
        Gravitational acceleration.
    mean_active_grain_size: float
        Mean grain size of the 'active' sediment parcels.
    frac_sand: float
        Fraction of the bed surface grain size composed of sand sized parcels.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import (
    ... _calculate_reference_shear_stress)
    >>> from numpy.testing import assert_almost_equal
    >>> assert_almost_equal(
    ...     _calculate_reference_shear_stress(1, 1, 1, 1, 0),
    ...     0.036,
    ...     decimal=2)
    >>> assert_almost_equal(
    ...     _calculate_reference_shear_stress(1000, 1.65, 9.8, 0.1, 0.9),
    ...     33.957,
    ...     decimal=2)

    """

    taursg = (
        fluid_density
        * R
        * g
        * mean_active_grain_size
        * (0.021 + 0.015 * np.exp(-20.0 * frac_sand))
    )

    if np.any(np.asarray(taursg < 0)):
        raise ValueError(
            "NetworkSedimentTransporter: Reference Shields stress is negative"
        )

    return taursg


def _calculate_parcel_volume_post_abrasion(
    starting_volume, travel_distance, abrasion_rate
):
    """Calculate parcel volumes after abrasion, according to Sternberg
    exponential abrasion.

    Parameters
    ----------
    starting_volume : float or array
        Starting volume of each parcel.
    travel_distance: float or array
        Travel distance for each parcel during this timestep, in ___.
    abrasion_rate: float or array
        Mean grain size of the 'active' sediment parcels.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_parcel_volume_post_abrasion
    >>> import pytest
    >>> _calculate_parcel_volume_post_abrasion(10,100,0.003)
    7.4081822068171785
    >>> _calculate_parcel_volume_post_abrasion(10,300,0.1)
    9.3576229688401746e-13
    >>> with pytest.raises(ValueError):
    ...     _calculate_parcel_volume_post_abrasion(10,300,-3)

    """

    volume = starting_volume * np.exp(travel_distance * (-abrasion_rate))

    if np.any(volume > starting_volume):
        raise ValueError("NST parcel volume *increases* due to abrasion")

    return volume


def _calculate_parcel_grain_diameter_post_abrasion(
    starting_diameter, pre_abrasion_volume, post_abrasion_volume
):
    """Calculate parcel grain diameters after abrasion, according to Sternberg
    exponential abrasion.

    Parameters
    ----------
    starting_diameter : float or array
        Starting volume of each parcel.
    pre_abrasion_volume: float or array
        Parcel volume before abrasion.
    post_abrasion_volume: float or array
        Parcel volume after abrasion.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_parcel_grain_diameter_post_abrasion
    >>> import numpy as np
    >>> from numpy.testing import assert_almost_equal

    If no abrasion happens, we should get the same value.

    >>> _calculate_parcel_grain_diameter_post_abrasion(10, 1, 1)
    10.0

    If some abrasion happens, test the value.

    >>> starting_diameter = 10
    >>> pre_abrasion_volume = 2
    >>> post_abrasion_volume = 1
    >>> expected_value = (
    ...     starting_diameter *
    ...     ( post_abrasion_volume / pre_abrasion_volume) ** (1. / 3.))
    >>> print(np.round(expected_value, decimals=3))
    7.937
    >>> assert_almost_equal(
    ...    _calculate_parcel_grain_diameter_post_abrasion(10, 2, 1),
    ...    expected_value)

    """

    abraded_grain_diameter = starting_diameter * (
        post_abrasion_volume / pre_abrasion_volume
    ) ** (1.0 / 3.0)

    return abraded_grain_diameter
