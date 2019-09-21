#!/usr/env/python

"""Landlab component that simulates xxxxxx

info about the component here

Fixes that need to happen:

    -- (9/16/19) Something is wrong with abrasion (or at least that's where
    it's obvious). See test comparing actual abrasion to calculated. Hmm..
    -- Check abrasion exponent units (per km or per m?)

    -- Channel width-- why is this anything other than a link attribute??

    -- JC: I found two items that I think should be changed in the _calc_transport_wilcock_crowe and I made these changes
            - frac_parcels was the inverse of what it should be so instead of a fraction it was a number >1
            - in unpacking W* we had a (1-frac_sand) this was fine when we were treating sand and gravel separately,
              but now that we are treating all parcels together, I no longer think this should be there, because if we are trying
              to move a sand parcel then this (1-frac_sand) does not make sense. I think this is  now equivalent to the original WC2003.
              Before it was equivalent to WC2003 as implemented in Cui TUGS formulation.
            - finally, I added the calculation of a parcel velocity instead of the travel time. I think this is
              better suited to the parcel centric spirit of the code. It was also needed to simplify move_parcel_downstream
              Plus, I think one day we will have a better way to parameterize parcel virtual velocity and this will then be
              easy to incorporate/update.

    -- Fix inelegant time indexing

.. codeauthor:: Jon Allison Katy

Created on Tu May 8, 2018
Last edit ---
"""

import numpy as np

# %% Import Libraries
from landlab import BAD_INDEX_VALUE, Component
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.components import FlowDirectorSteepest
from landlab.utils.decorators import use_file_name_or_kwds

_SUPPORTED_TRANSPORT_METHODS = ["WilcockCrowe"]

_OUT_OF_NETWORK = BAD_INDEX_VALUE - 1


_ACTIVE = 1
_INACTIVE = 0

class NetworkSedimentTransporter(Component):
    """Network bedload morphodynamic component.

    Landlab component designed to calculate _____.
    info info info

    **Usage:**
    Option 1 - Basic::
        NetworkSedimentTransporter(grid,
                             parcels,
                             transporter = asdfasdf,
                             discharge,
                             channel_geometry
                             )

    Examples
    ----------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorSteepest, FlowAccumulator, NetworkSedimentTransporter
    >>> from landlab.components.landslides import LandslideProbability
    >>> import numpy as np

    The NetworkSedimentTransporter moves parcels of sediment based on a given
    flow and a given sediment transport formulation. Thus we must first create
    parcels. The parcels must be a ``DataRecord`` with the following attributes

    - name
    - other thing
    -

    For example:

    >>> parcels = ...

    Next we must have a NetworkModelGrid on which the parcels are transported:

    >>> make nmg

    We are now ready to set up NetworkSedimentTransporter

    >>> nst = NetworkSedimentTransporter(grid, stuff)

    Finally, we un NetworkSedimentTransporter forward 10 timesteps of size 10
    time units.

    >>> for _ in range(10):
    ...     nst.run_one_step(10.)

    Look! Parcel 10 moved!

    >>>
    >>>

    """

    # component name
    _name = "NetworkSedimentTransporter"
    __version__ = "1.0"

    # component requires these values to do its calculation, get from driver
    _input_var_names = (
        "topographic__elevation",
        "channel_slope",
        "link_length",
        "channel_width",
        "flow_depth",
    )

    #  component creates these output values
    _output_var_names = ("thing", "thing")

    # units for each parameter and output
    _var_units = {
        "topographic__elevation": "m",
        "channel_slope": "m/m",
        "link_length": "m",
        "channel_width": "m",
        "flow_depth": "m",
    }

    # grid centering of each field and variable
    _var_mapping = {
        "topographic__elevation": "node",
        "channel_slope": "link",
        "link_length": "link",
        "channel_width": "link",
        "flow_depth": "link",
    }

    # short description of each field
    _var_doc = {
        "surface_water__depth": "Depth of streamflow on the surface",
        "topographic__elevation": "Topographic elevation at that node",
        "channel_slope": "Slope of the river channel through each reach",
        "link_length": "Length of each reach",
        "channel_width": "Flow width of the channel, assuming constant width",
        "flow_depth": "Depth of stream flow in each reach",
    }

    # Run Component
    #    @use_file_name_or_kwds
    #   Katy! We had to comment out ^ that line in order to get NST to instantiate. Back end changes needed.

    def __init__(
        self,
        grid,
        parcels,
        flow_director,
        flow_depth,
        bed_porosity=0.3,
        g=9.81,
        fluid_density=1000.,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
        **kwds
    ):
        """
        Parameters
        ----------
        grid: NetworkModelGrid
            A landlab network model grid in which links are stream channel
            segments.
        parcels: DataRecord
            A landlab DataRecord describing the characteristics and location of
            sediment "parcels".

            Either put more information about parcels here or put it above and
            reference it here.

        flow_director: FlowDirectorSteepest
            A landlab flow director. Currently, must be FlowDirectorSteepest.
        flow_depth: float, numpy array of shape (timesteps,links)
            Flow depth of water in channel at each link at each timestep. (m)
        bed_porosity: float, optional
            Proportion of void space between grains in the river channel bed.
            Default value is 0.3.
        g: float, optional
            Acceleration due to gravity. Default value is 9.81 (m/s^2)
        fluid_density: float, optional
            Density of the fluid (generally, water) in which sediment is
            moving. Default value is 1000 (kg/m^3)
        channel_width: float, optional
            DANGER DANGER-- Why don't we have this attached to the grid?
        transport_method: string
            Sediment transport equation option. Default (and currently only)
            option is "WilcockCrowe".
        """
        super(NetworkSedimentTransporter, self).__init__(grid, **kwds)

        self._grid = grid

        if not isinstance(grid, NetworkModelGrid):
            msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            raise ValueError(msg)

        self._parcels = parcels

        if not isinstance(parcels, DataRecord):
            msg = (
                "NetworkSedimentTransporter: parcels must be an instance"
                "of DataRecord"
            )
            raise ValueError(msg)

        self._num_parcels = self._parcels.dataset.element_id.size

        self.parcel_attributes = [
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
        self.fd = flow_director
        self.flow_depth = flow_depth
        self.bed_porosity = bed_porosity

        if not 0 <= self.bed_porosity < 1:
            msg = "NetworkSedimentTransporter: bed_porosity must be" "between 0 and 1"
            raise ValueError(msg)

        self.g = g
        self.fluid_density = fluid_density
        self._time_idx = 0
        self._time = 0.0

        if transport_method in _SUPPORTED_TRANSPORT_METHODS:
            self.transport_method = transport_method
        else:
            msg = "Transport Method not supported"
            raise ValueError(msg)
        # self.transport_method makes it a class variable, that can be accessed within any method within this class
        if self.transport_method == "WilcockCrowe":
            self.update_transport_time = self._calc_transport_wilcock_crowe

        self._width = self._grid.at_link[channel_width]

        if "channel_width" not in self._grid.at_link:
            msg = (
                "NetworkSedimentTransporter: channel_width must be assigned"
                "to the grid links"
            )
            raise ValueError(msg)

        if "topographic__elevation" not in self._grid.at_node:
            msg = (
                "NetworkSedimentTransporter: topographic__elevation must be "
                "assigned to the grid nodes"
            )
            raise ValueError(msg)

        if "bedrock__elevation" not in self._grid.at_node:
            msg = (
                "NetworkSedimentTransporter: topographic__elevation must be "
                "assigned to the grid nodes"
            )
            raise ValueError(msg)

        if "link_length" not in self._grid.at_link:
            msg = (
                "NetworkSedimentTransporter: link_length must be assigned"
                "to the grid links"
            )
            raise ValueError(msg)

        if "drainage_area" not in self._grid.at_link:
            msg = (
                "NetworkSedimentTransporter: channel_width must be assigned"
                "to the grid links"
            )
            raise ValueError(msg)

        # create field for channel slope if it doesnt exist yet.
        if "channel_slope" not in self._grid.at_link:
            self._channel_slope = self._grid.zeros(at="node")
            self._update_channel_slopes()
        else:
            self._channel_slope = self._grid.at_link["channel_slope"]

        if "time_arrival_in_link" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: time_arrival_in_link must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)

        if "starting_link" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: starting_link must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)

        if "abrasion_rate" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: abrasion_rate must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)

        if "density" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: density must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)

        if "active_layer" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: active_layer must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)

        if "location_in_link" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: location_in_link must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)

        if "D" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: D (grain size) must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)

        if "volume" not in self._parcels.dataset:
            msg = ("NetworkSedimentTransporter: volume must be"
                   "assigned to the parcels"
                   )
            raise ValueError(msg)


    @property
    def time(self):
        """Return current time."""
        return self._time

    def _create_new_parcel_time(self):
        """ If we are going to track parcels through time in DataRecord, we
        need to add a new time column to the parcels dataframe. This method simply
        copies over the attributes of the parcels from the former timestep.
        Attributes will be updated over the course of this step.
        """

        if self._time_idx != 0:

            self._parcels.add_record(time=[self._time])
            # ^ what's the best way to handle time?
            #            self._parcels.dataset['grid_element'].values[:,self._time_idx] = self._parcels.dataset[
            #                    'grid_element'].values[:,self._time_idx-1]
            #
            #            self._parcels.dataset['element_id'].values[:,self._time_idx] = self._parcels.dataset[
            #                    'element_id'].values[:,self._time_idx-1]

            self._parcels.ffill_grid_element_and_id()

            for at in self.parcel_attributes:
                self._parcels.dataset[at].values[
                    :, self._time_idx
                ] = self._parcels.dataset[at].values[:, self._time_idx - 1]

        self._find_now = self._parcels.dataset.time == self._time
        self._this_timesteps_parcels = np.zeros_like(
            self._parcels.dataset.element_id, dtype=bool
        )
        self._this_timesteps_parcels[:, -1] = True

        self._parcels_off_grid = (
            self._parcels.dataset.element_id[:, -1] == _OUT_OF_NETWORK
        )
        self._this_timesteps_parcels[self._parcels_off_grid, -1] = False

    def _update_channel_slopes(self):
        """text Can be simple-- this is what this does. 'private' functions can
        have very simple examples, explanations. Essentially note to yourself"""
        # Katy think this can be vectorized
        # Jon agrees, but is not sure yet how to do that
        for i in range(self._grid.number_of_links):

            upstream_node_id = self.fd.upstream_node_at_link()[i]
            downstream_node_id = self.fd.downstream_node_at_link()[i]

            self._channel_slope[i] = _recalculate_channel_slope(
                self._grid.at_node["topographic__elevation"][upstream_node_id],
                self._grid.at_node["topographic__elevation"][downstream_node_id],
                self._grid.length_of_link[i],
            )

    def _partition_active_and_storage_layers(self, **kwds):
        """For each parcel in the network, determines whether it is in the
        active or storage layer during this timestep, then updates node elevations
        """

        vol_tot = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._this_timesteps_parcels
        )
        vol_tot[np.isnan(vol_tot) == 1] = 0

        # Wong et al. (2007) approximation for active layer thickness.
        # NOTE: calculated using grain size and grain density calculated for
        # the active layer grains in each link at the **previous** timestep.
        # This circumvents the need for an iterative scheme to determine grain
        # size of the active layer before determining which grains are in the
        # active layer.

        if self._time_idx == 1:
            # In the first full timestep, we need to calc grain size & rho_sed.
            # Assume all parcels are in the active layer for the purposes of
            # grain size and mean sediment density calculations

            # FUTURE: make it possible to circumvent this if mean grain size
            # has already been calculated (e.g. during 'zeroing' runs)
            Linkarray = self._parcels.dataset.element_id[
                :, self._time_idx
                ].values  # current link of each parcel
            Darray = self._parcels.dataset.D[:, self._time_idx]
            Rhoarray = self._parcels.dataset.density.values
            Volarray = self._parcels.dataset.volume[:, self._time_idx].values

            d_mean_active = np.nan * np.zeros(self._grid.number_of_links)
            rhos_mean_active = np.nan * np.zeros(self._grid.number_of_links)

            for i in range(self._grid.number_of_links):
                d_i = Darray[Linkarray ==i]
                vol_i = Volarray[Linkarray ==i]
                rhos_i = Rhoarray[Linkarray ==i]
                vol_tot_i = np.sum(vol_i)

                d_mean_active[i] = np.sum(d_i * vol_i) / (
                    vol_tot_i
                )
                self.d_mean_active = d_mean_active

                rhos_mean_active[i] = np.sum(rhos_i * vol_i) / (
                    vol_tot_i
                )
                self.rhos_mean_active = rhos_mean_active

        tau = (self.fluid_density
               * self.g
               * self._grid.at_link["channel_slope"]
               * self.flow_depth[self._time_idx, :]
               )

        taustar = (tau
                   /((self.rhos_mean_active-self.fluid_density)
                   *self.g
                   *self.d_mean_active)
                   )

        self.active_layer_thickness = (0.515
                                  *self.d_mean_active
                                  *(3.09
                                    *(taustar-0.0549)**0.56)
                                  ) # in units of m

        capacity = (self._grid.at_link["channel_width"]
                    *self._grid.at_link["link_length"]
                    *self.active_layer_thickness
                    ) # in units of m^3

#       OLD capacity calculation
#       capacity = 2 * np.ones(
#            self._grid.number_of_links
#        )  # in units of m^

        for i in range(self._grid.number_of_links):

            if vol_tot[i] > 0:  # only do this check capacity if parcels are in link

                # First In Last Out.
                parcel_id_thislink = np.where(
                    self._parcels.dataset.element_id[:, self._time_idx] == i
                )[0]

                time_arrival_sort = np.flip(
                    np.argsort(
                        self._parcels.get_data(
                            time=[self._time],
                            item_id=parcel_id_thislink,
                            data_variable="time_arrival_in_link",
                        ),
                        0,
                    )
                )

                parcel_id_time_sorted = parcel_id_thislink[time_arrival_sort]

                cumvol = np.cumsum(
                    self._parcels.dataset.volume[parcel_id_time_sorted, self._time_idx]
                )

                idxinactive = np.where(cumvol > capacity[i])
                make_inactive = parcel_id_time_sorted[idxinactive]
                # idxbedabrade = np.where(cumvol < 2*capacity[i] and cumvol > capacity[i])
                # ^ syntax is wrong, but this is where we can identify the surface of the bed
                # for abrasion, I think we would abrade the particles in the active layer in active transport
                # and abrade the particles sitting on the bed. This line would identify those particles on
                # the bed that also need to abrade due to impacts from the sediment moving above.

                self._parcels.set_data(
                    time=[self._time],
                    item_id=parcel_id_thislink,
                    data_variable="active_layer",
                    new_value=_ACTIVE,
                )

                self._parcels.set_data(
                    time=[self._time],
                    item_id=make_inactive,
                    data_variable="active_layer",
                    new_value=_INACTIVE,
                )

        # Update Node Elevations

        # set active here. reference it below in wilcock crowe
        self._active_parcel_records = (
            self._parcels.dataset.active_layer == _ACTIVE
        ) * (self._this_timesteps_parcels)

        # print("active_parcel_records",self._active_parcel_records)

        vol_act = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._active_parcel_records
        )

        vol_act[np.isnan(vol_act) == 1] = 0

        self.vol_stor = (vol_tot - vol_act) / (1 - self.bed_porosity)

    # %%
    def _adjust_node_elevation(self):
        """Adjusts slope for each link based on parcel motions from last
        timestep and additions from this timestep.
        """

        number_of_contributors = np.sum(
            self.fd.flow_link_incoming_at_node() == 1, axis=1
        )
        downstream_link_id = self.fd.link_to_flow_receiving_node[
            self.fd.downstream_node_at_link()
        ]
        upstream_contributing_links_at_node = np.where(
            self.fd.flow_link_incoming_at_node() == 1, self._grid.links_at_node, -1
        )

        # Update the node topographic elevations depending on the quantity of stored sediment
        for n in range(self._grid.number_of_nodes):

            if number_of_contributors[n] > 0:  # we don't update head node elevations

                upstream_links = upstream_contributing_links_at_node[n]
                real_upstream_links = upstream_links[upstream_links != BAD_INDEX_VALUE]
                width_of_upstream_links = self._grid.at_link["channel_width"][
                    real_upstream_links
                ]
                length_of_upstream_links = self._grid.length_of_link[
                    real_upstream_links
                ]

                width_of_downstream_link = self._grid.at_link["channel_width"][
                    downstream_link_id
                ][n]
                length_of_downstream_link = self._grid.length_of_link[
                    downstream_link_id
                ][n]

                if (
                    downstream_link_id[n] == BAD_INDEX_VALUE
                ):  # I'm sure there's a better way to do this, but...
                    length_of_downstream_link = 0

                alluvium__depth = _calculate_alluvium_depth(
                        self.vol_stor[downstream_link_id][n],
                        width_of_upstream_links,
                        length_of_upstream_links,
                        width_of_downstream_link,
                        length_of_downstream_link, 
                        self.bed_porosity
                        )


                #                print("alluvium depth = ",alluvium__depth)
                #                print("Volume stored at n = ",n,"=",self.vol_stor[downstream_link_id][n])
                #                print("Denomenator",np.sum(width_of_upstream_links * length_of_upstream_links) + width_of_downstream_link * length_of_downstream_link)
                #
                self._grid.at_node["topographic__elevation"][n] = (
                    self._grid.at_node["bedrock__elevation"][n] + alluvium__depth
                )

    def _calc_transport_wilcock_crowe(self):
        """Method to determine the transport time for each parcel in the active
        layer using a sediment transport equation.

        Note: could have options here (e.g. Wilcock and Crowe, FLVB, MPM, etc)
        """
        # parcel attribute arrays from DataRecord

        Darray = self._parcels.dataset.D[:, self._time_idx]
        Activearray = self._parcels.dataset.active_layer[:, self._time_idx].values
        Rhoarray = self._parcels.dataset.density.values
        Volarray = self._parcels.dataset.volume[:, self._time_idx].values
        Linkarray = self._parcels.dataset.element_id[
            :, self._time_idx
        ].values  # link that the parcel is currently in

        R = (Rhoarray - self.fluid_density) / self.fluid_density

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
        self.Ttimearray = np.zeros(self._num_parcels)
        # ^ Ttimearray is the time to move through the entire length of a link
        self.pvelocity = np.zeros(self._num_parcels)
        # ^ pvelocity is the parcel virtual velocity = link length / link travel time

        # Calculate bed statistics for all of the links
        vol_tot = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._find_now
        )
        vol_tot[np.isnan(vol_tot) == 1] = 0

        vol_act = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._active_parcel_records
        )
        vol_act[np.isnan(vol_act) == 1] = 0

        # find active sand.
        findactivesand = (
            self._parcels.dataset.D < 0.002
        ) * self._active_parcel_records  # since find active already sets all prior timesteps to False, we can use D for all timesteps here.

        if np.any(findactivesand):
            # print("there's active sand!")
            vol_act_sand = self._parcels.calc_aggregate_value(
                np.sum, "volume", at="link", filter_array=findactivesand
            )
            vol_act_sand[np.isnan(vol_act_sand) == True] = 0
        else:
            vol_act_sand = np.zeros(self._grid.number_of_links)

        frac_sand = np.zeros_like(vol_act)
        frac_sand[vol_act != 0] = vol_act_sand[vol_act != 0] / vol_act[vol_act != 0]
        frac_sand[np.isnan(frac_sand) == True] = 0

        # Calc attributes for each link, map to parcel arrays
        for i in range(self._grid.number_of_links):

            active_here = np.where(np.logical_and(Linkarray == i, Activearray == 1))[0]
            d_act_i = Darray[active_here]
            vol_act_i = Volarray[active_here]
            rhos_act_i = Rhoarray[active_here]
            vol_act_tot_i = np.sum(vol_act_i)
            # ^ this behaves as expected. filterarray to create vol_tot above does not. --> FIXED?
            self.d_mean_active[i] = np.sum(d_act_i * vol_act_i) / (
                vol_act_tot_i
            )
            if vol_act_tot_i >0:
                self.rhos_mean_active[i] = np.sum(rhos_act_i * vol_act_i) / (
                    vol_act_tot_i
                )
            else:
                self.rhos_mean_active[i] = np.nan

            D_mean_activearray[Linkarray == i] = self.d_mean_active[i]
            frac_sand_array[Linkarray == i] = frac_sand[i]
            vol_act_array[Linkarray == i] = vol_act[i]
            Sarray[Linkarray == i] = self._grid.at_link["channel_slope"][i]
            Harray[Linkarray == i] = self.flow_depth[self._time_idx, i]
            Larray[Linkarray == i] = self._grid.at_link["link_length"][i]
            active_layer_thickness_array[Linkarray == i] = (
                    self.active_layer_thickness[i]
                    )

        Sarray = np.squeeze(Sarray)
        Harray = np.squeeze(Harray)
        Larray = np.squeeze(Larray)
        frac_sand_array = np.squeeze(frac_sand_array)

        # Wilcock and crowe claculate transport for all parcels (active and inactive)
        taursg = _calculate_reference_shear_stress(
                self.fluid_density,
                R,
                self.g,
                D_mean_activearray,
                frac_sand_array
                )

        #        print("d_mean_active = ", d_mean_active)
        #        print("taursg = ", taursg)

        # frac_parcel should be the fraction of parcel volume in the active layer volume
        # frac_parcel = vol_act_array / Volarray
        # ^ This is not a fraction
        # Instead I think it should be this but CHECK CHECK
        frac_parcel = np.nan * np.zeros_like(Volarray)
        frac_parcel[vol_act_array != 0] = (
            Volarray[vol_act_array != 0] / vol_act_array[vol_act_array != 0]
        )

        b = 0.67 / (1 + np.exp(1.5 - Darray / D_mean_activearray))

        tau = (self.fluid_density
               * self.g
               * Harray
               * Sarray
               )
        tau = np.atleast_1d(tau)

        taur = taursg * (Darray / D_mean_activearray) ** b
        tautaur = tau / taur
        tautaur_cplx = tautaur.astype(np.complex128)
        # ^ work around needed b/c np fails with non-integer powers of negative numbers

        W = 0.002 * np.power(tautaur_cplx.real, 7.5)
        W[tautaur >= 1.35] = 14 * np.power(
            (1 - (0.894 / np.sqrt(tautaur_cplx.real[tautaur >= 1.35]))), 4.5
        )
        W = W.real

        # compute parcel virtual velocity, m/s
        self.pvelocity[Activearray == 1] = (
            W[Activearray == 1]
            * (tau[Activearray == 1] ** (3 / 2))
            * frac_parcel[Activearray == 1]
            / (self.fluid_density ** (3 / 2))
            / self.g
            / R[Activearray == 1]
            / active_layer_thickness_array[Activearray == 1]
        )

        self.pvelocity[np.isnan(self.pvelocity)] = 0

        # Assign those things to the grid -- might be useful for plotting later...?
        self._grid.at_link["sediment_total_volume"] = vol_tot
        self._grid.at_link["sediment__active__volume"] = vol_act
        self._grid.at_link["sediment__active__sand_fraction"] = frac_sand

    def _move_parcel_downstream(self, dt):  # Jon
        """Method to update parcel location for each parcel in the active
        layer.
        """

        # we need to make sure we are pointing to the array rather than making copies
        current_link = self._parcels.dataset.element_id[
            :, self._time_idx
        ]  # same as Linkarray, this will be updated below
        location_in_link = self._parcels.dataset.location_in_link[
            :, self._time_idx
        ]  # updated below
        distance_to_travel_this_timestep = (
            self.pvelocity * dt
        )  # total distance traveled in dt at parcel virtual velocity
        # ^ movement in current and any DS links at this dt is at the same velocity as in the current link
        # ... perhaps modify in the future(?) or ensure this type of travel is kept to a minimum
        # ... or display warnings or create a log file when the parcel jumps far in the next DS link

        # print("distance traveled = ", distance_to_travel_this_timestep)

        #        if self._time_idx == 1:
        #            print("t", self._time_idx)

        for p in range(self._parcels.number_of_items):

            distance_to_exit_current_link = self._grid.at_link["link_length"][
                int(current_link[p])
            ] * (1 - location_in_link[p])

            # initial distance already within current link
            distance_within_current_link = self._grid.at_link["link_length"][
                int(current_link[p])
            ] * (location_in_link[p])

            running_travel_distance_in_dt = 0  # initialize to 0

            distance_left_to_travel = distance_to_travel_this_timestep[p]
            # if parcel in network at end of last timestep

            if self._parcels.dataset.element_id[p, self._time_idx] != _OUT_OF_NETWORK:
                # calc travel distances for all parcels on the network in this timestep

                # distance remaining before leaving current link

                while (
                    running_travel_distance_in_dt + distance_to_exit_current_link
                ) <= distance_to_travel_this_timestep[p]:
                    # distance_left_to_travel > 0:
                    # ^ loop through until you find the link the parcel will reside in after moving
                    # ... the total travel distance

                    # update running travel distance now that you know the parcel will move through the
                    # ... current link
                    running_travel_distance_in_dt = (
                        running_travel_distance_in_dt + distance_to_exit_current_link
                    )

                    # now in DS link so this is reset
                    distance_within_current_link = 0

                    # determine downstream link
                    downstream_link_id = self.fd.link_to_flow_receiving_node[
                        self.fd.downstream_node_at_link()[int(current_link[p])]
                    ]

                    # update current link to the next link DS
                    current_link[p] = downstream_link_id

                    if downstream_link_id == -1:  # parcel has exited the network
                        # (downstream_link_id == -1) and (distance_left_to_travel <= 0):  # parcel has exited the network
                        current_link[p] = _OUT_OF_NETWORK  # overwrite current link

                        # Keep parcel in data record but update its attributes so it is no longer accessed.
                        # Moving parcels into a separate exit array seems to be too computationally expensive.
                        # Probably worthwhile to update the following upon exit:
                        # parcels.dataset.element_id
                        # parcels.dataset.D
                        # parcels.dataset.volume
                        # and compute sub-dt time of exit

                        break  # break out of while loop

                    # ARRIVAL TIME in this link ("current_link") =
                    # (running_travel_distance_in_dt[p] / distance_to_travel_this_timestep[p]) * dt + "t" running time
                    # ^ DANGER DANGER ... if implemented make sure "t" running time + a fraction of dt
                    # ... correctly steps through time.

                    distance_to_exit_current_link = self._grid.at_link["link_length"][
                        int(current_link[p])
                    ]

                    distance_left_to_travel -= distance_to_exit_current_link

                # At this point, we have progressed to the link where the parcel will reside after dt
                distance_to_resting_in_link = (
                    distance_within_current_link  # zero if parcel in DS link
                    + distance_to_travel_this_timestep[p]
                    - running_travel_distance_in_dt  # zero if parcel in same link
                )

                # update location in current link
                if current_link[p] == _OUT_OF_NETWORK:
                    location_in_link[p] = np.nan

                else:
                    location_in_link[p] = (
                        distance_to_resting_in_link
                        / self._grid.at_link["link_length"][int(current_link[p])]
                    )

                # reduce D and volume due to abrasion
                vol = _calculate_parcel_volume_post_abrasion(
                        self._parcels.dataset.volume[p, self._time_idx],
                        distance_to_travel_this_timestep[p],
                        self._parcels.dataset.abrasion_rate[p]
                        )

                D = _calculate_parcel_grain_diameter_post_abrasion(
                        self._parcels.dataset.D[p, self._time_idx],
                        self._parcels.dataset.volume[p, self._time_idx],
                        vol
                        )
   
                # update parcel attributes
                self._parcels.dataset.location_in_link[
                    p, self._time_idx
                ] = location_in_link[p]
                self._parcels.dataset.element_id[p, self._time_idx] = current_link[p]
                self._parcels.dataset.active_layer[p, self._time_idx] = 1
                # ^ reset to 1 (active) to be recomputed/determined at next timestep

                # Jon -- I suggest we do this after the fact when plotting to reduce model runtime:
                # calculate the x and y value of each parcel at this time (based on squiggly shape file)
                # could also create a function that calculates the x and y value for all parcels at all time
                # that is just called once at the end of running the model.

                # self._parcels.dataset["x"] = x_value
                # self._parcels.dataset["y"] = y_value

                self._parcels.dataset.D[p, self._time_idx] = D
                self._parcels.dataset.volume[p, self._time_idx] = vol

    def run_one_step(self, dt):
        """Run NetworkSedimentTransporter forward in time.

        When the NetworkSedimentTransporter runs forward in time the following
        steps occur:

            1. A new set of records is created in the Parcels that cooreponds to the new time
            2. If parcels are remain on the network then:
                a. Active parcels are identifed based on entrainment critera.
                b. Effective bed slope is calculated based on inactive parcel volumes
                c. Transport rate is calculated...
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
            self._update_channel_slopes()  # I moved this down and commented out the second call to 'update channel slopes...'
            self._calc_transport_wilcock_crowe()
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
    >>> with pytest.raises(ValueError):
    ...     _recalculate_channel_slope(0., 10., 10.)

    """
    chan_slope = (z_up - z_down) / dx

    if chan_slope < 0.0:
        raise ValueError("NST Channel Slope Negative")

    if chan_slope < threshold:
        chan_slope = threshold

    return chan_slope


def _calculate_alluvium_depth(
                        stored_volume,
                        width_of_upstream_links,
                        length_of_upstream_links,
                        width_of_downstream_link,
                        length_of_downstream_link,
                        porosity
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
    >>> _calculate_alluvium_depth(100,np.array([0.5,1]),np.array([10,10]) 1, 10, 0.2)
    10.0
    >>>_calculate_alluvium_depth(24,np.array([0.1,3]),np.array([10,10]), 1, 1, 0.5)
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
        /(1-porosity)
    )
    # NOTE: Jon, porosity was left out in earlier version of the LL component,
    # but it seems it should be in here. Check me: is the eqn correct?

    if alluvium__depth < 0.0:
        raise ValueError("NST Alluvium Depth Negative")

    return alluvium__depth


def _calculate_reference_shear_stress(
                        fluid_density, 
                        R, 
                        g,
                        mean_active_grain_size,
                        frac_sand
                        ):
    """Calculate reference shields stress (taursg) using the sand content of
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
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter _calculate_reference_shear_stress
    >>> import pytest
    >>> _calculate_reference_shear_stress(1, 1, 1, 1, 0)
    0.036
    # Katy: I'm running in to a float representation issue here. The python 
    # result is 0.036000000000000004 , should be 0.036
    >>>_calculate_reference_shear_stress(1000, 1.65, 9.8, 0.1, 0.9)
    33.957000369403168

    """

    taursg = (
        fluid_density
        * R
        * g
        * mean_active_grain_size
        * (0.021 + 0.015 * np.exp(-20.0 * frac_sand))
    )

    if any(n < 0 for n in taursg):
        raise ValueError("NST reference Shields stress is negative")
        
    return taursg

def _calculate_parcel_volume_post_abrasion(
                        starting_volume,
                        travel_distance,
                        abrasion_rate
                        ):
    """Calculate parcel volumes after abrasion, according to Sternberg 
    exponential abrasion.

    Parameters
    ----------
    starting_volume : float
        Starting volume of each parcel.
    travel_distance: float
        Travel distance for each parcel during this timestep, in ___.
    abrasion_rate: float
        Mean grain size of the 'active' sediment parcels.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_parcel_volume_post_abrasion
    >>> import pytest
    >>> _calculate_parcel_volume_post_abrasion(10,100,0.003)
    7.4081822068171785
    >>>_calculate_parcel_volume_post_abrasion(10,300,0.1)
    9.3576229688401746e-13
    >>> with pytest.raises(ValueError):
    ...     _calculate_parcel_volume_post_abrasion(10,300,-3)

    """

    volume = starting_volume * np.exp(travel_distance * (-abrasion_rate))

    if volume > starting_volume:
        raise ValueError("NST parcel volume *increases* due to abrasion")

    return volume


def _calculate_parcel_grain_diameter_post_abrasion(
                        starting_diameter,
                        pre_abrasion_volume,
                        post_abrasion_volume
                        ):
    """Calculate parcel grain diameters after abrasion, according to Sternberg 
    exponential abrasion.

    Parameters
    ----------
    starting_diameter : float
        Starting volume of each parcel.
    pre_abrasion_volume: float
        Parcel volume before abrasion. 
    post_abrasion_volume: float
        Parcel volume after abrasion. 

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_parcel_grain_diameter_post_abrasion
    >>> import pytest
    >>> _calculate_parcel_grain_diameter_post_abrasion(10,1,1)
    10.0
    >>>_calculate_parcel_grain_diameter_post_abrasion(10,1,2)
    0.0001

    """

    abraded_grain_diameter = (starting_diameter
            * (post_abrasion_volume
            / pre_abrasion_volume
            ) ** (1 / 3)
            )

    return abraded_grain_diameter
