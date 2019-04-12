#!/usr/env/python

"""Landlab component that simulates xxxxxx

maybe it should be called "czuba_network_sediment_transporter"

info about the component here


Fixes that need to happen:

    -- Need to find better way to define filterarrays in time and item_id
	current option is very clunky. The one Nathan Lyons suggested doesn't work.

    -- What to do with parcels when they get to the last link? --> I didn't get to this.

    (x) tau/taur changes very subtly through time. It shouldn't. Track down this mystery.
            - storage volume stays the same.
            - channel slopes (!) change in time, even though parcels don't move to the next link
        --> Proposed solution at Ln 443. Check with Jon. Need to add 'bedrock__elevation' node attribute...

    -- Need to calculate distance a parcel travels in a timestep for abrasion

    !-- The abrasion exponent is applied to diameter, but doesn't impact parcel volume. Need to fix.

    -- Fix inelegant time indexing
    
    -- Looks to me that as part of run-one-step the element_id variable in parces is being changed from 
       An int to a float. I haven't tracked down why... but I think it should stay as an int. 

.. codeauthor:: Jon Allison Katy

Created on Tu May 8, 2018
Last edit ---
"""

import numpy as np

# %% Import Libraries
from landlab import BAD_INDEX_VALUE, Component
from landlab.utils.decorators import use_file_name_or_kwds
from landlab.grid.network import NetworkModelGrid
from landlab.data_record import DataRecord

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
                             channel_geometry,
                             active_layer_thickness)

    Examples
    ----------
    >>> from landlab import RasterModelGrid
    >>> from landlab.component import FlowDirectorSteepest, FlowAccumulator, NetworkSedimentTransporter
    >>> from landlab.components.landslides import LandslideProbability
    >>> import numpy as np

    Do setup of various sorts (grid, width, parcels, discharge, etc)

    Set up NetworkSedimentTransporter

    >>> nst = NetworkSedimentTransporter(grid, stuff)

    Run NetworkSedimentTransporter forward 10 timesteps of size 10 time units.

    >>> for _ in range(10):
    ...     nst.run_one_step(10.)

    Now lets double check we got the right answer

    We'd put code here that showed some results. Our goal here is not so much
    to test the code but to show how it is used and what it does. We would make
    addtional tests in a folder called test that would completely test the code.

    """

    # AP note: this is where the documentation ends and the component starts...
    # QUESTIONS:
    #   what about the variables that aren't on the grid? (e.g. active_layer_thickness, bed_material_porosity,etc?)
    #   what about the variables attached to the DataRecord, rather than the grid?

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
        active_layer_thickness,
        bed_porosity,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
        **kwds
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A raster grid.
        discharge : field name, optional
            Field name of discharge [units] at link.
        and more here...
        """
        super(NetworkSedimentTransporter, self).__init__(grid, **kwds)

        self._grid = grid

        if not isinstance(grid, NetworkModelGrid):
            msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            raise (ValueError, msg)

        self._parcels = parcels

        if not isinstance(parcels, DataRecord):
            msg = (
                "NetworkSedimentTransporter: parcels must be an instance"
                "of DataRecord"
            )
            raise (ValueError, msg)

        self._num_parcels = self._parcels.element_id.size

        self.parcel_attributes = [
            "time_arrival_in_link",
            "active_layer",
            "location_in_link",
            "D",
            "volume",
        ]

        # assert that the flow director is a component and is of type
        # FlowDirectorSteepest

        if not isinstance(flow_director, Component):
            msg = (
                "NetworkSedimentTransporter: value passed as flow_director "
                "is not a Landlab component."
            )
            raise (ValueError, msg)

        if flow_director.name != "FlowDirectorSteepest":
            msg = (
                "NetworkSedimentTransporter: flow_director must be "
                "FlowDirectorSteepest."
            )
            raise (ValueError, msg)

        # save reference to flow director
        self.fd = flow_director
        self.flow_depth = flow_depth
        self.bed_porosity = bed_porosity

        if not 0 <= self.bed_porosity < 1:
            msg = "NetworkSedimentTransporter: bed_porosity must be" "between 0 and 1"
            raise (ValueError, msg)

        self.active_layer_thickness = active_layer_thickness

        # NOTE: variable active_layer_thickness Wong et al 2007
        # "a predictor for active layer thickness that increases with both grain size and Shields number, i.e., equations (50), (17), (18), and (23);

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

        # save reference to discharge and width fields stored at-link on the
        # grid

        self._width = self._grid.at_link[channel_width]

        if "channel_width" not in self._grid.at_link:
            msg = (
                "NetworkSedimentTransporter: channel_width must be assigned"
                "to the grid links"
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
            self._update_channel_slopes()  # this would be a function that updates the channel slopes (if you use this more than once, make it a function)
        else:
            self._channel_slope = self._grid.at_link["channel_slope"]

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
            #            self._parcels['grid_element'].values[:,self._time_idx] = self._parcels[
            #                    'grid_element'].values[:,self._time_idx-1]
            #
            #            self._parcels['element_id'].values[:,self._time_idx] = self._parcels[
            #                    'element_id'].values[:,self._time_idx-1]

            self._parcels.ffill_grid_element_and_id()

            for at in self.parcel_attributes:
                self._parcels[at].values[:, self._time_idx] = self._parcels[at].values[
                    :, self._time_idx - 1
                ]
                
        self._find_now = self._parcels.time == self._time
        self._this_timesteps_parcels = np.zeros_like(self._parcels.element_id, dtype=bool)
        self._this_timesteps_parcels[:, -1] = True

    def _update_channel_slopes(self):
        """text Can be simple-- this is what this does. 'private' functions can
        have very simple examples, explanations. Essentially note to yourself"""
        # Katy think this can be vectorized
        # Jon agrees, but is not sure yet how to do that
        for l in range(self._grid.number_of_links):

            upstream_node_id = self.fd.upstream_node_at_link()[l]
            downstream_node_id = self.fd.downstream_node_at_link()[l]

            chan_slope = (
                self._grid.at_node["topographic__elevation"][upstream_node_id]
                - self._grid.at_node["topographic__elevation"][downstream_node_id]
            ) / self._grid.length_of_link[l]

            if chan_slope < 1e-4:
                chan_slope = 1e-4

            self._channel_slope[l] = chan_slope

        print("channel slopes = ", self._channel_slope)

    def _partition_active_and_storage_layers(
        self, **kwds
    ):  # Allison is working on this
        """For each parcel in the network, determines whether it is in the
        active or storage layer during this timestep, then updates node elevations
        """

        vol_tot = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._this_timesteps_parcels
        )

        # Katy made this change... I think you want this as number of parcels 
        # not the shape of element ID (which will be of size number of parcels 
        # times number of timesteps so far.)
        capacity = 2 * np.ones(self._num_parcels)  # REPLACE with real calculation for capacity

        # this is the old code that katy commented out.
#        capacity = 2 * np.ones(
#            np.size(self._parcels.element_id)
#        )  # REPLACE with real calculation for capacity

        for i in range(self._grid.number_of_links):

            if vol_tot[i] > 0:  # only do this check capacity if parcels are in link

                # First In Last Out.
                # parcel_id_thislink = np.where(self._parcels.DataFrame.element_id.values == i)[0]
                parcel_id_thislink = np.where(self._parcels.element_id[:, self._time_idx] == i)[0]

                # time_arrival_sort = np.flip(np.argsort(self._parcels.DataFrame.time_arrival_in_link.values[parcel_id_thislink]),0)
                # time_arrival_sort = np.flip(np.argsort(self._parcels['time_arrival_in_link'][parcel_id_thislink]),0)

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
                    self._parcels.volume[parcel_id_time_sorted, self._time_idx]
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
        self._active_parcel_records = self._parcels.active_layer == _ACTIVE * self._this_timesteps_parcels

        vol_act = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._active_parcel_records
        )
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

        # Update the node elevations depending on the quantity of stored sediment
        for l in range(self._grid.number_of_nodes):
            # ^ comment from Jon -- I usually don't like using l as an index because it looks too close to 1.

            if number_of_contributors[l] > 0:  # we don't update head node elevations

                upstream_links = upstream_contributing_links_at_node[l]
                real_upstream_links = upstream_links[upstream_links != BAD_INDEX_VALUE]
                width_of_upstream_links = self._grid.at_link["channel_width"][
                    real_upstream_links
                ]
                length_of_upstream_links = self._grid.length_of_link[
                    real_upstream_links
                ]

                width_of_downstream_link = self._grid.at_link["channel_width"][
                    downstream_link_id
                ][l]
                length_of_downstream_link = self._grid.length_of_link[
                    downstream_link_id
                ][l]

                if (
                    downstream_link_id[l] == BAD_INDEX_VALUE
                ):  # I'm sure there's a better way to do this, but...
                    length_of_downstream_link = 0

                # IMPROVE: deal with the downstream most link...
                # DANGER DANGER I think something is wrong here. Anytime there is volume in storage, the elevation changes, even if the storage volume hasn't changed.

                #                elev_change = (
                #                    2
                #                    * self.vol_stor[downstream_link_id][l]
                #                    / (
                #                        np.sum(width_of_upstream_links * length_of_upstream_links)
                #                        + width_of_downstream_link * length_of_downstream_link
                #                    )
                #                )
                #
                #                self._grid.at_node["topographic__elevation"][l] += elev_change
                #
                #
                # IDEA: Each node has a 'bedrock elevation' that is set at the start, then we can add the bed thickness to that.
                # perhaps ^ that is what you were modeling originally, but that's not how I had translated it originally.

                alluvium__depth = (
                    2
                    * self.vol_stor[downstream_link_id][l]
                    / (
                        np.sum(width_of_upstream_links * length_of_upstream_links)
                        + width_of_downstream_link * length_of_downstream_link
                    )
                )

                self._grid.at_node["topographic__elevation"][l] = (
                    self._grid.at_node["bedrock__elevation"][l] + alluvium__depth
                )

    def _calc_transport_wilcock_crowe(self):  # Allison
        """Method to determine the transport time for each parcel in the active
        layer using a sediment transport equation.

        Note: could have options here (e.g. Wilcock and Crowe, FLVB, MPM, etc)
        """
        # parcel attribute arrays from ItemCollector

        # another way of doing this --> check to see if this is copying. we don't want to be copying
        Darray = self._parcels.D[:, self._time_idx]

        #        Darray = np.array(parcels.DataFrame.D,copy=False) # this gives a copy, but we can set copy to false..?
        Activearray = self._parcels.active_layer[:, self._time_idx].values
        Rhoarray = self._parcels.density.values
        Volarray = self._parcels.volume[:, self._time_idx].values
        Linkarray = self._parcels.element_id[:, self._time_idx].values  # link that the parcel is currently in
        rho = self.fluid_density
        g = self.g
        R = (Rhoarray - rho) / rho

        # parcel attribute arrays to populate below
        frac_sand_array = np.zeros(self._num_parcels)
        vol_act_array = np.zeros(self._num_parcels)
        Sarray = np.zeros(self._num_parcels)
        Harray = np.zeros(self._num_parcels)
        Larray = np.zeros(self._num_parcels)
        d_mean_active = np.zeros(self._num_parcels)
        d_mean_active.fill(np.nan)
        self.Ttimearray = np.zeros(self._num_parcels)
        # ^ Ttimearray is the time to move through the entire length of a link

        # Calculate bed statistics for all of the links
        vol_tot = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._find_now
        )
        
        vol_act = self._parcels.calc_aggregate_value(
            np.sum, "volume", at="link", filter_array=self._active_parcel_records
        )

        # find active sand.
        findactivesand = (self._parcels.D < 0.002) * self._active_parcel_records # since find active already sets all prior timesteps to False, we can use D for all timesteps here. 

        if np.any(findactivesand):
            vol_act_sand = self._parcels.calc_aggregate_value(
                np.sum, "volume", at="link", filter_array=findactivesand
            )
            vol_act_sand[np.isnan(vol_act_sand) == True] = 0
        else:
            vol_act_sand = np.zeros(self._grid.number_of_links)

        frac_sand = vol_act_sand / vol_act
        frac_sand[np.isnan(frac_sand) == True] = 0

        # Calc attributes for each link, map to parcel arrays
        for i in range(self._grid.number_of_links):

            active_here = np.where(np.logical_and(Linkarray == i, Activearray == 1))[0]
            d_act_i = Darray[active_here]
            vol_act_i = Volarray[active_here]
            vol_act_tot_i = np.sum(vol_act_i)
            # ^ this behaves as expected. filterarray to create vol_tot above does not.
            d_mean_active[Linkarray == i] = np.sum(d_act_i * vol_act_i) / (
                vol_act_tot_i
            )

            frac_sand_array[Linkarray == i] = frac_sand[i]
            vol_act_array[Linkarray == i] = vol_act[i]
            Sarray[Linkarray == i] = self._grid.at_link["channel_slope"][i]
            Harray[Linkarray == i] = self.flow_depth[self._time_idx, i]
            Larray[Linkarray == i] = self._grid.at_link["link_length"][i]

        Sarray = np.squeeze(Sarray)
        Harray = np.squeeze(Harray)
        Larray = np.squeeze(Larray)
        frac_sand_array = np.squeeze(frac_sand_array)

        # Wilcock and crowe claculate transport for all parcels (active and inactive)
        taursg = (
            rho
            * R
            * g
            * d_mean_active
            * (0.021 + 0.015 * np.exp(-20.0 * frac_sand_array))
        )

        #        print("d_mean_active = ", d_mean_active)
        #        print("taursg = ", taursg)

        frac_parcel = vol_act_array / Volarray
        b = 0.67 / (1 + np.exp(1.5 - Darray / d_mean_active))
        tau = rho * g * Harray * Sarray
        taur = taursg * (Darray / d_mean_active) ** b

        tautaur = tau / taur
        print("tau / taur = ", tautaur)
        tautaur_cplx = tautaur.astype(np.complex128)
        # ^ work around needed b/c np fails with non-integer powers of negative numbers
        W = 0.002 * np.power(tautaur_cplx.real, 7.5)
        W[tautaur >= 1.35] = 14 * np.power(
            (1 - (0.894 / np.sqrt(tautaur_cplx.real[tautaur >= 1.35]))), 4.5
        )
        W = W.real

        # assign travel times only for active parcels
        self.Ttimearray[Activearray == 1] = (
            (rho ** (3 / 2))
            * g
            * R[Activearray == 1]
            * Larray[Activearray == 1]
            * self.active_layer_thickness
            / W[Activearray == 1]
            / (tau[Activearray == 1] ** (3 / 2))
            / (1 - frac_sand_array[Activearray == 1])
            / frac_parcel[Activearray == 1]
        )

        # self.Ttimearray = np.ones(np.shape(self._parcels["element_id"]))

        # Ttimearray[findactivesand==True] = rho**(3/2)*R[findactivesand==True]*g*Larray[findactivesand==True]*self.active_layer_thickness/W[findactivesand==True]/tau[findactivesand==True]**(3/2)/frac_sand_array[findactivesand==True]
        # ^ why?? if k = 1 ---> if it's sand...?  ASK JON about the logic here...

        # Assign those things to the grid -- might be useful for plotting later...?
        self._grid.at_link["sediment_total_volume"] = vol_tot
        self._grid.at_link["sediment__active__volume"] = vol_act
        self._grid.at_link["sediment__active__sand_fraction"] = frac_sand

    def _move_parcel_downstream(self, dt):  # Jon
        """Method to update parcel location for each parcel in the active
        layer.
        """
        # we need to make sure we are pointing to the array rather than making copies
        current_link = self._parcels.element_id[:, self._time_idx]  # same as Linkarray, this will be updated below
        location_in_link = self._parcels.location_in_link[:, self._time_idx]  # updated below
        distance_traveled = np.zeros(self._num_parcels)
        if self._time_idx == 1:
            print("t", self._time_idx)

        # However, a parcel is not always at the US end of a link, so need to determine
        # how much time it takes for that parcel to move out of the current link based on its
        # current location ...
        time_to_exit_current_link = self.Ttimearray * (1 - location_in_link)
        running_travel_time_in_dt = time_to_exit_current_link

        for p in range(self._parcels.number_of_items):
            # ^ loop through all parcels, this loop could probably be removed in future refinements
            # ... and compare to the timestep dt
            # loop through until you find the link the parcel will reside in after dt
            while running_travel_time_in_dt[p] <= dt:
                # determine downstream link
                current_link_of_parcel = self._parcels.element_id[p, self._time_idx].values

                downstream_link_id = self.fd.link_to_flow_receiving_node[
                    self.fd.downstream_node_at_link()[current_link_of_parcel]
                ]

                if downstream_link_id == -1:  # parcel has exited the network
                    downstream_link_id = _OUT_OF_NETWORK # Katy has added this potential approach. I also added it to an issue.

                    # I think we should then remove this parcel from the parcel item collector
                    # if so, we manipulate the exiting parcel here, but may want to note something about its exit
                    # such as output volume and output time into a separate outlet array
                    # Jon -- Agree a separate outlet accumulator array would be good.
                    # Similarly, if someone wanted to know the history of parcels passing through
                    # a particular link, we could create an accumulator array at that link and
                    # copy but not remove values at that location.

                    # ADD CODE FOR THIS HERE, right now these parcel will just cycle through not actually leaving the system

                    # This is why the code is breaking! These parcels need to be removed.

                    # self.Accumulator_Outlet.item_id = self._parcels.item_id[p]
                    # ^ need something like this but I don't know how to initialize or concatentate
                    # Then we need to remove the parcels from the self._parcel structure.

                    break  # break out of while loop

                current_link[p] = downstream_link_id
                location_in_link[p] = 0  # place parcel at upstream end of DS link
                # ARRIVAL TIME in this link ("current_link") is equal to "t" running time + "running_travel_time_in_dt"

                # movement in DS link is at the same velocity as in US link
                # perhaps modify in future or ensure this type of travel is kept to a minimum by
                # dt < travel time
                time_to_exit_current_link[p] = (
                    time_to_exit_current_link[p]
                    / self._grid.at_link["link_length"][
                        self._parcels.element_id[p, self._time_idx]
                    ]
                    * self._grid.at_link["link_length"][current_link[p]]
                )

                running_travel_time_in_dt[p] = (
                    running_travel_time_in_dt[p] + time_to_exit_current_link[p]
                )

                # TRACK RUNNING TRAVEL DISTANCE HERE SIMILAR TO RUNNING TRAVEL TIME

            time_in_link_before_dt = time_to_exit_current_link[p] - (
                running_travel_time_in_dt[p] - dt
            )
            # ^ if in same link, this equals dt

            # update location in current link
            location_in_link[p] = location_in_link[p] + (
                time_in_link_before_dt / time_to_exit_current_link[p]
            )

            # USE RUNNING TRAVEL DISTANCE TO UPDATE D AND VOL DUE TO ABRASION HERE

            distance_traveled[p] = 0  # NEED TO DEFINE

            # DANGER DANGER... abration_rate is now a diameter loss abration exponent
            D = (self._parcels.D[p, self._time_idx]) * (
                np.exp(distance_traveled[p] * (-self._parcels.abrasion_rate[p]))
            )
            # vol =  # DANGER DANGER... need to make total parcel volume as a function of change in grain size

            # update parcel attributes
            self._parcels.location_in_link[p, self._time_idx] = location_in_link[p]

            # calculate the x and y value of each parcel at this time (based on squiggly shape file)
            # could also create a function that calculates the x and y value for all parcels at all time
            # that is just called once at the end of running the model.
            # self._parcels["x"] = x_value
            # self._parcels["y"] = y)value

            self._parcels.element_id[p, self._time_idx] = current_link[p]
            self._parcels.active_layer[p, self._time_idx] = 1
            # ^ reset to 1 (active) to be recomputed/determined at next timestep
            # DANGER DANGER self._parcels["volume"][p,self._time_idx] = vol
            self._parcels.D[p, self._time_idx] = D

    # %%
    def run_one_step(self, dt):
        """stuff"""

        self._create_new_parcel_time()
        self._partition_active_and_storage_layers()
        self._adjust_node_elevation()
        self._update_channel_slopes()  # I moved this down and commented out the second call to 'update channel slopes...'
        self._calc_transport_wilcock_crowe()
        self._move_parcel_downstream(dt)

        self._time += (
            dt
        )  # DANGER DANGER. AP: should this happen here or before the prior six function calls?
        self._time_idx += 1
