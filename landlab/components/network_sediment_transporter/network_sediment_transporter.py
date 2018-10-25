
#!/usr/env/python

"""Landlab component that simulates xxxxxx

maybe it should be called "czuba_network_sediment_transporter"

info about the component here

.. codeauthor:: Jon Allison Katy

Created on Tu May 8, 2018
Last edit ---
"""

# %% Import Libraries
from landlab import Component
from landlab.utils.decorators import use_file_name_or_kwds
import numpy as np
from landlab import BAD_INDEX_VALUE

# %% Instantiate Object

class NetworkSedimentTransporter(Component):
    """Network bedload morphodynamic component.

    Landlab component designed to calculate _____.
    info info info

    **Usage:**
    Option 1 - Basic::
        NetworkSedimentTransporter(grid,
                             parcels,
                             parcel_reservoir,
                             transport_method = 'WilcockCrow',
                             transporter = asdfasdf
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
    _name = 'NetworkSedimentTransporter'
    __version__ = '1.0'

    # component requires these values to do its calculation, get from driver
    _input_var_names = (
        'topographic__elevation',    
        'channel_slope',
        'link_length',
        'channel_width',
        'flow_depth',
        'channel_discharge'
        )

    #  component creates these output values
    _output_var_names = (
        'thing',
        'thing',
        )

    # units for each parameter and output
    _var_units = {
        'topographic__elevation': 'm',
        'channel_slope': 'm/m',
        'link_length':'m',
        'channel_width':'m',
        'flow_depth': 'm',
        'channel_discharge': 'm3'
        }

    # grid centering of each field and variable
    _var_mapping = {
        'topographic__elevation': 'node',
        'channel_slope': 'link',
        'link_length':'link',
        'channel_width':'link',
        'flow_depth': 'link',
        'channel_discharge': 'link'        
        }

    # short description of each field
    _var_doc = {
        'surface_water__depth': 'Depth of streamflow on the surface',
        'topographic__elevation': 'Topographic elevation at that node',
        'channel_slope': 'Slope of the river channel through each reach',
        'link_length':'Length of each reach',
        'channel_width':'Flow width of the channel, assuming constant width',
        'flow_depth': 'Depth of stream flow in each reach',
        'channel_discharge': 'Stream channel discharge in each reach'         
        }

    # Run Component
    @use_file_name_or_kwds
    def __init__(self,
                 grid,
                 parcels,
                 parcel_reservoir,
                 flow_director,
                 flow_depth,
                 active_layer_thickness,
                 bed_porosity,
                 discharge='surface_water__discharge',
                 channel_width='channel_width',
                 transport_method = 'WilcockCrowe',
                 **kwds):
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
        self._parcels = parcels
        self._parcel_reservoir = parcel_reservoir


        # assert that the flow director is a component and is of type
        # FlowDirectorSteepest

        if not isinstance(flow_director, Component):
            msg = ("NetworkSedimentTransporter: value passed as flow_director "
                   "is not a Landlab component.")
            raise(ValueError, msg)

        if flow_director.__name__ is not "FlowDirectorSteepest":
            msg = ("NetworkSedimentTransporter: flow_director must be "
                   "FlowDirectorSteepest.")
            raise(ValueError, msg)

        # save reference to flow director
        self.fd = flow_director
        self.flow_depth = flow_depth
        self.bed_porosity = bed_porosity

        self.transport_method = transport_method # self.transport_method makes it a class variable, that can be accessed within any method within this class
        if self.transport_method =="WilcockCrowe":
            self.update_transport_time = self._calc_transport_wilcock_crowe

        # save reference to discharge and width fields stored at-link on the
        # grid
        self._discharge = self._grid.at_link[discharge]
        self._width = self._grid.at_link[channel_width]

        # create field for channel slope if it doesnt exist yet.
        if "channel_slope" not in self._grid.at_link():
            self._channel_slope = self._grid.zeros(at='node')
            self._update_channel_slopes() # this would be a function that updates the channel slopes (if you use this more than once, make it a function)
        else:
            self._channel_slope = self._grid.at_link["channel_slope"]


    def _update_channel_slopes(self):
        """text"""
        # Katy think this can be vectorized
        for l in range(self._grid.number_of_links):

            upstream_node_id = self.fd.upstream_node_at_link[l]
            downstream_node_id = self.fd.downstream_node_at_link[l]

            chan_slope = ((self._grid.at_node['topographic__elevation'][upstream_node_id]
                - self._grid.at_node['topographic__elevation'][downstream_node_id])
                / self._grid.length_of_link[l])

            if chan_slope < 1e-4:
                chan_slope = 1e-4

            self._channel_slope[l] = chan_slope

    def _partition_active_and_storage_layers(self, **kwds): # Allison is working on this
        """For each parcel in the network, determines whether it is in the
        active or storage layer during this timestep, then updates node elevations
        """
# %%
        vol_tot =  self._parcels.calc_aggregate_value(np.sum,'volume',at = 'link')

        capacity = 2* np.ones(np.size(self._parcels.DataFrame.element_id)) # REPLACE with real calculation for capacity

        for i in range(self._grid.number_of_links):

            if vol_tot[i]>0: #only do this check capacity if parcels are in link

                #First In Last Out
                parcel_id_thislink = np.where(self._parcels.DataFrame.element_id.values == i)[0]

                time_arrival_sort = np.flip(np.argsort(self._parcels.DataFrame.time_arrival_in_link.values[parcel_id_thislink]),0)
                parcel_id_time_sorted = parcel_id_thislink[time_arrival_sort]

                cumvol = np.cumsum(self._parcels.DataFrame.volume.values[parcel_id_time_sorted])

                idxinactive = np.where(cumvol>capacity[i])
                make_inactive = parcel_id_time_sorted[idxinactive]

                self._parcels.set_value(item_id = parcel_id_thislink,variable = 'active_layer', value = 1)
                self._parcels.set_value(item_id = make_inactive,variable = 'active_layer', value = 0)

        # Update Node Elevations
        findactive = self._parcels.DataFrame['active_layer']==1 # filter for only parcels in active layer
        vol_act = self._parcels.calc_aggregate_value(np.sum,
                                                  'volume',at = 'link',
                                                  filter_array = findactive)
        self.vol_stor = (vol_tot-vol_act)/(1-self.bed_porosity)
        # ^ Jon-- what's the rationale behind only calculating new node elevations
        # using the storage volume (rather than the full sediment volume)?
        # Jon response -- because during transport, that is the sediment defining the immobile
        # substrate that is setting the slope. Parcels that are actively transporting should not
        # affect the slope because they are moving.


# %%
    def _adjust_node_elevation(self, seed=0): # Allison is working on this as of July 23, 2018
        """Adjusts slope for each link based on parcel motions from last
        timestep and additions from this timestep.
        """
        number_of_contributors = np.sum(self.fd.flow__link_incoming_at_node == 1, axis=1)
        downstream_link_id = self.fd.link_to_flow_receiving_node[self.fd.downstream_node_at_link]
        upstream_contributing_links_at_node = np.where(self.fd.flow__link_incoming_at_node == 1, self._grid.links_at_node, -1)

        # Update the node elevations depending on the quantity of stored sediment
        for l in range(self._grid.number_of_nodes):

            if number_of_contributors[l] > 0: #we don't update head node elevations

                upstream_links = upstream_contributing_links_at_node[l]
                real_upstream_links = upstream_links[upstream_links != BAD_INDEX_VALUE]
                width_of_upstream_links = self._grid.at_link['channel_width'][real_upstream_links]
                length_of_upstream_links = self._grid.length_of_link[real_upstream_links]

                width_of_downstream_link = self._grid.at_link['channel_width'][downstream_link_id][l]
                length_of_downstream_link = self._grid.length_of_link[downstream_link_id][l]

                if downstream_link_id[l] == BAD_INDEX_VALUE: # I'm sure there's a better way to do this, but...
                    length_of_downstream_link = 0

                # IMPROVE: deal with the downstream most link...
                elev_change = 2*self.vol_stor[downstream_link_id][l]/(np.sum(width_of_upstream_links * length_of_upstream_links)
                    + width_of_downstream_link * length_of_downstream_link)

                self._grid.at_node['topographic__elevation'][l] += elev_change

        # Update channel slope
        self._update_channel_slopes()

# %%
    def _calc_transport_wilcock_crowe(self, flow_depth): # Allison
        """Method to determine the transport time for each parcel in the active
        layer using a sediment transport equation.

        Note: could have options here (e.g. Wilcock and Crowe, FLVB, MPM, etc)
        """

# %%
        # parcel attribute arrays from ItemCollector

        # another way of doing this --> check to see if this is copying. we don't want to be copying
        Darray = self._parcels.DataFrame.D.values

#        Darray = np.array(parcels.DataFrame.D,copy=False) # this gives a copy, but we can set copy to false..?
        Activearray = self._parcels.DataFrame.active_layer.values
        Rhoarray = self._parcels.DataFrame.density.values
        Volarray = self._parcels.DataFrame.volume.values
        Linkarray = self._parcels.DataFrame.element_id.values #link that the parcel is currently in
        rho = 1000 # kg/m3
        R = (Rhoarray-rho)/rho

        # parcel attribute arrays to populate below
        frac_sand_array = np.zeros(np.size(self._parcels.DataFrame.element_id))
        vol_act_array = np.zeros(np.size(self._parcels.DataFrame.element_id))
        Sarray = np.zeros(np.size(self._parcels.DataFrame.element_id))
        Harray = np.zeros(np.size(self._parcels.DataFrame.element_id))
        Larray = np.zeros(np.size(self._parcels.DataFrame.element_id))
        d_mean_active = np.zeros(np.size(self._parcels.DataFrame.element_id))
        d_mean_active.fill(np.nan)
        Ttimearray = np.zeros(np.size(self._parcels.DataFrame.element_id))

        # Calculate bed statistics for all of the links
        vol_tot =  self._parcels.calc_aggregate_value(np.sum,'volume',at = 'link')

        findactive = self._parcels.DataFrame['active_layer']==1 # filter for only parcels in active layer
        vol_act = self._parcels.calc_aggregate_value(np.sum,
                                                  'volume',at = 'link',
                                                  filter_array = findactive)

        findactivesand = np.logical_and(Darray<0.002,Activearray ==1)
        vol_act_sand = self._parcels.calc_aggregate_value(np.sum,
                                                'volume',at = 'link',
                                                filter_array = findactivesand)

        vol_act_sand[np.isnan(vol_act_sand)==True] = 0
        frac_sand = vol_act_sand/vol_act
        frac_sand[np.isnan(frac_sand)== True] = 0

        # Calc attributes for each link, map to parcel arrays
        for i in range(self._grid.number_of_links):

            active_here = np.where(np.logical_and(Linkarray == i,Activearray == 1))[0]

            d_act_i = Darray[active_here]
            vol_act_i = Volarray[active_here]
            d_mean_active[Linkarray == i] = np.sum(d_act_i * vol_act_i)/(vol_act[i])

            frac_sand_array[Linkarray == i] = frac_sand[i]
            vol_act_array[Linkarray == i] = vol_act[i]
            Sarray[Linkarray == i] = self._grid.at_link['channel_slope'][i]
            Harray[Linkarray == i] = self.flow_depth[i]
            Larray[Linkarray == i] = self._grid.at_link['link_length'][i]

        # Wilcock and crowe claculate transport for all parcels (active and inactive)
        taursg = rho * R * g * d_mean_active * (0.021 + 0.015*np.exp(-20.*frac_sand_array))
        frac_parcel = vol_act_array/Volarray;
        b = 0.67 / (1 + np.exp(1.5 - Darray/ d_mean_active))
        tau = rho * g * Harray * Sarray
        taur = taursg * (Darray / d_mean_active) **b
        tautaur = tau / taur
        tautaur_cplx = tautaur.astype(np.complex128)
        # ^ work around needed b/c np fails with non-integer powers of negative numbers
        W = 14 * np.power(( 1 - (0.894/np.sqrt(tautaur_cplx))), 4.5)
        W[tautaur_cplx<1.35] = 0.002 * np.power(tautaur[tautaur_cplx<1.35],7.5)
        W = W.real

        # assign travel times only where active==1
        Ttimearray[Activearray==1] = rho**(3/2)*R[Activearray==1]*g*Larray[Activearray==1]*active_layer_thickness/W[Activearray==1]/tau[Activearray==1]**(3/2)/(1-frac_sand_array[Activearray==1])/frac_parcel[Activearray==1]
        Ttimearray[findactivesand==True] = rho**(3/2)*R[findactivesand==True]*g*Larray[findactivesand==True]*active_layer_thickness/W[findactivesand==True]/tau[findactivesand==True]**(3/2)/frac_sand_array[findactivesand==True]
        # ^ why?? if k = 1 ---> if it's sand...?  ASK JON about the logic here...

        # Assign those things to the grid -- might be useful for plotting later...?
        self._grid.at_link['sediment_total_volume'] = vol_tot
        self._grid.at_link['sediment__active__volume'] = vol_act
        self._grid.at_link['sediment__active__sand_fraction'] = frac_sand


# %%
    def _move_parcel_downstream(self, i):    # Jon
        """Method to update parcel location for each parcel in the active
        layer.
        """

# %%
        # we need to make sure we are pointing to the array rather than making copies
        current_link = self._parcels.DataFrame.element_id.values # same as Linkarray, this will be updated below
        location_in_link = self._parcels.DataFrame.location_in_link.values # updated below
        # Ttimearray -- needs to be brought in here
        #^ Ttimearray is the time to move through the entire length of a link

        # However, a parcel is not always at the US end of a link, so need to determine
        # how much time it takes for that parcel to move out of the current link based on its
        # current location ...
        time_to_exit_current_link = Ttimearray*(1-location_in_link)
        running_travel_time_in_dt = time_to_exit_current_link


        for p in range(self._parcels.number_of_items):
            #^ loop through all parcels, this loop could probably be removed in future refinements

            # ... and compare to the timestep dt
            # loop through until you find the link the parcel will reside in after dt
            while running_travel_time_in_dt[p] <= dt :
                # determine downstream link
                downstream_link_id = fd.link_to_flow_receiving_node[fd.downstream_node_at_link[element_id[p]]]

                if downstream_link_id == -1 : #parcel has exited the network
                    # I think we should then remove this parcel from the parcel item collector
                    # if so, we manipulate the exiting parcel here, but may want to note something about its exit
                    # such as output volume and output time into a separate outlet array
                    
                    
                    # ADD CODE FOR THIS HERE, right now these parcel will just cycle through not actually leaving the system
                    
                    
                    break # break out of while loop

                current_link[p] = downstream_link_id
                location_in_link[p] = 0 # place parcel at upstream end of DS link
                # ARRIVAL TIME in this link ("current_link") is equal to "t" running time + "running_travel_time_in_dt"

                # movement in DS link is at the same velocity as in US link
                # perhaps modify in future or ensure this type of travel is kept to a minimum by
                # dt < travel time
                time_to_exit_current_link[p] = time_to_exit_current_link[p] / self._grid.at_link['link_length'][element_id[p]] * self._grid.at_link['link_length'][current_link[p]]
                running_travel_time_in_dt[p] = running_travel_time_in_dt[p] + time_to_exit_current_link[p]

                # TRACK RUNNING TRAVEL DISTANCE HERE SIMILAR TO RUNNING TRAVEL TIME

            time_in_link_before_dt = time_to_exit_current_link[p] - (running_travel_time_in_dt[p] - dt)
            #^ if in same link, this equals dt

            # update location in current link
            location_in_link[p] = location_in_link[p] + (time_in_link_before_dt / time_to_exit_current_link[p])

            # USE RUNNING TRAVEL DISTANCE TO UPDATE D AND VOL DUE TO ABRASION HERE
            
            distance_traveled = 0 # NEED TO DEFINE
            
            volume[p] = self._parcels.DataFrame.volume.values[p]*np.exp(distance_traveled * -self._parcels.DataFrame.abrasion_rate.values[p])

            D[p] = 2 * (volume[p] * 3 / (4 * np.pi))**(1/3)

            # update parcel attributes
            self._parcels.DataFrame.location_in_link.values[p] = location_in_link[p]
            self._parcels.DataFrame.element_id.values[p] = current_link[p]
            self._parcels.DataFrame.active_layer.values[p] = 1 # reset to 1 (active) to be recomputed/determined at next timestep
            self._parcels.DataFrame.volume.values[p] = volume[p]
            self._parcels.DataFrame.D.values[p] = D[p]


# %%
    def run_one_step(self,dt):
        """stuff"""

        self._partition_active_and_storage_layers(things)
        self._update_node_elevation(things)
        self._calc_transport_wilcock_crowe(things)
        self._move_parcel_downstream(things)