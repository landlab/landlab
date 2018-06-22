
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
import scipy.constants
import copy

# %% Instantiate Object


class NetworkSedimentTransporter(Component):
    """Network bedload morphodynamic component.

    Landlab component designed to calculate _____.
    info info info

    **Usage:**
    Option 1 - Uniform recharge::
        NetworkSedimentTransporter(grid, 
                             bed_parcels,
                             transport_method = 'WilcockCrow',
                             transporter = asdfasdf
                             discharge,
                             channel_geometry,
                             active_layer_thickness)

    Examples
    ----------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.landslides import LandslideProbability
    >>> import numpy as np

    Things here...


    """

# AP note: this is where the documentation ends and the component starts...

    # component name
    _name = 'NetworkSedimentTransporter'
    __version__ = '1.0'
    
    # component requires these values to do its calculation, get from driver
    _input_var_names = (
        'BMI things here',
        'soil__thickness',
        )

    #  component creates these output values
    _output_var_names = (
        'thing',
        'thing',
        )

    # units for each parameter and output
    _var_units = {
        'topographic__specific_contributing_area': 'm',
        'topographic__slope': 'tan theta',
        }

    # grid centering of each field and variable
    _var_mapping = {
        'topographic__specific_contributing_area': 'node',
        'topographic__slope': 'node',
        }

    # short description of each field
    _var_doc = {
        'topographic__specific_contributing_area':
            ('specific contributing (upslope area/cell face )' +
             ' that drains to node'),
        }

    # Run Component
    @use_file_name_or_kwds
    def __init__(self, grid, 
                 parcels,
                 discharge,
                 transport_method = 'WilcockCrowe',
                 channel_width,
                 flow_depth,
                 active_layer_thickness,
                 **kwds):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A raster grid.
        and more here...
        """
        
        self.transport_method = transport_method # self.transport_method makes it a class variable, that can be accessed within any method within this class
        if self.transport_method =="WilcockCrowe":
            self.update_transport_time = self._calc_transport_wilcock_crowe
        
        self.parcels = parcels
        self.grid = grid


        super(NetworkSedimentTransporter, self).__init__(grid, **kwds)
        
        
    def partition_active_and_storage_layers(self, **kwds): # Allison is working on this
        """For each parcel in the network, determines whether it is in the
        active or storage layer during this timestep, then updates node elevations
        """
# %%
        vol_tot =  parcels.calc_aggregate_value(np.sum,'volume',at = 'link')
        
        capacity = 3* np.ones(np.size(element_id)) # REPLACE with real calculation for capacity
        
        for i in range(grid.number_of_links):
            if vol_tot[i]>0: #only do this check capacity if parcels are in link
                
                #First In Last Out
                parcel_id_thislink = np.where(parcels.DataFrame.element_id.values == i)[0]
               
                time_arrival_sort = np.flip(np.argsort(parcels.DataFrame.time_arrival_in_link.values[parcel_id_thislink]),0)
                parcel_id_time_sorted = parcel_id_thislink[time_arrival_sort]
                
                cumvol = np.cumsum(parcels.DataFrame.volume.values[parcel_id_time_sorted])                

                idxinactive = np.where(cumvol>capacity[i])
                make_inactive = parcel_id_time_sorted[idxinactive]
                
                parcels.set_value(item_id = parcel_id_thislink,variable = 'active_layer', value = 1)
                parcels.set_value(item_id = make_inactive,variable = 'active_layer', value = 0)
        
        # Update Node Elevations
        findactive = parcels.DataFrame['active_layer']==1 # filter for only parcels in active layer        
        vol_act = parcels.calc_aggregate_value(np.sum,
                                                  'volume',at = 'link',
                                                  filter_array = findactive)
        vol_stor = (vol_tot-vol_act)/(1-Lp)
        # ^ Jon-- what's the rationale behind only calculating new node elevations 
        # using the storage volume (rather than the full sediment volume)?
        
        # in a for loop for now
# %%        
        # XXXXXXXXXXXXX Pseudo code...
        for l in range(grid.number_of_nodes):
            if somethingaboutheadnode[l]==0: #we don't update head node elevations
                
                elev change = 2*volstore(downstreamlink)/(np.sum(width_of_neighborLinks * length_of_neighbor_links))
                
                elevation(l) = elevation(l) + elev change
        
# %%
    def adjust_slope(self, seed=0):
        """Adjusts slope for each link based on parcel motions from last
        timestep and additions from this timestep.
        """
        #  COPIED FROM JON'S MATLAB CODE
            #%compute/update slope
            for i in LinkNum
                if i==OutletLinkID
                    Slope(i,1)=(Elev(t,i)-mnelev(i,1))./Length(i,1);
                else
                    Slope(i,1)=(Elev(t,i)-Elev(t,Connect(i,2)))./Length(i,1);
                
            
            Slope(Slope<1e-4)=1e-4;

# %%
    def _calc_transport_wilcock_crowe(self, H_foreachlink): # Allison
        """Method to determine the transport time for each parcel in the active
        layer using a sediment transport equation. 
        
        Note: could have options here (e.g. Wilcock and Crowe, FLVB, MPM, etc)
        """
        
        parcels = self._parcels
        grid = self._grid
        
# %%
        # parcel attribute arrays from ItemCollector
        
        # another way of doing this --> check to see if this is copying. we don't want to be copying
        Darray = parcels.DataFrame.D.values
        
#        Darray = np.array(parcels.DataFrame.D,copy=False) # this gives a copy, but we can set copy to false..?
        Activearray = parcels.DataFrame.active_layer.values
        Rhoarray = parcels.DataFrame.density.values
        Volarray = parcels.DataFrame.volume.values
        Linkarray = parcels.DataFrame.element_id.values #link that the parcel is currently in
        R = (Rhoarray-rho)/rho
        
        # parcel attribute arrays to populate below
        frac_sand_array = np.zeros(np.size(element_id))
        vol_act_array = np.zeros(np.size(element_id))
        Sarray = np.zeros(np.size(element_id)) 
        Harray = np.zeros(np.size(element_id)) 
        Larray = np.zeros(np.size(element_id))
        d_mean_active = np.zeros(np.size(element_id))
        d_mean_active.fill(np.nan)
        Ttimearray = np.zeros(np.size(element_id)) 
        
        # Calculate bed statistics for all of the links
        vol_tot =  parcels.calc_aggregate_value(np.sum,'volume',at = 'link')
        
        findactive = parcels.DataFrame['active_layer']==1 # filter for only parcels in active layer        
        vol_act = parcels.calc_aggregate_value(np.sum,
                                                  'volume',at = 'link',
                                                  filter_array = findactive)
              
        findactivesand = np.logical_and(Darray<0.002,active_layer ==1)        
        vol_act_sand = parcels.calc_aggregate_value(np.sum,
                                                'volume',at = 'link',
                                                filter_array = findactivesand)

        vol_act_sand[np.isnan(vol_act_sand)==True] = 0
        frac_sand = vol_act_sand/vol_act
        frac_sand[np.isnan(frac_sand)== True] = 0
        
        # Calc attributes for each link, map to parcel arrays
        for i in range(grid.number_of_links): 
            
            active_here = np.where(np.logical_and(Linkarray == i,active_layer == 1))[0]
            
            d_act_i = Darray[active_here]
            vol_act_i = Volarray[active_here]
            d_mean_active[Linkarray == i] = np.sum(d_act_i * vol_act_i)/(vol_act[i])
            
            frac_sand_array[Linkarray == i] = frac_sand[i]
            vol_act_array[Linkarray == i] = vol_act[i]
            Sarray[Linkarray == i] = grid.at_link['channel_slope'][i]
            Harray[Linkarray == i] = H[i]
            Larray[Linkarray == i] = grid.at_link['link_length'][i]

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
        Ttimearray[Activearray==1] = rho**(3/2)*R[Activearray==1]*g*Larray[Activearray==1]*theta/W[Activearray==1]/tau[Activearray==1]**(3/2)/(1-frac_sand_array[Activearray==1])/frac_parcel[Activearray==1]
        Ttimearray[findactivesand==True] = rho**(3/2)*R[findactivesand==True]*g*Larray[findactivesand==True]*theta/W[findactivesand==True]/tau[findactivesand==True]**(3/2)/frac_sand_array[findactivesand==True]
        # ^ why?? if k = 1 ---> if it's sand...?  ASK JON about the logic here...                            
                            
        #del i b tau taur tautaur tautaur_cplzx taursg W findactive findactivesand

        # Assign those things to the grid -- might be useful for plotting later...?
        grid.at_link['sediment_total_volume'] = vol_tot        
        grid.at_link['sediment__active__volume'] = vol_act
        grid.at_link['sediment__active__sand_fraction'] = frac_sand
        
        
# %%
    def move_parcel_downstream(self, i):    # Jon
        """Method to update parcel location for each parcel in the active 
        layer. 
        """
        for i in range(number of links):
            
            # if the current link is empty, go to the next link
            
            for p in range(number of parcels in that link at this time):
                
                # Move parcels downstream
                
                # COPIED FROM JON'S MATLAB CODE
                
                sm=[];#%time to move through current and downstream links
                csm=[];#%cumulative time to move through current and ds links
                sm=P_tt{t,i}(p)*(1-P_loc{t,i}(1,p));#%seconds, time to move out of current link
                csm=sm;#%seconds
                ii=i;
                
                while csm(end)<=dt #%loop until cumulative time toward outlet exceeds timestep
                    ii=Connect(ii,2); #%set ds link
                    if isnan(ii) #%if downstream is the outlet
                        #%leave outlet before dt ends
                        sm=cat(1,sm,NaN);
                        csm=cat(1,csm,sum(sm));
                        OutVol(t,1)=OutVol(t,1)+P_vol{t,i}(1,p);
                        break
    
                    sm=cat(1,sm,P_tt{t,i}(p)./Length(i).*Length(ii)); # travel at same velocity but in ds link
                    csm=cat(1,csm,sum(sm)); # cumulative time to move through all subsequent links
    
                # update parcel location
                
                if ~isnan(csm(end)) #%check to make sure parcel is still in the system
                    pi=Connect(i,length(csm)); #%parcel link
                    if length(csm)==1 #%still in same link
                        pl=(P_tt{t,i}(p)*P_loc{t,i}(1,p)+dt)/P_tt{t,i}(p);#%update location from current ploc
                    else #%moved to a ds link
                        pl=1-((csm(end)-dt)/P_tt{t,i}(p));#%update location from beginning of link
                        if pl<1 #%overshooting the next link, UPDATE IN FUTURE
                            pl=1;
                        #%note csm(end)-dt computes time remaining to move
                        #%through the rest of the link
                                  
                    P_loc{t+1,pi}=cat(2,P_loc{t+1,pi},pl);
    
                    P_idx{t+1,pi}=cat(2,P_idx{t+1,pi},P_idx{t,i}(1,p));
                    P_storage{t+1,pi}=cat(2,P_storage{t+1,pi},0);%activate
                    P_vol{t+1,pi}=cat(2,P_vol{t+1,pi},P_vol{t,i}(1,p));
                    P_d{t+1,pi}=cat(2,P_d{t+1,pi},P_d{t,i}(1,p));
                    P_tt{t+1,pi}=cat(2,P_tt{t+1,pi},0);%set to 0                