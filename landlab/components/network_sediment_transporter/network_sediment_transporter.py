
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
                 bed_parcels_item_collector,
                 discharge,
                 transport_method = 'WilcockCrow',
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

        self.transport_method = transport_method
        if self.transport_method =="WilcockCrow":
            self.update_transport_time = self._calc_transport_wilcock_crow
            
    def partition_active_and_storage_layers(self, **kwds):
        """For each parcel in the network, determines whether it is in the
        active or storage layer during this timestep
        """
        # COPIED FROM JON'S MATLAB CODE
        #cycle through each link
            for i=1:LinkNum
                if ~isempty(P_vol{t,i}) %only do this check capacity if parcels are in link
                    #First In Last Out
                    #compute cumulative volume in link beginning with last in
                    cvol=fliplr(cumsum(fliplr(P_vol{t,i})));
                    #%determine which parcels that were the first to enter are above capacity
                    exc=find(cvol>capacity(i,1),1,'last');
                    #%if parcels have been identified above capacity then
                    if ~isempty(exc)
                    #    %set their status to inactive, 1
                        P_storage{t,i}(1,1:exc)=1;
                    #    %set their travel time to 0
                        P_tt{t,i}(1,1:exc)=0;
                        
                    #    %UPDATE Elevations beginning here
                    #    %compute volume inactive
                        vstor=sum(P_vol{t,i}(1:exc))./(1-Lp);%m3 vol in storage, porosity Lp=0.4(?)
                    #    %length and upstream link lengths to compute new elev
                        Elev(t,i)=mxelevmod(i,1)+(2.*vstor)./sum(B(t,elevid)'.*Length(elevid,1));
                    #    %UPDATE Elevations complete
                        
                        clear cvol exc vstor usid elevid

            

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

    def _calc_transport_wilcock_crowe(self): # Allison
        """Method to determine the transport time for each parcel in the active
        layer using a sediment transport equation. 
        
        Note: could have options here (e.g. Wilcock and Crowe, FLVB, MPM, etc)
        """
# %%
        Darray = np.array(parcels.DataFrame.D)
        Activearray = np.array(parcels.DataFrame.active_layer)
        Rhoarray = np.array(parcels.DataFrame.density)
        Volarray = np.array(parcels.DataFrame.volume)
        Linkarray = np.array(parcels.DataFrame.element_id) #link that the parcel is currently in
        Ttimearray = np.zeros(np.size(element_id)) 
        Ttimearray.fill(np.nan)
        
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
        
        d_mean_active = np.zeros(np.size(element_id))
        d_mean_active.fill(np.nan)
        
        # Assign those things to the grid -- might be useful for plotting later...?
        
        grid.at_link['sediment_total_volume'] = vol_tot        
        grid.at_link['sediment__active__volume'] = vol_act
        grid.at_link['sediment__active__sand_fraction'] = frac_sand
        
# %%        
        for i in range(grid.number_of_links):
                 # xxx  if sediment_volume_active_layer > 0
                  # xxx    continue
                    
                # calculate the surface mean grain size
                
                active_here = np.where(np.logical_and(Linkarray == i,active_layer == 1))[0]
                
                d_act_i = Darray[active_here]
                vol_act_i = Volarray[active_here]
                
                d_mean_active[i] = np.sum(d_act_i * vol_act_i)/(vol_act[i])
                
                # Wilcock and crowe claculate transport for each active parcel in this link
                
                for j in range(0,np.size(active_here)): 
                    
                    R_j = (Rhoarray[j]-rho)/rho
                    
                    taursg = rho * R_j * g * d_mean_active[i] * (0.021 + 0.015*np.exp(-20.*frac_sand[i]))
                    
#                    actidxj=and(P_storage{t,i}==0,P_d{t,i}==Dpsd(k));
#                    actvolj=sum(P_vol{t,i}(actidxj));
                    
                    Fj = vol_act[j]/np.sum(vol_act_i);
                    
                    b = 0.67 / (1 + np.exp(1.5 - d_act_i[j] / d_mean_active[i]))
                    
                    tau = rho * g * H[i] * grid.at_link['channel_slope'][i]
                    
                    taurj = taursg * (d_act_i[j] / d_mean_active[i]) **b
                    
                    tautaurj = tau / taurj
                    
                    if tautaurj < 1.35:
                        
                        Wj = 0.002 * tautaurj **7.5
                        
                    else:
                        
                        Wj = 14.*( 1 - 0.894/np.sqrt(tautaurj))**4.5
                    
                    
                    if k==1 # UGH-- how to I index back to the original link number?
                        Ttimearray[j]=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./Fs;
                        
                    else
                        P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./(1-Fs)./Fj;
                    
                    
                    clear actidxj actvolj Fj b tau taurj tautaurj Wj
                
                clear actidx actvol actsandidx actsandvol Fs taursg


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