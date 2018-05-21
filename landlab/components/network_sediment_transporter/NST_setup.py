# -*- coding: utf-8 -*-
"""
This code outlines very basic inputs for the Network_bedload_morphodyn component

Created on Sun May 20 15:54:03 2018

@author: pfeif
"""
import numpy

from landlab.grid.network import NetworkModelGrid
from landlab.components import FlowAccumulator

from landlab.plot import graph


# %% Set the geometry using Network model grid (should be able to read in a shapefile here)

y_of_node = (0, 1, 2, 2, 3, 4, 4, 1.25)
x_of_node = (0, 0, 1, -0.5, -1, 0.5, -1.5, -1)

nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

graph.plot_graph(grid, at='node,link')

grid.at_node['topographic__elevation'] = [0., 1, 3, 2, 3, 4, 4, 2]

area = grid.add_ones('cell_area_at_node', at = 'node')

fa = FlowAccumulator(grid)

fa.run_one_step()

#%%
# each must also have
grid.at_link['Area'] = [,,,,,,]
grid.at_link['Slope'] = [,,,,,,]  
grid.at_link['Length'] = [1000,1000,1000,1000,1000,1000,1000]

# modify elevations so they are consistent with adjusted slopes


# %% initialize bed sediment (will become its own component)


parcels = SedimentParcels(grid,initialization_info_including_future_forcing)
    # parcels have 'time added' and 'starting location/link' as an attribute. 


populate item collection

each parcel has: ID, source, volume, D, link, location w/in link, Timestep arrival in link? age?, active layer 

        # COPIED FROM JON'S MATLAB CODE
            #compute capacity for all links
            capacity=Btmax.*Length.*theta;# m3 PROBLEM: we haven't given ourselves the info to determine Btmax yet...
            
            
            Dpsd=[2,2.83,5.66,16,45.3,90.5,181]'./1000;%m
            Fsspsd=[0.1,0.03,0.05,0.15,0.12,0.18,0.37]';
            
            VolIn=Btmax.*Length.*theta;%m3
            
            %subsurface
            VolInd=repmat(VolIn,1,gsclass).*repmat(Fsspsd',LinkNum,1);
            pidx=1;
            
            for i = 1:length(VolInd)
                for k = 1:7
                    %j=7-k+1;
                    j=k;
                    if VolInd(i,j) == 0
                        continue
                    np=ceil(VolInd(i,j)./10);
                    pvol=VolInd(i,j)./np;
                    P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
                    P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
                    pidx=pidx+np;
                    P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
                    P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
                    P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
                    P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
            
            %surface
            VolInd=repmat(VolIn,1,gsclass).*repmat(Fsfpsd',LinkNum,1);
            for i = 1:length(VolInd)
                for k = 1:7
                    %j=7-k+1;
                    j=k;
                    if VolInd(i,j) == 0
                        continue
            
                    np=ceil(VolInd(i,j)./10);
                    pvol=VolInd(i,j)./np;
                    P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
                    P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
                    pidx=pidx+np;
                    P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
                    P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
                    P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
                    P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
 


# Add parcels in at a given time --> attribute in the item collection

# %% Instantiate component(s)
dis = ExteralDischargeSetter(grid, filename, model='dhsvm')
        
#sq = SyntheticDischargeMaker(discharge,drainage_area) # OR read DHSVM. 
    # define a surface_water_discharge for each link and each timestep
                    
#sc = SyntheticChannelGeomMaker(hydraulic_geometry_scaling_rules,discharge)
    # 

nst = NetworkSedimentTransporter(so many inputs!)

# %% Run the component(s)

for t in range(timesteps):
    
    # move any sediment additions from forcing Item collector to bed item collector
    
    # sq.run_one_step  
    Qgage = 200. # for Methow, this is ~4*mean flow
    dt = 60*60*24;
    
    # sc.run_one_step   
    Bgage=30.906.*Qgage.^0.1215;
    Hgage=0.2703.*Qgage.^0.3447;
    Agage=4.5895e+9;#m2
    B=(repmat(Bgage,1,LinkNum)./(Agage.^0.5)).*repmat(usarea',timesteps,1).^0.5;
    H=(repmat(Hgage,1,LinkNum)./(Agage.^0.4)).*repmat(usarea',timesteps,1).^0.4;
    Btmax=max(B,[],1)';
    
    
    nst.run_one_step
    
    