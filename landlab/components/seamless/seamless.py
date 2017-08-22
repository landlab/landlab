#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 13:41:29 2017

@author: njlyons
"""

from copy import deepcopy
import Flowshed
import matplotlib.pyplot as plt
import numpy as np


class Seamless(object):
    """TODO: Description.
    """
    
    def __init__(self, grid, flowshed_ids, initial_species):
        """Initialize Seamless.
        """
        
        # Store grid and parameters.
        self._grid = deepcopy(grid)
        self.create_timestep_flowstep_data_structure(self._grid, flowshed_ids,
                                                     initial_species)
        
        # Initialize
        self.timestep = 0
        
    def run_one_step(self, dt, updated_grid):
                
        self.timestep += 1
        
        self.add_row_to_timestep_flowstep_data(self.timestep, updated_grid)
                
        for flowshed_id in self.flowshed_ids:
            
            # Get updated flowshed.
            con1 = (updated_grid.at_node['flowshed_id'] == flowshed_id)
            con2 = (updated_grid.node_is_core())
            updated_flowshed_ids = np.all([con1, con2], 0)
            
            # Get flowsheds that appeared in flowshed.
            updated_flowshed_with_old_ids = self._grid.at_node['flowshed_id'][updated_flowshed_ids]
            new_flowshed_ids = np.unique(updated_flowshed_with_old_ids)
            
            # Remove flowshed_id of self flowshed
            if any(flowshed_id == new_flowshed_ids):
                index_to_delete = np.where(new_flowshed_ids == flowshed_id)
                new_flowshed_ids = np.delete(new_flowshed_ids, index_to_delete)
                
            # Remove flowshed_id of null (-1) areas
            if any(-1 == new_flowshed_ids):
                index_to_delete = np.where(new_flowshed_ids == -1)
                new_flowshed_ids = np.delete(new_flowshed_ids, index_to_delete)
            
            # Add species from other flowsheds.
            self.timestep_flowshed_data[self.timestep][flowshed_id]['number_of_species'] += len(new_flowshed_ids)           
                        
            print(flowshed_id, new_flowshed_ids, self.timestep_flowshed_data[self.timestep][flowshed_id]['number_of_species'])
        
        # Update        
        
        self._grid = deepcopy(updated_grid)
        
    # timestep_flowshed_data management
        
    def create_timestep_flowstep_data_structure(self, grid, flowshed_ids, initial_species):
        
        # Convenience properties.
        self.flowshed_ids = flowshed_ids
        self.number_of_flowsheds = len(flowshed_ids)
        
        # Create structure with a row for timestep 0.
        time = 0
        self.timestep_flowshed_data = {}
        self.add_row_to_timestep_flowstep_data(time, grid)
        
        # Populate species count with initial values.
        for i,flowshed_id in enumerate(flowshed_ids):
            self.timestep_flowshed_data[time]['number_of_species'][flowshed_id] = initial_species[i]
            
    def add_row_to_timestep_flowstep_data(self, time, grid):
        
        self.timestep_flowshed_data[time] = np.zeros(self.number_of_flowsheds, dtype=[('flowshed_id', 'i4'),('number_of_species', 'i4'), ('stream_area', 'f4')])
        
        for i in range(self.number_of_flowsheds):
            # Populate known fields.
            self.timestep_flowshed_data[time]['flowshed_id'][i] = i
            self.timestep_flowshed_data[time]['stream_area'][i] = Flowshed.calculate_stream_area(grid, i)
            
            # Set unknown fields to -1.
            if time == 0:
                self.timestep_flowshed_data[time]['number_of_species'][i] = -1
            else:
                self.timestep_flowshed_data[time]['number_of_species'][i] = self.timestep_flowshed_data[time - 1]['number_of_species'][i]
            
    # Plotting
        
    def plot_number_of_species(self):
        
        # Get variables.
        timesteps = list(self.timestep_flowshed_data.keys())
        number_of_species = []
        for t in timesteps:
            number_of_species.append(sum(self.timestep_flowshed_data[t]['number_of_species']))
        
        plt.figure('Number of species')
        plt.plot(timesteps, number_of_species, 'k')
        plt.xlabel('Timestep')
        plt.ylabel('Number of species')
        
    def plot_area_versus_species(self):
        
        # Get variables.
        timesteps = list(self.timestep_flowshed_data.keys())
        area = []
        number_of_species = []
        for t in timesteps:
            area.append(self.timestep_flowshed_data[t]['stream_area'])
            number_of_species.append(self.timestep_flowshed_data[t]['number_of_species'])
        
        plt.figure('Number of species')
        plt.plot(area, number_of_species, 'k.')
        plt.xlabel('Stream area ($m^2$)')
        plt.ylabel('Number of species')
        
    def plot_area_versus_species_summation(self):
        
        # Get variables.
        timesteps = list(self.timestep_flowshed_data.keys())
        area = []
        number_of_species = []
        for t in timesteps:
            area.append(self.timestep_flowshed_data[t]['stream_area'])
            number_of_species.append(self.timestep_flowshed_data[t]['number_of_species'])
        
        plt.figure('Number of species')
        plt.plot(area, number_of_species, 'k.')
        plt.xlabel('Stream area ($m^2$)')
        plt.ylabel('Number of species')