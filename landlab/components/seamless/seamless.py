"""TODO: Description.
"""

from copy import deepcopy
from itertools import groupby
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from species import Species


class Seamless(object):
    """TODO: Description.
    """
    
    def __init__(self, grid, region_field, minimum_capture_area=0,
                 extinction_area_threshold=0):
        """Initialize Seamless.
        """
        
        # Store grid and set parameters.
        self._grid = grid
        self._grid_of_previous_timestep = deepcopy(grid)
        self._region_field = region_field
        self._minimum_capture_area = minimum_capture_area
        self._extinction_area_threshold = extinction_area_threshold
        self.timestep = -1
        self.time = {}
        self.species = []
            
    @property
    def region_ids(self):
        """Get the region ids of current previouse timestep."""
        return self._unique(self._grid.at_node[self._region_field])
    
    @property
    def previous_region_ids(self):
        """Get the region ids of the previouse timestep."""
        return self._unique(self._grid.at_node['region_at_previous_timestep'])
    
    @property
    def species_ids(self):
        """Get the species ids."""
        ids = []
        for species in self.species:
            ids.append(species.identifier)
        return ids
    
    def run_one_step(self, time):
                
        self.timestep += 1
        self.time[self.timestep] = time
        
        # Update regions.
        
        reg_update = {}
        nodes_above_area_threshold = (self._grid.at_node['drainage_area'] >= 
                                      self._minimum_capture_area)

        for region_id in self.region_ids: 
            
            region_grid = self._grid.at_node[self._region_field] == region_id

            reg_update[region_id] = {}
            
            # Mark region for extinction if area is below extinction parameter.
            A = self._grid.at_node['drainage_area'][region_grid]
            reg_update[region_id]['extinction'] = (max(A) < self._extinction_area_threshold)

            for previous_region_id in self.previous_region_ids:
                if previous_region_id not in reg_update.keys():
                    reg_update[previous_region_id] = {}
                reg_update[previous_region_id].setdefault('parent_region_ids', [])
                
                # Skip loop if comparing the same region.
                if previous_region_id == region_id:
                    continue
                
                # Identify nodes intersecting regions and those above the
                # drainage area threshold.
                previous_grid = self._grid.at_node['region_at_previous_timestep'] == previous_region_id
                intersecting_nodes = np.all([previous_grid, region_grid,
                                             nodes_above_area_threshold], 0)
                
                regions_intersected = any(intersecting_nodes)
                
                if regions_intersected:                    
                    reg_update[previous_region_id]['parent_region_ids'].append(region_id)

        self.grab_region_snapshot()
        
        pprint(reg_update)
        
        # Update species.
        
        new_species = []
        
        if self.timestep < 1:
            prior_timestep = 0
        else:
            prior_timestep = self.timestep - 1
        
        for species in self.species:
            if prior_timestep in species.phylogeny.keys():
                species_regions = species.phylogeny[prior_timestep]['regions']
                timestep_regions = []
                
                for region in species_regions:
                    parent_regions = reg_update[region]['parent_region_ids']

                    # Reassign region to species or add ew species based upon
                    # region history.
                    if reg_update[region_id]['extinction']:
                        continue
                    elif len(parent_regions) == 0:
                        timestep_regions.append(region)
                    else:
                        for pr in parent_regions:
                            new_species.append(Species(self.timestep, [pr], species.identifier))
                            new_species.append(Species(self.timestep, [region], species.identifier))
                    
                if len(timestep_regions) > 0:
                    species.record_timestep(self.timestep, timestep_regions)
                
        for species in new_species:
            self.add_species(species)
        
    # Species methods
    
    def add_species(self, species):
        self.species.append(species)
        
    def remove_species_from_region(self, timestep, species, region):
        updated_regions = deepcopy(species.phylogeny[timestep - 1]['regions'])
        updated_regions.remove(region)

        if len(updated_regions) > 0:
            species.record_timestep(timestep, updated_regions)
    
    def species_at_nodes_at_current_timestep(self, nodes):
        species = self.species_at_nodes(self.regions, nodes)
        return species
    
    def species_at_nodes_at_previous_timestep(self, nodes):
        species = self.species_at_nodes(self.regions_of_previous_timestep, nodes)
        return species
        
    def species_at_nodes(self, regions, nodes):
        species = []
        for region in self.previous_region_ids:
            if any(region.mask[0][nodes]):
                species.append(region.species)
        return species    
    
    def species_in_region(self, timestep, region):
        region_species = []
        for species in self.species:
            if species.exists_at_timestep(timestep):
                if region in species.phylogeny[timestep]['regions']:
                    region_species.append(species)
        return region_species
    
    def species_with_id(self, identifier):
        return self._object_with_identifier(self.species, identifier)
    
    def print_phylogeny(self):
        for species in self.species:
            print(species.identifier)
            pprint(species.phylogeny)
    
    # Region methods.
    
    def grab_region_snapshot(self):
        region_array = self._grid.at_node[self._region_field]
        self._grid.at_node['region_at_previous_timestep'] = deepcopy(region_array)
    
    # Plotting methods.
        
    def plot_number_of_species(self):
        plt.figure('Number of species')

        time_species_existed = np.array([])
        for species in self.species:
            time_species_existed = np.append(time_species_existed, 
                                             list(species.phylogeny.keys()))
            
        time_species_existed = time_species_existed.flatten()
                
        d = {}
        for i in range(int(max(time_species_existed)) + 1):
            d[i] = len(np.where(time_species_existed == i)[0])

        for key, value in d.items():
            plt.plot(key, value, 'k.')
        
        plt.xlabel('Timestep')
        plt.ylabel('Number of species')
        
    def create_tree(self):
        
        self.tree = {}
        
        steps = list(reversed(list(self.time.keys())))

        for step in steps:
            time = self.time[step]            

            # Get species that existed in step.
            self.tree[time] = {}
            step_species = []
            for species in self.species:
                if step in species.phylogeny.keys():
                    step_species.append(species)      
                    
            # Group species by parent.
            groups = groupby(step_species, lambda x: x.parent_species_id)
            for parent_id, group in groups:             
                
                group_species = [species for species in group]
                self.tree[time][parent_id] = group_species
    
    def plot_tree(self, x_multiplier=0.001):
        
        # Prepare figure.
        plt.figure('Phylogeny', facecolor='white')
        ax = plt.axes(frameon=False)
        
        y_interval = 1
       
        species_position = {}
        times = list(self.tree.keys())
        
        # Construct tree beginning at last time.
        for i,time in enumerate(times):
            step_tree = self.tree[time]
                        
            if time == max(times):
                previous_time = time
                next_time = times[i + 1]
            elif time == min(times):
                previous_time = times[i - 1]
                next_time = time
            else:
                previous_time = times[i - 1]
                next_time = times[i + 1]

            x_max = (time + (previous_time - time) * 0.5) * x_multiplier
            x_mid = time * x_multiplier
            x_min = (time - (time - next_time) * 0.5) * x_multiplier

            y = y_interval
            y_initial = y_interval
                  
            for parent_id, species_list in step_tree.items():
                group_only = False
                
                y_min = deepcopy(y)
                
                for species in species_list:
                    
                    if species.identifier in species_position.keys():
                        species_position[species.identifier]
                    else:
                        species_position[species.identifier] = y_initial
                        y_initial += y_interval
                    
                    y = species_position[species.identifier]
                    
                    existed_previous_step = (len(times) - i - 2) in species.phylogeny.keys()
                    if existed_previous_step:
                        # Draw line when a species continues across timesteps.
                        plt.plot([x_max, x_min], [y, y], 'k')
                        y_min += y_interval * 0.5
                        
                    else:   
                        # Draw the base of species line.
                        plt.plot([x_max, x_mid], [y, y], 'b')
                        
                        # Draw line that connects species at a timestep.
                        if time != 0:
                            y_max = y
                            plt.plot([x_mid, x_mid], [y_min, y_max], 'c')
                            group_only = True
                        
                    y += y_interval
                    y_max = y - y_interval
                    y_mid = (y_max - y_min) * 0.5 + y_min
                    
                    species_position[parent_id] = y_mid
                
                # Plot trunk
                if group_only:
                    y_mid = (y_max - y_min) * 0.5 + y_min
                    plt.plot([x_mid, x_min], [y_mid, y_mid], 'r')
        
        # Format figure.
        plt.xlim(xmin=0)
        plt.ylim(ymin=y_interval * 0.5)
        ax.get_xaxis().tick_bottom()
        ax.axes.get_yaxis().set_visible(False)
        xmin, xmax = ax.get_xaxis().get_view_interval()
        ymin, ymax = ax.get_yaxis().get_view_interval()
        ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=1.5))
        plt.xlabel('Time (ky)')
                
    def plot_area_versus_species(self):
        
        # Get variables.
        timesteps = list(self.region_data.keys())
        area = []
        number_of_species = []
        for t in timesteps:
            area.append(self.region_data[t]['stream_area'])
            number_of_species.append(self.region_data[t]['number_of_species'])
        
        plt.figure('Number of species')
        plt.plot(area, number_of_species, 'k.')
        plt.xlabel('Stream area ($m^2$)')
        plt.ylabel('Number of species')
        
    def plot_delta_area(self):
        
        # Get variables.
        timesteps = list(self.region_data.keys())[:-1]
        delta_area = []
        number_of_species = []
        for t in timesteps:
            delta_area.append(self.region_data[t + 1]['stream_area'] - self.region_data[t]['stream_area'])
            number_of_species.append(self.region_data[t+1]['number_of_species'])
        
        plt.figure('Number of species')
        plt.plot(delta_area, number_of_species, 'k.')
        plt.xlabel('Delta stream area ($m^2$)')
        plt.ylabel('Number of species')
    
    # Convenience methods.
        
    def _unique(self, values):
        unique = np.unique(values)
        unique_no_nulls = self._remove_value(-1, unique)
        return unique_no_nulls
    
    def _remove_value(self, value, array):
        if any(value == array):
            indices_to_delete = np.where(array == value)
            array = np.delete(array, indices_to_delete)
            
        return array     
    
    def _object_with_identifier(self, item_list, identifier):
        for item in item_list:
            if item.identifier == identifier:
                return item
                break
    