"""TODO: Description.
"""

from copy import deepcopy
from itertools import groupby
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from species import Species


class Seamless(object):
    """TODO: Description.
    """
    
    def __init__(self, grid, region_field, minimum_capture_area=0):
        """Initialize Seamless.
        """
        
        # Store grid and set parameters.
        self._grid = grid
        self._grid_of_previous_timestep = deepcopy(grid)
        self.region_field = region_field
        self._minimum_capture_area = minimum_capture_area
        self.timestep = 0
        self.species = []
            
    @property
    def region_ids(self):
        """Get the region ids of current previouse timestep."""
        return self._unique(self._grid.at_node[self.region_field])
    
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
    
    def run_one_step(self, dt):
                
        self.timestep += 1
        
        nodes_above_area_threshold = (self._grid.at_node['drainage_area'] >= 
                                      self._minimum_capture_area)
        
        for region in self.region_ids: 
            
            captured_regions = []         
            current_grid = self._grid.at_node[self.region_field] == region
            
            species_not_modified = self.species_in_region(self.timestep - 1, region)
            new_species = []
            
            for previous_region in self.previous_region_ids:
                
                # Skip loop if comparing the same region.
                if previous_region == region:
                    continue
                
                # Identify nodes intersecting regions and those above the
                # drainage area threshold.
                previous_grid = self._grid.at_node['region_at_previous_timestep'] == previous_region
                intersecting_nodes = np.all([previous_grid, current_grid,
                                             nodes_above_area_threshold], 0)
                
                regions_intersected = any(intersecting_nodes)
                                    
                if regions_intersected:                    
                    captured_regions.append(previous_region)
                    # Update species.
                    for species in self.species_in_region(self.timestep - 1, previous_region):
                        new_species.append(Species(self.timestep, [previous_region], species.identifier))
                        new_species.append(Species(self.timestep, [region], species.identifier))
                        self.remove_species_from_region(self.timestep, species, previous_region)
                        
                        if species.identifier in species_not_modified:
                            species_not_modified.remove(species.identifier)
                            
                    for species in new_species:
                        self.add_species(species)
            
            for species in species_not_modified:
                species.increment_timestep(self.timestep)
        
        self.grab_region_snapshot()
             
    # Species methods
    
    def add_species(self, species):
        self.species.append(species)
        
    def remove_species_from_region(self, timestep, species, region):
        species.phylogeny[timestep - 1]['regions'].remove(region)
        updated_regions = deepcopy(species.phylogeny[timestep - 1]['regions'])
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
    
    # Region methods.
    
    def grab_region_snapshot(self):
        region_array = self._grid.at_node[self.region_field]
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
        
    def plot_tree(self):
        plt.figure('Phylogeny', facecolor='white')
        ax = plt.axes(frameon=False)

        for step in range(10, -1, -1):
            # Get species that existed in step.
            step_species = []
            for species in self.species:
                if step in species.phylogeny.keys():
                    step_species.append(species)         
                    species.existed_previous_step = (step - 1) in species.phylogeny.keys()
            
            # Group species by parent.
            groups = groupby(step_species, lambda x: x.parent_species_id)
            tree = {}
            for parent_id, group in groups:  
                tree[parent_id] = [species for species in group]
            
            # Group by existance of a parent.
            tree2 = deepcopy(tree)
            
            for parent_id, species_list in tree.items():
                groups = groupby(species_list, lambda x: x.existed_previous_step)
                
                for parent_id, group in groups:  
                    tree2[parent_id] = deepcopy([species for species in group])
            print(555, tree.keys())
            print(666, tree2.keys())
            # Plot.
            x = 0
            y_max = step
            y_mid = y_max - 0.5
            y_min = y_max - 1

            for parent_id, species_list in tree2.items():
                if parent_id is not False and parent_id is not True:
                    # Draw vertical lines.
                    x_min = deepcopy(x)
                    for species in species_list:
                        if species.existed_previous_step:
                            plt.plot([x,x], [y_max, y_min], 'r')
                            x_min += 1
                        else:    
                            plt.plot([x,x], [y_max, y_mid], 'b')
                            
                            # Draw horizontal line.
                            x_max = x
                            plt.plot([x_min, x_max], [y_mid, y_mid], 'c')
                        x += 2
                    
                    x_max = x - 2
                    
                    # Plot trunk
                    if len(species_list) > 0:
                        x_mid = (x_max - x_min) * 0.5 + x_min
                        plt.plot([x_mid, x_mid], [y_mid, y_min], 'k')

#                print(step,parent_id, len(species_list))

        # Format figure.
        ax.get_yaxis().tick_left()
        ax.axes.get_xaxis().set_visible(False)
        xmin, xmax = ax.get_xaxis().get_view_interval()
        ymin, ymax = ax.get_yaxis().get_view_interval()
        ax.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=2))
        plt.ylabel('Timestep')
                
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
    