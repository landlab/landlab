"""TODO: Description.

Component written by Nathan Lyons beginning August 20, 2017.
"""

from copy import deepcopy
from landlab import Component
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from .species import Species


class CladeDiversifier(Component):
    """TODO: Description.
    
    SEAMLESS (Spatially Explicit Area Model of Landscape Evolution by SimulationS)
    
    The primary method of this class is :func:`run_one_step`.
    
    Construction:

        CladeDiversifier(grid, minimum_capture_area=0,
                         extinction_area_threshold=0)
        
    Parameters
    ----------
    grid : ModelGrid
        A grid.
    minimum_capture_area : float
        The minumum area that a region capture is recognized
    extinction_area_threshold : float
        The species in regions below this value will become extinct
    """
    
    _name = 'CladeDiversifier'

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}
    
    def __init__(self, grid, minimum_capture_area=0,
                 extinction_area_threshold=0, **kwds):
        """Initialize the CladeDiversifier.
    
        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        minimum_capture_area : float, optional (defaults to 0 m^2)
            The minumum area that a region capture is recognized
        extinction_area_threshold : float, optional (defaults to 0 m^2)
            The species in regions below this value will become extinct
        """
        
        # Store grid and set parameters.
        self._grid = grid
        self._minimum_capture_area = minimum_capture_area
        self._extinction_area_threshold = extinction_area_threshold
        self.timestep = -1
        self.at_timestep = {'time': [], 'captured_area': []}
        self.species = []
        
        
    
    @property
    def species_ids(self):
        """Get the species ids."""
        ids = []
        for species in self.species:
            ids.append(species.identifier)
        return ids
    
    @property
    def timesteps(self):
        number_of_timesteps = len(self.at_timestep['time'])
        return range(number_of_timesteps)
    
    def run_one_step(self, current_time):
        self.timestep += 1
        
        mg = self._grid
        A = mg.at_node['drainage_area']
        AE = self._extinction_area_threshold
        stream_nodes = A >= self._minimum_capture_area
        
        new_species = []
        
        if self.timestep == 0:
            prior_species = self.species
        else:
            prior_species = self.species_at_timestep(self.timestep - 1)
        
        captured_nodes = []
        
        for species in prior_species:
            prior_range = deepcopy(species.nodes)
            species.update_range(mg, self._minimum_capture_area)
            updated_range = species.nodes
            persistent_range = np.all([prior_range, updated_range, stream_nodes], 0)
            
            if len(np.where(persistent_range)[0]) == 0:
                persistent_range_max_A = 0
            else:
                persistent_range_max_A = max(A[persistent_range])
    
            new_range = np.all([prior_range, np.invert(updated_range)], 0)
            captured_nodes.extend(np.where(new_range)[0])
            
            barrier_created = any(persistent_range) and any(new_range)

            if barrier_created:
                # Add new species.
                parent_id = species.identifier

#                if max(A[new_range]) >= AE:
                species = Species(self.timestep, new_range, parent_id)
                new_species.append(species)
            
#                if persistent_range_max_A >= AE:
                species = Species(self.timestep, updated_range, parent_id)
                new_species.append(species)

            else:
                # Update existing species.
                if persistent_range_max_A <= AE:
                    # Extinction.
                    print(999)
                    species.nodes = prior_range  
                else:
                    # No changes so species persists.
                    species.timesteps_existed.append(self.timestep)
                    
#        print(len(self.species), self.species_at_timestep(self.timestep - 1), current_time)
        self.species.extend(new_species)
        
        self.at_timestep['time'].append(current_time)
        captured_area = len(np.unique(captured_nodes)) * mg.dx * mg.dy
        self.at_timestep['captured_area'].append(captured_area)
    
    # Species methods
    
    def add_new_species(self, timestep, nodes, parent_species_id=-1):
        species = Species(timestep, nodes, parent_species_id)
        self.species.append(species)
        return species
    
    def species_at_timestep(self, step):
        species = []
        for s in self.species:
            if step in s.timesteps_existed:
                species.append(s)
        return species
    
    def species_at_nodes_at_current_timestep(self, nodes):
        species = self.species_at_nodes(self.regions, nodes)
        return species
    
    def species_at_nodes_at_prior_timestep(self, nodes):
        species = self.species_at_nodes(self.prior_region_ids, nodes)
        return species
        
    def species_at_nodes(self, regions, nodes):
        species = []
        for region in self.prior_region_ids:
            if any(region.mask[0][nodes]):
                species.append(region.species)
        return species    
    
    def species_with_id(self, identifier):
        return self._object_with_identifier(self.species, identifier)
    
    def get_tree(self):
        
        tree = {}
        
        steps = list(reversed(self.timesteps))

        for step in steps:
            time = self.at_timestep['time'][step]           

            # Get species that existed in step.
            tree[time] = {}
            for species in self.species:
                if step in species.timesteps_existed:
                    key = species.parent_species_id
                        
                    if key not in tree[time].keys():
                        tree[time][key] = []
                    tree[time][key].append(species)
                
        return tree
    
    # Print methods.
    
    def print_phylogeny(self):
        for species in self.species:
            print('\n')
            print(species.identifier)
            pprint(species.phylogeny)
            
    def print_tree(self):
        tree = self.get_tree()
        pprint(tree)
    
    # Plotting methods.
        
    def plot_number_of_species(self):
        
        time_in_ky = np.multiply(self.at_timestep['time'], 1e-3)
        
        plt.figure('Number of species')

        time_species_existed = np.array([])
        for species in self.species:
            time_species_existed = np.append(time_species_existed, 
                                             species.timesteps_existed)
            
        time_species_existed = time_species_existed.flatten()
                
        d = {}
        for i in range(int(max(time_species_existed)) + 1):
            d[i] = len(np.where(time_species_existed == i)[0])

        for key, value in d.items():
            plt.plot(key, value, 'k.')
        
        plt.xlabel('Timestep')
        plt.ylabel('Number of species')
          
    def plot_tree(self, x_multiplier=0.001):
        
        tree = self.get_tree()
        
        # Prepare figure.
        plt.figure('Phylogeny', facecolor='white')
        ax = plt.axes(frameon=False)
        
        y_spacing = 1
       
        species_position = {}
        times = list(tree.keys())

        # Construct tree beginning at final time.
        for i, time in enumerate(times):

            if time == max(times):
                prior_time = time
                next_time = times[i + 1]
            elif time == min(times):
                prior_time = times[i - 1]
                next_time = time
            else:
                prior_time = times[i - 1]
                next_time = times[i + 1]

            x_max = (time + (prior_time - time) * 0.5) * x_multiplier
            x_min = (time - (time - next_time) * 0.5) * x_multiplier
            x_mid = np.mean([x_min, x_max])

            y_min = deepcopy(y_spacing)
                  
            for parent_id, species_list in tree[time].items():
                species_group = False

                y_species = deepcopy(y_min)

                for species in species_list:
                    
                    no_parent = species.parent_species_id == -1
                    
                    existed_prior_step = (len(times) - i - 2) in species.timesteps_existed

                    if species.identifier in species_position.keys():
                        y_species = species_position[species.identifier]

                    # Draw line when a species continues across timesteps.
                    if existed_prior_step or no_parent:
                        x = x_min
                    else:
                        x = x_mid
                        species_group = True
                    
                    plt.plot([x, x_max], [y_species, y_species], 'k')
                        
                    y_mid = np.mean([y_min, y_species])

                    if parent_id not in species_position.keys():
                        species_position[parent_id] = y_mid
                
                # Plot trunk
                if species_group:
                    
                    # Draw line that connects species at a timestep.
                    if time != 0:
                        plt.plot([x_mid, x_mid], [y_min, y_species], 'r')
                    
                    # Draw the base line of a species group.
                    plt.plot([x_min, x_mid], [y_mid, y_mid], 'c')
                    
                    species_position[parent_id] = y_mid
                                                            
                y_min = deepcopy(y_species + y_spacing)
                        
        # Format figure.
        plt.xlim(xmin=0, xmax=max(times) * x_multiplier)
        plt.ylim(ymin=y_spacing * 0.5)
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
            delta_area.append(self.region_data[t + 1]['stream_area'] - 
                              self.region_data[t]['stream_area'])
            number_of_species.append(self.region_data[t+1]['number_of_species'])
        
        plt.figure('Number of species')
        plt.plot(delta_area, number_of_species, 'k.')
        plt.xlabel('Delta stream area ($m^2$)')
        plt.ylabel('Number of species')
    
    def plot_captured_area(self):
        
        time_in_ky = np.multiply(self.at_timestep['time'], 1e-3)
        area_in_km = np.multiply(self.at_timestep['captured_area'], 1e-6)

        plt.figure('Captured area')
        plt.plot(time_in_ky, area_in_km, 'k')
        plt.xlabel('Time (ky)')
        plt.ylabel('Area captured ($km^2$)')
    
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
    