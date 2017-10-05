"""TODO: Description.

Component written by Nathan Lyons beginning August 20, 2017.
"""

from copy import deepcopy
from landlab import Component
from landlab.plot import imshow_grid
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import colorConverter
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
        A_min = self._minimum_capture_area
        A_e = self._extinction_area_threshold
        
        if self.timestep == 0:
            prior_species = self.species
        else:
            prior_species = self.species_at_timestep(self.timestep - 1)
        
        new_species = []
        captured_nodes = []
        
        for species in prior_species:
            child_species,captured = species.evolve(mg, self.timestep, A_min,
                                                    A_e)
            if len(child_species) > 0:
                new_species.extend(child_species)

            if len(captured) > 0:
                captured_nodes.extend(captured)
                
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
        
    def species_at_node(self, node):
        species_list = []
        for species in self.species:
            if species.nodes[node]:
                species_list.append(species)
        return species_list    
    
    def species_with_id(self, identifier):
        return self._object_with_identifier(self.species, identifier)
    
    def get_tree(self, selected_species=None):
        
        if selected_species == None:
            selected_species = self.species
        
        tree = {}
        
        steps = list(reversed(self.timesteps))
        
        for step in steps:
            time = self.at_timestep['time'][step]           

            # Get species that existed in step.
            tree[time] = {}
            for species in selected_species:
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
            
    def print_tree(self, selected_species=None):
        tree = self.get_tree(selected_species)
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
        number_of_timesteps = int(max(time_species_existed))
        print(number_of_timesteps)
        for i in range(number_of_timesteps):
            d[i] = len(np.where(time_species_existed == i)[0])

        for key, value in d.items():
            plt.plot(key, value, 'k.')
        
        plt.xlabel('Timestep')
        plt.ylabel('Number of species')
          
    def plot_tree(self, x_multiplier=0.001, selected_species=None):
        
        tree = self.get_tree(selected_species)
        
        # Prepare figure.
        plt.figure('Phylogeny')
        ax = plt.axes(frameon=False)
        
        y_spacing = 1
       
        species_position = {}
        times = list(tree.keys())
        print(times)
        timesteps = np.array(self.at_timestep['time'])

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

            timestep = np.where(timesteps == time)[0][0]
            prior_timestep = timestep-1
            
            x_max = (time + (prior_time - time) * 0.5) * x_multiplier
            x_min = (time - (time - next_time) * 0.5) * x_multiplier
            x_mid = np.mean([x_min, x_max])

            y_min = deepcopy(y_spacing)
                  
            for parent_id, species_list in tree[time].items():
                species_group = False

                y_species = deepcopy(y_min)

                for species in species_list:
                    
                    if species.identifier in species_position.keys():
                        y_species = species_position[species.identifier]
                    else:
                        y_species += y_spacing
                        
                    # Draw line when a species continues across timesteps.
                    existed_prior_step = prior_timestep in species.timesteps_existed
                    no_parent = species.parent_species_id == -1
                    if existed_prior_step or no_parent:
                        x = x_min
                    else:
                        x = x_mid
                        species_group = True
                    
                    plt.plot([x, x_max], [y_species, y_species], 'k')
                        
                y_mid = np.mean([y_min, y_species]) + y_spacing * 0.5
                species_position[parent_id] = y_mid
                
                # Plot trunk
                if species_group:
                    
                    # Draw line that connects species at a timestep.
                    if time != 0:
                        plt.plot([x_mid, x_mid], [y_min + y_spacing, y_species], 'r')
                    
                    # Draw the base line of a species group.
                    plt.plot([x_min, x_mid], [y_mid, y_mid], 'c')
                                                            
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
        plt.xlim(xmin=0, xmax=max(time_in_ky))
        plt.ylim(ymin=0)
    
    def plot_species_range(self, species):
        
        range_masks = []
        for s in species:
            range_masks.append(s.nodes)
        combined_range_mask = np.any(range_masks, 0)
        
#        range_mask = np.zeros(self.grid.number_of_nodes)
#        range_mask[species.nodes] = 1
              
        # generate the colors for your colormap
        c1 = colorConverter.to_rgba('white')
        c2 = colorConverter.to_rgba('cyan')
        
        cmap = colors.LinearSegmentedColormap.from_list('streamOverlay',[c1,c2], 2)
        cmap._init() # create the _lut array, with rgba values
        alphas = np.linspace(0, 1, cmap.N+3)
        cmap._lut[:,-1] = alphas
        
        plt.figure('Species range')
        
        imshow_grid(self.grid, 'topographic__elevation', cmap='gray')
        imshow_grid(self.grid, combined_range_mask, cmap=cmap,
                    allow_colorbar=False)
        
    
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
    