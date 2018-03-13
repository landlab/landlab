# BiotaEvolver plot functions.

from copy import deepcopy
from landlab.plot import imshow_grid
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np


def plot_area_versus_species(region_data):

    # Get variables.
    timesteps = list(region_data.keys())
    area = []
    number_of_species = []
    for t in timesteps:
        area.append(region_data[t]['stream_area'])
        number_of_species.append(region_data[t]['number_of_species'])

    plt.figure('Number of species')
    plt.plot(area, number_of_species, 'k.')
    plt.xlabel('Stream area ($m^2$)')
    plt.ylabel('Number of species')

def plot_delta_area(region_data):

    # Get variables.
    timesteps = list(region_data.keys())[:-1]
    delta_area = []
    number_of_species = []
    for t in timesteps:
        delta_area.append(region_data[t + 1]['stream_area'] - 
                          region_data[t]['stream_area'])
        number_of_species.append(region_data[t+1]['number_of_species'])

    plt.figure('Number of species')
    plt.plot(delta_area, number_of_species, 'k.')
    plt.xlabel('Delta stream area ($m^2$)')
    plt.ylabel('Number of species')

def plot_captured_area(at_timestep):

    time_in_ky = np.multiply(at_timestep['time'], 1e-3)
    area_in_km = np.multiply(at_timestep['captured_area'], 1e-6)

    plt.figure('Captured area')
    plt.plot(time_in_ky, area_in_km, 'k')
    plt.xlabel('Time (ky)')
    plt.ylabel('Area captured ($km^2$)')
    plt.xlim(xmin=0, xmax=max(time_in_ky))
    plt.ylim(ymin=0)

def plot_species_range(species, grid):

    range_masks = []
    for s in species:
        range_masks.append(s.range_mask)
    combined_range_mask = np.any(range_masks, 0)

#    range_mask = np.zeros(self.grid.number_of_nodes)
#    range_mask[species.range_mask] = 1

    # generate the colors for your colormap
    c1 = [1, 1, 1, 1]
    c2 = [0, 0, 1, 1]

    cmap = LinearSegmentedColormap.from_list('streamOverlay', [c1,c2], 2)
    cmap._init() # create the _lut array, with rgba values
    alphas = np.linspace(0, 1, cmap.N+3)
    cmap._lut[:,-1] = alphas

    plt.figure('Species range')

    imshow_grid(grid, 'topographic__elevation', cmap='gray')
    imshow_grid(grid, combined_range_mask, cmap=cmap,
                    allow_colorbar=False)

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

    timesteps = np.array(self.at_timestep['time'])

    # Construct tree beginning at final time.
    for i, time in enumerate(times):

        if time == max(times):
            later_time = time
            earlier_time = times[i + 1]
        elif time == min(times):
            later_time = times[i - 1]
            earlier_time = time
        else:
            later_time = times[i - 1]
            earlier_time = times[i + 1]

        timestep = np.where(timesteps == time)[0][0]
        prior_timestep = timestep - 1

        x_max = (time + (later_time - time) * 0.5) * x_multiplier
        x_min = (time - (time - earlier_time) * 0.5) * x_multiplier
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
                    plt.plot([x_mid, x_mid],
                             [y_min + y_spacing, y_species], 'r')

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
    ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black',
                         linewidth=1.5))
    plt.xlabel('Time (ky)')