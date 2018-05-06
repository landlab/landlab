# BiotaEvolver plot functions.

from copy import deepcopy
from landlab.plot import imshow_grid
from matplotlib.colors import LinearSegmentedColormap
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

def plot_number_of_species(record, species_data_frame, axes=None):

    # Prepare figure.
    if axes == None:
        fig = plt.figure('Number of species')
        axes = fig.add_axes(plt.axes())

    times = list(record.keys())
    count = []

    sdf = species_data_frame.species

    for t in times:
        species_t = sdf.object[record.time == record.time__last].tolist()
        count.append(len(species_data_frame.species_at_time(t)))

    times_in_ky = np.multiply(times, 1e-3)

    axes.plot(times_in_ky, count, 'k')

    axes.set_xlabel('Time (ky)')
    axes.set_ylabel('Number of species')


def plot_tree(species, times, axes=None, x_multiplier=0.001,
              x_axis_title='Time (ky)'):
    """Plot a phylogenetic tree.

    """
    # Add species identifier columns to a copy of the species DataFrame.

    species = deepcopy(species)
    species['parent_id'] = [np.nan] * len(species)
    species['species_number'] = [np.nan] * len(species)

    for i, row in species.iterrows():
        ps = row.object.parent_species
        species_num = row.object.identifier[1]
        species.loc[i, 'species_number'] = species_num
        if ps == -1:
            species.loc[i, 'parent_id'] = -1
        else:
            species.loc[i, 'parent_id'] = ps.identifier[1]

    # Sort inputs in reverse order for plotting.

    clades = species.index.get_level_values('clade').unique().tolist()
    clades.sort(reverse=True)

    times.sort(reverse=True)

    # Create matplotlib axes if axes were not inputted.

    if axes == None:
        plt.close('Phylogenetic tree')
        fig = plt.figure('Phylogenetic tree')
        axes = fig.add_axes(plt.axes())

    # Construct the tree by clade then by time.

    y_species_spacing = 0.5
    y_max = 0

    lines = []

    for clade in clades:
        clade_mask = species.index.get_level_values('clade') == clade
        species_clade = species.loc[clade_mask]

        # Store y-axis position of species to connect branches over time.
        y_species = {}

        for i, time in enumerate(times):
            # Get the temporal bounds represented by this time.
            if time == max(times):
                later_time = time
                earlier_time = times[i + 1]
            elif time == min(times):
                later_time = times[i - 1]
                earlier_time = time
            else:
                later_time = times[i - 1]
                earlier_time = times[i + 1]

            t_extinct = species_clade.time_disappeared.tolist()
            last_time = np.array(t_extinct)
            last_time[np.argwhere(np.isnan(t_extinct))] = max(times)

            # Get the species at time sorted by species identifier.
            time_mask = np.all([species_clade.time_appeared <= earlier_time,
                                last_time >= later_time], 0)
            species_time = species_clade.loc[time_mask]
            species_time = species_time.sort_values(['parent_id',
                                                     'species_number'],
                                                     ascending=False)

            # Extend the x values by half of time.
            x_max = (time + (later_time - time) * 0.5) * x_multiplier
            x_min = (time - (time - earlier_time) * 0.5) * x_multiplier
            x_max = later_time * x_multiplier
            x_min = earlier_time * x_multiplier

            y_time = []

            clade_branch = False

            for i, row in species_time.iterrows():

                # Get species object, s of current row.
                s = row.object

                no_parent = s.parent_species == -1

                if s in y_species.keys():
                    # Previously plotted.
                    y_next = y_species[s]
                else:
                    # First time plotting.
                    y_next = y_max + y_species_spacing

                if no_parent:
                    lines.extend(axes.plot([x_min, x_max], [y_next, y_next], 'm'))
                else:
                    lines.extend(axes.plot([x_min, x_max], [y_next, y_next], 'c'))

                y_species[s] = y_next

                existed_prior_time = row.time_appeared < earlier_time

                y_time.append(y_next)

                if y_next > y_max:
                    y_max = y_next

                if not existed_prior_time:
                    y_species[s.parent_species] = np.mean(y_time)
                    clade_branch = True

                # Plot species number.
                if time == max(times):
                    sid = s.identifier[1]
                    plt.text(x_max, y_next, sid, fontsize='small', va='center',
                             bbox=dict(boxstyle='square,pad=0.1', fc='w',
                                       ec='none'))

            # Plot vertical line that connects sibling species.
            if clade_branch:
                lines.extend(axes.plot([x_min, x_min], [min(y_time),
                                       max(y_time)], 'y'))

        # Plot clade id.
        if time == 0:
            plt.text(x_min, y_next, clade, fontsize='small', ha='right',
                     va='center', bbox=dict(boxstyle='square,pad=0.1', fc='w',
                                            ec='none'))

    # Format plot axes.

    axes.set_xlabel(x_axis_title)

    axes.set_yticks([])

    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)

    return lines
