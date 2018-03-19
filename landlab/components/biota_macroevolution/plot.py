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

def plot_number_of_species(record, axes=None):

    # Prepare figure.
    if axes == None:
        fig = plt.figure('Number of species')
        axes = fig.add_axes(plt.axes())

    times = list(record.keys())
    count = []
    for t in times:
        count.append(len(record[t]['species']))

    times_in_ky = np.multiply(times, 1e-3)

    axes.plot(times_in_ky, count, 'k')

    axes.set_xlabel('Time (ky)')
    axes.set_ylabel('Number of species')

def plot_tree(tree, x_multiplier=0.001, selected_species=None, axes=None):

    # Prepare figure.
    if axes == None:
        fig = plt.figure('Phylogenetic tree')
        axes = fig.add_axes(plt.axes())

    y_clade_spacing = 5
    y_species_spacing = 0.5
    y_next = 0
    y_max = 0

    species_position = {}
    times = list(tree.keys())

    lines = []

    # Plot the tree beginning at final time.
    clades = {}
    for t in times:
        for c in tree[t].keys():
            s = tree[t][c]
            clades.setdefault(c, [])
            clades[c] = list(set(clades[c]) | set(s))
    from pprint import pprint


    root_species = clades[-1]
    clades.pop(-1, None)
    pprint(root_species)
#    pid = np.array([s.identifier for s in root_species])

    for i, clade in enumerate(clades.keys()):
        for time in times:
            # Get the duration represented by this time.
            if time == max(times):
                later_time = time
                earlier_time = times[i + 1]
            elif time == min(times):
                later_time = times[i - 1]
                earlier_time = time
            else:
                later_time = times[i - 1]
                earlier_time = times[i + 1]

            # Extend the x values by half of time.
            x_max = (time + (later_time - time) * 0.5) * x_multiplier
            x_min = (time - (time - earlier_time) * 0.5) * x_multiplier
            x_mid = np.mean([x_min, x_max])
            root_parents = []
            for species in clades[clade]:
                existed_this_time = time in species.record.keys()
                existed_prior_time = earlier_time in species.record.keys()
                root_parents.append(clade)

                if existed_this_time:
                    if species.identifier in species_position.keys():
                        y_next = species_position[species.identifier]
                    else:
                        y_next += y_species_spacing

                    lines.extend(axes.plot([x_mid, x_max], [y_next, y_next], 'c'))

                    if existed_prior_time:
                        species_position[species.identifier] = y_next
                    else:
                        species_position[clade] = y_next

            for p in root_parents:
                existed_this_time = time in clade.record.keys()

                lines.extend(axes.plot([x_mid, x_max], [y_next, y_next], 'c'))

        y_next += y_clade_spacing

#    for i, time in enumerate(times):
#        # Get the duration represented by this time.
#        if time == max(times):
#            later_time = time
#            earlier_time = times[i + 1]
#        elif time == min(times):
#            later_time = times[i - 1]
#            earlier_time = time
#        else:
#            later_time = times[i - 1]
#            earlier_time = times[i + 1]
#
#        # Extend the x values by half of time.
#        x_max = (time + (later_time - time) * 0.5) * x_multiplier
#        x_min = (time - (time - earlier_time) * 0.5) * x_multiplier
#        x_mid = np.mean([x_min, x_max])
#
#
#        # Plot branches for each clade, where the clade includes the species
#        # that share a common parent species at the prior time.
#        for parent_id, species_list in tree[time].items():
#            clade_branch = False
#
#            y_clade_species = []
#
#            for species in species_list:
#                # Determine if species continues across timesteps (a continuing
#                # tree branch).
#
#                existed_prior_time = earlier_time in species.record.keys()
#                no_parent = parent_id == -1
#
#                if existed_prior_time or no_parent:
#                    x = x_min
#                else:
#                    clade_branch = True
#                    x = x_mid
#
#                if species.identifier in species_position.keys():
#                    # Clade trunk or branch.
#                    y_next = species_position[species.identifier]
#                    lines.extend(axes.plot([x, x_max], [y_next, y_next], 'k'))
#                elif no_parent:
#                    # Leaf node of clade with no branches.
#                    y_next = y_max + y_clade_spacing
#                    lines.extend(axes.plot([x, x_max], [y_next, y_next], 'k'))
#                else:
#                    if len(y_clade_species) > 0:
#                        # First leaf node of branched clade.
#                        y_next += y_clade_spacing * 0.5
#                    else:
#                        y_next += 1
#                    lines.extend(axes.plot([x, x_max], [y_next, y_next], 'k'))
#
##                if time == 0:
##                    print(species.identifier, y_next)
#
#                y_clade_species.append(y_next)
#
#                # Plot branches.
#
#                # Retain y position for other times of species.
#                if no_parent or not clade_branch:
#                    species_position[species.identifier] = y_next
#                else:
#                    species_position[parent_id] = np.mean(y_clade_species)
#
#                if y_next > y_max:
#                    y_max = y_next
#
#            if clade_branch:
#                # Draw the line that connects the clade species.
#                if time != 0:
#                    l = axes.plot([x_mid] * 2,
#                              [min(y_clade_species), max(y_clade_species)],
#                              'k')
#                    lines.extend(l)
#
#                # Draw the line that connects the above line with the trunk.
#                y_mid = np.mean(y_clade_species)
#                lines.extend(axes.plot([x_min, x_mid], [y_mid, y_mid], 'k'))
#
#            y_next = max(y_clade_species)

    return lines

    # Format plot.

    axes.set_xlabel('Time (ky)')
    axes.set_ylabel('Species')

    axes.set_xlim([0, max(times) * x_multiplier])
    axes.set_ylim(bottom=0)

    # Show only the x-axis.
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
