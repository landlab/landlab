# BiotaEvolver plot functions.

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


def plot_tree(tree, axes=None, x_multiplier=0.001, x_axis_title='Time (ky)'):
    """ Plot a phylogenetic tree.

    """
    times = list(tree.keys())

    # Get the species in each clade.
    clades = []
    for t in times:
        for p in tree[t].keys():
            for s in tree[t][p]:
                clades.append(s.identifier[0])
#                cid = -1
#                species = tree[t][-1]
#                for s in species:
#                    cid = s.identifier[0]
#                    clades.setdefault(cid, set())
#                    clades[cid].add(s)
#            else:
#                cid = p.identifier[0]
#                clades.setdefault(cid, set())
#                for s in tree[t][p]:
#                    clades[cid].add(s)
    clades = list(set(clades))
    clades.sort(reverse=True)
    print(clades)
    # Plot the tree beginning at final time.

    if axes == None:
        fig = plt.figure('Phylogenetic tree')
        axes = fig.add_axes(plt.axes())

    # Reverse time.
    times = times[::-1]

    y_clade_spacing = 0.5
    y_species_spacing = 0.5
    y_max = 0

    lines = []
    species_position = {}

    for cid in clades:
        y_max += y_clade_spacing

        for i, time in enumerate(times):
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

#            pid = []
#            unsort_species = list(clades[cid])
#            for p in unsort_species:
#                if p.parent_species == -1:
#                    pid.append(-1)
#                else:
#                    pid.append(p.identifier[1])
#            print(pid)
#            sorted_indices = np.argsort(pid)
#            print(np.array(unsort_species),[[sorted_indices]])
#            sorted_species = np.array(unsort_species)[[sorted_indices]]

            y_clade_species = []

            for time_clade in tree[time]:
                if time_clade == -1:
                    tid = -1
                else:
                    tid = time_clade.identifier[0]
                print(1,tid,cid)
                if tid is not cid:
                    print(2)
                    continue
                print(3)
                y_time_clade_species = []

                pid = []
                unsort_species = tree[time][time_clade]
                for p in unsort_species:
                    if p.parent_species == -1:
                        pid.append(-1)
                    else:
                        pid.append(p.identifier[1])
                print(pid)
                sorted_indices = np.argsort(pid)
                print(np.array(unsort_species),[[sorted_indices]])
                sorted_species = np.array(unsort_species)[[sorted_indices]]

                for species in sorted_species:
                    existed_this_time = time in species.record.keys()
                    existed_prior_time = earlier_time in species.record.keys()
                    existed_later_time = later_time in species.record.keys()
                    no_parent = species.parent_species == -1

                    if existed_this_time:
                        if species in species_position.keys():
                            # Previously plotted.
                            y_next = species_position[species]
                        else:
                            # First time plotting.
                            y_next = y_max + y_species_spacing

                        if no_parent:
                            lines.extend(axes.plot([x_min, x_max], [y_next, y_next], 'm'))
                        else:
                            lines.extend(axes.plot([x_min, x_max], [y_next, y_next], 'c'))

                        species_position[species] = y_next

                        if time == min(species.record.keys()):
                            y_clade_species.append(y_next)
                            y_time_clade_species.append(y_next)

                        if y_next > y_max:
                            y_max = y_next

    #                    if not existed_later_time or time == max(times):
                        sid = '{}.{}'.format(species.identifier[0],
                               species.identifier[1])
                        plt.text(x_max + (x_max * 0.01), y_next, sid,
                                 verticalalignment='center', fontsize='small')

                if not existed_prior_time:
                    species_position[species.parent_species] = np.mean(y_time_clade_species)

                if len(y_time_clade_species) > 0:
                    lines.extend(axes.plot([x_min, x_min],
                                           [min(y_time_clade_species),
                                            max(y_time_clade_species)], 'y'))

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

    # Format plot axes.

    axes.set_xlabel(x_axis_title)
    axes.set_ylabel('Species')

    axes.set_yticks([])

    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)

    return lines
