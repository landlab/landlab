# SpeciesEvolver plot functions.
from copy import deepcopy

from matplotlib.ticker import MaxNLocator
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from pandas import concat, DataFrame

from landlab.plot import imshow_grid


def plot_area_versus_species(area, number_of_species, axes=None):
    """ Plot the number of species as a function of zone area.

    Parameters
    ----------
    species_DataFrame : Pandas DataFrame
        The
    zone_DataFrame : Pandas DataFrame
        A
    time : float
        A
    axes : Matplotlob axes, optional
        A
    """
    # Prepare figure.
    if axes == None:
        title = 'Area versus number of species'
        plt.close(title)
        fig = plt.figure(title)
        axes = fig.add_axes(plt.axes())

    plt.plot(area, number_of_species, 'k.')
    plt.xlabel('Stream area ($m^2$)')
    plt.ylabel('Number of species')


def plot_number_of_species(record, species_DataFrame, axes=None,
                           x_multiplier=0.001, x_axis_title='Time (ky)',
                           max_time_to_plot=None):
    """Plot the number of species as a function of time.

    Parameters
    ----------
    record : Pandas DataFrame
        The
    species_DataFrame : Pandas DataFrame
        A
    axes : Matplotlob axes, optional
        A
    """
    # Prepare figure.
    if axes == None:
        title = 'Number of species'
        plt.close(title)
        fig = plt.figure(title)
        axes = fig.add_axes(plt.axes())

    count = []

    sdf = species_DataFrame

    for time in record.times:
        appeared_before_time = sdf.time_appeared <= time
        present_at_time = sdf.latest_time >= time
        extant_at_time = np.all([appeared_before_time, present_at_time], 0)

        species_time = sdf.object[extant_at_time].tolist()

        count.append(len(species_time))

    times = np.multiply(record.times, x_multiplier)

    max_time = times.max()
    min_count = min(count)
    max_count = max(count)

    if max_time_to_plot:
        i = np.array(record.times) <= max_time_to_plot
        times = times[i]
        count = np.array(count)[i]

    axes.plot(times, count, 'k')

    axes.set_xlim(left=0, right=max_time)
    axes.set_ylim(bottom=np.floor(min_count * 0.9),
                  top=np.ceil(max_count * 1.1))
    axes.set_xlabel(x_axis_title)
    axes.set_ylabel('Number of species')

    # Force y-axis label values to be integers.
    axes.yaxis.set_major_locator(MaxNLocator(integer=True))


def plot_zones(mg, zones):
    """

    """
    for z in zones:
        c = (z.plot_color[0], z.plot_color[1], z.plot_color[2], 0.7)
        cmap = colors.LinearSegmentedColormap.from_list('streamOverlay',
                                                        [(1,1,1,0), c], 2)
        imshow_grid(mg, z.mask, cmap=cmap, allow_colorbar=False)


def plot_tree(species_DataFrame, axes=None, x_multiplier=0.001,
         x_axis_title='Time (ky)', max_time_to_plot=None):
    """Plot a phylogenetic tree.

    Parameters
    ----------
    species_DataFrame : Pandas DataFrame
        The
    axes : Matplotlob axes, optional
        A
    x_multiplier : float, optional
        A
    x_axis_title : string, optional
        A
    max_time_to_plot : float, optional
        A
    """
    # Add the parent species number to a species DataFrame copy.

    sdf = deepcopy(species_DataFrame)

    for i, row in sdf.iterrows():
        ps = row.object.parent_species
        if ps == -1:
            sdf.loc[i, 'parent_number'] = -1
        else:
            sdf.loc[i, 'parent_number'] = ps.identifier[1]

    # Create matplotlib axes if none exists.

    if axes == None:
        title = 'Phylogenetic of species'
        plt.close(title)
        fig = plt.figure(title)
        axes = fig.add_axes(plt.axes())

    # Iterate clades in reverse alphabetical order for plotting.

    clades = sdf.clade.unique().tolist()
    clades.sort(reverse=True)

    # Initialize plotting parameters.
    duration = sdf.latest_time.max() - sdf.latest_time.min()
    label_offset = duration * x_multiplier * 0.05

    if max_time_to_plot == None:
        species_max_x = sdf.latest_time.max() * x_multiplier
    else:
        species_max_x = max_time_to_plot * x_multiplier

    y_species_spacing = 0.5
    y_max = 0

    # Add columns to position species on tree.
    sdf2 = DataFrame(columns=['number', 'x_min', 'x_max', 'y', 'subclades',
                              'children'])
    sdf = concat([sdf, sdf2])

    # Construct the tree by clade then by species.

    for clade in clades:
        # Get the species of the clade, sc.
        sc = sdf.loc[sdf.clade == clade]

        sc = sc.sort_values(['parent_number', 'latest_time', 'species'],
                            ascending=False)

        # Store y-axis position of sibling species to connect branches.
        y_sibling = {}

        sn_to_delete = []

        for i, row in sc.iterrows():
            # Get species number, sn and parent number, pn.
            sn = row.species
            pn = row.parent_number

            if row.time_appeared > max_time_to_plot:
                sn_to_delete.append(sn)
                continue

            # x position is set by

            x_min = row.time_appeared * x_multiplier

            if row.latest_time > max_time_to_plot:
                x_max = species_max_x
            else:
                x_max = row.latest_time * x_multiplier

            # y position is set by species history.

            if sn in y_sibling.keys():
                # Set y by the mean y children species.
                y_next = np.mean(y_sibling[sn])

            elif pn in y_sibling.keys():
                # Set y above a sibling of this species.
                y_next = max(y_sibling[pn]) + y_species_spacing
            else:
                # First species set in clade.
                y_next = y_max + y_species_spacing

            y_sibling.setdefault(pn, [])
            y_sibling[pn].append(y_next)

            if y_next > y_max:
                y_max = y_next

            cols = ['x_min', 'x_max', 'y', 'subclades']
            sc.loc[sn == sc.species, cols] = [x_min, x_max, y_next, pn]

        for n in sn_to_delete:
            sc.drop(sc[sc.species == n].index, inplace=True)

        if len(sc) == 0:
            continue

        # Plot tree.

        # Label clade.
        root = sc.loc[sc.species == 0]
        axes.text(root.x_min - label_offset, root.y, clade, fontsize='small',
                  ha='right', va='center')

        for i, row in sc.iterrows():
            # Plot tree edge.
            if row.x_max - row.x_min == 0:
                axes.plot(row.x_min, row.y, 'k.', markersize=3)
            else:
                axes.plot([row.x_min, row.x_max], 2 * [row.y], 'k')

            # Label species.
            axes.text(row.x_max + label_offset, row.y, str(row.species),
                      fontsize='small', va='center')

            # Plot inter timestep interpolation.
            if row.species in y_sibling.keys():
                y_smax = max(y_sibling[row.species])
                y_smin = min(y_sibling[row.species])
                y_smean = np.mean(y_sibling[row.species])

                x_max = sc.loc[sc.subclades == row.species].x_min.min()

                axes.plot(2 * [x_max], [y_smin, y_smax], color='0.8', zorder=0)
                axes.plot([row.x_max, x_max], 2 * [y_smean], color='0.8',
                          zorder=0)

#        print(clade, '---------')
#        print(sc)

    # Format plot axes.

    axes.set_xlim(left=0, right=sdf.latest_time.max() * x_multiplier)

    axes.set_xlabel(x_axis_title)

    axes.set_yticks([])

    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)
