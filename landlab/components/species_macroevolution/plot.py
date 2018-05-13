# SpeciesEvolver plot functions.

from copy import deepcopy
from landlab.plot import imshow_grid
from matplotlib.ticker import MaxNLocator
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame


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

    cmap = colors.LinearSegmentedColormap.from_list('streamOverlay', [c1, c2],
                                                    2)
    cmap._init() # create the _lut array, with rgba values
    alphas = np.linspace(0, 1, cmap.N + 3)
    cmap._lut[:,-1] = alphas

    plt.figure('Species range')

    imshow_grid(grid, 'topographic__elevation', cmap='gray')
    imshow_grid(grid, combined_range_mask, cmap=cmap,
                    allow_colorbar=False)


def number_of_species(record, species_DataFrame, axes=None, time_max=None):
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

    times = record.time.tolist()
    count = []

    sdf = species_DataFrame

    for time in times:
        appeared_before_time = sdf.time_appeared <= time
        present_at_time = sdf.latest_time >= time
        extant_at_time = np.all([appeared_before_time, present_at_time], 0)

        species_time = sdf.object[extant_at_time].tolist()

        count.append(len(species_time))

    times_in_ky = np.multiply(times, 1e-3)

    axes.plot(times_in_ky, count, 'k')

    axes.set_xlabel('Time (ky)')
    axes.set_ylabel('Number of species')

    # Force y-axis label values to be integers.
    axes.yaxis.set_major_locator(MaxNLocator(integer=True))


def imshow_zones(mg, zones):
    """

    """
    for z in zones:
        c = (z.plot_color[0], z.plot_color[1], z.plot_color[2], 0.7)
        cmap = colors.LinearSegmentedColormap.from_list('streamOverlay',
                                                        [(1,1,1,0), c], 2)
        imshow_grid(mg, z.mask, cmap=cmap, allow_colorbar=False)


def tree(species_DataFrame, axes=None, x_multiplier=0.001,
         x_axis_title='Time (ky)', time_min=None, time_max=None):
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
    time_min : float, optional
        A
    time_max : float, optional

    """
    # Add the parent species number to a species DataFrame copy.

    sdf = deepcopy(species_DataFrame)

    for i, row in sdf.iterrows():
        ps = row.object.parent_species
        if ps == -1:
            sdf.loc[i, 'parent_number'] = -1
        else:
            sdf.loc[i, 'parent_number'] = ps.identifier[1]

    # Create matplotlib axes if axes were not inputted.

    if axes == None:
        title = 'Phylogenetic of species'
        plt.close(title)
        fig = plt.figure(title)
        axes = fig.add_axes(plt.axes())

    # Iterate clades in reverse alphabetical order for plotting.

    clades = sdf.clade.unique().tolist()
    clades.sort(reverse=True)

    # Initialize plotting parameters.

    y_species_spacing = 0.5
    y_max = 0
    duration = sdf.latest_time.max() - sdf.time_appeared.min()
    label_offset = duration * x_multiplier * 0.01

    # Construct the tree by clade then by species.

    for clade in clades:
        species_clade = sdf.loc[sdf.clade == clade]
        species_clade = species_clade.sort_values(['parent_number',
                                                   'latest_time',
                                                   'species'], ascending=False)

        sc = DataFrame(columns=['number', 'x_min', 'x_max', 'y', 'subclades',
                                'children'])

        # Store y-axis position of sibling species to connect branches.
        y_sibling = {}

        for i, row in species_clade.iterrows():
            # Get species number, sn and parent number, pn.
            sn = row.species
            pn = row.parent_number

            # Get plot x, y of species.

            x_min = row.time_appeared * x_multiplier
            x_max = row.latest_time * x_multiplier

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

            data = {'number': sn, 'x_min': x_min, 'x_max': x_max, 'y': y_next,
                    'subclades': [pn]}
            sc = sc.append(DataFrame(data), ignore_index=True)

        # Plot tree.

        # Label clade.
        root = sc.loc[sc.number == 0]
        plt.text(root.x_min - label_offset, root.y, clade, fontsize='small',
                 ha='right', va='center')

        for i, row in sc.iterrows():
            # Plot tree edge.
            if row.x_max - row.x_min == 0:
                axes.plot(row.x_min, row.y, 'k.', markersize=3)
            else:
                axes.plot([row.x_min, row.x_max], 2 * [row.y], 'k')

            # Label species.
            plt.text(row.x_max + label_offset, row.y, row.number,
                     fontsize='small', va='center')

            # Plot inter timestep interpolation.
            if row.number in y_sibling.keys():
                y_smax = max(y_sibling[row.number])
                y_smin = min(y_sibling[row.number])
                y_smean = np.mean(y_sibling[row.number])

                x_max = sc.loc[sc.subclades == row.number].x_min.min()

                axes.plot(2 * [x_max], [y_smin, y_smax], color='0.8', zorder=0)
                axes.plot([row.x_max, x_max], 2 * [y_smean], color='0.8',
                          zorder=0)

        print(clade, '---------')
        print(sc)

    # Format plot axes.

    if time_min:
        axes.set_xlim(left=time_min * x_multiplier)
    if time_max:
        axes.set_xlim(right=time_max * x_multiplier)

    axes.set_xlabel(x_axis_title)

    axes.set_yticks([])

    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)





def plot_treex(species, times, axes=None, x_multiplier=0.001,
              x_axis_title='Time (ky)', max_plot_time=None):
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

    clades = species.clade.unique().tolist()
    clades.sort(reverse=True)

    times.sort(reverse=True)

    # Create matplotlib axes if axes were not inputted.

    if axes == None:
        title = 'Phylogenetic of species'
        plt.close(title)
        fig = plt.figure(title)
        axes = fig.add_axes(plt.axes())

    # Initialize plotting parameters.

    y_species_spacing = 0.5
    y_max = 0
    label_offset = (max(times) - min(times)) * x_multiplier * 0.01

    # Construct the tree by clade then by time.

    sc = DataFrame(columns=['number', 'x_min', 'x_max', 'y', 'subclades'])
    sc.subclades = [np.nan]

    for clade in clades:
        species_clade = species.loc[species.clade == clade]

        # Store y-axis position of species to connect branches over time.
        y_species = {}


        for i, time in enumerate(times):
            # Get the temporal bounds represented by this time.

            if time == max(times):
                time_later = max(times) + 1
                time_prior = times[i + 1]
            elif time == min(times):
                time_later = times[i - 1]
                time_prior = min(times) - 1
            else:
                time_later = times[i - 1]
                time_prior = times[i + 1]

            time_is_final = time == times[0]

            # Get the clade species extant at this time.

            appeared_before_time = species_clade.time_appeared <= time
            present_at_time = species_clade.latest_time >= time
            time_mask = np.all([appeared_before_time.tolist(),
                                present_at_time.tolist()], 0)
            species_time = species_clade.loc[time_mask]
            species_time = species_time.sort_values(['parent_id',
                                                     'species_number'],
                                                     ascending=False)

            # Set time position variables.

            x_time = time * x_multiplier
            x_time_prior = time_prior * x_multiplier
            y_time = []
            y_subclade = {}
            y_sibling = {}

            for i, row in species_time.iterrows():
                # Get species object, s of current row.
                s = row.object
                sid = s.identifier[1]

                if s in y_species.keys():
                    # Previously plotted.
                    y_next = y_species[s]
                elif s.parent_species in y_sibling.keys():
                    # Previously plotted.
                    y_next = max(y_sibling[s.parent_species]) + y_species_spacing
                else:
                    # First time plotting.
                    y_next = y_max + y_species_spacing

                y_species[s] = y_next

                y_sibling.setdefault(s.parent_species, [])
                y_sibling[s.parent_species].append(y_next)

                y_time.append(y_next)

                if y_next > y_max:
                    y_max = y_next

                exists_earlier_time = row.time_appeared <= time_prior
                exists_later_time = row.latest_time >= time_later

                if exists_earlier_time:
                    axes.plot([x_time_prior, x_time], 2 * [y_next], 'k')
                elif not exists_later_time:
                    axes.plot(x_time, y_next, 'k.', markersize=3)

                if not exists_earlier_time:
                    y_subclade.setdefault(s.parent_species, [])
                    y_subclade[s.parent_species].append(y_next)

                # Plot species number.
                if not exists_later_time or time_is_final:

                    plt.text(x_time + label_offset, y_next, sid,
                             fontsize='small', va='center')

                if s.parent_species == -1:
                    pid = -1
                else:
                    pid = s.parent_species.identifier[1]

                data = {'number': sid, 'x_min': x_time_prior,
                                'x_max': x_time, 'y': y_next,
                                'subclades': [pid]}
                sc = sc.append(DataFrame(data), ignore_index=True)

            # Plot lines that connects sibling species to parent species.
            for k in y_subclade.keys():
                if k not in species_time.object.tolist():
                    y_smax = max(y_subclade[k])
                    y_smin = min(y_subclade[k])
                    y_smean = np.mean(y_subclade[k])
    #            if clade_branch:
                    axes.plot(2 * [x_time], [y_smin, y_smax],
                              color='0.8', zorder=0)
                    axes.plot([x_time, x_time_prior], 2 * [y_smean],
                              color='0.8', zorder=0)

                    y_species[k] = y_smean

        # Plot clade id.
        if time == 0:
            plt.text(x_time_prior - label_offset, y_next, clade,
                     fontsize='small', ha='right', va='center')

    # Format plot axes.
    print(sc)
    if max_plot_time:
        axes.set_xlim(right=max_plot_time * x_multiplier)

    axes.set_xlabel(x_axis_title)

    axes.set_yticks([])

    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)

    return axes
