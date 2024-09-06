"""
Implements a series of functions to create and/or initialize the required
grain size distribution related fields

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager

"""

import numpy as np


def adds_2mm_to_gsd(gsd):
    """Adds the 2 mm fraction to the gsd. This allows faster calculations when
    sand fractions are required or have to be excluded. It is only added
    if the input gsd does not contain the specific 2 mm fraction

    Examples
    --------
    We will show different cases in which the gsd is and is not modified.

    >>> import numpy as np
    >>> from . import _initialize_gsd as initialize_gsd

    Case1: gsd will not be modified because we already know that all sizes are larger than 2 mm

    >>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]
    >>> gsd = initialize_gsd.adds_2mm_to_gsd(gsd)
    >>> gsd
    array([[ 32, 100, 100],
           [ 16,  25,  50],
           [  8,   0,   0]])

    So, in this case nothing happens. Now we will introduce a gsd that needs to include the
    2 mm fraction.

    Case2: gsd will  be modified because not all sizes are larger than 2 mm

    >>> gsd = [[32, 100, 100], [16, 25, 50], [1, 0, 0]]
    >>> gsd = initialize_gsd.adds_2mm_to_gsd(gsd)
    >>> gsd
    array([[ 32.  , 100.  , 100.  ],
           [ 16.  ,  25.  ,  50.  ],
           [  2.  ,   6.25,  12.5 ],
           [  1.  ,   0.  ,   0.  ]])

    """
    gsd = np.array(gsd)
    if np.any(gsd[:, 0] < 2) and not np.isin(2, gsd[:, 0]):
        gs_D = gsd[:, 0]  # Grain sizes
        gs_freq = gsd[:, 1:]  # Grain sizes frequency cumulative
        n_gsd_loc = gsd.shape[1] - 1  # Number of locations with different GSD
        id2mm = np.argmin(gs_D >= 2)  # row where 2 mm will be placed

        # interpolation to get values at 2mm
        i = id2mm - 1
        p_2 = np.log2(2)
        p_i = np.log2(gs_D[i])
        p_i_p1 = np.log2(gs_D[i + 1])
        p = (p_2 - p_i) / (p_i_p1 - p_i)

        f_2mm = np.zeros([1, n_gsd_loc + 1])
        f_2mm[0, 0] = 2
        for j in range(n_gsd_loc):
            f_i = gs_freq[i, j]
            f_i_p1 = gs_freq[i + 1, j]
            f_2mm[0, j + 1] = p * (f_i_p1 - f_i) + f_i

        gsd = np.vstack((gsd[:id2mm, :], f_2mm, gsd[id2mm:, :]))

    return gsd


def remove_sand_from_gsd(gsd, selected_eq):
    """When using Parker Eq. sand needs to be removed from gsd.
    This functions removes sand from gsd

    Examples
    --------
    We will show different cases in which the gsd is and is not modified.

    >>> import numpy as np
    >>> from . import _initialize_gsd as initialize_gsd

    Case1: gsd will not be modified because we already know that there is no sand

    >>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]
    >>> gsd = initialize_gsd.remove_sand_from_gsd(gsd, "Parker1990")
    >>> gsd
    array([[ 32, 100, 100],
           [ 16,  25,  50],
           [  8,   0,   0]])

    So, in this case nothing happens. Now we will introduce a gsd that needs to include the
    2 mm fraction.

    Case2: gsd will  be modified because not all sizes are larger than 2 mm. First we need to
    determine sand fractions, which is done by adding the 2 mm fraction

    >>> gsd = [[32, 100, 100], [16, 25, 50], [1, 0, 0]]
    >>> gsd = initialize_gsd.adds_2mm_to_gsd(gsd)
    >>> gsd
    array([[  32.  ,  100.  ,  100.  ],
           [  16.  ,   25.  ,   50.  ],
           [   2.  ,    6.25,   12.5 ],
           [   1.  ,    0.  ,    0.  ]])

    Now sand is removed

    >>> gsd = initialize_gsd.remove_sand_from_gsd(gsd, "Parker1990")
    >>> np.around(gsd, 2)
    array([[  32.  ,  100.  ,  100.  ],
           [  16.  ,   20.  ,   42.86],
           [   2.  ,    0.  ,    0.  ]])


    """
    gsd = np.array(gsd)
    if selected_eq == "Parker1990" and (np.any(gsd[:, 0] < 2)):
        id2mm = np.argmin(gsd[:, 0] >= 2)
        gsd_no_sand = gsd[:id2mm]
        for j in range(gsd.shape[1] - 1):
            corr_factor = gsd_no_sand[0, j + 1] - gsd_no_sand[-1, j + 1]
            gsd_no_sand[:, j + 1] = (
                (gsd_no_sand[:, j + 1] - gsd_no_sand[-1, j + 1]) * 100 / corr_factor
            )
        gsd = gsd_no_sand

    return gsd


def map_initial_bed_properties_to_nodes(gsd, gsd_loc_node):
    """Grain size distribution information is mapped from
    inputs to each nodes

    Examples
    --------
    In this case we have three locations (0, 1, and 2) which are specified in gsd_loc_node
    which has the same shape than the grid node

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from . import _initialize_gsd as initialize_gsd
    >>> from . import _initialize_fields as initialize

    >>> grid = RasterModelGrid((5, 5))

    >>> gsd_loc_node = [
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 2.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ... ]

    >>> gsd_loc_node = initialize.field_at_node(grid, gsd_loc_node)

    >>> gsd_loc_node
    array([0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 2, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1,
           1, 0])

    So, for example we can see location 0 at nodes 0, 4, 5 etc; location 1 at nodes
    1, 2, 3, etc; and location 2 at node 12.

    >>> gsd_loc_node[12]
    2

    The gsd for these three location is:

    >>> gsd = [
    ...     [32, 100, 100, 100],
    ...     [16, 25, 50, 60],
    ...     [8, 10, 20, 30],
    ...     [4, 5, 4, 2],
    ...     [1, 0, 0, 0],
    ... ]

    We need to go through the steps required in river bed dynamics. First, add the 2 mm fraction

    >>> gsd = initialize_gsd.adds_2mm_to_gsd(gsd)
    >>> gsd
    array([[ 32. , 100. , 100. , 100. ],
           [ 16. ,  25. ,  50. ,  60. ],
           [  8. ,  10. ,  20. ,  30. ],
           [  4. ,   5. ,   4. ,   2. ],
           [  2. ,   2.5,   2. ,   1. ],
           [  1. ,   0. ,   0. ,   0. ]])

    Second, remove sand (only when using Parker 1990)

    >>> gsd = initialize_gsd.remove_sand_from_gsd(gsd, "Parker1990")
    >>> np.around(gsd, 2)
    array([[ 32.  , 100.  , 100.  , 100.  ],
           [ 16.  ,  23.08,  48.98,  59.6 ],
           [  8.  ,   7.69,  18.37,  29.29],
           [  4.  ,   2.56,   2.04,   1.01],
           [  2.  ,   0.  ,   0.  ,   0.  ]])

    And now we can map to nodes

    >>> (
    ...     sand_fract,
    ...     gs_D_freq,
    ...     gs_D_equiv,
    ... ) = initialize_gsd.map_initial_bed_properties_to_nodes(gsd, gsd_loc_node)

    Sand fraction should be zero everywhere because we removed it by using Parker

    >>> sand_fract
    array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
           0., 0., 0., 0., 0., 0., 0., 0.])

    Let's check some of the gsd for selected nodes, representative of locations 1, 2, and 3.
    These nodes are: 5, 8, and 12, respectively

    >>> np.around(gs_D_freq[5], 2)
    array([0.77, 0.15, 0.05, 0.03])

    >>> np.around(gs_D_freq[8], 2)
    array([0.51, 0.31, 0.16, 0.02])

    >>> np.around(gs_D_freq[12], 2)
    array([0.4 , 0.3 , 0.28, 0.01])

    Another example, where the 2 mm fraction is not originally present and only sizes larger
    than 2 mm are present

    >>> gsd = [[32, 100, 100, 100], [16, 25, 50, 60], [8, 0, 0, 0]]
    >>> gsd = initialize_gsd.adds_2mm_to_gsd(gsd)
    >>> gsd
    array([[ 32, 100, 100, 100],
           [ 16,  25,  50,  60],
           [  8,   0,   0,   0]])

    >>> (
    ...     sand_fract,
    ...     gs_D_freq,
    ...     gs_D_equiv,
    ... ) = initialize_gsd.map_initial_bed_properties_to_nodes(gsd, gsd_loc_node)

    Sand fraction should be zero everywhere because there was none in the gsd

    >>> sand_fract
    array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
           0., 0., 0., 0., 0., 0., 0., 0.])

    Let's check some of the gsd for selected nodes, representative of locations 1, 2, and 3.
    These nodes are: 5, 8, and 12, respectively

    >>> np.around(gs_D_freq[5], 2)
    array([0.75, 0.25])

    >>> np.around(gs_D_freq[8], 2)
    array([0.5, 0.5])

    >>> np.around(gs_D_freq[12], 2)
    array([0.4, 0.6])

    """
    gsd = np.array(gsd)
    gsd_loc_node = np.array(gsd_loc_node, dtype=int)

    id2mm = np.where(gsd[:, 0] == 2)[0]

    # if id2mm is empty means that there is no sand
    if id2mm.size == 0:
        sand_fraction_0 = np.zeros([1, gsd.shape[1] - 1])
    else:
        sand_fraction_0 = gsd[id2mm, 1:]

    gs_freq = np.abs(
        -np.diff(gsd[:, 1:] / 100, axis=0)
    )  # Grain sizes frequency - Now in fraction
    gs_D_eq = (gsd[0:-1, 0] * gsd[1:, 0]) ** 0.5  # Equivalent grain sizes
    sand_fraction = np.zeros_like(
        gsd_loc_node, dtype=float
    )  # Sand fraction at each node

    # Bed grain sizes frequency in each node
    gs_D_equiv_freq = np.zeros([gsd_loc_node.size, gs_D_eq.shape[0]])

    for i in range(gs_freq.shape[1]):
        (id_gsd_loc) = np.where(gsd_loc_node == i)
        sand_fraction[id_gsd_loc] = sand_fraction_0[0, i] / 100
        gs_D_equiv_freq[id_gsd_loc, :] = gs_freq[:, i]

    return sand_fraction, gs_D_equiv_freq, gs_D_eq
