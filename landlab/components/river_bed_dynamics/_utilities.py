"""

This file contain a series of functions that helps with post-processing nodes
and links gsd that are results at different places in the river_bed_dynamics
landlab component. Examples of applications are given in each function.

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager

"""

import numpy as np
import pandas as pd


def vector_mapper(grid, vector):
    """Map vectors values from links onto nodes preserving the components.
    This is intended for graphic representation. We used a weighted average to account
    for direction and magnitude of each link.

    This method takes the vectors values on links and determines the
    vectors components in nodes.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import RiverBedDynamics
    >>> from . import _utilities as utilities
    >>> import matplotlib.pyplot as plt
    >>> from landlab import imshow_grid

    >>> grid = RasterModelGrid((5, 5))

    >>> u = np.full(grid.number_of_links, 0.0)

    >>> u[20] = 0.12
    >>> u[24] = 0.11
    >>> u[19] = -0.02
    >>> u[15] = -0.01

    >>> (u_vector, u_magnitude) = utilities.vector_mapper(grid, u)
    >>> u_x = u_vector[:, 0]
    >>> u_x.reshape(grid.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.  , -0.01,  0.05,  0.06,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])

    >>> u_y = u_vector[:, 1]
    >>> u_y.reshape(grid.shape)
    array([[ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ],
           [ 0.   ,  0.   , -0.005,  0.   ,  0.   ],
           [ 0.   ,  0.   ,  0.05 ,  0.   ,  0.   ],
           [ 0.   ,  0.   ,  0.055,  0.   ,  0.   ],
           [ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ]])

    To isually check the result use:
    imshow_grid(grid, u_magnitude,cmap='Blues',vmin=0,vmax=0.1)
    plt.quiver(grid.x_of_node,grid.y_of_node,u_x,u_y,scale = 0.75)

    """
    vector_x_r = vector[grid.links_at_node[:, 0]]
    vector_x_l = vector[grid.links_at_node[:, 2]]
    vector_y_t = vector[grid.links_at_node[:, 1]]
    vector_y_b = vector[grid.links_at_node[:, 3]]

    vector_x = 0.5 * (vector_x_r + vector_x_l)
    vector_y = 0.5 * (vector_y_t + vector_y_b)

    vector = np.transpose(np.vstack((vector_x, vector_y)))
    magnitude = np.sqrt(vector_x**2 + vector_y**2)

    return vector, magnitude


def map_gsd_from_link_to_node(self, location="bed_surf"):
    """Maps grain size distribution from links onto nodes.

    Given that all our calculations are conducted in links we implemented
    this function to display results in a raster or in nodes.
    default type is bedload, alternative type='bedload'

    Links with no bedload transport are ignored during mapping process

    Examples
    --------
    This is the same base example described extensively in river bed dynamics, so
    we removed comments that are already available in the main component

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import RiverBedDynamics
    >>> from . import _utilities as utilities

    >>> grid = RasterModelGrid((5, 5))

    >>> grid.at_node["topographic__elevation"] = [
    ...     [1.07, 1.06, 1.00, 1.06, 1.07],
    ...     [1.08, 1.07, 1.03, 1.07, 1.08],
    ...     [1.09, 1.08, 1.07, 1.08, 1.09],
    ...     [1.09, 1.09, 1.08, 1.09, 1.09],
    ...     [1.09, 1.09, 1.09, 1.09, 1.09],
    ... ]

    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])

    >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
    >>> grid.at_node["surface_water__velocity"] = np.full(grid.number_of_nodes, 0.25)
    >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.102)
    >>> grid.at_link["surface_water__velocity"] = np.full(grid.number_of_links, 0.25)

    >>> gsd_loc = [
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ... ]

    >>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]

    >>> rbd = RiverBedDynamics(
    ...     grid,
    ...     gsd=gsd,
    ...     bedload_equation="Parker1990",
    ...     bed_surf__gsd_loc_node=gsd_loc,
    ... )
    >>> rbd.run_one_step()

    All gsd calculations are done in links. Here we map the results from links to nodes
    For the bed surface gsd we have:

    >>> bed_surf__gsd_node = utilities.map_gsd_from_link_to_node(rbd)
    >>> np.round(bed_surf__gsd_node, 3)
    array([[ 0.656,  0.344],
           [ 0.562,  0.438],
           [ 0.531,  0.469],
           [ 0.562,  0.438],
           [ 0.656,  0.344],
           [ 0.688,  0.312],
           [ 0.531,  0.469],
           [ 0.5  ,  0.5  ],
           [ 0.531,  0.469],
           [ 0.688,  0.312],
           [ 0.688,  0.312],
           [ 0.531,  0.469],
           [ 0.5  ,  0.5  ],
           [ 0.531,  0.469],
           [ 0.688,  0.312],
           [ 0.688,  0.312],
           [ 0.531,  0.469],
           [ 0.5  ,  0.5  ],
           [ 0.531,  0.469],
           [ 0.688,  0.312],
           [ 0.656,  0.344],
           [ 0.562,  0.438],
           [ 0.531,  0.469],
           [ 0.562,  0.438],
           [ 0.656,  0.344]])

    And for the bed load gsd we have:

    >>> sed_transp__bedload_gsd_node = utilities.map_gsd_from_link_to_node(
    ...     rbd, location="bedload"
    ... )
    >>> np.round(sed_transp__bedload_gsd_node, 3)
    array([[ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.379,  0.621],
           [ 0.475,  0.525],
           [ 0.379,  0.621],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.282,  0.718],
           [ 0.33 ,  0.67 ],
           [ 0.282,  0.718],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.282,  0.718],
           [ 0.282,  0.718],
           [ 0.282,  0.718],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ],
           [ 0.   ,  0.   ]])

    A bed load gsd that shows [ 0.   ,  0.   ] means that there is zero transport at that node

    """
    gsd_nodes = np.zeros([self._grid.number_of_nodes, self._gsd.shape[0] - 1])

    if location == "bed_surf":
        gsd_l = self._bed_surf__gsd_link
    else:
        gsd_l = self._sed_transp__bedload_gsd_link

    # split into grain sizes
    for i in range(self._gsd.shape[0] - 1):
        for j in range(self._grid.number_of_nodes):
            gsd_nodes_0 = gsd_l[self._grid.links_at_node[j]][:, i]
            mask = gsd_nodes_0 > 0.0
            gsd_nodes_0 = gsd_nodes_0[mask]
            if gsd_nodes_0.size > 0:
                gsd_nodes[j, i] = np.mean(gsd_nodes_0)

    if location == "bed_surf":
        # Revert any changes to the fixed GSD nodes and fixed elevations nodes
        gsd_nodes[self._bed_surf__gsd_fix_node_id] = self._bed_surf__gsd_orig_node[
            self._bed_surf__gsd_fix_node_id
        ]

    return gsd_nodes


def format_gsd(self, bedload_gsd):
    """Gives a more friendly format for the bed surface or bed load GSD.
    Reads a bed load GSD, from links or nodes, and returns the GSD in
    cumulative percetage

    Examples
    --------
    This is the same base example described extensively in river bed dynamics, so
    we removed comments that are already available in the main component

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import RiverBedDynamics
    >>> from . import _utilities as utilities

    >>> grid = RasterModelGrid((5, 5))

    >>> grid.at_node["topographic__elevation"] = [
    ...     [1.07, 1.06, 1.00, 1.06, 1.07],
    ...     [1.08, 1.07, 1.03, 1.07, 1.08],
    ...     [1.09, 1.08, 1.07, 1.08, 1.09],
    ...     [1.09, 1.09, 1.08, 1.09, 1.09],
    ...     [1.09, 1.09, 1.09, 1.09, 1.09],
    ... ]

    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])

    >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
    >>> grid.at_node["surface_water__velocity"] = np.full(grid.number_of_nodes, 0.25)
    >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.102)
    >>> grid.at_link["surface_water__velocity"] = np.full(grid.number_of_links, 0.25)

    >>> gsd_loc = [
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ... ]

    >>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]

    >>> rbd = RiverBedDynamics(
    ...     grid,
    ...     gsd=gsd,
    ...     bedload_equation="Parker1990",
    ...     bed_surf__gsd_loc_node=gsd_loc,
    ... )
    >>> rbd.run_one_step()

    >>> utilities.format_gsd(rbd, rbd._sed_transp__bedload_gsd_link)
               8      16     32
    Link_0   0.0   0.000    0.0
    Link_1   0.0   0.000    0.0
    Link_2   0.0   0.000    0.0
    Link_3   0.0   0.000    0.0
    Link_4   0.0   0.000    0.0
    Link_5   0.0   0.000    0.0
    Link_6   0.0   0.000    0.0
    Link_7   0.0   0.000    0.0
    Link_8   0.0   0.000    0.0
    Link_9   0.0   0.000    0.0
    Link_10  0.0  52.498  100.0
    Link_11  0.0  52.498  100.0
    Link_12  0.0   0.000    0.0
    Link_13  0.0   0.000    0.0
    Link_14  0.0  71.774  100.0
    Link_15  0.0  52.498  100.0
    Link_16  0.0  71.774  100.0
    Link_17  0.0   0.000    0.0
    Link_18  0.0   0.000    0.0
    Link_19  0.0  71.774  100.0
    Link_20  0.0  71.774  100.0
    Link_21  0.0   0.000    0.0
    Link_22  0.0   0.000    0.0
    Link_23  0.0  71.774  100.0
    Link_24  0.0  71.774  100.0
    Link_25  0.0  71.774  100.0
    Link_26  0.0   0.000    0.0
    Link_27  0.0   0.000    0.0
    Link_28  0.0  71.774  100.0
    Link_29  0.0  71.774  100.0
    Link_30  0.0   0.000    0.0
    Link_31  0.0   0.000    0.0
    Link_32  0.0   0.000    0.0
    Link_33  0.0   0.000    0.0
    Link_34  0.0   0.000    0.0
    Link_35  0.0   0.000    0.0
    Link_36  0.0   0.000    0.0
    Link_37  0.0   0.000    0.0
    Link_38  0.0   0.000    0.0
    Link_39  0.0   0.000    0.0

    It also works for nodes. First we need to map bed load gsd from links to nodes

    >>> sed_transp__bedload_gsd_node = utilities.map_gsd_from_link_to_node(
    ...     rbd, location="bedload"
    ... )

    And now we can use the function

    >>> utilities.format_gsd(rbd, sed_transp__bedload_gsd_node)
               8      16     32
    Node_0   0.0   0.000    0.0
    Node_1   0.0   0.000    0.0
    Node_2   0.0   0.000    0.0
    Node_3   0.0   0.000    0.0
    Node_4   0.0   0.000    0.0
    Node_5   0.0   0.000    0.0
    Node_6   0.0  62.136  100.0
    Node_7   0.0  52.498  100.0
    Node_8   0.0  62.136  100.0
    Node_9   0.0   0.000    0.0
    Node_10  0.0   0.000    0.0
    Node_11  0.0  71.774  100.0
    Node_12  0.0  66.955  100.0
    Node_13  0.0  71.774  100.0
    Node_14  0.0   0.000    0.0
    Node_15  0.0   0.000    0.0
    Node_16  0.0  71.774  100.0
    Node_17  0.0  71.774  100.0
    Node_18  0.0  71.774  100.0
    Node_19  0.0   0.000    0.0
    Node_20  0.0   0.000    0.0
    Node_21  0.0   0.000    0.0
    Node_22  0.0   0.000    0.0
    Node_23  0.0   0.000    0.0
    Node_24  0.0   0.000    0.0

    """

    if bedload_gsd.shape[0] == self._grid.number_of_links:
        indexText = "Link_"
    else:
        indexText = "Node_"

    bedload_gsd = np.hstack((bedload_gsd, np.zeros([bedload_gsd.shape[0], 1])))
    bedload_gsd = np.cumsum(np.fliplr(bedload_gsd), axis=1) * 100
    bedload_gsd = np.around(bedload_gsd, 3)

    gs_D_ascending = np.sort(self._gsd[:, 0])
    columns = ["".join(item) for item in gs_D_ascending.astype(str)]

    df = pd.DataFrame(
        bedload_gsd,
        columns=columns,
        index=[indexText + str(i) for i in range(bedload_gsd.shape[0])],
    )
    return df


def get_available_fields():
    """Return a list of available fields and their units.

    Examples
    --------

    >>> from pprint import pprint
    >>> from . import _utilities as utilities
    >>> fields = utilities.get_available_fields()

    >>> pprint(fields)
    [('rbd._bed_subsurf__gsd_link', '[mm,%]'),
     ('rbd._bed_subsurf__gsd_node', '[mm,%]'),
     ('rbd._bed_surf__act_layer_thick_link', '[m]'),
     ('rbd._bed_surf__act_layer_thick_prev_time_link', '[m]'),
     ('rbd._bed_surf__elev_fix_node', '[m]'),
     ('rbd._bed_surf__geo_std_size_link', '[mm]'),
     ('rbd._bed_surf__geo_std_size_node', '[mm]'),
     ('rbd._bed_surf__geom_mean_size_link', '[mm]'),
     ('rbd._bed_surf__geom_mean_size_node', '[mm]'),
     ('rbd._bed_surf__gsd_fix_node', '[mm,%]'),
     ('rbd._bed_surf__gsd_link', '[mm,%]'),
     ('rbd._bed_surf__gsd_node', '[mm,%]'),
     ('rbd._bed_surf__gsd_orig_link', '[mm,%]'),
     ('rbd._bed_surf__gsd_orig_node', '[mm,%]'),
     ('rbd._bed_surf__median_size_link', '[mm]'),
     ('rbd._bed_surf__median_size_node', '[mm]'),
     ('rbd._bed_surf__sand_fract_link', '[-]'),
     ('rbd._bed_surf__sand_fract_node', '[-]'),
     ('rbd._bed_surf__thick_new_layer_link', '[m]'),
     ('rbd._sed_transp__bedload_gsd_fix_link', '[mm,%]'),
     ('rbd._sed_transp__bedload_gsd_link', '[mm,%]'),
     ('rbd._sed_transp__bedload_rate_fix_link', '[m^2/s]'),
     ('rbd._sed_transp__bedload_rate_link', '[m^2/s]'),
     ('rbd._sed_transp__net_bedload_node', '[m^2/s]'),
     ('rbd._surface_water__shear_stress_link', '[Pa]'),
     ('rbd._surface_water__velocity_prev_time_link', '[m/s]'),
     ('rbd._topogr__elev_orig_link', '[m]'),
     ('rbd._topogr__elev_orig_node', '[m]'),
     ('rbd._topogr__elev_subsurf_link', '[m]')]
    """

    # Define available fields and their units
    fields = sorted(
        [
            ("rbd._bed_subsurf__gsd_link", "[mm,%]"),
            ("rbd._bed_subsurf__gsd_node", "[mm,%]"),
            ("rbd._bed_surf__act_layer_thick_link", "[m]"),
            ("rbd._bed_surf__act_layer_thick_prev_time_link", "[m]"),
            ("rbd._bed_surf__elev_fix_node", "[m]"),
            ("rbd._bed_surf__geom_mean_size_link", "[mm]"),
            ("rbd._bed_surf__geom_mean_size_node", "[mm]"),
            ("rbd._bed_surf__geo_std_size_link", "[mm]"),
            ("rbd._bed_surf__geo_std_size_node", "[mm]"),
            ("rbd._bed_surf__gsd_fix_node", "[mm,%]"),
            ("rbd._bed_surf__gsd_link", "[mm,%]"),
            ("rbd._bed_surf__gsd_node", "[mm,%]"),
            ("rbd._bed_surf__gsd_orig_link", "[mm,%]"),
            ("rbd._bed_surf__gsd_orig_node", "[mm,%]"),
            ("rbd._bed_surf__median_size_link", "[mm]"),
            ("rbd._bed_surf__median_size_node", "[mm]"),
            ("rbd._bed_surf__sand_fract_link", "[-]"),
            ("rbd._bed_surf__sand_fract_node", "[-]"),
            ("rbd._bed_surf__thick_new_layer_link", "[m]"),
            ("rbd._sed_transp__bedload_gsd_fix_link", "[mm,%]"),
            ("rbd._sed_transp__bedload_gsd_link", "[mm,%]"),
            ("rbd._sed_transp__bedload_rate_link", "[m^2/s]"),
            ("rbd._sed_transp__net_bedload_node", "[m^2/s]"),
            ("rbd._sed_transp__bedload_rate_fix_link", "[m^2/s]"),
            ("rbd._surface_water__shear_stress_link", "[Pa]"),
            ("rbd._surface_water__velocity_prev_time_link", "[m/s]"),
            ("rbd._topogr__elev_orig_link", "[m]"),
            ("rbd._topogr__elev_orig_node", "[m]"),
            ("rbd._topogr__elev_subsurf_link", "[m]"),
        ]
    )

    return fields
