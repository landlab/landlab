"""
Implements a series of functions to create and/or initialize the required fields
where calculations are conducted.

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager

Examples
--------
In this case we provide examples of callings to all different functions.
These function is automatically called during the execution of the component
but for verification purposes we test them here.
This is the same base example described extensively in river bed dynamics, so
we removed comments that are already available in the main component

>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from . import _initialize_fields as initialize

>>> grid = RasterModelGrid((5, 5))

Case1a: Using default values

>>> bed_surf__gsd_loc_node = None  # Default value in river bed dynamics
>>> bed_surf__gsd_loc_node = initialize.field_at_node(grid, bed_surf__gsd_loc_node)
>>> bed_surf__gsd_loc_node
array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0])

Case1b: Setting bed_surf__gsd_loc_node values

>>> bed_surf__gsd_loc_node = [
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
... ]
>>> bed_surf__gsd_loc_node = initialize.field_at_node(grid, bed_surf__gsd_loc_node)
>>> bed_surf__gsd_loc_node
array([0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1,
       1, 0])

Case1c: Setting bed_surf__gsd_loc_node values with wrong size

>>> bed_surf__gsd_loc_node = [1, 2, 3]
>>> bed_surf__gsd_loc_node = initialize.field_at_node(grid, bed_surf__gsd_loc_node)
>>> bed_surf__gsd_loc_node
array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0])

Case2a: Using default values

>>> sed_transp__bedload_rate_fix_link = None  # Default value in river bed dynamics
>>> sed_transp__bedload_rate_fix_link = initialize.field_at_link(
...     grid, sed_transp__bedload_rate_fix_link
... )
>>> sed_transp__bedload_rate_fix_link
array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

Case2b: Setting sed_transp__bedload_gsd_fix_link as 6.7 at link 3, 6, and 10

>>> sed_transp__bedload_rate_fix_link = np.zeros([grid.number_of_links, 1])
>>> sed_transp__bedload_rate_fix_link[[3, 6, 10]] = 6.7
>>> sed_transp__bedload_rate_fix_link = initialize.field_at_link(
...     grid, sed_transp__bedload_rate_fix_link
... )
>>> sed_transp__bedload_rate_fix_link
array([ 0. ,  0. ,  0. ,  6.7,  0. ,  0. ,  6.7,  0. ,  0. ,  0. ,  6.7,
        0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,
        0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ])

Case2c: Setting sed_transp__bedload_gsd_fix_link with wrong size

>>> sed_transp__bedload_rate_fix_link = [1, 2, 3]
>>> sed_transp__bedload_rate_fix_link = initialize.field_at_link(
...     grid, sed_transp__bedload_rate_fix_link
... )
>>> sed_transp__bedload_rate_fix_link
array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

Case3a: Using default values

>>> sed_transp__bedload_gsd_fix_link = None
>>> gsd = np.array([[32, 100], [16, 25], [8, 0]])
>>> sed_transp__bedload_gsd_fix_link = initialize.gsd_at_link(
...     grid, sed_transp__bedload_gsd_fix_link, gsd
... )

For displaying purposes we only show the first 5 links, which are all zeros as
we did not specified any imposed gsd

>>> sed_transp__bedload_gsd_fix_link[:5, :]
array([[ 0.,  0.],
       [ 0.,  0.],
       [ 0.,  0.],
       [ 0.,  0.],
       [ 0.,  0.]])

Case3b: Setting sed_transp__bedload_gsd_fix_link at links 2 and 4 as
gsd = np.array([[32, 100], [16, 50], [8, 0]])

If we want to specify an imposed gsd to a link we can do this:
>>> sed_transp__bedload_gsd_fix_link = np.zeros(
...     [grid.number_of_links, gsd.shape[0] - 1]
... )
>>> sed_transp__bedload_gsd_fix_link[1] = [0.50, 0.50]
>>> sed_transp__bedload_gsd_fix_link[3] = [0.50, 0.50]
>>> sed_transp__bedload_gsd_fix_link = initialize.gsd_at_link(
...     grid, sed_transp__bedload_gsd_fix_link, gsd
... )

For displaying purposes we only show the first 5 links.

>>> sed_transp__bedload_gsd_fix_link[:5, :]
array([[ 0. ,  0. ],
       [ 0.5,  0.5],
       [ 0. ,  0. ],
       [ 0.5,  0.5],
       [ 0. ,  0. ]])

The rest of elements are [0., 0.] because they are not imposed. The importance
of this step is that in the component. It will know what links must preserve its gsd.

Case3c: Setting sed_transp__bedload_gsd_fix_link with wrong size

>>> sed_transp__bedload_gsd_fix_link = [0.50, 0.50]
>>> sed_transp__bedload_gsd_fix_link = initialize.gsd_at_link(
...     grid, sed_transp__bedload_gsd_fix_link, gsd
... )

For displaying purposes we only show the first 5 links, which are all zeros because
there was an error while trying to specify sed_transp__bedload_gsd_fix_link

>>> sed_transp__bedload_gsd_fix_link[:5, :]
array([[ 0.,  0.],
       [ 0.,  0.],
       [ 0.,  0.],
       [ 0.,  0.],
       [ 0.,  0.]])

Case4a: Using default values

>>> grid.at_link["surface_water__velocity"] = np.full(grid.number_of_links, 0.25)
>>> surface_water__velocity_prev_time_link = None

>>> surface_water__velocity_prev_time_link = initialize.velocity_at_link(
...     grid, surface_water__velocity_prev_time_link
... )
>>> surface_water__velocity_prev_time_link
array([ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25])

Case4b: Setting surface_water__velocity_prev_time_link with link 4 as 1.2 m/s
>>> surface_water__velocity_prev_time_link = np.full(grid.number_of_links, 0.50)
>>> surface_water__velocity_prev_time_link[4] = 1.2

>>> surface_water__velocity_prev_time_link = initialize.velocity_at_link(
...     grid, surface_water__velocity_prev_time_link
... )
>>> surface_water__velocity_prev_time_link
array([ 0.5,  0.5,  0.5,  0.5,  1.2,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
        0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
        0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
        0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5])

Case4c: Setting surface_water__velocity_prev_time_link with wrong size

>>> surface_water__velocity_prev_time_link = 1
>>> surface_water__velocity_prev_time_link = initialize.velocity_at_link(
...     grid, surface_water__velocity_prev_time_link
... )
>>> surface_water__velocity_prev_time_link
array([ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25])

"""

import numpy as np


def field_at_node(grid, field):
    n_nodes = grid.number_of_nodes

    if field is None:
        field = np.zeros(n_nodes, dtype=int)
    else:
        field = np.array(field, dtype=int)
        # Check that the input size is correct otherwise assigns zero
        if field.size > 0:
            if field.size == n_nodes:
                field = np.array(field, dtype=int).flatten()
            else:
                field = np.zeros(n_nodes, dtype=int)

    return field


def field_at_link(grid, field):
    n_links = grid.number_of_links

    if field is None:
        field = np.zeros(n_links, dtype=float)
    else:
        field = np.array(field, dtype=float)
        # Check that the input size is correct otherwise assigns zero
        if field.size > 0:
            if field.size == n_links:
                field = np.array(field, dtype=float).flatten()
            else:
                field = np.zeros(n_links, dtype=float)

    return field


def gsd_at_link(grid, field, gsd):
    n_links = grid.number_of_links

    if field is None:
        field = np.zeros((n_links, gsd.shape[0] - 1), dtype=float)
    else:
        field = np.array(field, dtype=float)
        # Check that the input size is correct otherwise assigns zero
        if field.shape[0] > 0:
            if not (
                (field.shape[0] == n_links) and (field.shape[1] == gsd.shape[0] - 1)
            ):
                field = np.zeros((n_links, gsd.shape[0] - 1), dtype=float)

    return field


def velocity_at_link(grid, field):
    n_links = grid.number_of_links

    if field is None:
        field = grid["link"]["surface_water__velocity"].copy()
    else:
        field = np.array(field, dtype=float)
        # Check that the input size is correct otherwise assigns zero
        if field.size > 0:
            if field.size == n_links:
                field = np.array(field, dtype=float).flatten()
            else:
                field = grid["link"]["surface_water__velocity"].copy()

    return field
