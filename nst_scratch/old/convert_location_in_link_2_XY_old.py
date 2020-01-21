# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 16:28:00 2019

This code converts location in a link to an X, Y


@author: Jon Czuba
"""
# %%
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

# from landlab.components import NetworkSedimentTransporter

# %%

# HELP HELP -- This all works if you specify a parcel_number and parcel_time
# we either need to make this into a function that returns parcel_x and parcel_y
# or loop through it in plot_network_parcels.py. I am not sure how to effectively
# do either.

# def plot_parcel_locations(nst_instance):
#    if not isinstance(nst_instance, NetworkSedimentTransporter):
#        raise ValueError("msg")

# identify the parcel to locate at a specified timestep
parcel_number = 6
parcel_time = 3

# find what link a given parcel is in at a specified time
parcel_link = parcels.element_id[parcel_number, parcel_time].values
# self._parcels.element_id[parcel_number,parcel_time].values
# parcels.element_id[:,parcel_time].values

# determine the location of that parcel in its link
parcel_loc = parcels.location_in_link[parcel_number, parcel_time].values

# get the X, Y vertices of the squiggly line for that link (loaded by import_shapefile.py)
link_x = x_of_polylines[int(parcel_link)]
link_y = y_of_polylines[int(parcel_link)]

# cumulative distance between squiggly-line link vertices [0 to link length]
link_dist = np.concatenate(
    [[0], np.cumsum(np.sqrt(np.diff(link_x) ** 2 + np.diff(link_y) ** 2))]
)

# cumulative relative distance between squiggly-line link verticies [0 to 1]
# divide cumulative distance by max distance to get vector of distances between
# 0 and 1
link_rel_dist = link_dist / np.max(link_dist)

# add condition if 0?

# determine which two points on squiggly line bound parcel location in that
# link
upper_loc_idx = np.array(np.where(np.greater(link_rel_dist - parcel_loc, 0) == True))

to_interp_link_loc = link_rel_dist[upper_loc_idx[0, 0] - 1 : upper_loc_idx[0, 0] + 1]
to_interp_link_x = link_x[upper_loc_idx[0, 0] - 1 : upper_loc_idx[0, 0] + 1]
to_interp_link_y = link_y[upper_loc_idx[0, 0] - 1 : upper_loc_idx[0, 0] + 1]

# check values are increasing
np.all(np.diff(to_interp_link_loc) > 0)

# interpolate the X,Y coordinates from the parcel location
parcel_x = np.interp(parcel_loc, to_interp_link_loc, to_interp_link_x)
parcel_y = np.interp(parcel_loc, to_interp_link_loc, to_interp_link_y)
