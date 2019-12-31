# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 16:28:00 2019

This code converts location in a link to an X, Y


@author: Jon Czuba
"""
import numpy as np


def locate_parcel_xy(grid, parcels, parcel_time, parcel_number, *args, **kwargs):

    _OUT_OF_NETWORK = grid.BAD_INDEX - 1

    if parcels.dataset.element_id[parcel_number, parcel_time] != _OUT_OF_NETWORK:

        # find what link a given parcel is in at a specified time
        parcel_link = parcels.dataset.element_id[parcel_number, parcel_time].values
        # self._parcels.element_id[parcel_number,parcel_time].values
        # parcels.element_id[:,parcel_time].values

        # determine the location of that parcel in its link
        parcel_loc = parcels.dataset.location_in_link[parcel_number, parcel_time].values

        # DANGER DANGER: This code assumes the verticies of links are ordered from upstream to downstream. Jon believes this will be the case for delineated river networks, so this line should not be necessary. For the test network, the links are ordered from downstream to upstream. The quickest and easiest fix is the following, but ideally we will not want this in the final code.
        parcel_loc = 0.9999 - parcel_loc

        # get the X, Y vertices of the squiggly line for that link (loaded by import_shapefile.py)

        # I am not sure these next 2 lines are necessary here, but I don't know how to do this differently and/or better.
        x_of_polylines = grid["link"]["x_of_polyline"]
        y_of_polylines = grid["link"]["y_of_polyline"]

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

        # determine which two points on squiggly line bound parcel location in that link
        upper_loc_idx = np.array(
            np.where(np.greater(link_rel_dist - parcel_loc, 0) == True)
        )

        to_interp_link_loc = link_rel_dist[
            upper_loc_idx[0, 0] - 1 : upper_loc_idx[0, 0] + 1
        ]
        to_interp_link_x = link_x[upper_loc_idx[0, 0] - 1 : upper_loc_idx[0, 0] + 1]
        to_interp_link_y = link_y[upper_loc_idx[0, 0] - 1 : upper_loc_idx[0, 0] + 1]

        # check values are increasing
        np.all(np.diff(to_interp_link_loc) > 0)

        # interpolate the X,Y coordinates from the parcel location
        parcel_x = np.interp(parcel_loc, to_interp_link_loc, to_interp_link_x)
        parcel_y = np.interp(parcel_loc, to_interp_link_loc, to_interp_link_y)

        # save data to a single variable. better would be to save this info as an element of parcels.dataset.X and ...Y
        XY = [parcel_x, parcel_y]

        # parcels.dataset["X"] = parcel_x
        # parcels.dataset["Y"] = parcel_y

    else:
        # if that parcel is no longer in the system do not try to compute X,Y and instead return NaN
        XY = [np.nan, np.nan]

    # return the X,Y values
    return XY
