"""
Created on Fri Oct 24 16:28:00 2019

This code converts location in a link to an X, Y


@author: Jon Czuba, Katy Barnhart
"""

import numpy as np


def locate_parcel_xy(grid, parcels, parcel_time_index, parcel_number):
    # determine the location of that parcel in its link
    parcel_loc = parcels.dataset.location_in_link[
        parcel_number, parcel_time_index
    ].values

    # parcels that end their timestep off network have the starting link id
    # recorded, and np.nan as the distance.

    if not np.isnan(parcel_loc):
        # get link id
        parcel_link = int(
            parcels.dataset.element_id[parcel_number, parcel_time_index].values
        )

        # DANGER DANGER: This code assumes the verticies of links are ordered
        # from upstream to downstream. This should be the case for delineated
        # river networks, so this line should not be necessary. A quick
        # work-around is the following, but ideally verticies should be flipped in GIS.
        # parcel_loc = 0.9999 - parcel_loc

        # get the X, Y vertices of the squiggly line for that link (loaded by
        # import_shapefile.py)

        # I am not sure these next 2 lines are necessary here, but I don't know
        # how to do this differently and/or better.
        if "x_of_polyline" in grid.at_link:
            link_x = grid["link"]["x_of_polyline"][parcel_link]
            link_y = grid["link"]["y_of_polyline"][parcel_link]
        else:
            flow_dir = grid.at_link["flow__link_direction"][parcel_link]
            head_node = grid.node_at_link_head[parcel_link]
            tail_node = grid.node_at_link_tail[parcel_link]
            if flow_dir == -1:
                link_x = [grid.x_of_node[head_node], grid.x_of_node[tail_node]]
                link_y = [grid.y_of_node[head_node], grid.y_of_node[tail_node]]
            elif flow_dir == 1:
                # 1 = with direction of from tail to head.
                link_x = [grid.x_of_node[tail_node], grid.x_of_node[head_node]]
                link_y = [grid.y_of_node[tail_node], grid.y_of_node[head_node]]
            else:
                raise ValueError(
                    "trying to plot on an inactive link. this should not happen."
                )
            # eventually need to use x_of_node, y_of_node, and nodes_at_link,
            # but the upstream to downstream ordering also matters.

        # cumulative distance between squiggly-line link vertices [0 to link length]
        link_dist = np.concatenate(
            [[0], np.cumsum(np.sqrt(np.diff(link_x) ** 2 + np.diff(link_y) ** 2))]
        )

        # cumulative relative distance between squiggly-line link verticies [0 to 1]
        # divide cumulative distance by max distance to get vector of distances between
        # 0 and 1
        link_rel_dist = link_dist / np.max(link_dist)

        # # determine which two points on squiggly line bound parcel location in that link
        # upper_loc_idx = np.argmax(link_rel_dist - parcel_loc > 0)
        #
        # to_interp_link_loc = link_rel_dist[upper_loc_idx - 1 : upper_loc_idx + 1]
        # to_interp_link_x = link_x[upper_loc_idx - 1 : upper_loc_idx + 1]
        # to_interp_link_y = link_y[upper_loc_idx - 1 : upper_loc_idx + 1]
        #
        # # check values are increasing
        # np.all(np.diff(to_interp_link_loc) > 0)

        # interpolate the X,Y coordinates from the parcel location
        parcel_x = np.interp(
            parcel_loc, link_rel_dist, link_x
        )  # , left=np.nan, right=np.nan)
        parcel_y = np.interp(parcel_loc, link_rel_dist, link_y)
        # assert np.isnan(parcel_x) == False

        # save data to a single variable. better would be to save this info as
        # an element of parcels.dataset.X and ...Y
        XY = [parcel_x, parcel_y]

        # parcels.dataset["X"] = parcel_x
        # parcels.dataset["Y"] = parcel_y

    else:
        # if that parcel is no longer in the system do not try to compute X,Y and
        # instead return NaN
        XY = [np.nan, np.nan]

    # return the X,Y values
    return XY
