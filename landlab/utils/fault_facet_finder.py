#! /usr/bin/env python

"""
This class is designed to provide functions to allow the automated
identification of planar facet surfaces above fault traces.
Module is SLOW (e.g., minutes+ per full analysis of a "large" data set). It is
only intended for model post-analysis or DEM analysis. Do not loop this class!!
This is part of the NSF funded project investigating fault scarp degradation,
Tucker, Hobley, McCoy.
"""
from __future__ import print_function

import sys

import numpy as np
import six
from pylab import figure, plot, show

from landlab.plot import imshow as gridshow

if six.PY3:

    def cmp(a, b):
        return (a > b) - (a < b)


class find_facets(object):
    """
    Note that this class assumes the grid does not change during the model
    run. Changes to data stored in the grid should (?) update automatically.

    If *fault_azimuth* is supplied, it should be -pi/2 < az <= pi/2 (i.e.,
    we don't consider fault dip, even if it's known).
    """

    def __init__(self, grid, elev_field="topographic__elevation", fault_azimuth=None):
        """
        Note that this class assumes the grid does not change during the model
        run. Changes to data stored in the grid should (?) update automatically.

        If *fault_azimuth* is supplied, it should be -pi/2 < az <= pi/2 (i.e.,
        we don't consider fault dip, even if it's known).
        """
        if not np.isclose(grid.dx, grid.dy):
            raise ValueError("row and column spacing must be the same")

        self.grid = grid
        self.elevs = self.grid.at_node[elev_field]
        self.az = fault_azimuth

    def analyse_fault_trace(self, fault_trace_node_ids):
        """
        This method takes the grid and an array listing the (contiguous) node
        ids of cells that contain a single fault segment trace of interest.

        It sets and returns the azimuth of the fault trace, az,
        -pi/2 < az <= pi/2.
        (i.e., smallest angle between north and the trace).
        """
        self.ft_trace_node_ids = fault_trace_node_ids
        x = self.grid.node_x[fault_trace_node_ids]
        y = self.grid.node_y[fault_trace_node_ids]
        (grad, offset) = np.polyfit(x, y, 1)
        angle = np.arctan(grad)
        if grad >= 0.0:
            self.az = np.pi / 2.0 - angle
        else:
            self.az = -np.pi / 2.0 - angle

        return self.az

    def set_slopes_aspects(self):
        self.slopes = self.grid.calc_slopes_of_nodes(elevs=self.elevs)
        self.aspect = self.grid.calc_aspect_of_node(elevs=self.elevs)
        print("Calculated and stored slopes and aspects...")

    def define_aspect_node_subset(self, angle_tolerance=5.0):
        """
        This method sets and returns a list of all nodes in the landscape which
        have an
        aspect within 5 degrees of perpendicular to the fault trace.

        It assumes self.az, the angle between north and the fault trace, has
        already been set, and also that self.slopes and self.aspect are also
        set.
        The returned boolean array is num_core_nodes long.
        *angle_tolerance* is the angle in degrees that the aspect must be within
        from the fault trace angle.
        NB: this version is too discriminating on the aspect restriction,
        presumably because we use only a single ft angle for what's really
        a 2d trace. Need to work with local aspect.
        """
        perp_to_az = np.pi / 2.0 + self.az
        five_degrees = np.pi / 180.0 * angle_tolerance  # note the silly naming
        perp_plus_five = (perp_to_az + five_degrees) % (2.0 * np.pi)
        perp_minus_five = (perp_to_az - five_degrees) % (2.0 * np.pi)
        other_dip_plus_five = (perp_to_az + np.pi + five_degrees) % (2.0 * np.pi)
        other_dip_minus_five = (perp_to_az + np.pi - five_degrees) % (2.0 * np.pi)

        # need to be careful near the 2pi->0 discontinuity...
        greater_condition = np.greater(self.aspect, perp_minus_five)
        lesser_condition = np.less(self.aspect, perp_plus_five)
        if (perp_to_az - five_degrees) < 0.0:
            condition_first_dip = np.logical_or(greater_condition, lesser_condition)
        else:
            condition_first_dip = np.logical_and(greater_condition, lesser_condition)
        greater_condition_2 = np.greater(self.aspect, other_dip_minus_five)
        lesser_condition_2 = np.less(self.aspect, other_dip_plus_five)
        if (perp_to_az + np.pi + five_degrees) // (
            2.0 * np.pi
        ):  # top condition exceeds 2pi
            condition_opposite_dip = np.logical_or(
                greater_condition_2, lesser_condition_2
            )
        else:
            condition_opposite_dip = np.logical_and(
                greater_condition_2, lesser_condition_2
            )

        self.aspect_close_nodes = np.logical_or(
            condition_first_dip, condition_opposite_dip
        )
        print("Calculated and stored nodes with aspects compatible with fault trace...")
        return self.aspect_close_nodes

    def define_aspect_node_subset_local(
        self, dist_tolerance=4.0, angle_tolerance=15.0, dip_dir="E"
    ):
        """
        """
        grid = self.grid
        try:
            print("using subset")
            # remember, steep_nodes is already core_nodes.size long
            subset = np.where(self.steep_nodes)[0]
        except NameError:
            print("using all nodes")
            subset = np.arange(grid.core_nodes.size)
        closest_ft_node = np.empty(subset.size, dtype=int)
        angle_to_ft = np.empty_like(closest_ft_node, dtype=float)
        new_angle_to_ft = np.empty_like(closest_ft_node, dtype=float)
        distance_to_ft = np.empty_like(closest_ft_node, dtype=float)
        distance_to_ft.fill(sys.float_info.max)
        new_distance_to_ft = np.empty_like(closest_ft_node, dtype=float)
        for i in self.ft_trace_node_ids:
            grid.calc_distances_of_nodes_to_point(
                (grid.node_x[i], grid.node_y[i]),
                node_subset=grid.core_nodes[subset],
                get_az="angles",
                out_distance=new_distance_to_ft,
                out_azimuth=new_angle_to_ft,
            )
            closer_nodes = new_distance_to_ft < distance_to_ft
            distance_to_ft[closer_nodes] = new_distance_to_ft[closer_nodes]
            angle_to_ft[closer_nodes] = new_angle_to_ft[closer_nodes]
            closest_ft_node[closer_nodes] = i
        self.closest_ft_node = -np.ones(grid.core_nodes.size)
        self.distance_to_ft = -np.ones(grid.core_nodes.size)
        self.angle_to_ft = -np.ones(grid.core_nodes.size)
        self.closest_ft_node[subset] = closest_ft_node
        self.distance_to_ft[subset] = distance_to_ft
        # angle_to_ft is actually the angle_from_ft! So we need to adjust.
        # a second problem is that pts downslope (opposite az) can also be on the line.
        # solution - take a dip_dir input...
        angle_to_ft = (angle_to_ft + np.pi) % (2.0 * np.pi)
        self.angle_to_ft[subset] = angle_to_ft
        # gridshow.imshow_grid_at_node(self.grid, self.distance_to_ft)
        # show()
        # gridshow.imshow_grid_at_node(self.grid, self.angle_to_ft)
        # show()
        # the relevant condition is now that the local aspect and angle to fault
        # are the same...
        # We need to bias the five degrees against distant points, as it's easier
        # to have similar angles in the far field. Rule should be in px - the
        # two angles should be within *angle_tol* px of each other at the ft
        # trace.
        divergence_at_ft = distance_to_ft * np.tan(
            (angle_to_ft - self.aspect[subset]) % np.pi
        )
        # might be *too* forgiving for close-in nodes
        condition = np.less(np.fabs(divergence_at_ft), grid.dx * dist_tolerance)
        # ...so add another tester; must be w/i 15 degrees of each other:
        diff_angles = np.min(
            [
                np.fabs(angle_to_ft - self.aspect[subset]),
                np.fabs(np.fabs(angle_to_ft - self.aspect[subset]) - 2.0 * np.pi),
            ],
            axis=0,
        )
        self.diff_angles = np.empty(grid.core_nodes.size, dtype=float)
        self.diff_angles.fill(sys.float_info.max)
        self.diff_angles[subset] = diff_angles
        # gridshow.imshow_grid_at_node(self.grid, self.angle_to_ft)
        # show()
        figure(6)
        gridshow.imshow_grid_at_node(
            self.grid, np.where(self.diff_angles < 100000.0, self.diff_angles, -1.0)
        )
        condition2 = np.less(diff_angles, angle_tolerance * np.pi / 180.0)
        condition = np.logical_and(condition, condition2)
        core_nodes_size_condition = np.zeros(grid.core_nodes.size, dtype=bool)
        core_nodes_size_condition[subset] = condition
        # gridshow.imshow_grid_at_node(self.grid, core_nodes_size_condition)
        # show()
        # core_nodes_size_condition = np.zeros(grid.core_nodes.size, dtype=bool)
        # core_nodes_size_condition[subset] = condition2
        # gridshow.imshow_grid_at_node(self.grid, core_nodes_size_condition)
        # show()
        self.aspect_close_nodes = core_nodes_size_condition
        print("Calculated and stored nodes with aspects compatible with fault trace...")
        return self.aspect_close_nodes

    def define_steep_nodes(self, threshold_in_degrees=5.0):
        """
        This method sets and returns a list of all nodes in the landscape which
        are "steep" and could be part of a facet.
        The critical hillslope angle is set by *threshold_in_degrees*, and
        defaults to 5.

        This assumes you have already called define_aspect_node_subset, in
        which self.slope is set.
        The returned boolean array is num_core_nodes long.
        """
        threshold_in_rads = threshold_in_degrees * np.pi / 180.0
        self.steep_nodes = np.greater(self.slopes, threshold_in_rads)
        print("Calculated and stored nodes with slopes exceeding slope threshold...")
        # gridshow.imshow_grid_at_node(self.grid, self.steep_nodes)
        # show()
        return self.steep_nodes

    def show_possible_nodes(self):
        """
        Once the subsets by aspect and slope have been set, call this function
        to see both the whole elevation map, and the subset of nodes that
        will be searched.
        """
        possible_core_nodes = np.logical_and(self.steep_nodes, self.aspect_close_nodes)
        figure(1)
        gridshow.imshow_grid_at_node(self.grid, self.elevs)
        figure(2)
        gridshow.imshow_grid_at_node(self.grid, self.slopes)
        figure(3)
        gridshow.imshow_grid_at_node(self.grid, self.aspect)
        figure(4)
        gridshow.imshow_grid_at_node(self.grid, possible_core_nodes)
        show()

    def find_coherent_facet_patches(self, tolerance=3.0, threshold_num_px=12):
        """
        This method searches the (already determined) possible pixels for
        patches with coherent slope angles, within a *tolerance* (in degrees).
        A patch is only recorded if it consists of at least *threshold_num_px*.

        The method records and returns:

        *  a ragged array of lists, where each list is the pixels comprising
           each facet patch, and
        *  a (num_patches, 2) array recording the mean slope and and its stdev
           for each patch.
        """
        self.possible_core_nodes = np.where(
            np.logical_and(self.steep_nodes, self.aspect_close_nodes)
        )[0]
        consistent_slope_patches = []
        for i in self.possible_core_nodes:
            nodes_in_patch = np.array([i])
            patch_size = 1
            mean_slope = self.slopes[nodes_in_patch]
            while 1:
                possible_neighbors = np.union1d(
                    self.grid.active_adjacent_nodes_at_node[nodes_in_patch].flat,
                    self.possible_core_nodes,
                )
                neighbor_slopes = self.slopes[possible_neighbors]
                low_tol_condition = np.greater(neighbor_slopes, mean_slope - tolerance)
                high_tol_condition = np.less(neighbor_slopes, mean_slope + tolerance)
                total_condition = np.logical_and(low_tol_condition, high_tol_condition)
                nodes_in_patch = possible_neighbors[total_condition]
                new_patch_size = nodes_in_patch.size
                if patch_size == new_patch_size:
                    break
                else:
                    patch_size = new_patch_size
                    mean_slope = np.mean(neighbor_slopes[total_condition])
            if new_patch_size > threshold_num_px:
                consistent_slope_patches.append(nodes_in_patch)
        # May also need a uniqueness test: a px should only appear in one list.
        # (i.e., all patches containing a given px should all be identical)
        self.consistent_slope_patches = consistent_slope_patches
        return consistent_slope_patches

    def find_slope_lines(self, tolerance=1.0):
        """
        This method attempts to find slope-consistent line profiles up facets,
        perpendicular to the fault.
        Assumes you used define_aspect_node_subset_local().
        """
        grid = self.grid
        self.possible_core_nodes = np.where(
            np.logical_and(self.steep_nodes, self.aspect_close_nodes)
        )[0]
        pcn = self.possible_core_nodes
        unique_starting_pts = np.unique(self.closest_ft_node[pcn])  # in real node nos
        # set up places to store the profile data:
        profile_ft_node_id = []
        profile_ft_node_x = []
        profile_ft_node_y = []
        profile_ft_node_z = []
        profile_ft_node_dist = []
        profile_x_facet_pts = []
        profile_z_facet_pts = []
        profile_S_facet_pts = []
        count = 0
        for i in unique_starting_pts:
            count += 1
            print("Running ", count, " of ", unique_starting_pts.size)
            # set the local angle of the ft trace:
            ft_pt_distances_to_node = self.grid.calc_distances_of_nodes_to_point(
                (grid.node_x[i], grid.node_y[i]), node_subset=self.ft_trace_node_ids
            )
            close_ft_nodes = np.less(ft_pt_distances_to_node, 5.0 * grid.dx)
            x = grid.node_x[self.ft_trace_node_ids[close_ft_nodes]]
            y = grid.node_y[self.ft_trace_node_ids[close_ft_nodes]]
            (grad, offset) = np.polyfit(x, y, 1)
            condition = np.equal(self.closest_ft_node[pcn], i)
            nodes_possible = pcn[condition]
            print(nodes_possible.size, " nodes")
            if nodes_possible.size > 10.0:
                # their_az = self.angle_to_ft[nodes_possible]
                # their_diff_angles = self.diff_angles[nodes_possible]
                their_elevs = self.elevs[grid.core_nodes][nodes_possible]
                # their_distances = self.distance_to_ft[nodes_possible]
                # need to make new distances so we remove the ambiguity of angle around the ft point (i.e., dists from a far-field pt on the ft normal)
                # now make a multiplier to make sure the reference point for
                # new distances is far from the actual pts:
                multiplier = 10.0 * np.ptp(grid.node_y[grid.core_nodes[nodes_possible]])
                # derive the position:
                x_ref = grid.node_x[i] + cmp(
                    grid.node_x[i],
                    np.mean(grid.node_x[grid.core_nodes[nodes_possible]]),
                ) * multiplier * abs(grad)
                y_ref = (
                    grid.node_y[i]
                    + cmp(
                        grid.node_y[i],
                        np.mean(grid.node_y[grid.core_nodes[nodes_possible]]),
                    )
                    * multiplier
                )
                # get new absolute distances
                dist_to_ft = self.grid.calc_distances_of_nodes_to_point(
                    (x_ref, y_ref), node_subset=np.array([i])
                )
                dists_along_profile = (
                    self.grid.calc_distances_of_nodes_to_point(
                        (x_ref, y_ref), node_subset=grid.core_nodes[nodes_possible]
                    )
                    - dist_to_ft
                )
                # note the ft is now the origin, but pts might be back-to-front (consistently, though)
                # sort the distances. Remove any pts that aren't in a "cluster".
                # We assume there will be one big "bunched" plane, then a load
                # of outliers
                dist_order = np.argsort(dists_along_profile)
                dist_diffs = np.diff(dists_along_profile[dist_order])
                print("dists along profile sorted: ", dists_along_profile[dist_order])
                print("dist diffs: ", dist_diffs)
                # max_diff = 3.*np.median(dist_diffs) #######this might need
                # attention if there's a heavy tail on the distances
                if grad < 1:
                    mod = np.sqrt(1.0 + grad ** 2.0)
                else:
                    mod = np.sqrt(1.0 + (1.0 / grad) ** 2.0)
                max_diff = 1.9 * mod * grid.dx
                locs_of_large_diffs = np.where(dist_diffs > max_diff)[0]
                # there should only be 1 place on the line where there's a cluster, i.e., a large pts_betw_of_max_diffs.
                # This is what we're seeking.
                # ...this can be empty quite easily
                pts_betw_large_diffs = np.diff(locs_of_large_diffs)
                # need to be careful here in case the where call gives an empty
                # array
                if locs_of_large_diffs.size > 1:
                    biggest_interval_loc = np.argmax(pts_betw_large_diffs)
                elif locs_of_large_diffs.size == 1:
                    # one side or the other must be bigger:
                    if 2.0 * locs_of_large_diffs[0] < dists_along_profile.size - 1:
                        locs_of_large_diffs = np.array(
                            [locs_of_large_diffs[0], (dists_along_profile.size - 1)]
                        )
                    else:
                        locs_of_large_diffs = np.array([0, locs_of_large_diffs[0]])
                    biggest_interval_loc = np.array([0])
                    # here we assume that the single large diff must be further
                    # from the ft than the plane
                else:
                    locs_of_large_diffs = np.array([0, (dists_along_profile.size - 1)])
                    biggest_interval_loc = np.array([0])
                    # ...all the pts in the line are one cluster
                # apply a test to ensure we only save "big" patches; a
                # threshold of 10 pts on the line
                try:
                    patch_size = pts_betw_large_diffs[biggest_interval_loc]
                except IndexError:  # pts_betw_large_diffs is empty
                    patch_size = locs_of_large_diffs[1] - locs_of_large_diffs[0]
                if patch_size > 10.0:
                    start_pt_of_cluster = locs_of_large_diffs[biggest_interval_loc] + 1
                    end_pt_of_cluster = (
                        locs_of_large_diffs[biggest_interval_loc + 1] + 1
                    )  # both referring to the sorted list
                    # both +1s are to account for required frame of ref changes - indices refer to where the big gaps start, not where they ends
                    # so:
                    dists_to_sorted_pts = dists_along_profile[dist_order][
                        start_pt_of_cluster:end_pt_of_cluster
                    ]
                    elevs_of_sorted_pts = their_elevs[dist_order][
                        start_pt_of_cluster:end_pt_of_cluster
                    ]
                    slopes_of_sorted_pts = self.slopes[nodes_possible][dist_order][
                        start_pt_of_cluster:end_pt_of_cluster
                    ]
                    profile_ft_node_id.append(i.copy())
                    profile_ft_node_x.append(grid.node_x[i].copy())
                    profile_ft_node_y.append(grid.node_y[i].copy())
                    profile_ft_node_z.append(self.elevs[i].copy())
                    profile_ft_node_dist.append(dist_to_ft.copy())
                    profile_x_facet_pts.append(dists_to_sorted_pts.copy())
                    profile_z_facet_pts.append(elevs_of_sorted_pts.copy())
                    profile_S_facet_pts.append(slopes_of_sorted_pts.copy())
                    figure(5)
                    plot(dists_to_sorted_pts, elevs_of_sorted_pts)
                    # dirty, but effective code!

        self.profile_ft_node_id = profile_ft_node_id
        self.profile_ft_node_x = profile_ft_node_x
        self.profile_ft_node_y = profile_ft_node_y
        self.profile_ft_node_z = profile_ft_node_z
        self.profile_ft_node_dist = profile_ft_node_dist
        self.profile_x_facet_pts = profile_x_facet_pts
        self.profile_z_facet_pts = profile_z_facet_pts
        self.profile_S_facet_pts = profile_S_facet_pts

    def fit_slopes_to_facet_lines(
        self, polynomial_degree=4, curvature_threshold=0.0004
    ):
        """
        Fits (linear) lines of best fit to extracted profiles, already stored as
        class properties.
        """
        avg_slopes_linear = []
        avg_slopes_poly = []
        curv_of_flattest_part_list = []
        slope_min_curv = []
        rsqd_list = []
        big_slope_small_curv = []
        elev_at_bssc = []
        for i in six.range(len(self.profile_x_facet_pts)):
            x = self.profile_x_facet_pts[i]
            z = self.profile_z_facet_pts[i]
            (grad, offset) = np.polyfit(x, z, 1)
            coeffs, residuals = np.polyfit(x, z, polynomial_degree, full=True)[:2]
            rsqd = 1.0 - residuals / (z.size * z.var())
            # differentiate the coeffs to get slope:
            diff_multiplier = np.arange(polynomial_degree + 1)[::-1]
            curv_multiplier = np.arange(polynomial_degree)[::-1]
            z_equ = np.poly1d(coeffs)
            S_equ = np.poly1d((coeffs * diff_multiplier)[:-1])
            curv_equ = np.poly1d(
                ((coeffs * diff_multiplier)[:-1] * curv_multiplier)[:-1]
            )
            S_at_each_pt = S_equ(x)
            curv_at_each_pt = curv_equ(x)
            avg_slopes_linear.append(abs(grad))
            avg_slopes_poly.append(np.amax(np.fabs(S_at_each_pt)))
            loc_of_flattest_part = np.argmin(np.fabs(curv_at_each_pt[2:-2])) + 2
            curv_of_flattest_part = curv_at_each_pt[loc_of_flattest_part]
            S_at_min_curve_untested = abs(S_at_each_pt[loc_of_flattest_part])
            small_curves = np.less(np.fabs(curv_at_each_pt[2:-2]), curvature_threshold)
            try:
                big_slope_small_curv.append(np.amax(S_at_each_pt[small_curves]))
                elev_at_bssc.append(z[np.argmax(S_at_each_pt[small_curves])])
            except ValueError:
                big_slope_small_curv.append(np.nan)
                elev_at_bssc.append(np.nan)
            slope_min_curv.append(S_at_min_curve_untested)
            curv_of_flattest_part_list.append(curv_of_flattest_part)
            rsqd_list.append(rsqd)
            # figure(8)
            # synthetic_z = grad*x + offset
            synthetic_z = z_equ(x)
            plot(x, z, "x")
            plot(x, synthetic_z, "-")
        self.avg_slopes_linear = np.array(avg_slopes_linear)
        self.avg_slopes_poly = np.array(avg_slopes_poly)
        self.curv_of_flattest_part = np.array(curv_of_flattest_part_list)
        self.slope_min_curv = np.array(slope_min_curv)
        self.big_slope_small_curv = np.array(big_slope_small_curv)
        self.elev_at_bssc = np.array(elev_at_bssc)
        self.rsqd = np.array(rsqd_list)
