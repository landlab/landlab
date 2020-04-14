#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Zone functions and class of SpeciesEvolver."""
from collections import OrderedDict
from enum import IntEnum, unique

import numpy as np
from pandas import Series


@unique
class Connection(IntEnum):
    """Zone connection type.

    Connection type represents the connectivity of zones in two time steps. It
    is described with the pattern, x-to-y where x and y are the descriptive
    counts (none, one, or many) of zone(s) at the earlier and later time step,
    respectively. For example, a connection type of one-to-many is a zone in
    the earlier time step that spatially overlaps multiple zones in the later
    time step.
    """

    NONE_TO_NONE = 0
    NONE_TO_ONE = 1
    ONE_TO_NONE = 2
    ONE_TO_ONE = 3
    ONE_TO_MANY = 4
    MANY_TO_ONE = 5
    MANY_TO_MANY = 6


def _update_zones(grid, prior_zones, new_zones, record):
    """Resolve the spatial connectivity of zones across two time steps.

    This method iterates over each zone of the prior time step to identify the
    zones of the current time step it spatially intersects. This method updates
    the zone attribute, ``successors``. Successor zones are the new zones
    existing at the current time step that are the continuation of prior time
    step zones referred to as the predecessor zones.

    The type of connection between predecessor and successor zones is described
    ``Connection``. The `_conn_type` property of zones are set by this method.

    In the `many` connections, a rule determines which of the prior zones
    persist as the zone in the current time step. The zone with the greatest
    area of intersection between the prior and current time steps persists to
    the current time step along with the others in `new_zones`.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab ModelGrid.
    prior_zones : Zone list
        The zones of the prior timestep.
    new_zones : Zone list
        The zones of the current timestep.
    record : Record
        The ZoneController record.

    Returns
    -------
    list of Zones
        The successor zones. The zones determined to exist at the current time
        step with an updated ``successors`` attribute.
    """
    # Stats are calculated for `record`.

    fragment_ct = 0
    capture_ct = 0
    area_captured = [0]

    # Get `successors`.

    if len(prior_zones) == 0 and len(new_zones) == 0:
        successors = []

    elif len(prior_zones) == 0 and len(new_zones) > 0:
        successors = new_zones
        conn_type = _determine_connection_type(0, 1)
        for n in new_zones:
            n._conn_type = conn_type
            n._successors = [n]

    elif len(prior_zones) > 0 and len(new_zones) == 0:
        successors = []
        conn_type = _determine_connection_type(1, 0)
        for p in prior_zones:
            p._conn_type = conn_type
            p._successors = []

    elif len(prior_zones) > 0 and len(new_zones) > 0:
        successors = []

        # Get index maps for the prior zones (`ps`) and new zones (`ns`).
        ps_index_map = _create_index_map(grid, prior_zones)
        ns_index_map = _create_index_map(grid, new_zones)

        replacements = OrderedDict()

        for i_p, p in enumerate(prior_zones):
            # Retain a copy of the prior zone mask to compare with the new zone
            # mask after the zone masks are updated.
            p_mask_copy = p.mask.copy()

            # Get the new zones that intersect (`i`) the prior zone.
            ns_i_p = _intersecting_zones(ps_index_map == i_p, ns_index_map, new_zones)
            ns_i_p_ct = len(ns_i_p)

            # Get the other prior zones that intersect the new zones.
            ns_mask = np.any([n.mask for n in ns_i_p], axis=0)
            ps_i_ns = _intersecting_zones(ns_mask, ps_index_map, prior_zones)
            ps_i_ns_ct = len(ps_i_ns)

            if ps_i_ns_ct == 0:
                ps_i_ns_ct = 1
                ps_i_ns = p

            # Get successors depending on connection type.

            conn_type = _determine_connection_type(ps_i_ns_ct, ns_i_p_ct)

            p_successors = _get_successors(
                p,
                conn_type,
                ps_i_ns,
                ns_i_p,
                prior_zones,
                ps_index_map,
                replacements,
                successors,
            )

            # Update statistics.

            if conn_type in [Connection.MANY_TO_ONE, Connection.MANY_TO_MANY]:
                # Copy successors to be `captured_zones`.
                captured_zones = list(p_successors)
                if p in captured_zones:
                    captured_zones.remove(p)
                capture_ct += len(captured_zones)

                for z in captured_zones:
                    captured_mask = np.all([~p_mask_copy, z.mask], 0)
                    area = grid.cell_area_at_node[captured_mask.flatten()].sum()
                    area_captured.append(area)

            elif conn_type in [Connection.ONE_TO_MANY, Connection.MANY_TO_MANY]:
                fragment_ct += ns_i_p_ct

            # Set connection.

            p._conn_type = conn_type
            p._successors = p_successors

            successors.extend(p_successors)

        for key, value in replacements.items():
            key._mask = value.mask

        # Get unique list of successors, preserving order.

        successors = Series(successors).drop_duplicates().tolist()

    # Update the record.

    record.increment_value("fragmentations", fragment_ct)
    record.increment_value("captures", capture_ct)
    record.increment_value("area_captured_sum", sum(area_captured))

    old_value = record.get_value("area_captured_max")
    old_value_is_nan = np.isnan(old_value)
    new_value_is_greater = max(area_captured) > old_value

    if old_value_is_nan or new_value_is_greater:
        record.set_value("area_captured_max", max(area_captured))

    return successors


def _create_index_map(grid, zones):
    i = np.where([z.mask for z in zones])

    index_map = np.zeros(grid.number_of_nodes, dtype=int)
    index_map[:] = -1
    index_map[i[1]] = i[0]

    return index_map


def _intersecting_zones(condition, zone_index_map, zone_list):
    mask = np.all([condition, zone_index_map > -1], 0)
    indices = np.unique(zone_index_map[np.argwhere(mask)])
    zones = [zone_list[i] for i in indices]
    return zones


def _determine_connection_type(prior_zone_count, new_zone_count):
    """Get the connection type based on the count of prior and new zones."""
    if prior_zone_count == 0:
        if new_zone_count == 1:
            return Connection.NONE_TO_ONE

    elif prior_zone_count == 1:
        if new_zone_count == 0:
            return Connection.ONE_TO_NONE
        elif new_zone_count == 1:
            return Connection.ONE_TO_ONE
        elif new_zone_count > 1:
            return Connection.ONE_TO_MANY

    elif prior_zone_count > 1:
        if new_zone_count == 1:
            return Connection.MANY_TO_ONE
        elif new_zone_count > 1:
            return Connection.MANY_TO_MANY


def _get_successors(
    p,
    conn_type,
    ps_i_ns,
    ns_i_p,
    prior_zones,
    ps_index_map,
    replacements,
    all_successors,
):
    if conn_type == Connection.ONE_TO_NONE:
        successors = []

    elif conn_type == Connection.ONE_TO_ONE:
        # The prior zone is set as the new zone because only the one new
        # and the one prior overlap.
        n = ns_i_p[0]
        replacements[p] = n
        successors = [p]

    elif conn_type in [Connection.ONE_TO_MANY, Connection.MANY_TO_MANY]:
        # Set the successors to the new zones that overlap p.
        # Although, replace the dominant n with p.

        dn = p._get_largest_intersection(ns_i_p, exclusions=list(replacements.values()))

        successors = []

        for i, n in enumerate(ns_i_p):
            dp = n._get_largest_intersection(
                ps_i_ns, exclusions=list(replacements.keys())
            )

            if n == dn and n in replacements.values():
                d = _get_replacement(replacements, n)
                successors.append(d)
            elif n == dn:
                replacements[dp] = n
                successors.append(dp)
            elif n in replacements.values():
                d = _get_replacement(replacements, n)
                successors.append(d)
            else:
                successors.append(dp)
                replacements[dp] = n

    elif conn_type == Connection.MANY_TO_ONE:
        # Set the successor to the prior zone that intersects n the most.
        n = ns_i_p[0]
        dp = n._get_largest_intersection(ps_i_ns)

        if p == dp and n in replacements.values():
            successors = [_get_replacement(replacements, n)]
        elif p == dp:
            successors = [p]
            replacements[p] = n
        elif n in replacements.values():
            successors = [_get_replacement(replacements, n)]
        else:
            successors = [dp]
            replacements[dp] = n

    return successors


def _get_replacement(replacements, new_zone):
    for key, value in replacements.items():
        if value == new_zone:
            return key


class Zone(object):
    """Zone object of SpeciesEvolver.

    The nodes and attributes of the spatial entities that taxa populate. This
    class is not intended to be managed directly.
    """

    def __init__(self, controller, mask):
        """Initialize a zone.

        Parameters
        ----------
        controller : ZoneController
            A SpeciesEvolver ZoneController.
        mask : ndarray
            The mask of the zone. True elements of this array correspond to the
            grid nodes of the zone.
        """
        self._controller = controller
        self._mask = mask.flatten()
        self._conn_type = None
        self._successors = []

    @property
    def mask(self):
        """The mask of the zone."""
        return self._mask

    @property
    def successors(self):
        """A list of zones connected to zone at the current time step."""
        return self._successors

    @property
    def area(self):
        """Zone area calculated as the sum of cell area at grid nodes."""
        area = self._controller._grid.cell_area_at_node[self._mask].sum()
        return area

    def _get_largest_intersection(self, zones, exclusions=[]):
        node_intersection_count = []
        for z in zones:
            if z in exclusions:
                node_intersection_count.append(-1)
            else:
                n = len(np.where(np.all([z.mask, self.mask], 0))[0])
                node_intersection_count.append(n)

        if all(x == -1 for x in node_intersection_count):
            return self

        return zones[np.argmax(node_intersection_count)]
