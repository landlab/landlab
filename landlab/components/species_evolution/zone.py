#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Zone functions and class of SpeciesEvolver."""
from collections import OrderedDict

import numpy as np


# Define connection types.

_NONE_TO_NONE = 'none_to_none'
_NONE_TO_ONE = 'none_to_one'
_ONE_TO_NONE = 'one-to-none'
_ONE_TO_ONE = 'one-to-one'
_ONE_TO_MANY = 'one-to-many'
_MANY_TO_ONE = 'many-to-one'
_MANY_TO_MANY = 'many-to-many'


def _update_zones(grid, prior_zones, new_zones, record):
    """Resolve the spatial connectivity of zones across two time steps.

    This method iterates over each zone of the prior time step to identify the
    zones of the current time step it spatially intersects. This method updates
    the zone attribute, ``successors``. Successor zones are the new zones
    existing at the current time step that are the continuation of prior time
    step zones referred to as the predecessor zones.

    The type of connection between predecessor and successor zones is described
    by the number of predecessors and successors that overlap each other. The
    connection type is represented by a string with the pattern, x-to-y where x
    and y are the descriptive counts (none, one, or many) of zone(s) at the
    prior and current time step, respectively. For example, a connection type
    of one-to-many is a zone in the prior time step that spatially overlaps
    multiple zones in the current time step.

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

    # Get `successors`, the zones of `time`.

    if len(prior_zones) == 0 and len(new_zones) == 0:
        successors = []

    elif len(prior_zones) > 0 and len(new_zones) == 0:
        successors = []
        conn_type = _determine_connection_type(1, 0)
        for p in prior_zones:
            p._conn_type = conn_type
            p._successors = []

    elif len(prior_zones) == 0 and len(new_zones) > 0:
        successors = new_zones
        conn_type = _determine_connection_type(0, 1)
        for n in new_zones:
            n._conn_type = conn_type
            n._successors = [n]

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
            ns_i_p = _intersecting_zones(
                ps_index_map == i_p, ns_index_map, new_zones
            )
            ns_i_p_ct = len(ns_i_p)

            # Get the other prior zones that intersect the new zones.
            ns_mask = np.any([n.mask for n in ns_i_p], axis=0)
            ps_i_ns = _intersecting_zones(ns_mask, ps_index_map, prior_zones)
            ps_i_ns_ct = len(ps_i_ns)

            # Get successors depending on connection type.

            conn_type = _determine_connection_type(ps_i_ns_ct, ns_i_p_ct)

            p_successors = _get_successors(
                p, conn_type, ps_i_ns, ns_i_p, prior_zones, ps_index_map,
                replacements, successors
            )

            # Update statistics.

            if conn_type in [_MANY_TO_ONE, _MANY_TO_MANY]:
                # Copy successors to be `captured_zones`.
                captured_zones = list(p_successors)
                if p in captured_zones:
                    captured_zones.remove(p)
                capture_ct += len(captured_zones)

                for z in captured_zones:
                    captured_mask = np.all([~p_mask_copy, z.mask], 0)
                    area = sum(grid.cell_area_at_node[captured_mask.flatten()])
                    area_captured.append(area)

            elif conn_type in [_ONE_TO_MANY, _MANY_TO_MANY]:
                fragment_ct += ns_i_p_ct

            # Set connection.

            p._conn_type = conn_type
            p._successors = p_successors

            successors.extend(p_successors)

        for key, value in replacements.items():
            key._mask = value.mask

        # Get unique list of successors, preserving order.

        ss = set()
        successors = [x for x in successors if not (x in ss or ss.add(x))]

    # Update the record.

    record.increment_value('fragmentation_count', fragment_ct)
    record.increment_value('capture_count', capture_ct)
    record.increment_value('area_captured_sum', sum(area_captured))

    old_value = record.get_value('area_captured_max')
    old_value_is_nan = np.isnan(old_value)
    new_value_is_greater = max(area_captured) > old_value

    if old_value_is_nan or new_value_is_greater:
        record.set_value('area_captured_max', max(area_captured))

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
            return _NONE_TO_ONE

    elif prior_zone_count == 1:
        if new_zone_count == 0:
            return _ONE_TO_NONE
        elif new_zone_count == 1:
            return _ONE_TO_ONE
        elif new_zone_count > 1:
            return _ONE_TO_MANY

    elif prior_zone_count > 1:
        if new_zone_count == 1:
            return _MANY_TO_ONE
        elif new_zone_count > 1:
            return _MANY_TO_MANY


def _get_successors(
    p, conn_type, ps_i_ns, ns_i_p, prior_zones, ps_index_map, replacements,
    all_successors
):

    if conn_type == _ONE_TO_ONE:
        # The prior zone is set as the new zone because only the one new
        # and the one prior overlap.
        n = ns_i_p[0]
        replacements[p] = n
        successors = [p]

    elif conn_type in [_ONE_TO_MANY, _MANY_TO_MANY]:
        # Set the successors to the new zones that overlap p.
        # Although, replace the dominant n with p.

        dn = p._get_largest_intersection(
            ns_i_p, exclusions=list(replacements.values())
        )

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

    elif conn_type == _MANY_TO_ONE:
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

    The nodes and attributes of the spatial entities that species populate.
    This class is not intended to be managed directly.
    """
    def __init__(self, mask):
        """Initialize a zone.

        Parameters
        ----------
        mask : ndarray
            The mask of the zone. True elements of this array correspond to the
            grid nodes of the zone.
        """
        self._mask = mask.flatten()
        self._species = set()
        self._conn_type = None
        self._successors = []

    @property
    def mask(self):
        """The mask of the zone."""
        return self._mask

    @property
    def species(self):
        """The set of species that inhabit the zone.

        The order that species were added to this attribute is not preserved.
        """
        return self._species

    @property
    def successors(self):
        """A list of zones connected to zone at the current time step."""
        return self._successors

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
