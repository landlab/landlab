#!/usr/bin/env python
"""Zone object of SpeciesEvolver."""
import numpy as np
from landlab.plot import imshow_grid

# Define path types.

_NONE_TO_NONE = 'none_to_none'
_NONE_TO_ONE = 'none_to_one'
_ONE_TO_NONE = 'one-to-none'
_ONE_TO_ONE = 'one-to-one'
_ONE_TO_MANY = 'one-to-many'
_MANY_TO_ONE = 'many-to-one'
_MANY_TO_MANY = 'many-to-many'


def _update_zones(grid, time, prior_zones, new_zones, record_add_on):
    """Resolve the spatial connectivity of zones across two timesteps.

    Paths represent the temporal connectivity of zones. The returned
    DataFrame describes all of the paths of the inputted zones. The origin
    of a path is a prior zone. The destinations of a path are the
    zones that spatially intersect the prior zone.

    Path type is determined by the number of zone connections between the
    prior and current timesteps. The path type refers to these connects and
    is designated by a string with the pattern, x-to-y where x and y are
    the descriptive counts (none, one, or many) of zone(s) at the prior and
    current timestep, respectively. For example, a path type of one-to-many
    is a zone in the prior timestep that was fragmented into many zones in
    the current timestep.

    In any of the `many` path type, a rule determines which of the prior
    zones persist as the zone in the current timestep. The zone with the
    greatest area of intersection between the prior and current timesteps
    persists to the current timestep along with the others in `new_zones`.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab ModelGrid.
    time : float
        The current simulation time.
    prior_zones : Zone list
        The zones of the prior timestep.
    new_patches : Zone list
        The zones of the current timestep.
    record_add_on : defaultdict
        A dictionary to pass values to the SpeciesEvolver record.

    Returns
    -------
    Zone list
        Zone connectivity across two timesteps in a DataFrame with the
        columns, time, origin, destinations, and path_type.
    """
    # Stats are calculated for `record_add_on`.
    fragment_ct = 0
    capture_ct = 0
    area_captured = [0]
    cell_area = grid.cellarea

    # Get `destinations`, the zones of `time`.

    if len(prior_zones) == 0 and len(new_zones) == 0:
        destinations = []

    elif len(prior_zones) > 0 and len(new_zones) == 0:
        destinations = []
        path_type = _determine_path_type(1, 0)
        for p in prior_zones:
            p.path[time] = (path_type, [])

    elif len(prior_zones) == 0 and len(new_zones) > 0:
        destinations = new_zones
        path_type = _determine_path_type(0, 1)
        for n in new_zones:
            n.path[time] = (path_type, [n])

    elif len(prior_zones) > 0 and len(new_zones) > 0:
        destinations = []
        
        # Get index maps for the prior zones (`ps`) and new zones (`ns`).
        ps_index_map = _create_index_map(grid, prior_zones)
        ns_index_map = _create_index_map(grid, new_zones)

        replacements = {}

        for i_p, p in enumerate(prior_zones):
            # Retain a copy of the prior zone mask to compare with the new zone
            # mask after the zone masks are updated.
            p_mask_copy = p.mask.copy()

            # Get the new zones that intersect (`i`) the prior zone.
            ns_i_p = _intersecting_zones(ps_index_map == i_p, ns_index_map,
                                         new_zones)
            ns_i_p_ct = len(ns_i_p)
            
            # Get the other prior zones that intersect the new zones.
            ns_mask = np.any([n.mask for n in ns_i_p], axis=0)
            ps_i_ns = _intersecting_zones(ns_mask, ps_index_map,
                                          prior_zones)
            ps_i_ns_ct = len(ps_i_ns)

            # Get destinations depending on path type.

            path_type = _determine_path_type(ps_i_ns_ct, ns_i_p_ct)

            p_destinations = _get_destinations(p, path_type, ps_i_ns, ns_i_p,
                                               prior_zones, ps_index_map,
                                               replacements, destinations)

            # Update statistics.

            if path_type in [_MANY_TO_ONE, _MANY_TO_MANY]:
                captured_zones = p_destinations.copy()
                if p in captured_zones:
                    captured_zones.remove(p)
                capture_ct += len(captured_zones)

                for z in captured_zones:
                    captured_nodes = np.all([~p_mask_copy, z.mask], 0)
                    number_of_captured_nodes = len(np.where(captured_nodes)[0])
                    area_captured.append(number_of_captured_nodes * cell_area)

            elif path_type in [_ONE_TO_MANY, _MANY_TO_MANY]:
                fragment_ct += len(ps_i_ns)

            # Update path.

            p.path[time] = (path_type, p_destinations)

            destinations.extend(p_destinations)

        for key, value in replacements.items():
            key.mask = value.mask

        destinations = list(set(destinations))

    # Update record add on.

    record_add_on['capture_count'] += capture_ct

    if max(area_captured) > record_add_on['area_captured_max']:
        record_add_on['area_captured_max'] = max(area_captured)
    record_add_on['area_captured_sum'] += sum(area_captured)
    record_add_on['fragmentation_count'] += fragment_ct
    
    return destinations


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


def _determine_path_type(prior_zone_count, new_zone_count):
    """Get the path type based on the count of prior and new zones."""
    if prior_zone_count == 0:
        if new_zone_count == 0:
            return _NONE_TO_NONE
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


def _get_destinations(p, path_type, ps_i_ns, ns_i_p, prior_zones,
                      ps_index_map, replacements, all_destinations):
    if path_type in [_NONE_TO_NONE, _ONE_TO_NONE]:
        destinations = []

    elif path_type == _ONE_TO_ONE:
        # The prior zone is set as the new zone because only the one new
        # and the one prior overlap.
        n = ns_i_p[0]
        replacements[p] = n
        destinations = [p]

    elif path_type in [_ONE_TO_MANY, _MANY_TO_MANY]:
        # Set the destinations to the new zones that overlap p.
        # Although, replace the dominant n with p.
        
        dn = p._get_largest_intersection(ns_i_p,
                                         exclusions=list(replacements.values()))

        destinations = []
        
        for i,n in enumerate(ns_i_p):
            dp = n._get_largest_intersection(ps_i_ns,
                                             exclusions=list(replacements.keys()))
            
            if n == dn and n in replacements.values():    
                d = _get_replacement(replacements, n)
                destinations.append(d)
            if n == dn:         
                replacements[dp] = n
                destinations.append(dp)  
            elif n in replacements.values():
                d = _get_replacement(replacements, n)
                destinations.append(d)
            else:
                destinations.append(dp)
                replacements[dp] = n

    elif path_type == _MANY_TO_ONE:
        # Set the destination to the prior zone that intersects n the most.
        n = ns_i_p[0]
        dp = n._get_largest_intersection(ps_i_ns)

        if p == dp and n in replacements.values():
            destinations = [_get_replacement(replacements, n)]
        elif p == dp:
            destinations = [p]
            replacements[p] = n
        elif n in replacements.values():
            destinations = [_get_replacement(replacements, n)]
        else:
            destinations = [dp]
            replacements[dp] = n

    return destinations

def _get_replacement(replacements, new_zone):
    for key, value in replacements.items():
        if value == new_zone:
            return key
    return None

class Zone(object):

    def __init__(self, mask):
        """Zone object of SpeciesEvolver.

        The nodes and attributes of the entities that species populate.

        This class is not intended to be managed directly.

        Parameters
        ----------
        mask : ndarray
            The mask of the zone. True elements of this array correspond to the
            grid nodes of the zone.
        """
        self.mask = mask.flatten()
        self._species = []
        self.path = {}

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
