"""Zone SpeciesEvolver object.
"""
from random import random

import numpy as np
from pandas import DataFrame


class Zone(object):
    """The nodes and attributes of the entities that species populate.

    A portion of the model domain. It is the model entity that recognizes the
    model grid changes relevant to a species.
    """

    subtype = 'base'

    # Define path types.
    NONE_TO_ONE = 'none-to-one'
    NONE_TO_MANY = 'none-to-many'
    ONE_TO_NONE = 'one-to-none'
    ONE_TO_ONE = 'one-to-one'
    ONE_TO_MANY = 'one-to-many'
    MANY_TO_NONE = 'many-to-none'
    MANY_TO_ONE = 'many-to-one'
    MANY_TO_MANY = 'many-to-many'

    def __init__(self, mask):
        """
        Parameters
        ----------
        mask : boolean ndarray
            The mask of the zone. True elements of this array correspond to the
            grid nodes of the zone.
        """
        self.mask = mask

        self.plot_color = (random(), random(), random(), 1)

    def __str__(self):
        return '<{}>'.format(self.__class__.__name__)

    @classmethod
    def _get_paths(cls, prior_zones, new_zones, prior_time, time, grid):
        """Get the data of connectivity paths across two timesteps.

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
        prior_zones : Zone list
            The zones of the prior timestep.
        new_patches : Zone list
            The zones of the current timestep.
        prior_time :

        time : float
            The current simulation time.
        grid : ModelGrid
            A Landlab ModelGrid.

        Returns
        -------
        paths : Pandas DataFrame
            Zone connectivity across two timesteps in a DataFrame with the
            columns, time, origin, destinations, and path_type.
        be_records_supplement : dictionary
            The items of this dictionary will become items in the
            SpeciesEvolver records for this time.
        """
        paths = DataFrame(columns=['time', 'origin', 'destinations',
                                   'path_type'])

        # Stack the masks for prior (p) and new (n) zones.
        ps = np.vstack(list(p.mask for p in prior_zones))
        ns = np.vstack(list(n.mask for n in new_zones))

        # Keep track of new zones replaced by prior zones. Zone in the
        # dictionary key will be replaced by their values after the prior
        # zones are processed.
        replacements = {}

        # New zones that do not intersect prior zones will be process after
        # the prior zones.
        new_overlap_not_found = np.array(new_zones)

        # The cumulation of captured grid area for zones is retained
        # in order to add it to be_records_supplement after all zones are
        # processed.
        number_of_captures = 0
        area_captured = [0]
        cell_area = grid.cellarea

        all_destinations = []

        for p in prior_zones:
            # Retain a copy of the prior zone mask to compare with the new zone
            # mask after the zone masks are updated.
            p_mask_copy = p.mask.copy()

            # Get the new zones that overlap the prior zone.
            p_in_ns = np.all([p.mask == ns, ns], 0)
            n_indices = np.argwhere(p_in_ns)
            n_indices = np.unique(n_indices[:, 0])
            n_overlaps_p = [new_zones[i] for i in n_indices]
            n_overlaps_p_count = len(n_indices)

            # Get the prior zone that is overlapped by the new zones.
            if n_overlaps_p_count > 0:
                n = new_zones[n_indices[0]]
                n_nodes = n.mask
                n_in_ps = np.all([n_nodes == ps, ps], 0)
                p_indices = np.argwhere(n_in_ps)
                p_indices = np.unique(p_indices[:, 0])
                p_overlaps_n = [prior_zones[i] for i in p_indices]
                p_overlaps_n_count = len(p_indices)
            else:
                p_overlaps_n_count = 1

            # The new zones that overlapped this p can be deleted. They
            # cannot have none-to-one path type because they intersect at
            # least this p.
            delete = np.where(new_overlap_not_found == n_overlaps_p)
            new_overlap_not_found = np.delete(new_overlap_not_found, delete)

            path_type = cls._determine_path_type(p_overlaps_n_count,
                                                 n_overlaps_p_count)

            # Determine path attributes depending upon the path type.

            if path_type == cls.ONE_TO_NONE:
                destinations = []

            elif path_type == cls.ONE_TO_ONE:
                # The prior zone is set as the new zone
                # because only the one new and the one prior overlap.
                p.mask = n.mask
                destinations = [p]

                replacements[n] = p

            elif path_type == cls.ONE_TO_MANY:
                # Set the destinations to the new zones that overlap p.
                # Although, replace the dominant n with p.
                dominant_n = p._get_largest_intersection(n_overlaps_p)
                p.mask = dominant_n.mask
                n_overlaps_p[n_overlaps_p.index(dominant_n)] = p
                destinations = n_overlaps_p

                replacements[dominant_n] = p

            elif path_type == cls.MANY_TO_ONE:
                # Set the destination to the prior zone that intersects n the
                # most.
                dominant_p = n._get_largest_intersection(p_overlaps_n)
                dominant_p.mask = n.mask
                destinations = [dominant_p]

                replacements[n] = dominant_p

            elif path_type == cls.MANY_TO_MANY:
                # Get the new zone that intersects p the most.
                dominant_n = p._get_largest_intersection(n_overlaps_p)
                dominant_n_nodes = dominant_n.mask
                dominant_n_in_ps = np.all([dominant_n_nodes == ps, ps], 0)

                # Get the prior zone that intersects the dominant n the most.
                p_indices = np.argwhere(dominant_n_in_ps)
                p_indices = np.unique(p_indices[:, 0])
                p_of_dominant_n = [prior_zones[i] for i in p_indices]
                dominant_p = dominant_n._get_largest_intersection(p_of_dominant_n)

                # Set the destination to the n that overlaps p. Although, the
                # p that greatest intersects the dominant n replaces the
                # dominant n in the destination list.
                dominant_p.mask = dominant_n.mask
                n_overlaps_p[n_overlaps_p.index(dominant_n)] = dominant_p
                destinations = n_overlaps_p

                replacements[dominant_n] = dominant_p

            # Update the capture statistics.
            if path_type in [cls.MANY_TO_ONE]:
                number_of_captures += len(n_overlaps_p)

                for n_i in n_overlaps_p:
                    captured_nodes = np.all([~p_mask_copy, n_i.mask], 0)
                    number_of_captured_nodes = len(np.where(captured_nodes)[0])
                    area_captured.append(number_of_captured_nodes * cell_area)

            paths.loc[len(paths)] = {'time': time, 'origin': p,
                      'destinations': destinations,
                      'path_type': path_type}

            all_destinations.extend(destinations)

        # Handle new zones that did not overlap prior zones.
        for n in new_overlap_not_found:
            path_type = cls._determine_path_type(p_overlaps_n_count,
                                                 n_overlaps_p_count)
            paths.loc[len(paths)] = {'time': time, 'origin': np.nan,
                      'destinations': [n], 'path_type': path_type}

        # Update zones that were replaced using the dominant zones above.
        for zone in replacements.keys():
            for index, row in paths.iterrows():
                if zone in row.destinations:
                    i = np.where(zone == np.array(row.destinations))[0][0]
                    paths.loc[index, 'destinations'][i] = replacements[zone]

        # Construct output.
        add_on = {'number_of_captures': number_of_captures,
                  'area_captured_max': max(area_captured),
                  'area_captured_sum': sum(area_captured)}

        output = {'paths': paths, 'species_evolver_records_add_on': add_on}

        return output

    @classmethod
    def _determine_path_type(cls, prior_zone_count, new_zone_count):
        # Set path type.

        if prior_zone_count == 0:
            if new_zone_count == 1:
                return cls.NONE_TO_ONE
            elif new_zone_count > 1:
                return cls.NONE_TO_MANY

        elif prior_zone_count == 1:
            if new_zone_count == 0:
                return cls.ONE_TO_NONE
            elif new_zone_count == 1:
                return cls.ONE_TO_ONE
            elif new_zone_count > 1:
                return cls.ONE_TO_MANY

        elif prior_zone_count > 1:
            if new_zone_count == 0:
                return cls.MANY_TO_NONE
            elif new_zone_count == 1:
                return cls.MANY_TO_ONE
            elif new_zone_count > 1:
                return cls.MANY_TO_MANY

    def _get_largest_intersection(self, zones):
        n_zone_overlap_nodes = []
        for z in zones:
            n = len(np.where(np.all([z.mask, self.mask], 0))[0])
            n_zone_overlap_nodes.append(n)

        return zones[np.argmax(n_zone_overlap_nodes)]

    @staticmethod
    def get_zones(grid, field_name='zone_id'):
        """Get zones using a grid field.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        field_name: string
            The name of the grid field.
        """
        # Get the unique field values.

        values = np.unique(grid.at_node[field_name])
        values = values[~np.isnan(values)]

        # Create a zone for each field value.

        zones = [None] * len(values)

        for i, value in enumerate(values):
            mask = grid.at_node[field_name] == value
            zones[i] = Zone(mask)

        return zones
