"""Simulate the macroevolution processes of species.

Landlab component that populates a model domain with species and evolves the
species over time.

Component written by Nathan Lyons beginning August 20, 2017.
"""

from collections import OrderedDict
from itertools import product
from landlab import Component
from landlab.components.biota_macroevolution import Record
from landlab.core.messages import warning_message
import numpy as np
from pandas import DataFrame
from string import ascii_uppercase


class BiotaEvolver(Component):
    """Simulate the macroevolution processes of species.

    Component attributes species and zones are Pandas DataFrames that contain
    all of the species and zones, respectively.

    This component is adapted from SEAMLESS (Spatially Explicit Area Model of
    Landscape Evolution by SimulationS). See Albert et al., 2017, Systematic
    Biology 66.

    The primary method of this class is :func:`run_one_step`.

    Construction:

        BiotaEvolver(grid)
    """

    _name = 'BiotaEvolver'

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}

    def __init__(self, grid, **kwds):
        """Initialize BiotaEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        """
        Component.__init__(self, grid, **kwds)

        # Set DataFrames.
        self.record = Record()
        self.species = DataFrame(columns=['clade', 'species', 'time_appeared',
                                 'time_disappeared', 'latest_time', 'subtype', 'object'])
        self.zones = DataFrame(columns=['time_appeared', 'time_disappeared', 'latest_time',
                                        'subtype', 'object'])
        self.zone_paths = DataFrame(columns=['time', 'origin', 'destinations',
                                             'path_type'])

        # Track the clade names (keys) and the max species number (values).
        self._species_ids = {}

    # Update methods

    def run_one_step(self, time, zones_at_time, **kwds):
        """Run the macroevolution processes for a single timestep.

        Data describing the connectivity of zones over time is stored in the
        zone_paths DataFrame.

        Parameters
        ----------
        time : float
            The time in the simulation.
        zones_at_time : Zone list
            A list of the BiotaEvolver Zone objects at time.
        """
        if len(self.species) == 0:
            warning_message('No species exist. Introduce species to '
                            'BiotaEvolver.')
        else:
            self.record.append_entry(time)

            paths = self._get_zone_paths(time, zones_at_time, **kwds)

            child_species, extinct_species = self._update_species(time, paths,
                                                                  **kwds)

            # Flatten and get unique zone path destinations.
            destinations = paths.destinations.tolist()
            updated_zones = list(set(sum(destinations, [])))

            # Update DataFrames.
            self._set_disappearance_time_for_species(extinct_species, time)
            self._add_species_to_DataFrame(child_species, time)
            self._update_zones_DataFrame(updated_zones, time)
            self.zone_paths = self.zone_paths.append(paths, ignore_index=True)

    def _get_zone_paths(self, time, new_zones, **kwds):
        prior_time = self.record.get_time_prior_to_time(time)
        prior_zones = self.zones_at_time(prior_time)
        zone_types = set([type(p) for p in new_zones])

        paths = DataFrame(columns=['time', 'origin', 'destinations',
                                   'path_type'])

        for zt in zone_types:
            priors_with_type = list(filter(lambda z: isinstance(z, zt),
                                           prior_zones))
            new_with_type = list(filter(lambda z: isinstance(z, zt),
                                        new_zones))
            output = zt._get_paths(priors_with_type, new_with_type, time,
                                   self._grid, **kwds)

            # Parse path output.
            paths = paths.append(output['paths'], ignore_index=True)

            if 'biota_evolver_record_add_on' in output.keys():
                add_on = output['biota_evolver_record_add_on']
                for key in add_on.keys():
                    if key not in self.record.columns:
                        self.record[key] = np.NaN * len(self.record)

                    i = self.record.time == time
                    self.record.loc[i, key] = add_on[key]

        return paths

    def _update_species(self, time, zone_paths, **kwds):
        # Process only the species extant at the prior time.
        prior_time = self.record.get_time_prior_to_time(time)
        extant_species = self.species_at_time(prior_time)

        # Get the species that persist in `time` given the outcome of the
        # macroevolution processes of the species.

        child_species = []
        extinct_species = []

        species_types = set([type(s) for s in extant_species])

        for st in species_types:
            s_with_type = list(filter(lambda s: isinstance(s, st),
                                      extant_species))

            output = st._evolve_type(prior_time, time, s_with_type, zone_paths)

            child_species.extend(output['child_species'])
            extinct_species.extend(output['extinct_species'])

            # Set id for child species.
            for cs in output['child_species']:
                clade = cs.parent_species.clade
                cs._identifier = self._get_unused_species_id(clade)

        return child_species, extinct_species

    # Update DataFrame methods

    def _set_disappearance_time_for_species(self, species, time):
        mask = np.in1d(self.species.object.values, species)
        self.species.time_disappeared[mask] = [time] * len(species)
        self.species.latest_time[mask] = [time] * len(species)

    def _add_species_to_DataFrame(self, new_species, time):
        clade, num = ([s.identifier[0] for s in new_species],
                      [s.identifier[1] for s in new_species])
        t = [time] * len(new_species)
        st = [(s.subtype) for s in new_species]
        data = OrderedDict({'clade': clade, 'species': num, 'time_appeared': t,
                            'time_disappeared': np.nan, 'latest_time': t, 'subtype': st,
                            'object': new_species})

        new_species_df = DataFrame(data)
        self.species = self.species.append(new_species_df, ignore_index=True)

    def _update_zones_DataFrame(self, zones_at_time, time):
        z = self.zones
        z_updated, z_new = self._object_set_difference(zones_at_time, z)

        if z_updated:
            # Update the latest time value of the updated zones.
            z_updated_mask = np.in1d(self.zones.object.values, z_updated)
            z.loc[z_updated_mask, 'latest_time'] = [time] * len(z_updated)

        if z_new:
            # Add the new zones to the dataframe.
            t = [time] * len(z_new)
            st = [z.subtype for z in z_new]
            data = OrderedDict({'time_appeared': t, 'latest_time': t,
                                'subtype': st, 'object': z_new})
            new_zones = DataFrame(data)
            self.zones = z.append(new_zones, ignore_index=True)

    def _object_set_difference(self, objects, dataframe):
        objs = set(objects)
        df_objs = set(dataframe.object.values)

        # Identify the objects already in the dataframe.
        updated_objects = list(objs.intersection(df_objs))

        # Identify the objects that are new.
        new_objects = list(objs - df_objs)

        return updated_objects, new_objects

    # Species methods

    def introduce_species(self, species, time):
        """Add a species to BiotaEvolver.

        Parameters
        ----------
        species : BiotaEvolver Species
            The species to introduce.
        time : float
            The time in the simulation.
        """
        cn = self._get_unused_clade_name()
        sid = self._get_unused_species_id(cn)
        species._identifier = sid

        species_zones = species.record.loc[time, 'zones']

        self._add_species_to_DataFrame([species], time)
        self._update_zones_DataFrame(species_zones, time)
        if time not in self.record.times:
            self.record.loc[len(self.record), 'time'] = time

    def _get_unused_clade_name(self):
        alphabet = list(ascii_uppercase)
        used_ids = list(self._species_ids.keys())
        potential_clade_name = np.setdiff1d(alphabet, used_ids)

        size = 1
        while len(potential_clade_name) == 0:
            a = product(ascii_uppercase, repeat=size)
            a = [''.join(s) for s in a]
            potential_clade_name = np.setdiff1d(a, used_ids)
            size += 1

        clade_name = potential_clade_name[0]

        self._species_ids[clade_name] = -1

        return clade_name

    def _get_unused_species_id(self, clade_name):
        self._species_ids[clade_name] += 1
        return (clade_name, self._species_ids[clade_name])

    # Query methods

    def species_at_time(self, time):
        """Get the species extant at a time.

        Parameters
        ----------
        time : float
            The time in the simulation.

        Returns
        -------
        species : BiotaEvolver Species list
            The species that exist at the inputted time.
        """
        s = self.species

        appeared_before_time = s.time_appeared <= time
        disappeared_after_time = np.logical_or(s.time_disappeared > time,
                                               s.time_disappeared.isnull())

        extant_at_time = np.all([appeared_before_time, disappeared_after_time],
                                0)

        species = s.object[extant_at_time].tolist()

        return species

    def zones_at_time(self, time):
        """Get the zones extant at a time.

        Parameters
        ----------
        time : float
            The time in the simulation.

        Returns
        -------
        zones : BiotaEvolver Zone list
            The zones that exist at the inputted time.
        """
        con1 = time >= self.zones.time_appeared
        con2 = time <= self.zones.latest_time
        is_at_time = np.all([con1, con2], axis=0)

        zones = self.zones.object[is_at_time].tolist()

        return zones

    def species_with_identifier(self, identifier_element,
                                return_objects=False):
        """Get species with an identifier element.

        A singular species is returned when *identifier_element* is a tuple
        with elements that match a species in the *species* DataFrame. The
        first element is the clade name and the second element is the species
        number.

        The species of a clade are returned when *identifier_element* is
        a string that matches a clade name in the *species* DataFrame.

        The species that share a species number are returned when
        *identifier_element* is an integer that matches species number in the
        *species* DataFrame.

        By default, the species with *identifier_element* will be returned in a
        DataFrame. Alternatiely, a list of Species objects can be returned by
        setting *return_objects* to ``True``. A singular species is returned
        when *identifier_element* is a tuple. Otherwise, multiple species are
        returned.

        Parameters
        ----------
        identifier_element : tuple, str, or int
            The identifier element of the species to return.
        return_objects : boolean, optional
            ``True`` returns species as BiotaEvolver objects. ``False``, the
            default, returns a DataFrame.

        Returns
        -------
        Pandas DataFrame, BiotaEvolver Species, or BiotaEvolver Species list
            The species with identifiers that matched *identifier_element*. The
            type of the return object is set by *return_objects*.
        """
        s_in = self.species

        element_type = type(identifier_element)

        if element_type == tuple:
            # Get a singular species using a clade name and number.
            clade = identifier_element[0]
            num = identifier_element[1]

            if not np.all([len(identifier_element) == 2,
                           isinstance(clade, str), isinstance(num, int)], 0):
                raise TypeError('*identifier_element* when it is a tuple must '
                                'have a length of 2. The first element must '
                                'be a string, and the second must be an '
                                'integer.')

            s_out = self.species.loc[np.all([s_in.clade == clade,
                                             s_in.species == num], 0)]

        elif element_type == str:
            # Get the species of a clade.
            s_out = s_in.loc[s_in.clade == identifier_element]

        elif element_type == int:
            # Get the species that match a number.
            s_out = self.species.loc[s_in.species == identifier_element]

        else:
            raise TypeError('*identifier_element* must be a tuple, string, or '
                            'integer.')

        if len(s_out) == 0:
            return None
        elif return_objects:
            s_out_list = s_out.object.tolist()
            if element_type == tuple:
                return s_out_list[0]
            else:
                return s_out_list
        else:
            return s_out
