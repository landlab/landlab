"""Simulate the macroevolution processes of species.

Landlab component that populates a model domain with species and evolves the
species over time.

Component written by Nathan Lyons beginning August 2017.
"""

from collections import OrderedDict
from itertools import product
from landlab import Component
from landlab.components.species_evolution import RecordCollection
from landlab.core.messages import warning_message
import numpy as np
from pandas import DataFrame
from string import ascii_uppercase


class SpeciesEvolver(Component):
    """Simulate the macroevolution processes of species.

    Component attributes species and zones are Pandas DataFrames that contain
    all of the species and zones, respectively.

    This component is adapted from SEAMLESS (Spatially Explicit Area Model of
    Landscape Evolution by SimulationS). See Albert et al., 2017, Systematic
    Biology 66.

    The primary method of this class is :func:`run_one_step`.

    Construction:

        SpeciesEvolver(grid)
    """

    _name = 'SpeciesEvolver'

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}

    def __init__(self, grid, **kwds):
        """Initialize SpeciesEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        """
        Component.__init__(self, grid, **kwds)

        # Set DataFrames.
        self.records = RecordCollection()
        self.species = DataFrame(columns=['clade', 'species', 'time_appeared',
                                          'latest_time', 'subtype', 'object'])
        self.zones = DataFrame(columns=['time_appeared', 'latest_time',
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
            A list of the SpeciesEvolver Zone objects at time.
        """
        if len(self.species) == 0:
            print(warning_message('No species exist. Introduce species to '
                                  'SpeciesEvolver.'))
        elif len(self.records.DataFrame) == 0:
            print(warning_message('Species must be introduced at a time prior'
                                  ' to the ``run_one_step`` time.'))
        else:
            self.records.insert_time(time)

            paths = self._get_zone_paths(time, zones_at_time, **kwds)
            survivors = self._get_surviving_species(time, paths, **kwds)

            # Flatten and get unique zone path destinations.
            destinations = paths.destinations.tolist()
            updated_zones = list(set(sum(destinations, [])))

            # Update DataFrames.
            self._update_species_DataFrame(survivors, time)
            self._update_zones_DataFrame(updated_zones, time)
            self.zone_paths = self.zone_paths.append(paths, ignore_index=True)

    def _get_zone_paths(self, time, new_zones, **kwds):
        prior_time = self.records.get_time_prior_to_time(time)
        prior_zones = self.zones_at_time(prior_time, return_objects=True)
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

            if 'species_evolver_records_add_on' in output.keys():
                add_on = output['species_evolver_records_add_on']
                self.records.modify_time(time, add_on)

        return paths

    def _get_surviving_species(self, time, zone_paths, **kwds):
        # Process only the species extant at the prior time.
        prior_time = self.records.get_time_prior_to_time(time)
        extant_species = self.species_at_time(prior_time, return_objects=True)

        # Get the species that persist in `time` given the outcome of the
        # macroevolution processes of the species.

        surviving_species = []

        species_types = set([type(s) for s in extant_species])

        for st in species_types:
            s_with_type = list(filter(lambda s: isinstance(s, st),
                                      extant_species))

            output = st._evolve_type(prior_time, time, s_with_type, zone_paths)

            surviving_species.extend(output['surviving_parent_species'])
            surviving_species.extend(output['child_species'])

            # Set id for child species.
            for cs in output['child_species']:
                clade = cs.parent_species.clade
                cs._identifier = self._get_unused_species_id(clade)

        return surviving_species

    # Update DataFrame methods

    def _update_species_DataFrame(self, species_at_time, time):
        sdf = self.species
        s_updated, s_new = self._object_set_difference(species_at_time, sdf)

        if s_updated:
            # Update the latest time value of the updated species.
            s_updated_mask = np.in1d(sdf.object.values, s_updated)
            sdf.latest_time[s_updated_mask] = [time] * len(s_updated)

        if s_new:
            clade = [s.identifier[0] for s in s_new]
            s_number = [s.identifier[1] for s in s_new]
            t = [time] * len(s_new)
            st = [(s.subtype) for s in s_new]
            data = OrderedDict({'clade': clade, 'species': s_number,
                                'time_appeared': t, 'latest_time': t,
                                'subtype': st, 'object': s_new})
            new_species = DataFrame(data)
            self.species = sdf.append(new_species, ignore_index=True)

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

    def introduce_species(self, species):
        """Add a species to SpeciesEvolver.

        Parameters
        ----------
        species : SpeciesEvolver Species
            The species to introduce.
        """
        cn = self._get_unused_clade_name()
        sid = self._get_unused_species_id(cn)
        species._identifier = sid

        time = species.records.time__latest

        species_zones = species.records.get_value(time, 'zones')

        self._update_species_DataFrame([species], time)
        self._update_zones_DataFrame(species_zones, time)
#        if time not in self.records.times:
#            self.records.DataFrame.loc[len(self.records), 'time'] = time

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

    def species_at_time(self, time, return_objects=False):
        """Get the species extant at a time.

        Parameters
        ----------
        time : float
            The time in the simulation.
        return_objects : boolean, optional
            ``True`` returns species as SpeciesEvolver objects. ``False``, the
            default, returns a DataFrame.

        Returns
        -------
        species : SpeciesEvolver Species list
            The species that exist at the inputted time.
        """
        return self._object_at_time(self.species, time, return_objects)

    def zones_at_time(self, time, return_objects=False):
        """Get the zones extant at a time.

        Parameters
        ----------
        time : float
            The time in the simulation.
        return_objects : boolean, optional
            ``True`` returns zones as SpeciesEvolver objects. ``False``, the
            default, returns a DataFrame.

        Returns
        -------
        zones : SpeciesEvolver Zone list
            The zones that exist at the inputted time.
        """
        return self._object_at_time(self.zones, time, return_objects)

    def _object_at_time(self, df, time, return_objects):
        appeared_before_time = df.time_appeared <= time
        present_at_time = df.latest_time >= time
        extant_at_time = np.all([appeared_before_time, present_at_time], 0)

        df_time = df.loc[extant_at_time]

        if return_objects:
            return df_time.object.tolist()
        else:
            return df_time

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
        DataFrame. Alternatively, a list of Species objects can be returned by
        setting *return_objects* to ``True``. A singular species is returned
        when *identifier_element* is a tuple. Otherwise, multiple species are
        returned.

        Parameters
        ----------
        identifier_element : tuple, str, or int
            The identifier element of the species to return.
        return_objects : boolean, optional
            ``True`` returns species as SpeciesEvolver objects. ``False``, the
            default, returns a DataFrame.

        Returns
        -------
        DataFrame, SpeciesEvolver Species, or SpeciesEvolver Species list
            The species with identifiers that matched *identifier_element*. The
            type of the return object is set by *return_objects*.
        """
        sdf = self.species

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

            s_out = self.species.loc[np.all([sdf.clade == clade,
                                             sdf.species == num], 0)]

        elif element_type == str:
            # Get the species of a clade.
            s_out = sdf.loc[sdf.clade == identifier_element]

        elif element_type == int:
            # Get the species that match a number.
            s_out = self.species.loc[sdf.species == identifier_element]

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

    def get_area_species_data(self):
        """Get the area and number of species of each zone.

        """
        time = self.records.time__latest

        species_time = self.species_at_time(time, return_objects=True)
        zones_time = self.zones_at_time(time, return_objects=True)

        data = {'area': [], 'number_of_species': []}

        for z in zones_time:
            mask = z.mask
#            data['area'].append(np.sum(mask) * cell_area)
            if len(np.where(mask)[0]) == 0:
                data['area'].append(0)
            else:
                data['area'].append(self._grid.at_node['drainage_area'][mask].max())

            z_species_count = 0

            for s in species_time:
                zones_species = s.records.get_value(time, 'zones')
                if z in zones_species:
                    z_species_count += 1

            data['number_of_species'].append(z_species_count)

        return data
