"""Simulate the macroevolution processes of species.

Landlab component that populates a model domain with species and evolves the
species over time.

Component written by Nathan Lyons beginning August 2017.
"""
from collections import OrderedDict
from itertools import product
from string import ascii_uppercase

import numpy as np
import pandas as pd

from landlab import Component
from landlab.core.messages import warning_message
from landlab.data_record import DataRecord
from .zone import Zone


class SpeciesEvolver(Component):
    """Simulate the macroevolution processes of species.

    This component evolves the species introduced to a model grid. The
    evolution rules are programmed in `Species` objects. The geographic range
    of species is controlled by `Zone` objects. Component attributes `species`
    and `zones` are Pandas DataFrames that contain all of the species and zones
    of a model, respectively.

    A general workflow:

    1.  The component is instantiated.
    2.  Species are introduced using the :func:`introduce_species` method.
    3.  The primary method of this class is :func:`run_one_step`. The temporal
        connectivity of zones is identified by this method and stored in
        `zone_paths`. The Species method :func:`evolve` is called through
        :func:`run_one_step`.

        Species evolution rules take into account
    4.

    Species are automatically assigned identifiers in the order that they are
    introduced and created by parent species. Clades are lettered from A to Z
    then AA to AZ and so forth. Clade members are numbered in the order of
    appearance.

    This component is adapted from SEAMLESS (Spatially Explicit Area Model of
    Landscape Evolution by SimulationS). See Albert et al., 2017, Systematic
    Biology 66.

    Attributes
    ----------
    records : RecordCollection
        A Landlab RecordCollection.
    species : DataFrame
        A Pandas DataFrame containing species data.
    zones : DataFrame
        A Pandas DataFrame containing zone data.
    zone_paths : DataFrame
        A Pandas DataFrame containing zone path data.
    """
    _name = 'SpeciesEvolver'

    _input_var_names = ()

    _output_var_names = ()

    _var_units = {}

    _var_mapping = {}

    _var_doc = {}

    def __init__(self, grid, initial_time=0):
        """Instantiate SpeciesEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        initial_time :

        Examples
        --------
        >>> from landlab.components import (FastscapeEroder, FlowRouter,
                                            LinearDiffuser, SpeciesEvolver)
        >>> from landlab import RasterModelGrid
        >>>

        """
        Component.__init__(self, grid)

        # Create DataFrames.

        self.dataRecord = DataRecord(grid, time=[initial_time])
        self.species = pd.DataFrame(columns=['clade', 'species',
                                             'time_appeared', 'latest_time',
                                             'subtype', 'object'])
        self.zones = pd.DataFrame(columns=['time_appeared', 'latest_time',
                                           'subtype', 'object'])
        self.zone_paths = pd.DataFrame(columns=['time', 'origin',
                                                'destinations', 'path_type'])

        # Track the clade names (keys) and the max species number (values).

        self._species_ids = {}

        # Set initial zones.

        zones_at_time = Zone.get_zones(self._grid)
        self._update_zones_DataFrame(initial_time, zones_at_time)

    # Update methods

    def run_one_step(self, dt):
        """Run macroevolution processes for a single timestep.

        Data describing the connectivity of zones over time is stored in the
        zone_paths DataFrame.

        Parameters
        ----------
        dt : float
            The model time step.
        """
        # Update time.

        time = self.dataRecord.latest_time + dt
        self.dataRecord.add_record(time=[time])

        # Process zones and species.

        if len(self.species) == 0:
            print(warning_message('No species exist. Introduce species to '
                                  'SpeciesEvolver.'))
            return

        zones_at_time = Zone.get_zones(self._grid)

        paths = self._get_zone_paths(time, zones_at_time)
        survivors = self._get_surviving_species(time, paths)

        # Flatten and get unique zone path destinations.

        destinations = paths.destinations.tolist()
        updated_zones = list(set(sum(destinations, [])))

        # Update DataFrames.

        self._update_species_DataFrame(time, survivors)
        self._update_zones_DataFrame(time, updated_zones)
        self.zone_paths = self.zone_paths.append(paths, ignore_index=True)

    def _get_zone_paths(self, time, new_zones):
        prior_time = self.dataRecord.prior_time
        prior_zones = self.zones_at_time(prior_time, return_objects=True)

        zone_types = set([type(p) for p in new_zones])

        paths = pd.DataFrame(columns=['time', 'origin', 'destinations',
                                      'path_type'])

        for zt in zone_types:
            priors_with_type = list(filter(lambda z: isinstance(z, zt),
                                           prior_zones))

            new_with_type = list(filter(lambda z: isinstance(z, zt),
                                        new_zones))

            output = zt._get_paths(priors_with_type, new_with_type, prior_time,
                                   time, self._grid)

            # Parse path output.
            paths = paths.append(output['paths'], ignore_index=True)

            if 'species_evolver_records_add_on' in output.keys():
                add_on = output['species_evolver_records_add_on']
                for key, value in add_on.items():
                    self.dataRecord.add_record(time=[time],
                                               new_record={key: (['time'],
                                                                 [value])})

        return paths

    def _get_surviving_species(self, time, zone_paths):
        # Process only the species extant at the prior time.
        prior_time = self.dataRecord.prior_time
        extant_species = self.species_at_time(prior_time, return_objects=True)

        # Get the species that persist in `time` given the outcome of the
        # macroevolution processes of the species.

        surviving_species = []

        species_types = set([type(s) for s in extant_species])

        for st in species_types:
            s_with_type = list(filter(lambda s: isinstance(s, st),
                                      extant_species))

            output = st.evolve_type(s_with_type, zone_paths)

            surviving_species.extend(output['surviving_parent_species'])
            surviving_species.extend(output['child_species'])

            # Set id for child species.
            for cs in output['child_species']:
                clade = cs.parent_species.clade
                cs._identifier = self._get_unused_species_id(clade)

            if 'species_evolver_records_add_on' in output.keys():
                add_on = output['species_evolver_records_add_on']
                for key, value in add_on.items():
                    self.dataRecord.add_record(time=[time],
                                               new_record={key: (['time'],
                                                                 [value])})

        return surviving_species

    # Update DataFrame methods

    def _update_species_DataFrame(self, time, species_at_time):
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
            new_species = pd.DataFrame(data)
            self.species = sdf.append(new_species, ignore_index=True)

    def _update_zones_DataFrame(self, time, zones_at_time):
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
            new_zones = pd.DataFrame(data)
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

        The species is introduced at the latest time in the dataRecord.

        Parameters
        ----------
        species : SpeciesEvolver Species
            The species to introduce.

        Examples
        --------
        """
        if species in self.species.object.values:
            msg = 'The species object, {} was already introduced.'
            raise Exception(msg.format(species))

        # Set species identifier.

        if species.identifier == None:
            if species.parent_species == None:
                cn = self._get_unused_clade_name()
            else:
                cn = species.parent_species.clade

            sid = self._get_unused_species_id(cn)
            species._identifier = sid

        # Update dataframes.

        time = self.dataRecord.latest_time

        self._update_species_DataFrame(time, [species])
#        self._update_zones_DataFrame(time, species.zones)

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

        self._species_ids[clade_name] = None

        return clade_name

    def _get_unused_species_id(self, clade_name):
        if self._species_ids[clade_name] == None:
            self._species_ids[clade_name] = 0
        else:
            self._species_ids[clade_name] += 1
        return (clade_name, self._species_ids[clade_name])

    # Query methods

    def species_at_time(self, time, return_objects=False):
        """Get the species that exist at a time.

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

        Examples
        --------

        """
        return self._objects_at_time(time, self.species, return_objects)

    def zones_at_time(self, time, return_objects=False):
        """Get the zones that exist at a time.

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

        Examples
        --------


        """
        return self._objects_at_time(time, self.zones, return_objects)

    def _objects_at_time(self, time, df, return_objects):
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
        """Get species using identifiers.

        A singular species is returned when *identifier_element* is a tuple
        with elements that match species in the *species* DataFrame. The
        first element of the tuple is the clade name and the second element is
        the species number.

        The species of a clade are returned when *identifier_element* is
        a string that matches a clade name in the *species* DataFrame.

        The species that share a species number are returned when
        *identifier_element* is an integer that matches species number in the
        *species* DataFrame.

        By default, the species with *identifier_element* will be returned in a
        DataFrame. Alternatively, a list of Species objects can be returned by
        setting *return_objects* to ``True``. A singular species is returned
        when *identifier_element* is a tuple. Otherwise, the species object(s)
        are returned in a list.

        `None` is returned if no species have an identifier that matches
        *identifier_element*.

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

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import SpeciesEvolver
        >>> from landlab.components.species_evolution import Species, Zone
        >>> import numpy as np
        >>> mg = RasterModelGrid((3, 5))
        >>> zone_id = np.array([np.nan, np.nan, np.nan, np.nan, np.nan,
        ...                     np.nan,      1,      2,      3, np.nan,
        ...                     np.nan, np.nan, np.nan, np.nan, np.nan])
        >>> zone_field = mg.add_field('node', 'zone_id', zone_id)
        >>> zones = Zone.get_zones(mg)
        >>> se = SpeciesEvolver(mg)

        Instantiate and introduce species.

        >>> species1 = Species(zones[0])
        >>> species2 = Species(zones[1])
        >>> species3 = Species(zones[2], parent_species=species2)
        >>> se.introduce_species(species1)
        >>> se.introduce_species(species2)
        >>> se.introduce_species(species3)

        Get all the species introduced in a dataframe.

        >>> se.species
          clade species time_appeared latest_time subtype         object
        0     A       0             0           0    base  <Species A.0>
        1     B       0             0           0    base  <Species B.0>
        2     B       1             0           0    base  <Species B.1>

        Get the species, B.0.

        >>> se.species_with_identifier(('B', 0))
          clade species time_appeared latest_time subtype         object
        1     B       0             0           0    base  <Species B.0>

        Get all of the species in clade, B.

        >>> se.species_with_identifier('B')
          clade species time_appeared latest_time subtype         object
        1     B       0             0           0    base  <Species B.0>
        2     B       1             0           0    base  <Species B.1>

        Get all of the species with the same number, 0, despite the clade.

        >>> se.species_with_identifier(0)
          clade species time_appeared latest_time subtype         object
        0     A       0             0           0    base  <Species A.0>
        1     B       0             0           0    base  <Species B.0>

        Get the species, B.0 as an object rather than dataframe.

        >>> species_obj = se.species_with_identifier(('B', 0),
                                                     return_objects=True)
        >>> species_obj.identifier
        ('B', 0)
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
