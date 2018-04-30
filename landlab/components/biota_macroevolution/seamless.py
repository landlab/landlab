"""Simulate the macroevolution processes of species.

Landlab component that tracks the species in a model domain.

Component written by Nathan Lyons beginning August 20, 2017.
"""

from collections import OrderedDict
from landlab import Component
from landlab.components.biota_macroevolution import Record
from landlab.core.messages import warning_message
import numpy as np
import pandas as pd
from pprint import pprint
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

    def __init__(self, grid, map_vars=None, **kwds):
        """Initialize BiotaEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        """
        Component.__init__(self, grid, map_vars=None, **kwds)

        # Create DataFrames.
        columns = ['time_appeared', 'latest_time', 'time_disappeared', 'subtype', 'object']
        index = pd.MultiIndex(levels=[[],[]], labels=[[],[]],
                      names=['clade', 'number'])
        self.species = pd.DataFrame(columns=columns, index=index)
        self.zones = pd.DataFrame(columns=columns)

        self.record = Record()
        self.record['time'] = np.nan

        self.zone_paths = pd.DataFrame(columns=['time', 'origin',
                                                'destinations', 'path_type'])

        # Keep track of the clades (keys) and the max species number (values).
        self._clades = {}

    # Update methods

    def run_one_step(self, time, zones_at_time, **kwargs):
        """Run the macroevolution processes for a single timestep.

        Parameters
        ----------
        time : float
            The time in the simulation.
        zones_at_time : Zone list
            A list of the BiotaEvolver Zone objects at time.
        """
        self.evolve_species(time, zones_at_time, **kwargs)

    def evolve_species(self, time, zones_at_time, **kwargs):
        """Run the macroevolution processes of all extant species.

        Data describing the connectivity of zones over time is stored in the
        zone_paths DataFrame.

        Parameters
        ----------
        time : float
            The time in the simulation.
        zones_at_time : Zone list
            A list of BiotaEvolver Zone objects.
        """
        self._update_record(time)

        if len(self.species) == 0:
            warning_message('No species exist. Introduce species to '
                            'BiotaEvolver.')
        else:
            prior_time = self.record.time__latest

            paths = self._get_zone_paths(prior_time, time, zones_at_time,
                                         **kwargs)

            time_species = self._get_surviving_species(prior_time, time, paths)

            # Flatten and get unique zone path destinations.
            destinations = paths.destinations.tolist()
            updated_zones = list(set(sum(destinations, [])))
            self.zone_paths = self.zone_paths.append(paths, ignore_index=True)

            # Update DataFrames.
            self._update_species_data_frame(time_species, time)
            self._update_zones_data_frame(updated_zones, time)

    def _get_zone_paths(self, prior_time, time, new_zones, **kwargs):
        prior_zones = self.zones_at_time(prior_time)
        zone_types = set([type(p) for p in new_zones])

        paths = pd.DataFrame(columns=['time', 'origin', 'destinations',
                                      'path_type'])

        for zt in zone_types:
            priors_with_type = list(filter(lambda z: isinstance(z, zt),
                                           prior_zones))
            new_with_type = list(filter(lambda z: isinstance(z, zt),
                                        new_zones))
            output = zt._get_paths(priors_with_type, new_with_type, time,
                                   self._grid, **kwargs)

            # Parse path output.
            paths = paths.append(output['paths'], ignore_index=True)

            if 'biota_evolver_record_add_on' in output.keys():
                for key in output['biota_evolver_record_add_on'].keys():
                    if key not in self.record.columns:
                        self.record[key] = np.NaN * len(self.record)

                    i = self.record.time == time
                    self.record.loc[i, key] = output[1][key]

        return paths

    def _get_surviving_species(self, prior_time, time, zone_paths):
        # Process only the species extant at the prior time.
        extant_species = self.species_at_time(prior_time)

        origins = zone_paths.origin.tolist()

        # Get the species that persist in `next_time` given the outcome of the
        # macroevolution processes of the species.

        surviving_species = []

        for es in extant_species:
            # Get paths that include the zone origin of this species.
            species_zones = es.record.loc[prior_time, 'zones']
            indices = np.where(np.isin(origins, species_zones))[0]

            if len(indices) > 0:
                es_paths = zone_paths.loc[indices]

                output = es.run_macroevolution_processes(time, es_paths)
                if output['species_persists']:
                    surviving_species.extend(es)
                if output['child_species']:
                    surviving_species.extend(output['child_species'])

                    # Set id for child species.
                    for cs in output['child_species']:
                        cs._identifier = self._get_unused_species_id(es.clade)

        return surviving_species

    # Update data frame methods

    def _update_species_data_frame(self, species_at_time, time):
        s_updated, s_new = self._object_set_difference(species_at_time,
                                                       self.species)

        if s_updated:
            # Update the latest time value of the updated species.
            s_updated_mask = np.in1d(self.species.object.values, s_updated)
            self.species.latest_time[s_updated_mask] = [time] * len(s_updated)

        if s_new:
            ids = [s.identifier for s in s_new]
            index = pd.MultiIndex.from_tuples(ids, names=['clade', 'number'])
            t = [time] * len(s_new)
            st = [(s.subtype) for s in s_new]
            data = OrderedDict({'time_appeared': t, 'latest_time': t,
                                'subtype': st, 'object': s_new})
            new_species_df = pd.DataFrame(data, index=index)
            self.species = self.species.append(new_species_df).unstack().stack()

    def _update_zones_data_frame(self, zones_at_time, time):
        z_updated, z_new = self._object_set_difference(zones_at_time,
                                                       self.zones)

        if z_updated:
            # Update the latest time value of the updated zones.
            z_updated_mask = np.in1d(self.zones.object.values, z_updated)
            self.zones.latest_time[z_updated_mask] = [time] * len(z_updated)

        if z_new:
            # Add the new zones to the dataframe.
            t = [time] * len(z_new)
            st = [z.subtype for z in z_new]
            data = OrderedDict({'time_appeared': t, 'latest_time': t,
                                'subtype': st, 'object': z_new})
            new_zones = pd.DataFrame(data)
            self.zones = self.zones.append(new_zones, ignore_index=True)
            self.zones.index.name = 'index'

    def _update_record(self, time):
        data = OrderedDict({'time': [time]})
        new_record = pd.DataFrame(data)
        self.record = self.record.append(new_record, ignore_index=True)

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
        cid = self._get_unused_clade_id()
        sid = self._get_unused_species_id(cid)
        species._identifier = sid

        species_zones = species.record.loc[time, 'zones']

        self._update_species_data_frame([species], time)
        self._update_zones_data_frame(species_zones, time)
        if time not in self.record.times:
            self.record.loc[len(self.record), 'time'] = time

    def _get_unused_clade_id(self):
        used_ids = list(self._clades.keys())

        alphabet = list(ascii_uppercase)
        clade_id = np.setdiff1d(alphabet, used_ids)

        duplicator = alphabet
        while len(clade_id) == 0:
           duplicator = np.core.defchararray.add(alphabet, duplicator)
           clade_id = np.setdiff1d(duplicator, used_ids)

        self._clades[clade_id[0]] = 0

        return clade_id[0]

    def _get_unused_species_id(self, clade_id):
        self._clades[clade_id] += 1
        return (clade_id, self._clades[clade_id])

    #TODO
    def get_tree(self):
        """Get phylogenetic tree.
        """
        tree = {}

        clades = list(set(self.species.index.get_level_values('clade').tolist()))
        clades.sort()

        for c in clades:
            clade_species = self.species_of_clade(c)
            for s in clade_species:
                if s.parent_species == -1:
                    parent_id = -1
                else:
                    parent_id = s.parent_species
                tree.setdefault(parent_id, [])

                sid = s.identifier
                tree[parent_id].append(s)

        return tree

    #TODO
    def print_tree(self):
        tree = self.get_tree()

        readable_tree = {}
        for time in tree.keys():
            for p in tree[time].keys():
                if p == -1:
                    pid = p
                else:
                    pid = p.identifier
                species_ids = [s.identifier for s in tree[time][p]]
                readable_tree.setdefault(time, {})
                readable_tree[time][pid] = species_ids

        pprint(readable_tree)

    def calculate_extinctions_per_million_species_per_year(self, time):
        """Get the quantity, extinctions per million species per year (E/MSY).

        E/MSY is calculated for the species that existed at or before time.

        Parameters
        ----------
        time : float
            The final time included in the calculation.

        Returns
        -------
        e_per_msy : float
            Extinctions per million species per year.
        """
        # Get all of the species that existed at or before time.
        species = self.species.loc[self.species.time_appeared <= time]

        # Cumulate the extant years of species.
        species_years = sum(species.latest_time - species.time_appeared)

        number_of_all_species = len(species)
        number_of_extant_species_at_time = len(self.species_at_time(time))
        number_of_extinctions = (number_of_all_species -
                                 number_of_extant_species_at_time)

        e_per_msy = number_of_extinctions / species_years * 1e6

        return e_per_msy

    # Convenience methods.

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
        species = self.species.object[self.species.latest_time == time].tolist()
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

    def species_of_clade(self, clade):
        """Get the species of a clade.

        """
        clade_species = self.species.index.get_level_values('clade') == clade
        return self.species.object[clade_species].tolist()

    def species_with_identifier(self, clade, species_number):
        clade_species = self.species.index.get_level_values('clade') == clade
        species = self.species.index.get_level_values('number') == species_number
        clade_and_species = np.all([clade_species, species], 0)
        return self.species.object[clade_and_species].tolist()[0]
