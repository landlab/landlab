#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Evolve species in a landscape.

Component written by Nathan Lyons beginning August 2017.
"""
from collections import OrderedDict
from itertools import count, product
from string import ascii_uppercase

import numpy as np
from pandas import DataFrame, isnull

from landlab import Component
from landlab.components.species_evolution.record import Record


class SpeciesEvolver(Component):
    """Evolve species in a landscape.

    This component along with ``Species`` objects can be used in models to
    develop organismal lineages as they develop in response to landscape
    change. The history of lineages (phylogeny) are tracked through species.
    Evolution processes are conducted by Species objects. The species type,
    ``ZoneSpecies`` is distributed with SpeciesEvolver. Species are intended
    to be subclassed for unique species behavior, attributes, and model
    approaches. The component stores species objects and calls the species
    evolutionary process methods.

    The general workflow to use this component in a model is

    1. Instantiate the component.
    2. Instantiate species.
    3. Introduce species using the SpeciesEvolver ``introduce_species`` method.
    4. Advance the model in time using the ``run_one_step`` method.

    Species can be introduced at model onset and later time steps. Species can
    be ``ZoneSpecies`` and other ``Species`` subclass types, and multiple types
    may be introduced to the same SpeciesEvolver instance.

    ``ZoneSpecies`` evolves through the processes of dispersal, speciation,
    and extinction. The spatial aspect of these processes is determined using
    ``Zone`` objects. A ``ZoneController`` is used to create and manage
    zones as well as efficiently create multiple ZoneSpecies. See the
    documentation of ``ZoneController`` for more information. SpeciesEvolver
    knows nothing about zones and their controller, meaning the concept of
    zones are not required for other species types.

    Model time and other variables can be accessed with the class attribute,
    ``record_data_frame``. Rime is tracked to construct phylogeny. The unit of
    time is not considered within the component and can be thought of as in
    years or whatever unit is needed. Time is advanced using the `dt` parameter
    of the ``run_one_step`` method.

    The species extant at the current model time are evolved during the
    ``run_one_step`` method. The persistence or extinction of a species is
    determined by the species. Creation of child species is also determined by
    species. Species metadata can be accessed with the attribute,
    ``species_data_frame``.

    Species are automatically assigned identifiers. The identifier of a species
    is a two element tuple. The first element is the clade id. Clades are
    lettered from A to Z then AA to AZ and so forth as more clades are created.
    The second identifier element designates the clade members numbered in the
    order of appearance. For example, the first species introduced is A.0 and
    if that species speciates, the first child species is A.1.

    The development of this component was inspired by SEAMLESS (Spatially
    Explicit Area Model of Landscape Evolution by SimulationS). See Albert et
    al., 2017, Systematic Biology 66.

    Examples
    --------
    The evolution of a lowland species lineage in response to mountain range
    formation is simulated using ZoneSpecies managed by ZoneController.
    Mountain range formation is forced without processes for simplicity in this
    example.

    Import modules used in the following examples.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import SpeciesEvolver
    >>> from landlab.components.species_evolution import ZoneController

    Create a model grid with mountain scale resolution. The elevation is
    equally low throughout the grid at model onset.

    >>> mg = RasterModelGrid((3, 7), 1000)
    >>> z = mg.add_ones('topographic__elevation', at='node')
    >>> z.reshape(mg.shape)
    array([[ 1.,  1.,  1.,  1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  1.]])

    Instantiate the component with the grid.

    >>> se = SpeciesEvolver(mg)

    ZoneController requires a function that returns a mask of the total
    extent of species habitat. The mask is a boolean array where `True` values
    represent nodes that satisfy habitat conditions. Zone objects of species
    are not created here. The mask only masks the extent where species can
    exist. This function returns `True` where elevation is below 100, which is
    where the simulated lowland species of this model can inhabit.

    >>> def zone_func(grid):
    ...     return grid.at_node['topographic__elevation'] < 100

    Instantiate ZoneController with the grid and zone function. The
    initial zones are created at instantiation. In this example, one zone is
    created because all nodes of the zone mask are adjacent to each other.

    >>> zc = ZoneController(mg, zone_func)
    >>> len(zc.zones) == 1
    True

    Additional examples of controller usage are provided in ZoneController.

    The ``mask`` of the zone is True where the conditions of the zone function
    are met. All nodes of the grid are included because the elevation of each
    node is below 100.

    >>> zc.zones[0].mask
    array([ True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True], dtype=bool)

    Introduce a species to the zone. The attribute, ``species_data_frame``
    indicates only the one species exists because we indicated to populate each
    one with one species, and only one zone exists.

    >>> species = zc.populate_zones_uniformly(1)
    >>> se.introduce_species(species)
    >>> se.species_data_frame
      clade  number  time_appeared  latest_time
    0     A       0              0            0

    Force a change in the zone mask to demonstrate component functionality.
    Here we begin a new time step where topography is uplifted by 200 that
    forms a ridge trending north-south in the center of the grid.

    >>> z[[3, 10, 17]] = 200
    >>> z.reshape(mg.shape)
    array([[   1.,    1.,    1.,  200.,    1.,    1.,    1.],
           [   1.,    1.,    1.,  200.,    1.,    1.,    1.],
           [   1.,    1.,    1.,  200.,    1.,    1.,    1.]])

    The current elevation, the elevation following uplift, is represented here.

    - - - ^ - - -       elevation:  - 1
    - - - ^ - - -                   ^ 200
    - - - ^ - - -

    The updated zone mask is below.

    . . . x . . .       key:    . node in zone mask
    . . . x . . .               x node outside of zone mask
    . . . x . . .

    Run a step of both the ZoneController and SpeciesEvolver. Two zones exist
    following this time step.

    >>> delta_time = 1000
    >>> zc.run_one_step(delta_time)
    >>> se.run_one_step(delta_time)
    >>> len(zc.zones) == 2
    True

    A additional zone was created because the zone mask was not continuous.

    . . . ^ * * *       key:    . a zone
    . . . ^ * * *               * another zone
    . . . ^ * * *               ^ mountain range

    The split of the initial zone triggered speciation. Species A.0 became
    extinct as it speciated to child species, A.1 and A.2.

    >>> se.species_data_frame
      clade  number  time_appeared  latest_time
    0     A       0              0            0
    1     A       1           1000         1000
    2     A       2           1000         1000

    The phylogeny of the species is represented below.
    ::

                   ┌─ A.1
        A.0 ───────┤
                   └─ A.2
            _________
            0    1000
              time
    """
    _name = 'SpeciesEvolver'

    _info = {
        "species__richness": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The number of species at each node",
        },
    }

    def __init__(self, grid, initial_time=0):
        """Instantiate SpeciesEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        initial_time : float, int, optional
            The initial time. The unit of time is unspecified within the
            component. The default is 0.
        """
        super(SpeciesEvolver, self).__init__(grid)

        # Create data structures.

        self._record = Record(initial_time)
        self._record.set_value('species_count', 0)

        self._species = OrderedDict([('clade', []),
                                     ('number', []),
                                     ('time_appeared', []),
                                     ('latest_time', []),
                                     ('object', [])])

        # Create clade name generator.

        def _get_next_clade_name():
            for size in count(1):
                for s in product(ascii_uppercase, repeat=size):
                    yield ''.join(s)

        self._clade_generator = _get_next_clade_name()

        # Create a species richness field.

        _ = grid.add_zeros(
            'species__richness', at='node', dtype=int, clobber=True
        )

    @property
    def record_data_frame(self):
        """A Pandas DataFrame of SpeciesEvolver variables over time."""
        return self._record.data_frame

    @property
    def species_data_frame(self):
        """A Pandas DataFrame of species metadata."""
        cols = list(self._species.keys())
        cols.remove('object')
        sort_cols = ['clade', 'number']
        return DataFrame(
            self._species,
            columns=cols).sort_values(by=sort_cols).reset_index(drop=True)

    def run_one_step(self, dt):
        """Update the species for a single time step.

        This method advances the model time in the component record, calls the
        evolve functions of species extant at the current time, and updates the
        variables in the record.

        Parameters
        ----------
        dt : float
            The model time step duration. Time in the record is advanced by
            the value of this parameter.
        """
        self._record.advance_time(dt)

        # Process species.

        extant_species = self.filter_species(time=self._record.prior_time)

        for es in extant_species:
            es._evolve_stage_1(dt, self._record)

        survivors = []

        for es in extant_species:
            es_persists, child_species = es._evolve_stage_2(dt, self._record)

            if es_persists:
                survivors.append(es)

            if len(child_species) > 0:
                survivors.extend(child_species)

        self._update_species_data(survivors)

        self._grid.at_node['species__richness'] = (
            self._get_species_richness_map()
        )

    def introduce_species(self, species):
        """Add species to SpeciesEvolver.

        The species are introduced at the latest time in the record. Each
        species is assigned an identifier and is included in
        ``species_data_structure``.

        Parameters
        ----------
        species : Species or list of Species
            The species to introduce.

        Examples
        --------
        ZoneSpecies are used to demonstrate this method.

        Import modules used in the following examples.

        >>> from landlab import RasterModelGrid
        >>> from landlab.components import SpeciesEvolver
        >>> from landlab.components.species_evolution import ZoneController

        Create a model grid with flat topography.

        >>> mg = RasterModelGrid((3, 7), 1000)
        >>> z = mg.add_ones('topographic__elevation', at='node')

        Instantiate SpeciesEvolver and a ZoneController. Instantiate the
        latter with a function that masks the low elevation zone extent. Only
        one zone is created.

        >>> se = SpeciesEvolver(mg)
        >>> def zone_func(grid):
        ...     return grid.at_node['topographic__elevation'] < 100
        >>> zc = ZoneController(mg, zone_func)
        >>> len(zc.zones) == 1
        True

        Introduce a species to the one zone.

        >>> introduced_species = zc.populate_zones_uniformly(1)
        >>> se.introduce_species(introduced_species)
        """
        if not isinstance(species, list):
            species = [species]

        self._update_species_data(species)

    def _update_species_data(self, species_at_time):
        """Update the species data structure and set identifiers, if needed.

        This method sets identifiers and species metadata for the newly
        introduced species. For the previously introduced, this method updates
        the 'latest_time` value of the species metadata.

        Parameters
        ----------
        species_at_time : list of species
            The species at the current model time.
        """
        time = self._record.latest_time

        s_at_time = species_at_time
        s_recorded = self._species['object']

        # Identify species previously introduced.
        s_introduced = [s for s in s_at_time if s in s_recorded]

        # Identify new species.
        s_new = [s for s in s_at_time if s not in s_recorded]

        # Update previously introduced species.
        for s in s_introduced:
            idx = self._species['object'].index(s)
            self._species['latest_time'][idx] = time

        # Set the data of new species.
        if s_new:
            self._set_species_identifiers(s_new)

            clade = [s.identifier[0] for s in s_new]
            s_number = [s.identifier[1] for s in s_new]

            t = [time] * len(s_new)

            self._species['clade'].extend(clade)
            self._species['number'].extend(s_number)
            self._species['time_appeared'].extend(t)
            self._species['latest_time'].extend(t)
            self._species['object'].extend(s_new)

        # Update species count.
        self._record.increment_value('species_count', len(species_at_time))

    def _set_species_identifiers(self, species):
        """Set identifiers of species.

        Parameters
        ----------
        species : a list of Species
            The species to set identifiers.
        """
        # Get existing identifiers.

        clades = np.array(self._species['clade'])
        nums = np.array(self._species['number'])

        for s in species:
            # Set species identifier.

            if s.parent_species is None:
                clade = next(self._clade_generator)
                num = 0
            else:
                clade = s.parent_species.identifier[0]
                num = nums[np.where(clades == clade)[0]].max() + 1

            s._identifier = (clade, num)

            # Update existing identifiers.

            clades = np.append(clades, clade)
            nums = np.append(nums, num)

    def filter_species(
        self, species_subset=np.nan, time=np.nan, clade=np.nan, number=np.nan
    ):
        """Get species objects.

        This method returns all of the species introduced to the component when
        no optional parameters are included. Optionally, the species in
        ``species_subset`` are filtered by the other parameters of this method.

        The parameter, ``time`` filters species by the time the species is
        known to exist as indicated in ``species_data_frame``. An exception is
        raised when ``time`` is not in the component record. The parameter,
        ``clade`` and ``number`` filter species by their identifier elements.
        An empty list is returned if no species have an identifier that matches
        these parameters. Multiple parameters can be combined to limit results.

        Parameters
        ----------
        species_subset : list of Species, optional
            The species to filter. By default, all species introduced to the
            component are included.
        time : float, int, optional
            The model time.  By default, all species extant at all times in the
            component record can be returned.
        clade : string, optional
            Species of this clade will be returned. By default, species of all
            clades can be returned.
        number : int, optional
            Species with this identifier number will be returned. By default,
            species with all numbers can be returned.

        Returns
        -------
        species : a list of Species
            The SpeciesEvolver species that pass through the filter. An
            empty list is returned if no species pass through the filter. The
            list is sorted first by clade, then by species number.

        Examples
        --------
        ZoneSpecies are used to demonstrate this method.

        Import modules used in the following examples.

        >>> from landlab import RasterModelGrid
        >>> from landlab.components import SpeciesEvolver
        >>> from landlab.components.species_evolution import ZoneController

        Create a model grid with flat topography.

        >>> mg = RasterModelGrid((3, 7), 1000)
        >>> z = mg.add_ones('topographic__elevation', at='node')

        Instantiate SpeciesEvolver and a ZoneController. Instantiate the
        latter with a function that masks the low elevation zone extent.
        Only one zone is created.

        >>> se = SpeciesEvolver(mg)
        >>> def zone_func(grid):
        ...     return grid.at_node['topographic__elevation'] < 100
        >>> zc = ZoneController(mg, zone_func)
        >>> len(zc.zones) == 1
        True

        Introduce two species to the one zone.

        >>> introduced_species = zc.populate_zones_uniformly(2)
        >>> se.introduce_species(introduced_species)

        Increment the model time, force mountain range formation to demonstrate
        this method, and then increment model time again.

        >>> zc.run_one_step(1000)
        >>> se.run_one_step(1000)

        >>> z[[3, 10, 17]] = 200

        >>> zc.run_one_step(1000)
        >>> se.run_one_step(1000)

        Display metadata of all the species.

        >>> se.species_data_frame
          clade  number  time_appeared  latest_time
        0     A       0              0         1000
        1     A       1           2000         2000
        2     A       2           2000         2000
        3     B       0              0         1000
        4     B       1           2000         2000
        5     B       2           2000         2000

        All species objects are returned when no parameters are inputted.

        >>> filtered_species = se.filter_species()
        >>> [s.identifier for s in filtered_species]
        [('A', 0), ('A', 1), ('A', 2), ('B', 0), ('B', 1), ('B', 2)]

        The species extant at a time are returned when ``time`` is specified.
        Here we get the species extant at time 0. Species A.0 and B.0 are the
        species extant at this time.

        >>> filtered_species = se.filter_species(time=0)
        >>> [s.identifier for s in filtered_species]
        [('A', 0), ('B', 0)]

        Get the species, B.2. The one species with this identifier is returned.

        >>> filtered_species = se.filter_species(clade='B', number=2)
        >>> [s.identifier for s in filtered_species]
        [('B', 2)]

        Get all of the species of clade, B.

        >>> filtered_species = se.filter_species(clade='B')
        >>> [s.identifier for s in filtered_species]
        [('B', 0), ('B', 1), ('B', 2)]

        Get all of the species with the same number, 0, despite the clade.

        >>> filtered_species = se.filter_species(number=0)
        >>> [s.identifier for s in filtered_species]
        [('A', 0), ('B', 0)]

        An empty list is return when no species match a valid value type of
        ``identifier element.``

        >>> se.filter_species(clade='C')
        []

        Parameters, ``time`` and species identifier elements can be combined.
        >>> filtered_species = se.filter_species(time=2000, clade='B')
        >>> [s.identifier for s in filtered_species]
        [('B', 1), ('B', 2)]

        By default, this method filters all species introduced to
        SpeciesEvolver. The initial collection of filtered objects can be
        limited using ``species_subset``, which makes possible iterative use of
        this method.

        >>> filtered_species_1 = se.filter_species(time=2000)
        >>> [s.identifier for s in filtered_species_1]
        [('A', 1), ('A', 2), ('B', 1), ('B', 2)]

        >>> filtered_species_2 = se.filter_species(number=1)
        >>> [s.identifier for s in filtered_species_2]
        [('A', 1), ('B', 1)]
        """
        # Handle inputted species.

        if type(species_subset) is list:
            all_objects = self._species['object']
            idx = [i for i, e in enumerate(all_objects) if e in species_subset]

            input_species = OrderedDict()
            for k, v in self._species.items():
                input_species[k] = np.array(v)[idx].tolist()
        elif isnull(species_subset):
            input_species = self._species

        # Handle time.

        if isnull(time):
            # Get the indices of all species.
            idx_time = np.ones(len(input_species['number']), dtype=bool)
        else:
            idx_time = self._mask_species_by_time(input_species, time)

        # Handle clade.

        if isnull(clade):
            idx_clade = np.ones(len(input_species['number']), dtype=bool)
        else:
            idx_clade = np.array(input_species['clade']) == clade

        # Handle number.

        if isnull(number):
            idx_number = np.ones(len(input_species['number']), dtype=bool)
        else:
            idx_number = np.array(input_species['number']) == number

        # Get the species objects.

        idx = np.all([idx_time, idx_clade, idx_number], 0)
        species = np.array(input_species['object'])[idx].tolist()

        species.sort(key=lambda s: (s.identifier[0], s.identifier[1]))

        return species

    def _mask_species_by_time(self, species, time):
        """Get a mask of ``species`` that exist at ``time``.

        Parameters
        ----------
        species : a list of Species
            The species to query.
        time : float, int
            The model time of ``species`` existence to mask.

        Returns
        -------
        ndarray
            A mask of ``species`` extant at ``time``.
        """
        if time not in self._record.times:
            raise ValueError('The time, {} is not in the record.'.format(time))

        # Get indices of species within time bounds.

        t_appeared = np.array(species['time_appeared'])
        t_latest = np.array(species['latest_time'])

        t_prior = t_appeared <= time
        t_later = t_latest >= time
        species_mask = np.all([t_prior, t_later], 0)

        return species_mask

    def _get_species_richness_map(self):
        """Get a map of the number of species."""
        species = self.filter_species(time=self._record.latest_time)

        if species:
            masks = np.stack([s.range_mask for s in species])
            richness_mask = sum(masks).astype(int)
        else:
            richness_mask = np.zeros(self._grid.number_of_nodes, dtype=int)

        return richness_mask
