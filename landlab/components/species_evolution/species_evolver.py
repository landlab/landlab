#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Evolve life in a landscape.

Life evolves alongside landscapes by biotic and abiotic processes under complex
dynamics at Earth's surface. Researchers who wish to explore these dynamics can
use this component as a tool for them to build landscape-life evolution models.
Landlab components, including SpeciesEvolver are designed to work with a shared
model grid. Researchers can build novel models using plug-and-play surface
process components to evolve the grid's landscape alongside the life tracked by
SpeciesEvolver. The simulated life evolves following customizable processes.

Component written by Nathan Lyons beginning August 2017.
"""
from collections import OrderedDict

import numpy as np
from pandas import DataFrame, isnull

from landlab import Component

from .record import Record


class SpeciesEvolver(Component):
    """Evolve life in a landscape.

    This component tracks ``Taxon`` objects as they evolve in a landscape. The
    component calls the evolutionary process methods of tracked ``Taxon``
    objects. ``Taxon`` are intended to be subclassed for unique behavior,
    attributes, and model approaches, including different implementations of
    evolutionary processes.

    The general workflow to use this component in a model is

    1. Instantiate the component.
    2. Instantiate taxa.
    3. Introduce taxa to SpeciesEvolver using the ``track_taxon`` method.
    4. Advance the component instance in time using ``run_one_step`` method.

    Taxa can be introduced at model onset and later time steps. Multiple types
    can be tracked by the same SpeciesEvolver instance.

    The taxon type, ``ZoneTaxon`` is distributed with SpeciesEvolver. The
    spatial aspect of ``ZoneTaxon`` macroevolutionary processes is determined
    using ``Zone`` objects. A ``ZoneController`` is used to create and manage
    zones as well as efficiently create multiple ZoneTaxon objects. See the
    documentation of ``ZoneController`` and ``ZoneTaxon`` for more information.
    SpeciesEvolver knows nothing about zones and their controller, meaning the
    concept of zones are not required for other taxon types.

    Model time and other variables can be viewed with the class attribute,
    ``record_data_frame``. Time is recorded to track the history of taxa
    lineages. The unit of time is not considered within the component other
    than the record, and can be thought of as in years or whatever unit is
    needed. Time is advanced with the ``dt`` parameter of the ``run_one_step``
    method.

    The geographic ranges of the taxa at the current model time are evaluated
    during the ``run_one_step`` method. Each taxon object determines if it
    persists or becomes extinct, and if it creates child ``Taxon`` objects.
    Taxa metadata can be viewed with the attribute, ``taxa_data_frame``.

    Taxa are automatically assigned unique identifiers, ``uid``. Identifiers
    are used to reference and retrieve taxon objects. Identifiers are assigned
    in the order taxa are introduced to SpeciesEvolver.

    Examples
    --------
    The evolution of a lowland taxa lineage in response to mountain range
    formation is simulated using ZoneTaxon managed by ZoneController. Mountain
    range formation is forced without processes for simplicity in this example.

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

    Instantiate the component with the grid as the first parameter.

    >>> se = SpeciesEvolver(mg)

    ZoneController requires a function that returns a mask of the total extent
    of taxa habitat. The mask is a boolean array where `True` values represent
    nodes that satisfy habitat conditions. Zone objects are not created here.
    The mask only maps the extent where taxa can exist. This function returns
    `True` where elevation is below 100, which is where the simulated lowland
    taxa of this model can inhabit.

    >>> def zone_func(grid):
    ...     return grid.at_node['topographic__elevation'] < 100

    Instantiate ZoneController with the grid and zone function. The initial
    zones are created at controller instantiation. In this example, one zone is
    created because all nodes of the zone mask are adjacent to each other.

    >>> zc = ZoneController(mg, zone_func)
    >>> len(zc.zones) == 1
    True

    Additional examples of controller usage are provided in ``ZoneController``
    documentation.

    The ``mask`` of the zone is True where the conditions of the zone function
    are met. All nodes of the grid are included because the elevation of each
    node is below 100. The ``zone`` attribute of ``ZoneController`` returns a
    list of the zones that currently exist in the model. Below we return the
    mask of the single zone by indexing this list.

    >>> zc.zones[0].mask
    array([ True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True,  True,  True,  True,  True,  True,  True,
            True,  True,  True], dtype=bool)

    Populate a taxon to the zone. The attribute, ``taxa_data_frame`` indicates
    only the one taxon exists because we populated each zone with one taxon,
    and only the one zone exists.

    >>> taxon = zc.populate_zones_uniformly(1)
    >>> se.track_taxa(taxon)
    >>> se.taxa_data_frame
         appeared  latest_time  extant
    uid
    0           0            0    True

    Force a change in the zone mask to demonstrate component functionality.
    Here we begin a new time step where topography is uplifted by 200 that forms
    a ridge trending north-south in the center of the grid.

    >>> z[[3, 10, 17]] = 200
    >>> z.reshape(mg.shape)
    array([[   1.,    1.,    1.,  200.,    1.,    1.,    1.],
           [   1.,    1.,    1.,  200.,    1.,    1.,    1.],
           [   1.,    1.,    1.,  200.,    1.,    1.,    1.]])

    The current elevation, the elevation following uplift, is represented here.
    ::

        - - - ^ - - -       elevation:  - 1
        - - - ^ - - -                   ^ 200
        - - - ^ - - -

    The updated zone mask is below.
    ::

        . . . x . . .       key:    . node in zone mask
        . . . x . . .               x node outside of zone mask
        . . . x . . .

    Run a step of both the ZoneController and SpeciesEvolver. Both are run to
    keep time in sync between the ``ZoneController``and ``SpeciesEvolver``
    instances.
    >>> delta_time = 1000
    >>> zc.run_one_step(delta_time)
    >>> se.run_one_step(delta_time)

    Two zones exist following this time step.

    >>> len(zc.zones) == 2
    True

    An additional zone was created because the zone mask was not continuous.
    ::

        . . . ^ * * *       key:    . a zone
        . . . ^ * * *               * another zone
        . . . ^ * * *               ^ mountain range

    The split of the initial zone triggered speciation. Taxon 0 became extinct
    as it speciated to child taxa 1 and 2.

    >>> se.taxa_data_frame
         appeared  latest_time  extant
    uid
    0           0         1000   False
    1        1000         1000    True
    2        1000         1000    True

    The phylogenetic tree of the simulated taxa is represented below. The
    number at the line tips are the taxa identifiers.
    ::

                   ┌─ 1
          0 ───────┤
                   └─ 2
            _________
            0    1000
              time

    The split of the initial zone into two zones at time 1000 triggered taxon 0
    to evolve into two child taxon objects. Taxon 1 occupies a zone on one side
    of the mountain range, and taxon 2 occupies a zone on the other side. This
    outcome is the result of the evolutionary processes programmed within
    ``ZoneTaxon`` as well as the parameters used in this example (default
    values were used because optional parameters were not set). Different
    behavior can be achieved by subclassing ``ZoneTaxon`` or ``Taxon``.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Lyons, N.J., Albert, J.S., Gasparini, N.M. (in review). SpeciesEvolver: A
    Landlab component to evolve life in simulated landscapes. Journal of Open
    Source Software.

    **Additional References**

    Albert, J.S., Schoolmaster Jr, D.R., Tagliacollo, V., Duke-Sylvester, S.M.
    (2016). Barrier displacement on a neutral landscape: Toward a theory of
    continental biogeography. Systematic Biology 66(2), 167–182.

    Lyons, N.J., Val, P., Albert, J.S., Willenbring, J.K., Gasparini, N.M., in
    review. Topographic controls on divide migration, stream capture, and
    diversification in riverine life. Earth Surface Dynamics.

    """

    _name = "SpeciesEvolver"

    _info = {
        "taxa__richness": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The number of taxa at each node",
        }
    }

    _cite_as = """@article{lyons2020species,
        author = {Lyons, N.J. and Albert, J.S. and Gasparini, N.M.},
        title = {SpeciesEvolver: A Landlab component to evolve life in simulated landscapes},
        journal = {Journal of Open Source Software},
        }"""

    def __init__(self, grid, initial_time=0):
        """Instantiate SpeciesEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        initial_time : float, int, optional
            The initial time. The unit of time is not considered within the
            component, except in the record. The default is 0.
        """
        super(SpeciesEvolver, self).__init__(grid)

        # Create data structures.

        self._record = Record(initial_time)
        self._record.set_value("taxa", 0)

        self._taxa = OrderedDict(
            [
                ("uid", []),
                ("appeared", []),
                ("latest_time", []),
                ("extant", []),
                ("object", []),
            ]
        )

        # Create a taxa richness field.

        _ = grid.add_zeros("taxa__richness", at="node", dtype=int, clobber=True)

    @property
    def record_data_frame(self):
        """A Pandas DataFrame of SpeciesEvolver variables over time.

        Each row is data of a model time step. The time of the step is recorded
        in the `time` column. `taxa` is the count of taxa extant at a time.
        Additional columns can be added and updated by SpeciesEvolver objects
        during the component ``run_one_step`` method. For example, ``ZoneTaxon``
        add the columns, 'speciations', 'extinctions', and 'pseudoextinctions'
        that are the counts of these variables for this type. Other types may
        also increment these values.

        The DataFrame is created from a dictionary associated with a
        SpeciesEvolver ``Record`` object. nan values in Pandas DataFrame force
        the column to become float values even when data are integers. The
        original value type is retained in the ``Record`` object.
        """
        return self._record.data_frame

    @property
    def taxa_data_frame(self):
        """A Pandas DataFrame of taxa metadata.

        Each row is the metadata of a taxon. The column, 'uid' is the unique
        identifier. The column, `appeared` is the first model time that the
        taxon was tracked by the component. The column, `latest_time` is the
        latest model time the taxon was extant. The column, `extant` indicates
        if the taxon is extant or not (extinct) at the latest (current) time in
        the model.

        The DataFrame is created from a data structure within the component.
        """
        cols = list(self._taxa.keys())
        cols.remove("uid")
        cols.remove("object")
        df = DataFrame(self._taxa, columns=cols, index=self._taxa["uid"])
        df.index.name = "uid"

        return df

    def run_one_step(self, dt):
        """Update the taxa for a single time step.

        This method advances the model time in the component record, calls the
        evolve functions of taxa extant at the current time, and updates the
        variables in the record.

        Parameters
        ----------
        dt : float
            The model time step duration. Time in the record is advanced by
            the value of this parameter.
        """
        self._record.advance_time(dt)

        # Create a dictionary of the taxa to update, `taxa_evolving`.

        extant_taxa = np.array(self._taxa["object"])[self._taxa["extant"]]
        taxa_evolving = OrderedDict.fromkeys(extant_taxa, True)

        # Run stages of taxa evolution.

        stage = 0

        while any(taxa_evolving.values()):
            # Run stage for the evolving taxa.

            evolving_list = filter(taxa_evolving.get, taxa_evolving)
            stage_children = []

            for taxon in evolving_list:
                taxon_is_evolving = taxon._evolve(dt, stage, self._record)

                if len(taxon.children) > 0:
                    stage_children.extend(taxon.children)

                taxa_evolving[taxon] = taxon_is_evolving and taxon.extant

            children_dict = OrderedDict(
                zip(stage_children, [c.extant for c in stage_children])
            )
            children_dict.update(taxa_evolving)
            taxa_evolving = children_dict

            stage += 1

        self._update_taxa_data(taxa_evolving.keys())

    def track_taxa(self, taxa):
        """Add taxa to be tracked over time by SpeciesEvolver.

        The taxon/taxa are introduced at the latest time in the record and
        also tracked during following model times. Each taxon is assigned an
        identifier and then can be viewed in ``taxa_data_frame``.

        Parameters
        ----------
        taxa : Taxon or list of Taxon
            The taxa to introduce.

        Examples
        --------
        ZoneTaxon are used to demonstrate this method.

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

        Track the taxon of the one zone.

        >>> taxon = zc.populate_zones_uniformly(1)
        >>> se.track_taxa(taxon)

        The one taxon is now tracked by SpeciesEvolver as indicated by the taxa
        DataFrame.

        >>> se.taxa_data_frame
             appeared  latest_time  extant
        uid
        0           0            0    True
        """
        if not isinstance(taxa, list):
            taxa = [taxa]

        self._update_taxa_data(taxa)

    def _update_taxa_data(self, taxa_at_time):
        """Update the taxa data structure, set identifiers, and taxa statistics.

        This method sets identifiers and metadata for the newly introduced
        taxa. For the previously introduced, this method updates the
        'latest_time` value of the taxa metadata.

        Parameters
        ----------
        taxa_at_time : list of Taxon
            The taxa at the current model time.
        """
        time = self._record.latest_time

        t_recorded = self._taxa["object"]

        # Update previously introduced taxa.

        t_introduced = [taxon for taxon in taxa_at_time if taxon in t_recorded]

        for taxon in t_introduced:
            # Update taxon data.

            idx = self._taxa["object"].index(taxon)
            self._taxa["latest_time"][idx] = time
            self._taxa["extant"][idx] = taxon.extant

        # Set the data of new taxa.

        t_new = [taxon for taxon in taxa_at_time if taxon not in t_recorded]

        for taxon in t_new:
            # Set identifiers.

            if self._taxa["uid"]:
                taxon._uid = max(self._taxa["uid"]) + 1
            else:
                taxon._uid = 0

            # Append taxon data.

            self._taxa["uid"].append(taxon.uid)
            self._taxa["appeared"].append(time)
            self._taxa["latest_time"].append(time)
            self._taxa["extant"].append(taxon.extant)
            self._taxa["object"].append(taxon)

        # Update taxa stats.

        self._record.set_value("taxa", sum(self._taxa["extant"]))

        self._grid.at_node["taxa__richness"] = self._get_taxa_richness_map()

    def get_taxon_objects(
        self, uid=np.nan, time=np.nan, extant_at_latest_time=np.nan, ancestor=np.nan
    ):
        """Get taxon objects filtered by parameters.

        This method returns all taxon objects tracked by the component when no
        optional parameters are included. The objects returned can be limited
        using one or more parameters. An exception is raised when ``time`` is
        not in the component record. An empty list is returned if no taxon
        match the parameters.

        Parameters
        ----------
        uid : int, optional
            The taxon with this identifier will be returned. An object is
            returned in a list if a taxon with this identifier exists. By
            default, taxa with any identifier can be returned.
        time : float, int, optional
            Limit the taxa returned to those that existed at this time as
            listed in ``taxa_data_frame``. By default, taxa at all of the times
            listed in the component record can be returned.
        extant_at_latest_time : boolean, optional
            Limit the taxa returned to the extant state of taxa at the latest
            time in the record. By default, taxa with any extant state can be
            returned.
        ancestor : int, optional
            Limit the taxa returned to those descending from the taxon
            designated as the ancestor. The ancestor is designated using its
            ``uid``. By default, taxa with any or no ancestors are returned.

        Returns
        -------
        taxa : a list of Taxon
            The Taxon objects that pass through the filter. The list is sorted
            by ``uid``. An empty list is returned if no taxa pass through the
            filter.

        Examples
        --------
        ZoneTaxon are used to demonstrate this method.

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

        Introduce two taxa to the one zone.

        >>> taxa = zc.populate_zones_uniformly(2)
        >>> se.track_taxa(taxa)

        Increment the model time, force mountain range formation to demonstrate
        this method, and then advance model time again.

        >>> z[[3, 10, 17]] = 200
        >>> zc.run_one_step(1000)
        >>> se.run_one_step(1000)

        Remove the westernmost zone. Advance model time once again.

        >>> z[mg.node_x < 3000] = 200
        >>> zc.run_one_step(1000)
        >>> se.run_one_step(1000)

        Display taxa metadata.

        >>> se.taxa_data_frame
             appeared  latest_time  extant
        uid
        0           0         1000   False
        1           0         1000   False
        2        1000         2000   False
        3        1000         2000    True
        4        1000         2000   False
        5        1000         2000    True

        All taxon objects are returned when no parameters are inputted.

        >>> taxa = se.get_taxon_objects()
        >>> [t.uid for t in taxa]
        [0, 1, 2, 3, 4, 5]

        Get the taxon with the identifier, 4.

        >>> se.get_taxon_objects(uid=4)
        [<ZoneTaxon, uid=4>]

        The species extant at a time are returned when ``time`` is specified.
        Here we get the species extant at time 0. Species 0 and 1 are the
        species extant at this time.

        >>> se.get_taxon_objects(time=0)
        [<ZoneTaxon, uid=0>, <ZoneTaxon, uid=1>]

        Get the taxa that descended from species 0.

        >>> se.get_taxon_objects(ancestor=0)
        [<ZoneTaxon, uid=2>, <ZoneTaxon, uid=3>]

        Get the taxa that descended from species 0 and that are currently
        extant. The latest time of taxa 2 and 3 is equal to the latest time in
        the record, altough only one of these taxa are extant.

        >>> se.get_taxon_objects(extant_at_latest_time=True, ancestor=0)
        [<ZoneTaxon, uid=3>]

        Note that combining `extant_at_latest_time` with `time` does not return
        the taxa extant at the inputted time. Rather the taxa present at `time`
        that are still extant are returned. An empty list is returned if no
        taxa match the criteria.

        >>> se.get_taxon_objects(extant_at_latest_time=True, time=0)
        []

        An empty list is also return when no taxa match a valid value of
        ``uid.``

        >>> se.get_taxon_objects(uid=11)
        []
        """
        # Handle ancestor.

        if isnull(ancestor):
            taxa = self._taxa
        elif ancestor in self._taxa["uid"]:
            idx_number = self._taxa["uid"].index(ancestor)
            taxon = self._taxa["object"][idx_number]

            descendants = []
            stack = [taxon]

            while stack:
                if taxon.children:
                    descendants.extend(taxon.children)
                    stack.extend(taxon.children)

                stack.remove(taxon)

                if stack:
                    taxon = stack[0]

            descendants = list(set(descendants))
            taxa = self._subset_taxa_data_structure(descendants)
        else:
            return []

        # Handle identifier.

        if isnull(uid):
            idx_id = np.ones(len(taxa["uid"]), dtype=bool)
        else:
            idx_id = np.array(taxa["uid"]) == uid

        # Handle time.

        if isnull(time):
            idx_time = np.ones(len(taxa["uid"]), dtype=bool)
        else:
            idx_time = self._mask_taxa_by_time(taxa, time)

        # Handle extant state.

        if isnull(extant_at_latest_time):
            idx_ext = np.ones(len(taxa["uid"]), dtype=bool)
        else:
            idx_ext = np.array(taxa["extant"]) == extant_at_latest_time

        # Get the Taxon list.

        idx = np.all([idx_time, idx_id, idx_ext], 0)
        taxa = np.array(taxa["object"])[idx].tolist()
        taxa.sort(key=lambda taxon: taxon.uid)

        return taxa

    def _subset_taxa_data_structure(self, subset):
        all_objects = self._taxa["object"]
        idx = [i for i, e in enumerate(all_objects) if e in subset]

        structure = OrderedDict()
        for k, v in self._taxa.items():
            structure[k] = np.array(v)[idx].tolist()

        return structure

    def _mask_taxa_by_time(self, taxa, time):
        """Get a mask of ``taxa`` that exist at ``time``.

        Parameters
        ----------
        taxa : a list of Taxon
            The taxa to query.
        time : float, int
            The model time of ``taxa`` existence to mask.

        Returns
        -------
        ndarray
            A mask of ``taxa`` extant at ``time``.
        """
        if time not in self._record.times:
            raise ValueError("The time, {} is not in the record.".format(time))

        # Create a mask of taxa within time bounds.

        t_appeared = np.array(taxa["appeared"])
        t_latest = np.array(taxa["latest_time"])

        t_prior = t_appeared <= time
        t_later = t_latest >= time
        taxa_mask = np.all([t_prior, t_later], 0)

        return taxa_mask

    def _get_taxa_richness_map(self):
        """Get a map of the number of taxa."""
        taxa = self.get_taxon_objects(extant_at_latest_time=True)

        if taxa:
            masks = np.stack([taxon.range_mask for taxon in taxa])
            richness_mask = sum(masks).astype(int)
        else:
            richness_mask = np.zeros(self._grid.number_of_nodes, dtype=int)

        return richness_mask
