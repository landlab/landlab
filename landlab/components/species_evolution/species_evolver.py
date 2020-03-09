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
from pandas import DataFrame

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

    Taxa are automatically assigned unique taxon identifiers, ``tid``.
    Identifiers are used to reference and retrieve taxon objects. Identifiers
    are assigned in the order taxa are introduced to SpeciesEvolver.

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
          pid       type  t_first  t_final
    tid
    0    <NA>  ZoneTaxon        0     <NA>

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
          pid       type  t_first  t_final
    tid
    0    <NA>  ZoneTaxon        0     <NA>
    1       0  ZoneTaxon     1000     <NA>

    The phylogenetic tree of the simulated taxa is represented below. The
    number at the line tips are the taxa identifiers.
    ::

          0 ──────┬── 0
                  │
                  └── 1
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

    Lyons, N.J., Albert, J.S., Gasparini, N.M. (2020). SpeciesEvolver: A
    Landlab component to evolve life in simulated landscapes. Journal of Open
    Source Software 5(46), 2066, https://doi.org/10.21105/joss.02066

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
        year = {2020},
        journal = {Journal of Open Source Software},
        volume = {5},
        number = {46},
        doi = {10.21105/joss.02066},
        url = {https://doi.org/10.21105/joss.02066}
        }"""

    def __init__(self, grid, initial_time=0):
        """Instantiate SpeciesEvolver.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        initial_time : float, int, optional
            The initial time. The unit of time is not considered within the
            component, with the exception that time is logged in the record.
            The default value of this parameter is 0.
        """
        super(SpeciesEvolver, self).__init__(grid)

        # Create data structures.

        self._record = Record(initial_time)
        self._record.set_value("taxa", 0)

        self._taxa_data = OrderedDict(
            [
                ("tid", []),
                ("pid", []),
                ("type", []),
                ("t_first", []),
                ("t_final", [])
            ]
        )

        self._taxon_objs = []

        # Create a taxa richness field.

        _ = grid.add_zeros("taxa__richness", at="node", dtype=int, clobber=True)

    @property
    def record_data_frame(self):
        """A Pandas DataFrame of SpeciesEvolver variables over time.

        Each row is data of a model time step. The time of the step is recorded
        in the `time` column. `taxa` is the count of taxa extant at a time.
        Additional columns can be added and updated by SpeciesEvolver objects
        during the component ``run_one_step`` method. See documention of Taxon
        objects for an explanation of these columns.

        The DataFrame is created from a dictionary associated with a
        SpeciesEvolver ``Record`` object. nan values in Pandas DataFrame force
        the column to become float values even when data are integers. The
        original value type is retained in the ``Record`` object.
        """
        return self._record.data_frame

    @property
    def taxa_data_frame(self):
        """A Pandas DataFrame of taxa metadata.

        Each row is the metadata of a taxon. The column, 'tid' is the taxon
        identifier. The column, `appeared` is the first model time that the
        taxon was tracked by the component. The column, `latest_time` is the
        latest model time the taxon was extant. The column, `extant` indicates
        if the taxon is extant or not (extinct) at the latest (current) time in
        the model.

        The DataFrame is created from a data structure within the component.
        """
        data = self._taxa_data
        cols = list(data.keys())
        cols.remove("tid")
        df = DataFrame(data, columns=cols, index=data["tid"])
        df.index.name = "tid"

        # Change column number type because pandas makes a column float if it
        # has any nan values.
        df["pid"] = df["pid"].astype("Int64")
        if all(isinstance(item, int) for item in data['t_final'] if not np.isnan(item)):
            df["t_final"] = df["t_final"].astype("Int64")

        return df

    def run_one_step(self, dt):
        """Update the taxa for a single time step.

        This method advances the model time in the component record, calls the
        evolve method of taxa extant at the current time, and updates the
        variables in the record and taxa dataframes.

        Parameters
        ----------
        dt : float
            The model time step duration. Time in the record is advanced by the
            value of this parameter.
        """
        record = self._record
        record.advance_time(dt)

        # Create a dictionary of the taxa to update at the current model time.
        # Keys are objects of extant taxa. Values are booleans indicating if
        # stages remain for respective taxon.

        time_dict = OrderedDict.fromkeys(self._taxon_objs, True)

        # Iteratively call taxa ``_evolve`` method until all stages of all taxa
        # have run.

        stage = 0

        while any(time_dict.values()):
            # Run evolution stage.

            stage_dict = OrderedDict([])
            evolving_taxa = filter(time_dict.get, time_dict)

            for taxon in evolving_taxa:
                # Run evolution stage of taxon with remaining stages.
                stages_remain, taxon_children = taxon._evolve(dt, stage, record)

                if taxon_children:
                    stage_dict.update(
                        OrderedDict.fromkeys(taxon_children, stages_remain)
                    )

                stage_dict[taxon] = stages_remain and taxon.extant

            time_dict.update(stage_dict)

            stage += 1

        self._update_taxa_data(time_dict.keys())

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
              pid       type  t_first  t_final
        tid
        0    <NA>  ZoneTaxon        0     <NA>
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
        data = self._taxa_data

        t_recorded = self._taxon_objs
        t_introduced = [taxon for taxon in taxa_at_time if taxon in t_recorded]
        t_new = [taxon for taxon in taxa_at_time if taxon not in t_recorded]

        # Update previously introduced taxa.

        for taxon in t_introduced:
            # Update taxon data.

            if not taxon.extant:
                idx = data["tid"].index(taxon.tid)
                data["t_final"][idx] = time
                self._taxon_objs.remove(taxon)

        # Set the data of new taxa.

        for taxon in t_new:
            # Set identifiers.

            if data["tid"]:
                taxon._tid = max(data["tid"]) + 1
            else:
                taxon._tid = 0

            # Append taxon data.

            data["tid"].append(taxon.tid)
            if taxon.parent != None:
                data["pid"].append(taxon.parent.tid)
            else:
                data["pid"].append(np.nan)
            data["type"].append(type(taxon).__name__)
            data["t_first"].append(time)
            if taxon.extant:
                data["t_final"].append(np.nan)
                self._taxon_objs.append(taxon)
            else:
                data["t_final"].append(time)

        # Update taxa stats.

        self._record.set_value("taxa", len(self._taxon_objs))

        self._grid.at_node["taxa__richness"] = self._get_taxa_richness_map()


    def get_taxon_objects(
        self, tids=np.nan, min_time=np.nan, max_time=np.nan, ancestor=np.nan
    ):
        """Get extant taxon objects filtered by parameters.

        This method returns all taxon objects tracked by the component when no
        optional parameters are included. The objects returned can be limited
        using one or more parameters.

        Parameters
        ----------
        tids : list of int, optional
            The taxa with these identifiers will be returned. A list is
            returned even if only one object is contained within the list. By
            default, when `tids` is not specified, extant taxa with any
            identifier can be returned.
        min_time : float, int, optional
            Limit the taxa returned to those extant after this time. By
            default, extant taxa at all of the times listed in the component
            record can be returned. Must be less than ``max_time`` if both
            specified.
        max_time : float, int, optional
            Limit the taxa returned to those extant prior to this time. By
            default, extant taxa at all of the times listed in the component
            record can be returned. Must be greater than ``min_time`` if both
            specified.
        ancestor : int, optional
            Limit the taxa returned to those descending from the taxon
            designated as the ancestor. The ancestor is designated using its
            ``tid``. By default, taxa with any or no ancestors are returned.

        Returns
        -------
        taxa : a list of Taxon
            The Taxon objects that pass through the filter. The list is sorted
            by ``tid``. An empty list is returned if no taxa pass through the
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

        Force mountain range formation to demonstrate this method, and then
        advance model time.

        >>> z[mg.x_of_node == 5000] = 200
        >>> zc.run_one_step(1000)
        >>> se.run_one_step(1000)

        Form a second mountain range, and advance the model time again.

        >>> z[mg.x_of_node == 3000] = 200
        >>> zc.run_one_step(1000)
        >>> se.run_one_step(1000)

        Uplift the east, forming a high elevation plataue, removing the
        easternmost zone. Advance model time once again.

        >>> z[mg.x_of_node < 3000] = 200
        >>> zc.run_one_step(1000)
        >>> se.run_one_step(1000)

        Display taxa metadata.

        >>> se.taxa_data_frame
              pid       type  t_first  t_final
        tid
        0    <NA>  ZoneTaxon        0     <NA>
        1    <NA>  ZoneTaxon        0     <NA>
        2       0  ZoneTaxon     1000     3000
        3       1  ZoneTaxon     1000     3000
        4       0  ZoneTaxon     2000     <NA>
        5       1  ZoneTaxon     2000     <NA>

        Objects of all extant taxon are returned when no parameters are
        inputted.

        >>> se.get_taxon_objects()
        [<ZoneTaxon, tid=0>,
         <ZoneTaxon, tid=1>,
         <ZoneTaxon, tid=4>,
         <ZoneTaxon, tid=5>]

        Get the taxon objects with identifiers, 0 and 1.

        >>> se.get_taxon_objects(tids=4, 5])
        [<ZoneTaxon, tid=0>, <ZoneTaxon, tid=1>]

        The species extant at a time are returned when ``time`` is specified.
        Here we get the species extant at time 0. Species 0 and 1 are the
        species extant at this time.

        >>> se.get_taxon_objects(time=1000)
        [<ZoneTaxon, tid=0>, <ZoneTaxon, tid=1>]

        Get the taxa that descended from species 0.

        >>> se.get_taxon_objects(ancestor=0)
        [<ZoneTaxon, tid=2>, <ZoneTaxon, tid=3>]

        Get the taxa that descended from species 0 and that are currently
        extant. The latest time of taxa 2 and 3 is equal to the latest time in
        the record, altough only one of these taxa are extant.

        >>> se.get_taxon_objects(time=1000, ancestor=0)
        [<ZoneTaxon, tid=3>]

        Note that combining `extant_at_latest_time` with `time` does not return
        the taxa extant at the inputted time. Rather the taxa present at `time`
        that are still extant are returned. An empty list is returned if no
        taxa match the criteria.

        >>> se.get_taxon_objects(extant_at_latest_time=True, time=0)
        []

        An empty list is also return when no taxa match a valid value of
        ``tid.``

        >>> se.get_taxon_objects(tid=11)
        []
        """
        # Create `results` where each element is a list tids output by a
        # parameter query of taxa.

        extant_tids = [taxon.tid for taxon in self._taxon_objs]
        results = [extant_tids]

        data = self._taxa_data

        # Query by identifiers.

        if isinstance(tids, list):
            results.append(tids)

        # Query by time.

        if not np.isnan(min_time) or not np.isnan(max_time):
            if min_time > max_time:
                raise ValueError("`min_time` cannot be greater than `max_time`")

            if not np.isnan(min_time):
                earliest_time = min_time
            else:
                earliest_time = self._record.earliest_time

            if not np.isnan(max_time):
                latest_time = max_time
            else:
                latest_time = self._record.latest_time

            t_first = np.array(data["t_first"])
            t_latest = np.nan_to_num(data["t_final"], nan=self._record.latest_time)

            mask = np.all([t_first <= latest_time, t_latest >= earliest_time], 0)
            results.append(np.array(data["tid"])[mask].tolist())

        # Query by ancestor.

        if not np.isnan(ancestor) and ancestor in data["tid"]:
            df = self.taxa_data_frame
            df['pid'] = df['pid'].fillna(-1)

            taxon = ancestor

            descendants = []
            stack = [taxon]

            while stack:
                children = df.index[df['pid'] == taxon].tolist()

                if children:
                    descendants.extend(children)
                    stack.extend(children)

                stack.remove(taxon)

                if stack:
                    taxon = stack[0]

            descendants = list(set(descendants))
            results.append(descendants)

        # Get the Taxon objects that match all parameter query results.

        matched_tids = list(set.intersection(*map(set, results)))
        taxa = [taxon for taxon in self._taxon_objs if taxon.tid in matched_tids]
        taxa.sort(key=lambda taxon: taxon.tid)

        return taxa

    def _get_taxa_richness_map(self):
        """Get a map of the number of taxa."""
        taxa = self._taxon_objs

        if taxa:
            masks = np.stack([taxon.range_mask for taxon in taxa])
            richness_mask = masks.sum(axis=0).astype(int)
        else:
            richness_mask = np.zeros(self._grid.number_of_nodes, dtype=int)

        return richness_mask
