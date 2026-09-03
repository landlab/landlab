#!/usr/bin/env python3
from collections.abc import Iterable
from collections.abc import Mapping
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import ArrayLike
from numpy.typing import NDArray
from requireit import raise_as
from requireit import require_array
from requireit import require_contains
from requireit import require_dtype
from requireit import require_one_of
from requireit import require_shape
from requireit import require_sorted


class DataRecord:
    """Data structure to store variables in time and/or space dimensions.

    This class uses a xarray Dataset to store variables. This datastructure is
    located at the property ``dataset``. The DataRecord expands xarray Dataset
    with additional attributes and functions, including the ability to
    aggregate values on Landlab grid elements.

    DataRecord uses the concept of an "item", a physical thing that is located
    on the grid and has some properties, stored as variables. DataRecord tracks
    variables through time, variables associated with a Landlab grid across
    items, or both.

    Thus data variables can vary along one or both of the following dimensions:

        - time (model time)
        - item_id: variables can characterize a set of items (each identified
          by an individual id) that reside on the grid.

    If an item or set of items is defined, each item must be defined by the
    grid element and the element id at which it resides, e.g.:

        grid_element = 'node'
        element_id = 9.

    When items are defined, each item is given a unique id and the underlying
    Dataset uses a dimension "item_id". **Items are assigned ids beginning with
    0 followed by consecutively increasing integers.**

    Examples:

        - the variable 'mean_elevation' characterizes the grid and varies with
          time,
        - the variable 'clast__rock_type' characterizes a set of items (clasts)
          and varies with item_id,
        - the variable 'clast__size' can vary with both time and item_id

    In the above case, `grid_element` and `element_id` are default data
    variables (in addition to any user-specified variables).

    For each item, `element_id` must be less than the number of this item's
    grid_element that exist on the grid or be one of the dummy element values.
    For example, if the grid has 100 links, and no dummy link values are
    indicated, then, no item can live at link 100 or link -3 because only links
    0 to 99 exist in this example.

    Anything that the DataRecord keeps track of is considered a "record",
    whether it uses one or both of the two standard dimensions (**time** and
    **item_id**).

    DataRecord provides two method to assist with adding new records. The
    method ``add_item`` should be used when no new variables are being added.
    The method ``add_record`` should be used when new variables are being
    added or when a variable is only tracked over the **time** dimension.
    """

    _name = "DataRecord"

    def __init__(
        self,
        grid,
        dummy_elements=None,
        time: ArrayLike | None = None,
        items=None,
        data_vars=None,
        attrs=None,
    ):
        """
        Parameters
        ----------
        grid : ModelGrid
        dummy_elements : dict
            Dictionary indicating valid values for dummy grid elements. For
            example, if you need an "exit" off of a grid with  100 links, you
            could indicate `dummy_elements = {"link": [9999]}`
            to set a link id of 9999 as a dummy link. Multiple dummy elements
            are possible and we recommend using values larger than the number
            of grid elements for the dummy values.
        time : list or 1-D array of float or int (optional)
            The initial time(s) to add to the record. A time dimension is not
            created if the value is 'None' (default).
        items : dict (optional)
            Generic items that live on grid elements. No item is created if the
            value is 'None' (default). Otherwise, dictionary describes the
            position of generic items on the grid. The structure is:

            .. code-block:: python

                {"grid_element": [grid_element], "element_id": [element_id]}

            where:

                - [grid_element] is a str or number-of-items-long array
                  containing strings of the grid element(s) on which the items
                  live. Valid locations depend on the grid type. If provided as a
                  string it is assumed that all items live on the same type of
                  grid element.
                - [element_id] is an array of integers identifying the grid
                  element ID on which each item resides.

            An example argument would be:

            .. code-block:: python

                {
                    "grid_element": numpy.array(["node"], ["node"], ["link"]),
                    "element_id": numpy.array([1], [5], [1]),
                }

        data_vars : dict (optional)
            Dictionary of the data variables to be recorded. The structure is:

            .. code-block:: python

                {
                    "variable_name_1": (["dimensions"], variable_data_1),
                    "variable_name_2": (["dimensions"], variable_data_2),
                }

            where:

                - 'variable_name...' is a string of the variable name (label)
                - ['dimensions'] is the dimension(s) over which the variable
                  exists: can be ['time'], ['item_id'] or ['item_id', 'time'].
                - variable_data is an array containing the data, its size must
                  match that of the variable dimension(s).
        attrs : dict (optional)
            Dictionary of global attributes on the DataRecord (metadata).
            Example: {'time_units' : 'y'}

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3, 3))

        Example of a DataRecord with time as the only dimension:

        >>> dr1 = DataRecord(
        ...     grid,
        ...     time=[0.0],
        ...     data_vars={"mean_elevation": (["time"], np.array([100]))},
        ...     attrs={"time_units": "y"},
        ... )

        DataRecord builds off of xarray Dataset, a multi-dimensional, in
        memory, array  database. Dataset implements the mapping interface with
        keys given by variable names and values given by DataArray objects for
        each variable name.

        A DataRecord can have dimensions 'time' and/or 'item_id'.

        The xarray Dataset is stored in the public attribute ``dataset``.

        Coordinates are one dimensional arrays used for label-based indexing.
        DataRecord inherits all the methods and attributes from
        ``xarray.Dataset``.

        >>> dr1.dataset.to_dataframe()
              mean_elevation
        time
        0.0              100
        >>> dr1.dataset.time.values
        array([0.])
        >>> dr1.variable_names
        ['mean_elevation']
        >>> dr1.dataset["mean_elevation"].values
        array([100])
        >>> list(dr1.dataset.attrs.items())
        [('time_units', 'y')]

        >>> list(dr1.dataset.attrs.items())
        [('time_units', 'y')]

        Example of a DataRecord with item_id as the only dimension:

        >>> my_items2 = {
        ...     "grid_element": np.array(("node", "link"), dtype=str),
        ...     "element_id": np.array([1, 3]),
        ... }
        >>> dr2 = DataRecord(grid, items=my_items2)

        Note that both arrays (grid_element and element_id) have 1 dimension
        as they only vary along the dimension 'item_id'.

        >>> dr2.dataset.to_dataframe()[["grid_element", "element_id"]]
                grid_element  element_id
        item_id
        0               node           1
        1               link           3

        Example of a DataRecord with dimensions time and item_id:

        >>> my_items3 = {
        ...     "grid_element": np.array([["node"], ["link"]]),
        ...     "element_id": np.array([[1], [3]]),
        ... }
        >>> dr3 = DataRecord(grid, time=[0.0], items=my_items3)

        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.

        >>> dr3.dataset.to_dataframe()[["grid_element", "element_id"]]
                     grid_element  element_id
        item_id time
        0       0.0          node           1
        1       0.0          link           3

        """
        with_time = False
        if time is not None:
            with_time = True
            time = self._norm_time(time)

        # set attributes, if any
        attrs = attrs or {}
        if not isinstance(attrs, Mapping):
            raise TypeError("Attributes (attrs) passed to DataRecord must be a mapping")

        # save a reference to the grid
        self._grid = grid

        # depending on the grid type, permitted locations for items vary
        self._permitted_locations = self._grid.groups

        # save dummy elements reference
        # check dummies and reformat into {"node": [0, 1, 2]}
        self._dummy_elements = dummy_elements or {}
        for at in self._permitted_locations:
            size = self._grid[at].size
            for item in self._dummy_elements.get(at, ()):
                if 0 <= item < size:
                    raise ValueError(
                        f"dummy id {at} {item} must be outside valid range [0, {size})"
                    )

        # set initial items, if any
        if items is not None:
            if not isinstance(items, Mapping):
                raise TypeError("items must be mapping")
            with raise_as(TypeError):
                require_contains(
                    items, required=("grid_element", "element_id"), name="items"
                )

            grid_elements = items["grid_element"]
            element_ids = items["element_id"]

            n_items = len(element_ids)

            if with_time:
                shape = (n_items, len(time))
            else:
                shape = (n_items,)

            # check that grid_element and element_id exist on the grid and
            # have valid format:
            grid_elements, element_ids = self._norm_elements_and_ids(
                grid_elements, element_ids
            )

            grid_elements, element_ids = self._broadcast_elements_and_ids(
                grid_elements, element_ids, shape=shape
            )

            # check that element IDs do not exceed number of elements
            # on the grid:
            self._check_element_ids(grid_elements, element_ids)

            dims = ("item_id",)
            coords = {"item_id": np.arange(n_items, dtype=np.intp)}

            if with_time:
                dims += ("time",)
                coords["time"] = time

            data_vars_dict = {
                "grid_element": (dims, grid_elements),
                "element_id": (dims, element_ids, {"dtype": np.intp}),
            }

        else:
            # no items, initial dictionary of variables is empty:
            data_vars_dict = {}
            coords = {"time": time} if with_time else {}

        # set variables, if any
        if data_vars is not None:
            if not isinstance(data_vars, Mapping):
                raise TypeError("data_vars must be mapping")

            valid_dims = {("item_id",)}
            if with_time:
                valid_dims |= {("time",), ("time", "item_id"), ("item_id", "time")}

            for key in data_vars:
                dims = tuple(data_vars[key][0])
                with raise_as(ValueError):
                    require_one_of(dims, allowed=valid_dims)

            # create complete dictionary of variables
            # (= initial data_vars_dict + additional user-defined data_vars):
            data_vars_dict.update(data_vars)

        # create an xarray Dataset:
        self._dataset = xr.Dataset(data_vars=data_vars_dict, coords=coords, attrs=attrs)

    def _norm_time(self, time: ArrayLike) -> NDArray[np.floating]:
        time = np.atleast_1d(time)
        with raise_as(TypeError):
            require_dtype(time, dtype=float, allow_cast=True, name="time")
        with raise_as(ValueError):
            require_shape(time, shape=("n_times",), name="time")
            require_sorted(time, strict=True, name="time")

        return np.asarray(time, dtype=float)

    def _norm_elements_and_ids(
        self,
        elements: ArrayLike,
        ids: ArrayLike,
    ) -> tuple[NDArray[np.str_], NDArray[np.intp]]:
        return (
            norm_grid_element(elements, allowed=self._permitted_locations),
            norm_element_id(ids),
        )

    def _broadcast_elements_and_ids(
        self,
        elements: NDArray[np.str_],
        ids: NDArray[np.intp],
        *,
        shape: tuple[int, ...],
    ) -> tuple[NDArray[np.str_], NDArray[np.intp]]:
        return (
            np.broadcast_to(elements, shape).copy(),
            np.broadcast_to(ids, shape).copy(),
        )

    def _check_element_ids(
        self,
        elements: NDArray[np.str_],
        ids: NDArray[np.integer],
    ) -> None:
        """Check that element_id values are valid."""
        require_dtype(ids, dtype=np.integer, name="ids")

        for at in self._permitted_locations:
            selected = ids[elements == at]
            if selected.size == 0:
                continue

            size = self._grid[at].size
            dummy_values = self._dummy_elements.get(at, ())

            invalid = (selected < 0) | (selected >= size)
            if dummy_values:
                invalid &= ~np.isin(selected, dummy_values)

            if np.any(invalid):
                raise ValueError(f"invalid ids for grid element {at!r}.")

    def add_record(
        self,
        time: ArrayLike | None = None,
        item_id: ArrayLike | None = None,
        new_item_loc: dict[str, Any] | None = None,
        new_record: dict[str, Any] | None = None,
    ) -> None:
        """Add data to the DataRecord.

        Append new values to existing variables in the DataRecord, optionally
        at specified ``time`` and/or ``item_id`` coordinates. It can also
        update the location of existing items.

        Parameters
        ----------
        time : array_like of float or int, optional
            Time(s) at which to add data. Must be one-dimensional and strictly
            increasing. Required if the DataRecord has a ``time`` dimension and
            ``new_item_loc`` is provided. Must be ``None`` for time-independent
            DataRecords.
        item_id : array_like of int, optional
            IDs of existing items to update. Must be one-dimensional. All values
            must already exist in the DataRecord. Use :meth:`add_item` to create
            new items.
        new_item_loc : dict, optional
            Mapping with keys ``"grid_element"`` and ``"element_id"`` specifying
            updated locations for items. Both ``time`` and ``item_id`` must be
            provided when using this argument.
        new_record : dict, optional
            Mapping of variable names to values to add. Values should be compatible
            with the dimensions implied by ``time`` and/or ``item_id``.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3, 3))

        Example of a DataRecord with dimensions time and item_id:

        >>> items = {
        ...     "grid_element": [["node"], ["link"]],
        ...     "element_id": [[1], [3]],
        ... }

        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.

        >>> dr = DataRecord(grid, time=0.0, items=items)

        Records relating to pre-existing items can be added to the DataRecord
        using the method 'add_record':

        >>> dr.add_record(
        ...     time=2.0,
        ...     item_id=0,
        ...     new_item_loc={
        ...         "grid_element": [["node"]],
        ...         "element_id": [[6]],
        ...     },
        ...     new_record={"item_size": (["item_id", "time"], [[0.2]])},
        ... )
        >>> dr.dataset["element_id"].values
        array([[ 1.,  6.],
               [ 3., nan]])
        >>> dr.get_data([2.0], [0], "item_size")
        array([0.2])

        The 'add_record' method can also be used to add a non item-related
        record:

        >>> dr.add_record(time=50.0, new_record={"mean_elev": (("time",), [110])})
        >>> dr.dataset["mean_elev"].to_dataframe()
              mean_elev
        time
        0.0         NaN
        2.0         NaN
        50.0      110.0
        """
        if time is not None and "time" not in self._dataset:
            raise KeyError("this DataRecord is time-independent; time must be None.")
        if item_id is not None and "item_id" not in self._dataset:
            raise KeyError("This DataRecord does not have an item_id dimension.")
        if new_item_loc is not None and item_id is None:
            raise ValueError("new_item_loc requires item_id.")
        if new_item_loc is not None and time is None:
            raise ValueError(
                "new_item_loc requires time; use set_data() to change item locations"
                " instead."
            )

        coords_to_add = {}
        new_data_vars = {}

        if item_id is not None:
            item_id = np.atleast_1d(item_id)

            require_dtype(item_id, dtype=np.intp, allow_cast=True)
            require_shape(item_id, shape=("n_items",))

            item_id = np.asarray(item_id, dtype=np.intp)

            if not np.all(np.isin(item_id, self._dataset["item_id"].values)):
                raise ValueError(
                    "all item_id values must already exist in this DataRecord;"
                    " use add_item() to add new items."
                )

            coords_to_add["item_id"] = item_id

        if time is not None:
            coords_to_add["time"] = self._norm_time(time)

        if new_item_loc is not None:
            with raise_as(KeyError):
                require_contains(
                    new_item_loc,
                    required=("grid_element", "element_id"),
                    name="new_item_loc",
                )

            new_grid_element = new_item_loc["grid_element"]
            new_element_id = new_item_loc["element_id"]

            new_grid_element, new_element_id = self._norm_elements_and_ids(
                new_grid_element, new_element_id
            )
            self._check_element_ids(new_grid_element, new_element_id)

            dims = tuple(dim for dim in ("item_id", "time") if dim in coords_to_add)
            new_data_vars = {
                "grid_element": (dims, new_grid_element),
                "element_id": (dims, new_element_id),
            }

        if new_record is not None:
            new_data_vars.update(new_record)

        ds_to_add = xr.Dataset(data_vars=new_data_vars, coords=coords_to_add)

        self._dataset = xr.merge((self._dataset, ds_to_add), compat="no_conflicts")

    def add_item(
        self,
        time: ArrayLike | None = None,
        new_item: dict[str, Any] | None = None,
        new_item_spec: dict[str, Any] | None = None,
    ) -> None:
        """Add new items to a DataRecord.

        Parameters
        ----------
        time : array_like of number, optional
            Time coordinate(s) for the new items. Required if this DataRecord has
            a ``time`` dimension and must be ``None`` otherwise.
        new_item : mapping
            Mapping that defines the new items. Must contain the keys
            ``"grid_element"`` and ``"element_id"``.
        new_item_spec : mapping, optional
            Additional variables to add for the new items, in xarray-style
            ``{name: (dims, values)}`` form.

        Notes
        -----
        New items are assigned consecutive ``item_id`` values. For time-dependent
        records, ``grid_element`` and ``element_id`` are broadcast over the new
        ``item_id`` and ``time`` coordinates as needed before merging into the
        underlying dataset.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord

        Create a time-independent DataRecord with two items:

        >>> grid = RasterModelGrid((3, 3))
        >>> items = {
        ...     "grid_element": ["node", "node"],
        ...     "element_id": [1, 3],
        ... }
        >>> dr = DataRecord(grid, items=items)

        Add two new items:

        >>> dr.add_item(
        ...     new_item={
        ...         "grid_element": ["node", "node"],
        ...         "element_id": [4, 5],
        ...     }
        ... )
        >>> dr.number_of_items
        4

        Create a time-dependent DataRecord with one time coordinate:

        >>> items = {
        ...     "grid_element": np.array([["node"], ["link"]]),
        ...     "element_id": np.array([[1], [3]]),
        ... }
        >>> dr = DataRecord(grid, time=[0.0], items=items)

        Add two new items at a new time:

        >>> dr.add_item(
        ...     time=[1.0],
        ...     new_item={
        ...         "grid_element": [["node"], ["node"]],
        ...         "element_id": [[4], [5]],
        ...     },
        ...     new_item_spec={
        ...         "size": (["item_id", "time"], [[10], [5]]),
        ...     },
        ... )
        >>> dr.number_of_items
        4

        The time coordinate is extended:

        >>> dr.time_coordinates
        [0.0, 1.0]

        Values for a new variable are defined only for the new items at the new time:

        >>> dr.dataset["size"][:, 1].values
        array([nan, nan, 10.,  5.])
        """
        has_time = "time" in self._dataset

        if time is None and has_time:
            raise ValueError("time is required for a time-dependent DataRecord")
        if time is not None and not has_time:
            raise KeyError("time must be None for a time-independent DataRecord")
        if not isinstance(new_item, Mapping):
            raise TypeError("new_item must be a mapping")

        time = self._norm_time(time) if has_time else None

        new_item_spec = {} if new_item_spec is None else dict(new_item_spec)

        with raise_as(KeyError):
            require_contains(
                new_item, required=("grid_element", "element_id"), name="new_item"
            )
        elements = new_item["grid_element"]
        ids = new_item["element_id"]

        elements, ids = self._norm_elements_and_ids(elements, ids)
        self._check_element_ids(elements, ids)

        n_existing = len(self._dataset["item_id"])
        if n_existing == 0:
            new_item_ids = np.arange(len(ids), dtype=np.intp)
        else:
            start = self._dataset["item_id"][-1] + 1
            new_item_ids = np.arange(start, start + len(ids), dtype=np.intp)

        coords_to_add = {"item_id": new_item_ids}
        dims = ("item_id",)
        shape = (len(ids),)
        if has_time:
            coords_to_add["time"] = time
            dims += ("time",)
            shape += (len(time),)

        elements, ids = self._broadcast_elements_and_ids(elements, ids, shape=shape)

        data_vars_dict = {"grid_element": (dims, elements), "element_id": (dims, ids)}

        # other variables:
        data_vars_dict.update(new_item_spec)

        # Dataset of new record:
        ds_to_add = xr.Dataset(data_vars=data_vars_dict, coords=coords_to_add)

        # Merge new record and original dataset:
        self._dataset = xr.merge((self._dataset, ds_to_add), compat="no_conflicts")

    def get_data(self, time=None, item_id=None, data_variable=None):
        """Return values of a variable at a given time and/or for selected items.

        Parameters
        ----------
        time : array_like of float or int, optional
            Time coordinate to select. Must be length 1 if provided.
        item_id : array_like of int, optional
            Indices of items to select.
        data_variable : str
            Name of the variable to retrieve.

        Returns
        -------
        ndarray
            Values of *data_variable* after applying any provided selections. The
            shape depends on the dimensions of the variable and the selections
            applied.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3, 3))

        >>> items = {
        ...     "grid_element": "node",
        ...     "element_id": np.array([[1], [3], [3], [7]]),
        ... }
        >>> data = {
        ...     "item_size": (
        ...         ["item_id", "time"],
        ...         np.array([[0.3], [0.4], [0.8], [0.4]]),
        ...     )
        ... }
        >>> dr = DataRecord(grid, time=[50.0], items=items, data_vars=data)

        Select by time and item_id:
        >>> dr.get_data(time=[50.0], item_id=[2], data_variable="element_id")
        array([3])

        Select all items at a time:
        >>> dr.get_data(time=[50.0], data_variable="item_size")
        array([0.3, 0.4, 0.8, 0.4])

        Select items without time:
        >>> dr.get_data(item_id=[1, 2], data_variable="grid_element")
        array([['node'], ['node']], dtype='<U6')
        """
        if data_variable is None:
            raise ValueError("data_variable must not be None")

        with raise_as(KeyError):
            require_contains(self._dataset, required=(data_variable,), name="dataset")

        required = ()
        if item_id is not None:
            required += ("item_id",)
        if time is not None:
            required += ("time",)

        with raise_as(KeyError):
            require_contains(self._dataset.coords, required=required, name="dataset")

        data = self._dataset[data_variable]

        selectors = {}
        if time is not None:
            time = require_array(np.atleast_1d(time), dtype=np.number, shape=(1,))
            selectors["time"] = _find_time(time[0], self._dataset["time"])

        if item_id is not None:
            selectors["item_id"] = _norm_item_id(item_id, self._dataset["item_id"])

        if selectors:
            data = data.isel(**selectors)

        return data.values

    def set_data(self, time=None, item_id=None, data_variable=None, new_value=np.nan):
        """Set a variable value at a model time and/or an item to a new value.

        The value of only one variable can be changed at a time using this
        method.

        Parameters
        ----------
        time : list or 1-D array of float or int
            The time coordinate of the record to set.
        item_id : list or 1-D array of int
            The item id of the record to set.
        data_variable : string
            The label of the variable to set.
        new_value : list or 1-D array
            The new value to give to the variable data.

        Returns
        -------
        DataRecord with updated data.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3, 3))

        Example of a DataRecord with dimensions time and item_id:

        >>> my_items4 = {
        ...     "grid_element": "node",
        ...     "element_id": np.array([[1], [3], [3], [7]]),
        ... }

        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.

        >>> my_data4 = {
        ...     "item_size": (
        ...         ["item_id", "time"],
        ...         np.array([[0.3], [0.4], [0.8], [0.4]]),
        ...     )
        ... }
        >>> dr4 = DataRecord(grid, time=[50.0], items=my_items4, data_vars=my_data4)
        >>> dr4.dataset["item_size"].values
        array([[0.3],
               [0.4],
               [0.8],
               [0.4]])
        >>> dr4.set_data([50.0], [2], "item_size", [0.5])
        >>> dr4.dataset["item_size"].values
        array([[0.3],
               [0.4],
               [0.5],
               [0.4]])
        """
        if data_variable not in self.variable_names:
            raise KeyError(
                "the variable '{}' is not in the " "DataRecord".format(data_variable)
            )

        # If record to be changed is 'grid_element' or 'element_id',
        # check that provided grid_element is valid and that new
        # grid_element+element_id combination exist on the grid and
        # have valid format:
        if data_variable in ("grid_element", "element_id"):
            if data_variable == "grid_element":
                assoc_grid_element = new_value
                assoc_element_id = self.get_data(time, item_id, "element_id")[0]
            if data_variable == "element_id":
                if not isinstance(new_value, int):
                    raise ValueError(
                        "You have passed a non-integer "
                        "element_id to DataRecord, this is not "
                        "permitted"
                    )
                if new_value < 0:
                    raise ValueError(
                        "You have passed an element id below "
                        "zero. This is not permitted"
                    )
                assoc_element_id = new_value
                assoc_grid_element = self.get_data(time, item_id, "grid_element")[0]

            _ = norm_grid_element(assoc_grid_element, allowed=self._permitted_locations)

            if assoc_element_id >= self._grid[assoc_grid_element].size:
                raise ValueError(
                    f"The location {assoc_grid_element} {assoc_element_id}"
                    " does not exist on this grid"
                )

        if time is None:
            self._dataset[data_variable].values[item_id] = new_value
        else:
            try:
                len(time)
            except TypeError as exc:
                raise TypeError("time must be a list or a 1-d array") from exc
            try:
                # check that time coordinate already exists
                time_index = np.where(self._dataset.time.values == time)[0][0]
            except IndexError as exc:
                raise IndexError(
                    "The time you passed is not currently"
                    " in the DataRecord, you must change the value"
                    " you pass or first create the new time "
                    " coordinate using the add_record method"
                ) from exc

            if item_id is None:
                self._dataset[data_variable].values[time_index] = new_value
            else:
                try:
                    len(item_id)
                except TypeError as exc:
                    raise TypeError("item_id must be a list or a 1-d array") from exc
                try:
                    self._dataset["item_id"]
                    self._dataset[data_variable].values[item_id, time_index] = new_value
                except KeyError as exc:
                    raise KeyError("This DataRecord does not hold items") from exc

    def calc_aggregate_value(
        self,
        func,
        data_variable,
        at="node",
        filter_array=None,
        fill_value=np.nan,
        args=(),
        **kwargs,
    ):
        """Apply a function to a variable aggregated at grid elements.

        Parameters
        ----------
        func : function
            Function to apply to be aggregated.
        data_variable : str
            Name of variable on which to apply the function.
        at : str, optional
            Name of grid element at which to apply the function.
            Default is "node".
        filter_array: boolean array with dimensions matching that of the
            DataRecord (optional)
            Array to filter the DataRecord before aggregation.
        fill_value: float
            Fill value for array. Default is np.nan.
        args : tuple (optional)
            Additional positional arguments to pass to the function.
        **kwargs : key value pairs (optional)
            Additional keyword arguments to pass to func.

        Returns
        -------
        out : ndarray
            Array of size number-of-grid_elements (grid_elements is the group
            passed as 'at' argument).

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.data_record import DataRecord
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3))
        >>> element_id = [0, 0, 0, 0, 1, 2, 3, 4, 5, 9999]
        >>> volumes = [4, 5, 1, 2, 3, 4, 5, 6, 7, 1234]
        >>> ages = [10, 11, 12, 13, 14, 15, 16, 8, 10, 3456]
        >>> grid_element = "node"
        >>> data = {"ages": ages, "volumes": volumes}
        >>> dr = DataRecord(
        ...     grid,
        ...     dummy_elements={"node": [9999]},
        ...     items={"grid_element": "node", "element_id": np.array(element_id)},
        ...     data_vars={
        ...         "ages": (["item_id"], np.array(ages)),
        ...         "volumes": (["item_id"], np.array(volumes)),
        ...     },
        ... )
        >>> s = dr.calc_aggregate_value(func=xr.Dataset.sum, data_variable="ages")
        >>> s
        array([46., 14., 15., 16.,  8., 10., nan, nan, nan])
        >>> len(s) == grid.number_of_nodes
        True

        If you want to first filter the DataRecord and then aggregate, first
        create a filter array with dimensions matching that of the DataRecord
        and has `True` for entries that should be retained and False for
        entries that should be ignored.

        For example, if we wanted to aggregate volume for items with an age
        greater than 10 we would to the following:

        >>> f = dr.dataset["ages"] > 10.0
        >>> v_f = dr.calc_aggregate_value(
        ...     func=xr.Dataset.sum, data_variable="volumes", filter_array=f
        ... )
        >>> v_f
        array([ 8.,  3.,  4.,  5., nan, nan, nan, nan, nan])

        If we wanted the value for elements with no volume to be zero instead
        of np.nan we could use the keyword argument ``fill_value``.

        >>> f = dr.dataset["ages"] > 10.0
        >>> v_f = dr.calc_aggregate_value(
        ...     func=xr.Dataset.sum,
        ...     data_variable="volumes",
        ...     filter_array=f,
        ...     fill_value=0.0,
        ... )
        >>> v_f
        array([8., 3., 4., 5., 0., 0., 0., 0., 0.])

        An array of ``fill_value`` is returned when ``filter_array`` is all
        ``False`` (np.nan is the default value).

        >>> f = dr.dataset["ages"] > 4000.0
        >>> v_f = dr.calc_aggregate_value(
        ...     func=xr.Dataset.sum, data_variable="volumes", filter_array=f
        ... )
        >>> v_f
        array([nan, nan, nan, nan, nan, nan, nan, nan, nan])

        Other values can be specified for ``fill_value``.

        >>> f = dr.dataset["ages"] > 4000.0
        >>> v_f = dr.calc_aggregate_value(
        ...     func=xr.Dataset.sum,
        ...     data_variable="volumes",
        ...     filter_array=f,
        ...     fill_value=0.0,
        ... )
        >>> v_f
        array([0., 0., 0., 0., 0., 0., 0., 0., 0.])
        """
        filter_at = self._dataset["grid_element"] == at

        filter_valid_element = (self._dataset["element_id"] >= 0) * (
            self._dataset["element_id"] < self._grid[at].size
        )

        if filter_array is None:
            my_filter = filter_at * filter_valid_element
        else:
            my_filter = filter_at * filter_valid_element * filter_array

        if np.any(my_filter):
            # Filter DataRecord with my_filter and groupby element_id:
            filtered = self._dataset.where(my_filter).groupby("element_id")

            # Calculate values
            vals = filtered.map(func, *args, **kwargs)  # .reduce

            # Create a nan array that we will fill with the results of the sum
            # this should be the size of the number of elements, even if there are
            # no items living at some grid elements.
            out = fill_value * np.ones(self._grid[at].size)

            # put the values of the specified variable into the correct location
            # of the out array.
            out[vals.element_id.values.astype(int)] = vals[data_variable]

            return out
        else:
            return np.repeat(fill_value, self._grid[at].size)

    def ffill_grid_element_and_id(self):
        """Fill NaN values of the fields 'grid_element' and 'element_id'.

        Fields are filled by propagating values forward in time.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3, 3))

        Example of a DataRecord with dimensions time and item_id:

        >>> my_items3 = {
        ...     "grid_element": np.array([["node"], ["link"]]),
        ...     "element_id": np.array([[1], [3]]),
        ... }

        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.

        >>> dr3 = DataRecord(grid, time=[0.0], items=my_items3)

        Records relating to pre-existing items can be added to the DataRecord
        using the method 'add_record':

        >>> dr3.add_record(
        ...     time=[2.0, 3.0],
        ...     new_record={"mean_elevation": (["time"], np.array([200.0, 250.0]))},
        ... )

        Adding this data record created two new time coordinates. The
        grid_element and element_id of the items has been filled with 'nan'
        for these time coordinates.

        >>> dr3.dataset["grid_element"].values
        array([['node', nan, nan],
               ['link', nan, nan]], dtype=object)
        >>> dr3.dataset["element_id"].values
        array([[ 1., nan, nan],
               [ 3., nan, nan]])

        To fill these values with the last valid value, use the method
        ffill_grid_element_and_id:

        >>> dr3.ffill_grid_element_and_id()
        >>> dr3.dataset["grid_element"].values
        array([['node', 'node', 'node'],
               ['link', 'link', 'link']], dtype=object)
        >>> dr3.dataset["element_id"].values
        array([[1., 1., 1.],
               [3., 3., 3.]])

        In some applications, there may be no prior valid value. Under these
        circumstances, those values will stay as NaN. That is, this only
        forward fills, and does not backfill.

        >>> my_items3 = {
        ...     "grid_element": np.array([["node"], ["link"]]),
        ...     "element_id": np.array([[1], [3]]),
        ... }
        >>> dr3 = DataRecord(grid, time=[0.0], items=my_items3)
        >>> dr3.dataset["element_id"].values
        array([[1], [3]])
        >>> dr3.dataset["grid_element"].values
        array([['node'], ['link']], dtype='<U6')

        Next add some new items at a new time.

        >>> dr3.add_item(
        ...     time=[1.0],
        ...     new_item={
        ...         "grid_element": np.array([["node"], ["node"]]),
        ...         "element_id": np.array([[4], [4]]),
        ...     },
        ...     new_item_spec={"size": (["item_id", "time"], [[10], [5]])},
        ... )

        Two items have been added at a new timestep 1.0:

        >>> dr3.number_of_items
        4
        >>> dr3.time_coordinates
        [0.0, 1.0]
        >>> dr3.dataset["element_id"].values
        array([[ 1., nan],
               [ 3., nan],
               [nan,  4.],
               [nan,  4.]])
        >>> dr3.dataset["grid_element"].values
        array([['node', nan],
               ['link', nan],
               [nan, 'node'],
               [nan, 'node']], dtype=object)

        We expect that the NaN's to the left of the 4.s will stay NaN. And they
        do.

        >>> dr3.ffill_grid_element_and_id()
        >>> dr3.dataset["element_id"].values
        array([[ 1.,  1.],
               [ 3.,  3.],
               [nan,  4.],
               [nan,  4.]])
        >>> dr3.dataset["grid_element"].values
        array([['node', 'node'],
               ['link', 'link'],
               [nan, 'node'],
               [nan, 'node']], dtype=object)

        Finally, if we add a new time, we see that we need to fill in the
        full time column.

        >>> dr3.add_record(time=[2])
        >>> dr3.dataset["element_id"].values
        array([[ 1.,  1., nan],
               [ 3.,  3., nan],
               [nan,  4., nan],
               [nan,  4., nan]])
        >>> dr3.dataset["grid_element"].values
        array([['node', 'node', nan],
               ['link', 'link', nan],
               [nan, 'node', nan],
               [nan, 'node', nan]], dtype=object)

        And that forward filling fills everything as expected.

        >>> dr3.ffill_grid_element_and_id()
        >>> dr3.dataset["element_id"].values
        array([[ 1.,  1.,  1.],
               [ 3.,  3.,  3.],
               [nan,  4.,  4.],
               [nan,  4.,  4.]])

        >>> dr3.dataset["grid_element"].values
        array([['node', 'node', 'node'],
               ['link', 'link', 'link'],
               [nan, 'node', 'node'],
               [nan, 'node', 'node']], dtype=object)
        """

        ei = self._dataset["element_id"].values

        for i in range(ei.shape[0]):
            for j in range(1, ei.shape[1]):
                if np.isnan(ei[i, j]):
                    ei[i, j] = ei[i, j - 1]

        self._dataset["element_id"] = (["item_id", "time"], ei)

        ge = self._dataset["grid_element"].values
        for i in range(ge.shape[0]):
            for j in range(1, ge.shape[1]):
                if ge[i, j] not in self._permitted_locations:
                    ge[i, j] = ge[i, j - 1]
        self._dataset["grid_element"] = (["item_id", "time"], ge)

    @property
    def dataset(self):
        """The xarray Dataset that serves as the core datastructure."""
        return self._dataset

    @property
    def variable_names(self):
        """Return the name(s) of the data variable(s) in the record as a
        list."""
        _keys = []
        for key in self._dataset.to_dataframe().keys():
            _keys.append(key)
        return _keys

    @property
    def number_of_items(self):
        """Return the number of items in the DataRecord."""
        return len(self._dataset.item_id)

    @property
    def item_coordinates(self):
        """Return a list of the item_id coordinates in the DataRecord."""
        return self._dataset.item_id.values.tolist()

    @property
    def number_of_timesteps(self):
        """Return the number of time steps in the DataRecord."""
        return len(self._dataset.time)

    @property
    def time_coordinates(self):
        """Return a list of the time coordinates in the DataRecord."""
        return self._dataset.time.values.tolist()

    @property
    def earliest_time(self):
        """Return the earliest time coordinate in the DataRecord."""
        return min(self._dataset.time.values)

    @property
    def latest_time(self):
        """Return the latest time coordinate in the DataRecord."""
        return max(self._dataset.time.values)

    @property
    def prior_time(self):
        """Return the penultimate time coordinate in the DataRecord."""
        if self.number_of_timesteps < 2:
            return np.nan
        else:
            return sorted(self.time_coordinates)[-2]


def norm_element_id(ids: ArrayLike) -> NDArray[np.intp]:
    """Normalize element IDs to a 1-D array of ``np.intp``.

    Parameters
    ----------
    ids : array_like
        Element IDs.

    Returns
    -------
    ndarray of intp, shape (n_ids,)
        Normalized element IDs with dtype ``np.intp``.

    Raises
    ------
    ValueError
        If ``ids`` cannot be safely interpreted as integers.

    Examples
    --------
    >>> norm_element_id([0, 1, 2])
    array([0, 1, 2])

    >>> norm_element_id(5)
    array([5])

    >>> ids = norm_element_id(np.array([1, 2], dtype=np.int32))
    >>> ids.dtype == np.intp
    True
    """
    with raise_as(ValueError):
        ids = require_dtype(np.asarray(np.atleast_1d(ids)), dtype=np.integer)
    return ids.astype(np.intp, copy=False)


def norm_grid_element(
    elements: ArrayLike,
    *,
    allowed: Iterable[str] = (),
) -> NDArray[np.str_]:
    """Normalize grid element names to a NumPy string array.

    Parameters
    ----------
    elements : array_like of str or str
        Grid element name(s). A scalar string is treated as a single-element array.
    allowed : iterable of str, optional
        Valid grid element names. If provided, all values in ``elements`` must
        be contained in ``allowed``.

    Returns
    -------
    ndarray of str, shape (n_elements,)
        Normalized grid element names.

    Raises
    ------
    ValueError
        If any value in ``elements`` is not in ``allowed``.

    Examples
    --------
    >>> norm_grid_element("node")
    array(['node'], dtype='<U4')

    >>> norm_grid_element(["node", "link"])
    array(['node', 'link'], dtype='<U4')

    >>> norm_grid_element(["node", "link"], allowed=["node", "link", "patch"])
    array(['node', 'link'], dtype='<U5')
    """
    allowed = set(allowed)

    ge_dtype = f"<U{max(len(x) for x in allowed)}" if allowed else str
    elements = np.asarray(
        [elements] if isinstance(elements, str) else elements,
        dtype=ge_dtype,
    )

    if allowed and not set(np.unique(elements)) <= allowed:
        raise ValueError("elements must contain valid locations for this grid type.")

    return elements


def _find_time(time: float | int, times: ArrayLike) -> int:
    """Find the index of a time value in an array of times.

    Parameters
    ----------
    time : float or int
        Time value to locate.
    times : array_like of float or int
        Array of times to search.

    Returns
    -------
    int
        Index of the matching time in ``times``.

    Notes
    -----
    Floating-point times are matched with ``numpy.isclose``. Integer times are
    matched with exact equality.

    Examples
    --------
    >>> _find_time(2, [0, 1, 2, 3])
    2

    >>> _find_time(2.0, [0.0, 1.0, 2.0, 3.0])
    2
    """
    times = np.asarray(times)

    require_shape(times, shape=("n_times",), name="times")
    with raise_as(TypeError):
        require_dtype(times, dtype=np.number, name="times")

    time_is_float = np.issubdtype(np.asarray(time).dtype, np.floating)
    times_is_float = np.issubdtype(times.dtype, np.floating)

    if time_is_float or times_is_float:
        matches = np.isclose(times, time)
    else:
        matches = times == time

    n_matches = np.count_nonzero(matches)
    if n_matches == 0:
        raise IndexError("time does not match any of the provided times")
    if n_matches > 1:
        raise IndexError("time matches more than one of the provided times")

    return int(np.flatnonzero(matches)[0])


def _norm_item_id(item_id: ArrayLike, values: ArrayLike) -> int:
    item_id = np.asarray(np.atleast_1d(item_id))

    require_shape(item_id, shape=("n_items",), name="item_id")
    with raise_as(TypeError):
        require_dtype(item_id, dtype=np.integer, name="item_id")

    try:
        values[item_id]
    except IndexError as exc:
        raise IndexError(
            "The item_id you passed does not exist in this DataRecord"
        ) from exc

    return item_id
