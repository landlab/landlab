#!/usr/bin/env python3

import numpy as np
import xarray as xr


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
        time=None,
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

        # save a reference to the grid
        self._grid = grid

        # depending on the grid type, permitted locations for items vary
        self._permitted_locations = self._grid.groups

        # save dummy elements reference
        # check dummies and reformat into {"node": [0, 1, 2]}
        self._dummy_elements = dummy_elements or {}
        for at in self._permitted_locations:
            for item in self._dummy_elements.get(at, []):
                if (item < self._grid[at].size) and (item >= 0):
                    raise ValueError(f"Dummy id {at} {item} invalid")

        # set initial time coordinates, if any
        if isinstance(time, (list, np.ndarray)):
            self._times = np.array(time)
            self._number_of_times = len(self._times)
        elif time is not None:
            raise TypeError("time must be a list or numpy array")

        # set initial items, if any
        if items is not None:
            try:
                items.keys()
            except AttributeError as exc:
                # items is not a dict
                raise TypeError(
                    "You must provide an `items` dictionary "
                    "(see documentation for required format)"
                ) from exc
            try:
                _grid_elements, _element_ids = (
                    items["grid_element"],
                    items["element_id"],
                )
            except KeyError as exc:
                # grid_element and/or element_id not provided
                raise TypeError(
                    "You must provide an `items` dictionary,"
                    "(see documentation for required format)"
                ) from exc

            self._number_of_items = len(_element_ids)
            if len(_grid_elements) != self._number_of_items:
                if isinstance(_grid_elements, str):
                    pass
                else:
                    raise ValueError(
                        "The number of grid_element passed "
                        "to DataRecord must be 1 or equal "
                        "to the number of element_id."
                    )

            # check that grid_element and element_id exist on the grid and
            # have valid format:
            _grid_elements, _element_ids = self._check_grid_element_and_id(
                _grid_elements, _element_ids
            )

            # check that element IDs do not exceed number of elements
            # on the grid:
            self._check_element_id_values(_grid_elements, _element_ids)

            # create coordinates for the dimension 'item_id':
            self._item_ids = np.array(range(self._number_of_items), dtype=int)

            # create initial dictionaries of variables:
            if time is not None:
                data_vars_dict = {
                    "grid_element": (["item_id", "time"], _grid_elements),
                    "element_id": (
                        ["item_id", "time"],
                        _element_ids,
                        {"dtype": int},
                    ),
                }
                coords = {"time": self._times, "item_id": self._item_ids}
            else:
                # no time
                data_vars_dict = {
                    "grid_element": (["item_id"], _grid_elements),
                    "element_id": (["item_id"], _element_ids, {"dtype": int}),
                }
                coords = {"item_id": self._item_ids}

        else:
            # no items, initial dictionary of variables is empty:
            data_vars_dict = {}
            if time is not None:
                coords = {"time": self._times}
            else:  # no item and no time = no dimension
                coords = {}

        # set variables, if any
        if data_vars is not None:
            try:
                # check format (dict)
                data_vars.keys()
            except AttributeError as exc:
                raise TypeError(
                    "Data variables (data_vars) passed to "
                    "DataRecord must be a dictionary (see "
                    "documentation for valid structure)"
                ) from exc
            for key in data_vars.keys():
                # check dict structure and dims:
                if data_vars[key][0] not in (
                    ["time"],
                    ["item_id"],
                    ["time", "item_id"],
                    ["item_id", "time"],
                ):
                    raise ValueError(
                        "Data variable dimensions must be " "time and/or item_id"
                    )

            # create complete dictionary of variables
            # (= initial data_vars_dict + additional user-defined data_vars):
            data_vars_dict.update(data_vars)

        # set attributes, if any
        attrs = attrs or {}
        if not isinstance(attrs, dict):
            raise TypeError(
                "Attributes (attrs) passed to DataRecord" "must be a dictionary"
            )

        # create an xarray Dataset:
        self._dataset = xr.Dataset(data_vars=data_vars_dict, coords=coords, attrs=attrs)

    def _check_grid_element_and_id(self, grid_element, element_id):
        """Check the location and size of grid_element and element_id."""
        if isinstance(grid_element, str):
            # create list of grid_element for all items
            ge_name = grid_element
            if hasattr(self, "_number_of_times"):
                # if time
                grid_element = np.array(
                    np.empty(
                        (self._number_of_items, self._number_of_times), dtype=object
                    )
                )

                if element_id.shape != grid_element.shape:
                    element_id = np.broadcast_to(element_id, grid_element.shape)

            else:
                # no time
                grid_element = np.array(
                    np.empty((self._number_of_items,), dtype=object)
                )
            grid_element.fill(ge_name)

        # verify all grid elements are valid.
        for loc in grid_element.flatten():
            if loc not in self._permitted_locations:
                raise ValueError(
                    "One or more of the grid elements"
                    " provided is/are not permitted location"
                    " for this grid type"
                )

        return grid_element, element_id

    def _check_element_id_values(self, grid_element, element_id):
        """Check that element_id values are valid."""
        for at in self._permitted_locations:
            max_size = self._grid[at].size

            # this needs to work with 2d arrays (rows, col = np.where (so grid
            # element always needs to be at least 2d.))
            ind = np.nonzero(grid_element == at)
            selected_elements = element_id[ind]

            if selected_elements.size > 0:
                dummy_values = self._dummy_elements.get(at, [])
                index_values = np.arange(0, max_size)
                valid_values = np.concatenate((dummy_values, index_values))

                valid_elements = np.isin(selected_elements, valid_values)

                if not np.all(valid_elements):
                    raise ValueError("Invalid element_ids provided.")

        if not np.issubdtype(element_id.dtype, np.integer):
            raise ValueError(
                "You have passed a non-integer element_id to "
                "DataRecord, this is not permitted"
            )

    def add_record(self, time=None, item_id=None, new_item_loc=None, new_record=None):
        """Add a new record to the DataRecord.

        Unlike add_item, this method can support adding records that include
        new variables to the DataRecord. It can also support adding records
        that do not include time.

        Parameters
        ----------
        time : list or 1-D array of float or int
            Time step at which the record is to be added.
        item_id : list or 1-D array of int (optional)
            ID of the item to which the new record relates.
        new_item_loc: dict (optional)
            Dictionary of the new item location. If the new record is a change
            in the item location (grid_element and/or element_id), this field
            must be provided as:

            .. code-block:: python

                {"grid_element": [grid_element], "element_id": [element_id]}

            Both must be provided even if only one is being changed.

        new_record : dict
            Dictionary containing the new record. Structure should be:
            {'variable_name_1' : (['dimensions'], variable_data_1)}
            with:

                - 'variable_name_1' : name of the (potentially new) variable
                - ['dimensions'] : dimension(s) along which the new record
                  varies; can be ['time'], ['item_id] or ['item_id', 'time']
                - variable_data_1 : new data array, size must match the
                  variable dimension(s)

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
        ...     time=[2.0],
        ...     item_id=[0],
        ...     new_item_loc={
        ...         "grid_element": np.array([["node"]]),
        ...         "element_id": np.array([[6]]),
        ...     },
        ...     new_record={"item_size": (["item_id", "time"], np.array([[0.2]]))},
        ... )
        >>> dr3.dataset["element_id"].values
        array([[ 1.,  6.],
               [ 3., nan]])
        >>> dr3.get_data([2.0], [0], "item_size")
        array([0.2])

        The 'add_record' method can also be used to add a non item-related
        record:

        >>> dr3.add_record(time=[50.0], new_record={"mean_elev": (["time"], [110])})
        >>> dr3.dataset["mean_elev"].to_dataframe()
              mean_elev
        time
        0.0         NaN
        2.0         NaN
        50.0      110.0
        """
        if time is not None:
            try:
                # check that time is a dim of the DataRecord
                self._dataset["time"]
            except KeyError as exc:
                raise KeyError("This DataRecord does not record time") from exc

            if not isinstance(time, (list, np.ndarray)):
                # check input type
                raise TypeError(
                    "You have passed a time that is"
                    " not permitted, must be list or array"
                )
            else:
                if item_id is not None:
                    try:
                        # check that DataRecord holds items
                        self._dataset["item_id"]
                    except KeyError as exc:
                        raise KeyError("This DataRecord does not hold items") from exc
                    try:
                        # check that item_id is list or array
                        len(item_id)
                    except TypeError as exc:
                        raise TypeError("item_id must be a list or a 1D array") from exc
                    if not all(i in self._dataset["item_id"].values for i in item_id):
                        # check that item_id already exist
                        raise ValueError(
                            "One or more of the value(s) you "
                            "passed as item_id is/are not "
                            "currently in the DataRecord. Change"
                            " the input values create a new item"
                            "using the method add_item"
                        )
                    coords_to_add = {"time": np.array(time), "item_id": item_id}

                    # if item location is changed by this new record, check
                    # that both grid_element and element_id are provided:
                    if new_item_loc is not None:
                        try:
                            new_grid_element = new_item_loc["grid_element"]
                            new_element_id = new_item_loc["element_id"]
                        except KeyError as exc:
                            raise KeyError(
                                "You must provide a "
                                "new_item_loc dictionary with both "
                                "grid_element and element_id"
                            ) from exc
                        # check that grid_element and element_id exist
                        # on the grid and have valid format:
                        (
                            new_grid_element,
                            new_element_id,
                        ) = self._check_grid_element_and_id(
                            new_grid_element, new_element_id
                        )

                        # check that element IDs do not exceed number
                        # of elements on this grid:
                        self._check_element_id_values(new_grid_element, new_element_id)

                        _new_data_vars = {
                            "grid_element": (["item_id", "time"], new_grid_element),
                            "element_id": (["item_id", "time"], new_element_id),
                        }
                    else:
                        # new_item_loc is `None`
                        _new_data_vars = {}
                else:
                    # no item
                    coords_to_add = {"time": np.array(time)}
                    _new_data_vars = {}

        else:
            # no time
            if item_id is not None:
                if not all(i in self._dataset["item_id"].values for i in item_id):
                    # check that item_id already exist
                    raise ValueError(
                        "One or more of the value(s) you "
                        "passed as item_id is/are not "
                        "currently in the DataRecord. Change"
                        " the input values create a new item"
                        "using the method add_item"
                    )

                coords_to_add = {"item_id": np.array(item_id)}
                _new_data_vars = {}

                # no time so if item location needs to be changed,
                # user should use set_data
                if new_item_loc is not None:
                    raise ValueError(
                        "Use the method set_data to change the "
                        "location of an item in this DataRecord"
                    )
            else:
                # no item
                _new_data_vars = {}
                coords_to_add = {}

        if new_record is not None:
            # add new_record to dict of variables to add
            _new_data_vars.update(new_record)

        # create dataset of new record
        ds_to_add = xr.Dataset(data_vars=_new_data_vars, coords=coords_to_add)

        # merge new record and original dataset
        self._dataset = xr.merge((self._dataset, ds_to_add), compat="no_conflicts")

    def add_item(self, time=None, new_item=None, new_item_spec=None):
        """Add new item(s) to the current DataRecord.

        Parameters
        ----------
        time : list or 1-D array of float or int
            Time step at which the items are to be added.
        new_item : dict
            Structure is:

            .. code-block:: python

                {"grid_element": [grid_element], "element_id": [element_id]}

            where:

                - [grid_element] is str or number-of-items long array
                  containing strings of the grid element(s) on which the items
                  live. Valid locations depend on the grid type. If provided as
                  a string it is assumed that all items live on the same type of
                  grid element.
                - [element_id] is an array of integers identifying the grid
                  element ID on which each item resides.

            An example argument would be:

            .. code-block:: python

                {
                    "grid_element": numpy.array([["node"], ["node"], ["link"]]),
                    "element_id": numpy.array([[1], [5], [1]]),
                }

        new_item_spec : dict (optional)
            Dictionary containing any data variables (other than
            'grid_element' and 'element_id') relating to the new item(s) to be
            added. Structure is:

            .. code-block:: python

                {"variable_name_1": (["dimensions"], variable_data_1)}

            where:

                - 'variable_name_1' : name of the (potentially new) variable
                - ['dimensions'] : dimension(s) along which the new record
                  varies; can be ['time'], ['item_id] or ['item_id', 'time']
                - variable_data_1 : new data array, size must match the
                  variable dimension(s)

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

        Items can be added to a DataRecord that already holds similar items,
        using the method 'add_item':

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

        If a data variable is also added with the new items ('size' in this
        example), the values for this variable are filled with 'nan' for the
        pre-existing items:

        >>> dr3.dataset["size"][:, 1].values
        array([nan, nan, 10.,  5.])

        The previous line calls the values of the variable 'size', for all
        items, at time=1; the first two items don't have a value for the
        variable 'size'.
        """
        if time is None and "time" in self._dataset["grid_element"].coords:
            raise ValueError(
                "The items previously defined in this DataRecord"
                ' have dimensions "time" and "item_id", '
                'you must provide a "time" for the new item(s)'
            )

        if not isinstance(new_item, dict):
            raise TypeError(
                "You must provide an new_item dictionary "
                "(see documentation for required format)"
            )

        try:
            # check that dict contains correct entries
            _grid_elements, _element_ids = (
                new_item["grid_element"],
                new_item["element_id"],
            )

        except KeyError as exc:
            raise KeyError(
                "You must provide a new_item dictionary "
                "(see documentation for required format)"
            ) from exc

        number_of_new_items = len(new_item["element_id"])
        # first id of new item = last item in existing datarecord+1
        new_first_item_id = self._dataset["item_id"][-1].values + 1
        new_item_ids = np.array(
            range(new_first_item_id, new_first_item_id + number_of_new_items)
        )

        if time is not None:
            try:
                self._dataset["time"]
            except KeyError as exc:
                raise KeyError("This DataRecord does not record time") from exc
            if not isinstance(time, (list, np.ndarray)):
                raise TypeError(
                    "You have passed a time that is not "
                    "permitted, must be list or a 1-D array"
                )
            else:
                coords_to_add = {
                    "time": np.array(time),
                    "item_id": np.array(new_item_ids),
                }
                # check that grid_element and element_id exist
                # on the grid and have valid format
                _grid_elements, _element_ids = self._check_grid_element_and_id(
                    _grid_elements, _element_ids
                )

                # check that element IDs do not exceed number
                # of elements on this grid
                self._check_element_id_values(_grid_elements, _element_ids)

                data_vars_dict = {
                    "grid_element": (["item_id", "time"], _grid_elements),
                    "element_id": (["item_id", "time"], _element_ids),
                }

        else:
            # no time
            coords_to_add = {"item_id": np.array(new_item_ids)}
            # check that grid_element and element_id exist on
            # the grid and have valid format:
            _grid_elements, _element_ids = self._check_grid_element_and_id(
                _grid_elements, _element_ids
            )
            # check that element IDs do not exceed number of
            # elements on this grid
            self._check_element_id_values(_grid_elements, _element_ids)

            data_vars_dict = {
                "grid_element": (["item_id"], _grid_elements),
                "element_id": (["item_id"], _element_ids),
            }

        # other variables:
        if new_item_spec is not None:
            data_vars_dict.update(new_item_spec)

        # Dataset of new record:
        ds_to_add = xr.Dataset(data_vars=data_vars_dict, coords=coords_to_add)

        # Merge new record and original dataset:
        self._dataset = xr.merge((self._dataset, ds_to_add), compat="no_conflicts")

    def get_data(self, time=None, item_id=None, data_variable=None):
        """Get the value of a variable at a model time and/or for an item.

        Parameters
        ----------
        time : list or 1-D array of float or int (optional)
            The time coordinate of the record to get.
        item_id : list or 1-D array of int (optional)
            The item id of the record to get.
        data_variable : string
            The label of the variable to get.

        Returns
        -------
        object
            The value of *variable* at *time* and/or for *item_id*. The type of
            the returned object is dependent on the type of the variable value.

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
        >>> dr4.get_data([50.0], [2], "element_id")
        array([3])
        >>> dr4.get_data(time=[50.0], data_variable="item_size")
        array([0.3, 0.4, 0.8, 0.4])
        >>> dr4.get_data(item_id=[1, 2], data_variable="grid_element")
        array([['node'],
               ['node']], dtype=object)
        """
        try:
            self._dataset[data_variable]
        except KeyError as exc:
            raise KeyError(
                f"the variable {data_variable!r} is not in the DataRecord"
            ) from exc
        if time is None:
            if item_id is None:
                return self._dataset[data_variable].values
            else:
                try:
                    self._dataset["item_id"]
                except KeyError as exc:
                    raise KeyError("This DataRecord does not hold items") from exc
                try:
                    len(item_id)
                except TypeError as exc:
                    raise TypeError("item_id must be a list or a 1-D array") from exc
                try:
                    self._dataset["item_id"].values[item_id]
                except IndexError as exc:
                    raise IndexError(
                        "The item_id you passed does not exist " "in this DataRecord"
                    ) from exc

                return self._dataset.isel(item_id=item_id)[data_variable].values

        else:  # time is not None
            try:
                self._dataset["time"]
            except KeyError as exc:
                raise KeyError("This DataRecord does not record time") from exc
            try:
                len(time)
            except TypeError as exc:
                raise TypeError("time must be a list or a 1-D array") from exc
            try:
                time_index = int(self.time_coordinates.index(time[0]))
            except ValueError as exc:
                raise IndexError(
                    "The time you passed is not currently"
                    " in the DataRecord, you must change the value"
                    " you pass or first create the new time "
                    " coordinate using the add_record method"
                ) from exc
            if item_id is None:
                return self._dataset.isel(time=time_index)[data_variable].values
            else:
                try:
                    self._dataset["item_id"]
                except KeyError as exc:
                    raise KeyError("This DataRecord does not hold items") from exc
                try:
                    len(item_id)
                except TypeError as exc:
                    raise TypeError("item_id must be a list or a 1-D array") from exc
                try:
                    self._dataset["item_id"].values[item_id]
                except IndexError as exc:
                    raise IndexError(
                        "The item_id you passed does not exist " "in this DataRecord"
                    ) from exc
                return self._dataset.isel(time=time_index, item_id=item_id)[
                    data_variable
                ].values

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
            self._check_grid_element_and_id(assoc_grid_element, assoc_element_id)
            if assoc_element_id >= self._grid[assoc_grid_element].size:
                raise ValueError(
                    "The location "
                    + assoc_grid_element
                    + " "
                    + str(assoc_element_id)
                    + " does not exist on this grid"
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
        array([['node'],
               ['link']],
              dtype='<U4')

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
