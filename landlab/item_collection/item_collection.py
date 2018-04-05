#! /usr/bin/env python
"""

"""
import numpy as np
from six import string_types

from pandas import DataFrame

from landlab.field import GroupError

_LOCATIONS = {'node': 'number_of_nodes',
              'patch': 'number_of_patches',
              'link': 'number_of_links',
              'corner': 'number_of_corners',
              'face': 'number_of_faces',
              'cell': 'number_of_cells'}

_FILL_VALUE = np.nan

class ItemCollection(object):
    """
    Datastructure to hold generic items that live on grid elements.

    This is a base class to contain the majority of core ItemCollection
    functionality. It inherits from the Pandas Dataframe.

    It is useful to define a few terms used in this documentation.

    1. A **grid_element** is a geometric part of the model grid. Links, nodes,
    cells, and patches are all examples of grid elements. A specific instance of
    a model grid will have a fixed number of each of these grid elements.

    2. An **item** is a generic thing that, at minimum, has a location on the
    grid defined by a grid element name, and a grid element id.

    3. A **variable** is an attribute of an **item** that is not related to its
    location on the grid.

    4. Each **item** has a **value** for each **variable**.

    ItemCollection requires that all data live on a grid element, with an
    element id that is less than the number of those elements. For example if
    you have a grid with only 100 links, no item can live at link 100 or
    link -3.

    Methods
    -------
    add_variable
    add_item
    get_value
    set_value
    get_items_on_grid_element
    calc_aggregate_value

    """

    def __init__(self, grid, data=None, grid_element=None, element_id=None):
        """
        Parameters
        ----------
        grid : ModelGrid
        data : dictionary
            A group of number-of-items long arrays. All arrays must be the same
            length. The dictionary keys will become the column names of the
            ItemCollection.
        grid_element : str or number-of-items long array
            The type of grid element each element lives on. The element type must
            be consistent with the type of grid provided (e.g. only nodes and links
            are valid if grid is of type NetworkModelGrid). If provided as a string
            it is assumed that all items live on the same type of grid element.
            Valid locations are: node, link, cell, patch, corner, face
        element_id : number-of-items long array
            The grid element id where each item resides.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.item_collection import ItemCollection
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid(3,3)
        >>> element_id = [0, 0, 1, 1, 2, 3, 5]
        >>> volume = [1, 2, 3, 4, 5, 6, 7]
        >>> age = [10, 11, 12, 13, 14, 15, 16]
        >>> rock = ['sand', 'sand', 'silt', 'sand', 'silt', 'clay', 'clay']
        >>> data = {'age': age,
        ...         'volume': volume,
        ...         'rock': rock}
        >>> ic = ItemCollection(grid,
        ...                     data = data,
        ...                     grid_element ='node',
        ...                     element_id = element_id)

        The ItemCollectionData is stored in a Pandas DataFrame. Pandas
        Dataframes can store all different types of variables (e.g. float, int,
        string).

        >>> print(ic.DataFrame)
          grid_element  element_id  age  rock  volume
        0         node           0   10  sand       1
        1         node           0   11  sand       2
        2         node           1   12  silt       3
        3         node           1   13  sand       4
        4         node           2   14  silt       5
        5         node           3   15  clay       6
        6         node           5   16  clay       7

        The column on the left is the item id (or Pandas Index), grid_element
        and element_id area always the first two variables in the dataframe,
        and the remaining variables are provided in the order given by
        `np.sort`.

        Grid element can be any of the valid grid elements for the grid type
        provided.

        >>> data = {'age': age,
        ...         'volume': volume}
        >>> grid_element = ['node', 'link', 'node', 'link', 'node', 'link', 'node']
        >>> ic = ItemCollection(grid,
        ...                     data = data,
        ...                     grid_element = grid_element,
        ...                     element_id = element_id)

        To add another variable, use the `add_variable` function.

        >>> density = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        >>> ic.add_variable('density', density)

        To add one or more additional items, or set of items, use the
        `add_item` function.

        >>> new_element_id = [5, 8, 8]
        >>> new_volume = [1, 5, 9]
        >>> new_age = [13, 14, 15]
        >>> new_density = [0.4, 0.5, 0.7]
        >>> data = {'age': new_age,
        ...         'volume': new_volume,
        ...         'density': new_density}
        >>> ic.add_item(data = data,
        ...              grid_element = 'node',
        ...              element_id = new_element_id)
        >>> print(ic.DataFrame)
          grid_element  element_id  age  density  volume
        0         node           0   10      0.1       1
        1         link           0   11      0.2       2
        2         node           1   12      0.3       3
        3         link           1   13      0.4       4
        4         node           2   14      0.5       5
        5         link           3   15      0.6       6
        6         node           5   16      0.7       7
        7         node           5   13      0.4       1
        8         node           8   14      0.5       5
        9         node           8   15      0.7       9

        To get the value of an item, use the `get_value` function. If you
        request one item, a single value is returned.

        >>> val = ic.get_value(item_id = 0, variable = 'age')
        >>> print(val)
        10

        You can also request multiple items. You can get multiple items for the
        same variable.

        >>> val = ic.get_value(item_id=[0,6], variable='age')
        >>> print(val)
        0    10
        6    16
        Name: age, dtype: int64
        >>> type(val)
        <class 'pandas.core.series.Series'>

        Note here that a Pandas object has been returned instead of a numpy
        array. It has information about the variable and the item id. To get
        the values only, do the following:

        >>> val.values
        array([10, 16])

        You can also get multiple variables for the same item.

        >>> val = ic.get_value(item_id = 0, variable = ['age', 'volume'])
        >>> print(val)
        age       10
        volume     1
        Name: 0, dtype: object

        Finally, a combination of is also possible. For example if we wanted
        the age and volume of of items 3 and 5, we would use:

        >>> val = ic.get_value(item_id = [3, 5], variable = ['age', 'volume'])
        >>> print(val)
           age  volume
        3   13       4
        5   15       6

        Note here that you are getting both the age and volume of both items
        3 and 5. They will be returned in the order specified with the keyword
        arguments to`get_value`.

        To change the value of an item, use the `set_value` function.

        >>> ic.set_value(item_id = 0, variable = 'age', value = 20.)
        >>> print(ic.DataFrame)
          grid_element  element_id   age  density  volume
        0         node           0  20.0      0.1       1
        1         link           0  11.0      0.2       2
        2         node           1  12.0      0.3       3
        3         link           1  13.0      0.4       4
        4         node           2  14.0      0.5       5
        5         link           3  15.0      0.6       6
        6         node           5  16.0      0.7       7
        7         node           5  13.0      0.4       1
        8         node           8  14.0      0.5       5
        9         node           8  15.0      0.7       9

        To get all items that are at a particular grid element, use the
        `get_items_on_grid_element` function.

        >>> items = ic.get_items_on_grid_element(at='node', element_id=1)
        >>> print(items)
          grid_element  element_id   age  density  volume
        2         node           1  12.0      0.3       3

        You can calculate aggregated summary statistics using the
        `calc_aggregate_value` function.

        >>> s = ic.calc_aggregate_value(np.sum, 'age', at='node')
        >>> print(s)
        [ 20.  12.  14.  nan  nan  29.  nan  nan  29.]

        Here the function passed can be any function that returns a single
        value. If aggregated at node, the result will be a numpy array with
        size `grid.number_of_nodes`. Nodes that have no items at them will have
        a detault fill value of `np.nan`

        """
        # save a reference to the grid
        self._grid = grid

        # Get the locations that are permitted on this grid.
        permitted_locations = []
        for loc in _LOCATIONS:
            try:
                self._grid.keys(loc)
                permitted_locations.append(loc)
            except GroupError:
                pass

        self.permitted_locations = permitted_locations

        self._check_sizes(data)

        grid_element = self._check_grid_element_and_id(grid_element, element_id)

        # data column names =
        self.variable_names = list(np.sort(list(data.keys())))

        # add grid element and element ID to data frame
        data['grid_element'] = grid_element
        data['element_id'] = element_id

        # initialized the PD dataframe now that we've done checks
        self.DataFrame = DataFrame(data)

        # order DataFrame
        self._column_order = ['grid_element', 'element_id'] + self.variable_names
        self.DataFrame = self.DataFrame[self._column_order]

        # check that element IDs do not exceed number of elements on this grid
        self._check_element_id_values()

    def _check_sizes(self, data):
        """Check the sizes of data provided"""
        # get the number of elements in the dataset:
        num_items = []
        for dat in data.keys():
            # check the size of all parts of data, either length 1 or
            # length num_items
            ni = len(data[dat])
            num_items.append(ni)

        if np.all(num_items[0] == np.array(num_items)):
            self.number_of_items = num_items[0]
        else:
            raise ValueError(('Data passed to ItemCollection must be '
                              ' the same length.'))

    def _check_element_id_values(self):
        """Check that element_id values are valid."""
        for at in self.permitted_locations:

            max_size = self._grid[at].size

            selected_elements = self.DataFrame.loc[self.DataFrame['grid_element'] == at, 'element_id']

            if selected_elements.size > 0:
                if max(selected_elements) >= max_size:
                    raise ValueError(('An item residing at ' + at + ' has an '
                                      'element_id larger than the size of this '
                                      'part of the grid.'))
                less_than_zero = selected_elements < 0
                if any(less_than_zero):
                    raise ValueError(('An item residing at ' + at + ' has '
                                      'an element id below zero. This is not '
                                      'permitted.'))
        dtype = self.DataFrame['element_id'].dtype
        if dtype != int:
            raise ValueError(('You have passed a non integer element id. to '
                             'ItemCollection, this is not permitted.'))

    def _check_grid_element_and_id(self, grid_element, element_id):
        """Check that grid_element and element_id are the right size."""
        # make sure that grid element is of a permitted type and the correct size.
        if isinstance(grid_element, string_types):
            if grid_element in self.permitted_locations:
                pass
            else:
                raise ValueError(('Location index provided: ' + grid_element +
                                  ' is not a permitted location for this grid '
                                  'type.'))
            ge_name = grid_element
            grid_element = np.empty((self.number_of_items, ), dtype=object)
            grid_element.fill(ge_name)

        else:
            for loc in grid_element:
                if loc in self.permitted_locations:
                    pass
                else:
                    raise ValueError(('Location index provided: ' + loc + ' is not'
                                     ' a permitted location for this grid type.'))

        if len(grid_element) != self.number_of_items:
            raise ValueError(('grid_element passed to ItemCollection must be '
                              ' the same length as the data or 1.'))

        if len(element_id) != self.number_of_items:
            raise ValueError(('element_id passed to ItemCollection must be '
                              ' the same length as the data.'))

        return grid_element

    def add_variable(self, variable, values):
        """Add a new variable to the ItemCollection.

        Parameters
        ----------
        variable : str
            The name of the new variable.
        values : array-like
            The values for the new variable. The number of values must match the
            number of items.
        """
        if isinstance(variable, string_types):
            if np.array(values).size == self.DataFrame.shape[0]:
                pass
            else:
                raise ValueError(('Values passed to add_variable must have '
                                  'the same size as the current ItemCollection.'))

            if variable in self.DataFrame.columns.values:
                raise ValueError(('Variable name ' + variable + ' passed to add_variable '
                                  'already exists. This is not permitted.'))
        else:
            raise ValueError(('Variable name passed to add_variable must be of '
                              'type string.'))

        # assign variable to values
        self.DataFrame[variable] = values

        # reassign variable names and re-order dataframe
        self.variable_names = list(np.sort(list(self.DataFrame)[2:]))
        self._column_order = ['grid_element', 'element_id'] + self.variable_names
        self.DataFrame = self.DataFrame[self._column_order]

    def add_item(self, data=None, grid_element=None, element_id=None):
        """Add new items to the ItemCollection

        Parameters
        ----------
        data : dictionary
            A group of number-of-items long arrays. All arrays must be the same
            length. The dictionary keys will become the column names of the
            ItemCollection. No new variables are permitted (use 'add_variable')
            and all existing variables must be included.
        grid_element : str or number-of-items long array
            The type of grid element each element lives on. The element type must
            be consistent with the type of grid provided (e.g. only nodes and links
            are valid if grid is of type NetworkModelGrid). If provided as a string
            it is assumed that all items live on the same type of grid element.
            Valid locations are: node, link, cell, patch, corner, face
        element_id : number-of-items long array
            The grid element id where each item resides.

        """

        self._check_sizes(data)

        grid_element = self._check_grid_element_and_id(grid_element, element_id)

        # add grid element and element ID to data frame
        data['grid_element'] = grid_element
        data['element_id'] = element_id

        new_data = DataFrame(data)

        old_columns = self.DataFrame.columns.values
        new_columns = new_data.columns.values

        for colname in new_columns:
            if colname not in old_columns:
                raise ValueError(('A data associated with a new column name is '
                                  'being passed to ItemCollection using the '
                                  'method add_item. You must first add this '
                                  'variable using the method add_variable.'))
        for colname in old_columns:
            if colname not in new_columns:
                raise ValueError(('New items are being added to an '
                                  'ItemCollection that do not include already '
                                  'existing variables. When adding items to '
                                  'ItemCollection you must pass values for all '
                                  'existing variables.'))

        # append new data frame, ingoring its current index (which just adds)
        # additional indicies.
        self.DataFrame = self.DataFrame.append(new_data, ignore_index=True)

        # enforce column order
        self.DataFrame = self.DataFrame[self._column_order]

        # check that element IDs do not exceed number of elements on this grid
        self._check_element_id_values()

    def get_value(self, item_id=None, variable=None):
        """Get the value of an item.

        Parameters
        ----------
        item_id : int or iterable
            The id or ids of the items to query.
        variable : str or iterable of string_types
            The name or names of the variables to query.

        Returns
        -------
        Dataframe
        """
        # if either are provided as str/int, convert to list
        if isinstance(item_id, int):
            item_id = list([item_id])
        if isinstance(variable, string_types):
            variable = list([variable])
        # calculate the number returned
        number_returned = len(item_id) * len(variable)
        # return based on number returned.
        if number_returned == 1:
            return self.DataFrame.loc[item_id[0], variable[0]]
        elif len(item_id) == 1:
            return self.DataFrame.loc[item_id[0], variable]
        elif len(variable) == 1:
            return self.DataFrame.loc[item_id, variable[0]]
        else:
            return self.DataFrame.loc[item_id, variable]

    def set_value(self, item_id=None, variable=None, value=None):
        """Set the value of an item.

        Parameters
        ----------
        item_id : int
            The id of the item.
        variable : str
            The variable name.
        value :
            The value to set for ``variable`` and ``item_id``.

        """
        self.DataFrame.loc[item_id, variable] = value

    def get_items_on_grid_element(self, at=None, element_id = 0):
        """Get all items on a grid element with a particular element_id.

        Parameters
        ----------
        at : str
            Name of grid element type.
        element_id : int
            Element id.

        Returns
        -------
        Dataframe
        """

        vals = self.DataFrame.loc[(self.DataFrame['grid_element'] == at) &
                                  (self.DataFrame['element_id'] == element_id)]
        return vals

    def calc_aggregate_value(self, func, var, at='node', fill_value=_FILL_VALUE, args=(), **kwargs):
        """Apply a function to a variable aggregated at grid elements.

        Parameters
        ----------
        func : function
            Function to apply to aggregated
        var : str
            Column name of variable to sum
        at : str, optional
            Name of grid element at which to sum. Default is "node".
        fill_value : str, float, int, optional
            Value to use for grid element locations where no items are located.
            Default is np.nan.
        args : tuple, optional
            Additional positional arguments to pass to function in after the
            array/series
        **kwargs : key value pairs, optional
            Additional keyword arguments to pass to func.

        Returns
        -------
        out : ndarray
            Array of size (num_grid_elements,) where num_grid_elements is the
            number of elements at the location specified with the `at`
            keyword argument.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.item_collection import ItemCollection
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid(3,3)
        >>> element_id = [0, 0, 0, 0, 1, 2, 3, 4, 5]
        >>> volumes = [4, 5, 1, 2, 3, 4, 5, 6, 7]
        >>> ages = [10, 11, 12, 13, 14, 15, 16, 8, 10]
        >>> grid_element = 'node'
        >>> data = {'ages': ages,
        ...         'volumes':volumes}
        >>> ic = ItemCollection(grid,
        ...                     data = data,
        ...                     grid_element ='node',
        ...                     element_id = element_id)

        We can calculate aggregate values by passing the function we want to
        use to `calc_aggregate_value`.

        >>> s = ic.calc_aggregate_value(np.sum, 'ages')
        >>> print(s)
        [ 46.  14.  15.  16.   8.  10.  nan  nan  nan]
        >>> len(s) == grid.number_of_nodes
        True

        You can even pass functions that require additional positional
        arguments or keyword arguments. For example, in order to get the 25th
        percentile, we we do the following.

        >>> s = ic.calc_aggregate_value(np.percentile, 'ages', q=25)
        >>> print(s)
        [ 10.75  14.    15.    16.     8.    10.      nan    nan    nan]


        """
        # select those items located on the correct type of element,
        # group by element_id and sum.
        grouped = self.DataFrame.loc[self.DataFrame['grid_element'] == at].groupby('element_id')
        vals = grouped.agg(func, *args, **kwargs)

        # create a nan array that we will fill with the results of the sum
        # this should be the size of the number of elements, even if there are
        # no items living at some grid elements.
        out = fill_value * np.ones(self._grid[at].size)

        # put the values of the specified variable in to the correct location
        # of the out array.
        out[vals.index] = vals[var]

        return out
