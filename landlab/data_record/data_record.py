
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from six import string_types

from xarray import Dataset


class DataRecord(Dataset):
    """
    Datastructure to hold variables relating to the grid or to generic items
    that live on grid elements.

    This is a base class to contain the majority of core Datarecord
    functionality. It inherits from the xarray Dataset.

    Data variables can vary along one or both of the following dimensions:
        - time (model time)
        - item_id: variables can characterize a set of items (each identified
            by an individual id) that reside on the grid.

    Examples:
        - the variable 'mean_elevation' characterizes the grid and
            varies with time,
        - the variable 'clast__rock_type' characterizes a set of items (clasts)
            and varies with item_id,
        - the variable 'clast__size' can vary with both time and item_id

    If an item or set of items is defined, each item must be defined by the
    grid element and the element id at which it resides: e.g.,
    grid_element = 'node', element_id = 9. In this case, 'grid_element'
    and 'element_id' are default data variables (in addition to any
    user-specified variables).
    For each item, the element_id must be less than the number of this item's
    grid_element that exist on the grid: e.g, if the grid has 100 links, no
    item can live at link 100 or link -3.
    """
    _name = "DataRecord"

    _cite_as = """"""

    _input_var_names = ("",)

    _output_var_names = ("",)

    _var_units = {"": ""}

    _var_mapping = {"": ""}

    _var_doc = {"": ""}

    def __init__(self,
                 grid,
                 time=None,
                 items=None,
                 data_vars=None,
                 attrs=None,
                 compat='broadcast_equals'):
        """
        Parameters
        ----------
        grid : ModelGrid
        time : list or array (optional)
            Time steps to be recorded.
            - if length=1: corresponds to first recorded timestep
            - if length=3: corresponds to all recorded timesteps,
                            required format is [first_timestep,
                                                last_step,
                                                timestep_size]

        items : dict (optional)
            If None: no item is created
            Otherwise, dictionnary describing the position of generic items on
            the grid.
            Structure is:
                {'grid_element' : [grid_element],
                 'element_id' : [element_id]}
            With:
                - [grid_element] a str or number-of-items-long array containing
                strings of the grid element(s) on which the items live. Valid
                locations depend on the grid type. If provided as a string
                it is assumed that all items live on the same type of grid
                element.
                - [element_id] is an array of integers identifying the grid
                element ID on which each item resides.
            An example argument would be:
                {'grid_element' : numpy.array(['node'], ['node'], ['link']),
                 'element_id' :   numpy.array([1],      [5],      [1]     ])}

        data_vars : dict (optional)
            Dictionary of the data variables to be recorded.
            Structure is:
                {'variable_name_1' : (['dimensions'], variable_data_1),
                 'variable_name_2' : (['dimensions'], variable_data_2)}
                - 'variable_name' is a string of the variable name (label)
                - ['dimensions'] is the dimension(s) over which the variable
                exists: can be ['time'], ['item_id'] or ['item_id', 'time'].
                - variable_data is an array containing the data, its size must
                match that of the variable dimension(s).

        attrs : dict (optional)
            Dictionary of global attributes on the Datarecord (metadata).
            Example: {'time_units' : 'y'}

        compat: str (optional)
            String indicating how to compare variables of the same name
            for potential conflicts:
                - 'broadcast_equals': all values must be equal when variables
                are broadcast against each other to ensure common dimensions.
                - 'equals': all values and dimensions must be the same.
                - 'identical': all values, dimensions and attributes must be
                the same.
                Default is 'broadcast_equals'.


        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3,3))

        Example of a Datarecord with time as the only dimension:
        >>> dr1=DataRecord(grid, time=[0.],
        ...                data_vars={'mean_elevation' : (['time'],
        ...                                               np.array([100]))},
        ...                attrs={'time_units' : 'y'})

        Datarecord is a xarray.Dataset, a multi-dimensional, in memory, array
        database. Dataset implements the mapping interface with keys given by
        variable names and values given by DataArray objects for each variable
        name.
        A Datarecord can have dimensions 'time' and/or 'item_id'.
        Coordinates are one dimensional arrays used for label-based indexing.

        >>> dr1
        <xarray.DataRecord>
        Dimensions:         (time: 1)
        Coordinates:
          * time            (time) float64 0.0
        Data variables:
            mean_elevation  (time) int64 100
        Attributes:
            time_units:  y

        Datarecord inherits all the methods and attributes from xarray.Dataset.
        >>> dr1['mean_elevation'].values
        array([100])
        >>> dr1.dims
        Frozen(SortedKeysDict({'time': 1}))
        >>> dr1.attrs['time_units']
        'y'
        >>> dr1.to_dataframe()
              mean_elevation
        time
        0.0              100

        Example of a Datarecord with item_id as the only dimension:
        >>> my_items2 = {'grid_element': np.array(('node', 'link'), dtype=str),
        ...              'element_id': np.array([1, 3])}
        Note that both arrays have 1 dimension as they only vary along
        the dimension 'item_id'.
        >>> dr2=DataRecord(grid,
        ...                items=my_items2)
        >>> dr2.to_dataframe()
                    grid_element  element_id
            item_id
            0               node           1
            1               link           3

        Example of a Datarecord with dimensions time and item_id:
        >>> my_items3 = {'grid_element':np.array([['node'], ['link']]),
        ...              'element_id': np.array([[1],[3]])}
        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.
        >>> dr3=DataRecord(grid,
        ...                time=[0.],
        ...                items=my_items3)
        >>> dr3
        <xarray.DataRecord>
        Dimensions:       (item_id: 2, time: 1)
        Coordinates:
          * time          (time) float64 0.0
          * item_id       (item_id) int64 0 1
        Data variables:
            grid_element  (item_id, time) <U4 'node' 'link'
            element_id    (item_id, time) int64 1 3
        """

        # save a reference to the grid
        self._grid = grid

        # depending on the grid type, permitted locations (for items) vary:
        self.permitted_locations = self._grid.groups

        # TIME
        if time is not None:
            # create coordinates for the dimension 'time':
            if len(time) == 1:
                self._times = np.array(time)
                self.number_of_times = 1
            elif isinstance(time, (list, np.ndarray)):
                self._times = np.array(range(time[0], time[1]+time[2], time[2]))
                # time[1]+time[2] so that last timestep is not excluded
                self.number_of_times = len(self._times)
            else:
                raise TypeError(('time must be a list or an array'))

        # ITEMS
        if items is not None:
            if isinstance(items, dict):
                try:
                    self.grid_elements, self.element_ids = (
                            items['grid_element'], items['element_id'])
                except KeyError:
                    print('You must provide an ''items'' dictionnary,'
                          '(see documentation for required format)')
                self.number_of_items = len(self.element_ids)

                # check that grid_element and element_id exist on the grid and
                # have valid format:
                self.grid_elements = self._check_grid_element_and_id(
                        self.grid_elements, self.element_ids)

                # check that element IDs do not exceed number of elements
                # on this grid:
                self._check_element_id_values(self.grid_elements,
                                              self.element_ids)

                # create coordinates for the dimension 'item_id':
                self.item_ids = np.array(range(self.number_of_items))

                # create initial dictionaries of variables and coordinates:
                if time is not None:
                    # create coordinates for the dimension 'item_id':
                    self.item_ids = np.array(range(self.number_of_items))
                    data_vars_dict = {
                            'grid_element' : (
                                    ['item_id', 'time'], self.grid_elements),
                             'element_id' : (
                                     ['item_id', 'time'], self.element_ids)}
                    coords = {'time' : self._times,
                              'item_id' : self.item_ids}

                else: # no time
                    self.item_ids = np.array(range(self.number_of_items))
                    data_vars_dict = {
                            'grid_element' : (
                                    ['item_id'], self.grid_elements),
                             'element_id' : (
                                     ['item_id'], self.element_ids)}
                    coords = {'item_id' : self.item_ids}

            else:
                raise TypeError(('You must provide an ''items'' dictionary '
                                 '(see documentation for required format)'))
        else:
            # no items, initial dictionary of variables is empty:
            data_vars_dict = {}
            if time is not None:
                coords = {'time' : self._times}
            else: # no item and no time = no dimension
                coords = {}


        # VARIABLES
        if data_vars is not None:

            # TO DO: CHECK FORMAT OF THE PASSED DATA_VARS DICTIONARY:
            # - must be dictionnary,
            # - must have format: {'data_var_name' : ([dim(s)], data)}
            # - [dim(s)] must be ['time'], ['item_id'] or ['time', 'item_id']
            # - data must have size of its dimension(s): no need to check that,
            # this is done internally when instantiating the Dataset, raises
            # ValueError with problematic variables if conflict : TO TEST

            # Checks (to test):
            try:
                data_vars.keys() # check if dict
                for key in data_vars.keys(): # check dict format(?) and dims:
                    if ((data_vars[key][0] in (['time'],
                                          ['item_id'],
                                          ['time', 'item_id'],
                                          ['item_id', 'time'])) == False):
                        print('Data variable dimensions must be time and/or'
                              'item_id')

            except AttributeError:
                raise TypeError(('Data variables (data_vars) passed to'
                                 ' DataRecord must be a dictionary (see '
                                 'documentation for valid structure)'))

            # create complete dictionary of variables
            # (= initial data_vars_dict + additional user-defined data_vars):
            data_vars_dict.update(data_vars)

        # ATTRIBUTES
        if attrs is not None:
            try:
                attrs.keys()
            except AttributeError:
                raise TypeError(('Attributes (attrs) passed to DataRecord'
                                'must be a dictionary'))

        # initialize the xarray Dataset:
        # self.Dataset = Dataset(data_vars_dict, coords, attrs, compat)
        super().__init__(data_vars=data_vars_dict,
                                         coords=coords,
                                         attrs=attrs,
                                         compat=compat)

# Modified from ItemCollection:
    def _check_grid_element_and_id(self, grid_element, element_id, flag=None):
        """Check that grid_element and element_id are permitted and the
        right size."""
        if isinstance(grid_element, string_types):
            # all items are on same type of grid_element
            if grid_element in self.permitted_locations:
                pass
            else:
                raise ValueError(('Location provided: ' + grid_element +
                                  ' is not a permitted location for this grid '
                                  'type.'))

            #create list of grid_element for all items:
            ge_name = grid_element
            try: # if time
                grid_element = np.array(
                    np.empty((self.number_of_items,
                              self.number_of_times), dtype=object))
            except AttributeError: # no time
                grid_element = np.array(
                    np.empty((self.number_of_items, ), dtype=object))
            grid_element.fill(ge_name)

        else: # each item on different grid element
            for loc in grid_element:
                if loc[0] in self.permitted_locations or (
                        loc in self.permitted_locations): #depending on dims
                    pass
                else:
                    raise ValueError((
                            'One or more of the grid elements provided is/are'
                            ' not permitted location for this grid type.'))

        if flag==1: # testing values passed to add_record: don't check length
            pass
        else:
            if len(grid_element) != self.number_of_items:
                raise ValueError(('grid_element passed to DataRecord must be '
                                  ' the same length as the number of items or '
                                  '1.'))

            if len(element_id) != self.number_of_items:
                raise ValueError(('element_id passed to DataRecord must be '
                                  ' the same length as the number of items.'))

        return grid_element

    def _check_element_id_values(self, grid_element, element_id):
        """Check that element_id values are valid."""
        for at in self.permitted_locations:
            max_size = self._grid[at].size
            selected_elements_ind = [i for i in range(
                    len(grid_element)) if grid_element[i] == at]
            selected_elements = element_id[selected_elements_ind]

            if selected_elements.size > 0:
                if max(selected_elements) >= max_size:
                    raise ValueError(('An item residing at ' + at + ' has an '
                                      'element_id larger than the number of'
                                      + at + 'on the grid.'))
                less_than_zero = selected_elements < 0
                if any(less_than_zero):
                    raise ValueError(('An item residing at ' + at + ' has '
                                      'an element id below zero. This is not '
                                      'permitted.'))
        dtype = element_id.dtype
        if dtype != int:
            raise ValueError(('You have passed a non-integer element_id to '
                             'DataRecord, this is not permitted.'))

    def add_record(self,
                   model__time = None,
                   item_id=None,
                   new_item_loc=None,
                   new_record=None,
                   **kwargs):
        """ Add a time-related record to an existing DataRecord.

        Parameters
        ----------
        model__time : list or array of size 1
            Time step at which the record is to be added.
        item_id : int (optional)
            ID of the item to which the new record relates.
        new_item_loc: dict (optional)
            Dictionary of the new item location. If the new record is a change
            in the item location (grid_element and/or element_id), this field
            must be provided as:
                {'grid_element' : [grid_element],
                 'element_id' : [element_id]}
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
        >>> grid = RasterModelGrid((3,3))

        Example of a Datarecord with dimensions time and item_id:
        >>> my_items3 = {'grid_element': np.array([['node'], ['link']]),
        ...              'element_id': np.array([[1],[3]])}
        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.
        >>> dr3=DataRecord(grid,
        ...                time=[0.],
        ...                items=my_items3)

        Records relating to pre-existing items can be added to the Datarecord
        using the method 'add_record':
        >>> dr3.add_record(model__time=[2.0],
        ...                item_id=[0],
        ...                new_item_loc={'grid_element' : np.array([['node']]),
        ...                             'element_id' : np.array([[6]])},
        ...                new_record={'item_size':(
        ...                           ['item_id', 'time'], np.array([[0.2]]))})
        >>> dr3['element_id'].values
        array([[  1.,   6.],
               [  3.,  nan]])
        >>> dr3.get_data(2.0,0,'item_size')
        0.2

        The 'add_record' method can also be used to add a non item-related
        record:
        >>> dr3.add_record(model__time=[50.0],
        ...                new_record={'mean_elev': (['time'], [110])})
        >>> dr3.to_dataframe()
                     grid_element  element_id  item_size  mean_elev
        item_id time
        0       0.0          node         1.0        NaN        NaN
                2.0          node         6.0        0.2        NaN
                50.0          NaN         NaN        NaN      110.0
        1       0.0          link         3.0        NaN        NaN
                2.0           NaN         NaN        NaN        NaN
                50.0          NaN         NaN        NaN      110.0
        """

        if model__time is not None:
            if isinstance(model__time, (list, np.array)) == False:
                raise TypeError(('You have passed a model__time that is not'
                             'permitted, must be list or array.'))
            else:
                if item_id is not None:
                    # check that item_id already exist:
                    for i in item_id:
                        if i in self['item_id'].values:
                            pass
                        else:
                            raise ValueError ('There is no item with item_id' +
                                              i +', modify the value(s) you'
                                              'pass as ''item_id'' or create'
                                              'a new item using the method'
                                              'add_item.')



                    coords_to_add = {'time' : np.array(model__time),
                                     'item_id' : item_id}
                    # check that new_record has grid_element and element_id
                    # and that they are valid?

                    # if item location is changed by this new record, check
                    # that both grid_element and element_id are provided:
                    if new_item_loc is not None:
                        try:
                            new_grid_element = new_item_loc['grid_element']
                            new_element_id = new_item_loc['element_id']

                            # check that grid_element and element_id exist on
                            # the grid and have valid format:
                            new_grid_element = (
                                    self._check_grid_element_and_id(
                                    new_grid_element,
                                    new_element_id,
                                    flag=1))

                            # check that element IDs do not exceed number of
                            # elements on this grid:
                            self._check_element_id_values(new_grid_element,
                                                          new_element_id)

                            _new_data_vars = {
                                'grid_element' : (
                                        ['item_id', 'time'],
                                        new_grid_element),
                                 'element_id' : (
                                         ['item_id', 'time'],
                                         new_element_id)}

                        except KeyError:
                            raise KeyError(('You must provide an new_item_loc '
                                  'dictionnary with both grid_element and '
                                  'element_id'))

                else: # no item
                    coords_to_add = {'time' : np.array(model__time)}
                    _new_data_vars = {}

        else: # no time
            if item_id is not None:
                coords_to_add = {'item_id' : np.array(item_id)}
                _new_data_vars = {}

                # if item location is changed by this new record, check
                    # that both grid_element and element_id are provided:
                if new_item_loc is not None:
                    try:
                        new_grid_element, new_element_id = (
                                new_item_loc['grid_element'],
                                new_item_loc['element_id'])

                        # check that grid_element and element_id exist on
                        # the grid and have valid format:
                        new_grid_element = (
                                self._check_grid_element_and_id(
                                new_grid_element,
                                new_element_id,
                                flag=1))

                        # check that element IDs do not exceed number of
                        # elements on this grid:
                        self._check_element_id_values(new_grid_element,
                                                      new_element_id)

                        _new_data_vars = {
                            'grid_element' : (
                                    ['item_id'],
                                    new_grid_element),
                             'element_id' : (
                                     ['item_id'],
                                     new_element_id)}

                    except KeyError:
                        raise KeyError(('You must provide an new_item_loc '
                              'dictionnary with both grid_element and '
                              'element_id'))

        # add new_record to dict of variables to add
        _new_data_vars.update(new_record)

        # Dataset of new record:
        ds_to_add = Dataset(data_vars=_new_data_vars,
                            coords=coords_to_add,
                            **kwargs)
        # Merge new record and original dataset:
        #new_ds = xr.concat((self, ds_to_add), compat='equals')
        self.merge(ds_to_add, inplace='True', compat='no_conflicts')


    def add_item(self,
                  model__time = None,
                  new_item=None,
                  new_item_spec=None,
                  **kwargs):

        """ Add new item(s) to the current Datarecord.

        model__time : list or array of size 1
            Time step at which the items are to be added.
        new_item : dict
            Structure is:
                {'grid_element' : [grid_element],
                 'element_id' : [element_id]}
            with:
                - [grid_element] is str or number-of-items long array
                containing strings of the grid element(s) on which the items
                live. Valid locations depend on the grid type. If provided as
                a string it is assumed that all items live on the same type of
                grid element.
                - [element_id] is an array of integers identifying the grid
                element ID on which each item resides.
            An example argument would be:
                {'grid_element' : [['node'], ['node'], ['link']],
                 'element_id' :   [[1],      [5],      [1]     ]}
        new_item_spec : dict (optional)
            Dictionary containing any data variables (other than
            'grid_element' and 'element_id') relating to the new item(s) to be
            added.


        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3,3))

        Example of a Datarecord with dimensions time and item_id:
        >>> my_items3 = {'grid_element':np.array([['node'], ['link']]),
        ...              'element_id': np.array([[1],[3]])}
        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.
        >>> dr3=DataRecord(grid,
        ...                time=[0.],
        ...                items=my_items3)

        Items can be added to a Datarecord that already holds similar items,
        using the method 'add_item':
        >>> dr3.add_item(model__time=[1.0],
        ...              new_item={'grid_element' : np.array(
        ...                                              [['node'], ['node']]),
        ...                        'element_id' : np.array([[4],[4]])},
        ...              new_item_spec={'size': (
        ...                              ['item_id', 'time'], [[10],[5]])})

        Two items have been added at a new timestep 1.0:
        >>> dr3.coords
        Coordinates:
          * item_id  (item_id) int64 0 1 2 3
          * time     (time) float64 0.0 1.0

        If a data variable is also added with the new items ('size' in this
        example), the values for this variable are filled with 'nan' for the
        pre-existing items:
        >>> dr3['size'][:,1].values
        array([ nan,  nan,  10.,   5.])

        The previous line calls the values of the variable 'size', for all
        items, at time=1; the first two items don't have a value for the
        variable 'size'.

        """

        if model__time is None and 'time' in self['grid_element'].coords:
            raise ValueError(('The items previously defined in this Datarecord'
                              ' have dimensions "time" and "item_id", '
                              'please provide a "model__time" for the new'
                              'item(s)'))

        if isinstance(new_item, dict):
                try:
                    self.grid_elements, self.element_ids = (
                            new_item['grid_element'], new_item['element_id'])

                    number_of_new_items = len(new_item['element_id'])
                    # first id of new item = last item in existing datarecord+1
                    new_first_item_id = self['item_id'][-1].values + 1
                    new_item_ids = np.array(range(new_first_item_id,
                                     new_first_item_id+number_of_new_items))

                    if model__time is not None:
                            if isinstance(
                                    model__time, (list, np.array)) == False:
                                raise TypeError(('You have passed a '
                                                 'model__time that is not '
                                                 'permitted, must be list or '
                                                 'array.'))

                            else:
                                coords_to_add = {
                                        'time' : np.array(model__time),
                                        'item_id' : np.array(new_item_ids)}
                                # check that grid_element and element_id exist
                                # on the grid and have valid format:
                                self.grid_elements = (
                                        self._check_grid_element_and_id(
                                                self.grid_elements,
                                                self.element_ids))

                                # check that element IDs do not exceed number
                                # of elements on this grid
                                self._check_element_id_values(
                                        self.grid_elements, self.element_ids)

                                data_vars_dict = {'grid_element' : (
                                                ['item_id', 'time'],
                                                self.grid_elements),
                                                  'element_id' : (
                                                ['item_id', 'time'],
                                                self.element_ids)}

                    else: # model time is None
                        coords_to_add = {'item_id' : np.array(new_item_ids)}
                        # check that grid_element and element_id exist on
                        # the grid and have valid format:
                        self.grid_elements = self._check_grid_element_and_id(
                                                self.grid_elements,
                                                self.element_ids)
                        # check that element IDs do not exceed number of
                        # elements on this grid
                        self._check_element_id_values(
                                self.grid_elements, self.element_ids)

                        data_vars_dict = {'grid_element' : (
                                ['item_id'], self.grid_elements),
                             'element_id' : (['item_id'], self.element_ids)}

                except KeyError: # items dict does not contain correct entries
                    print('You must provide an ''items'' dictionnary'
                          '(see documentation for required format)')
        else: # new_item is not a dict
            print('You must provide an "items" dictionnary'
                          ' (see documentation for required format)')

        # other variables:
        if new_item_spec is not None:
            data_vars_dict.update(new_item_spec)

        # Dataset of new record:
        ds_to_add = Dataset(data_vars=data_vars_dict,
                            coords=coords_to_add)

        # Merge new record and original dataset:
        #new_ds = xr.concat((self, ds_to_add), compat='equals')
        self.merge(ds_to_add, inplace='True', compat='no_conflicts')

    def get_data(self, model__time=None, item_id=None, data_variable=None):
        """Get the value of a variable in the record at a model time
        and/or for an item.

        Parameters
        ----------
        model__time : integer or float (optional)
            The time coordinate of the record to get.
        item_id : integer (optional)
            The item id of the record to set.
        data_variable : string
            The label of the variable to get.

        Returns
        -------
        object
            The value of *variable* at *model__time* and/or for *item_id.
            The type of the returned object is dependent on the type of
            the variable value.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3,3))

        Example of a Datarecord with dimensions time and item_id:
        >>> my_items4 = {'grid_element' : 'node',
        ...              'element_id': np.array([[1],[3],[3],[7]])}
        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.
        >>> my_data4 = {'item_size': (['item_id', 'time'], np.array(
        ...                                 [[0.3],[0.4],[0.8],[0.4]]))}
        >>> dr4=DataRecord(grid,
        ...                time=[50.],
        ...                items=my_items4,
        ...                data_vars=my_data4)
        >>> dr4.get_data(50.,2,'element_id')
        3
        >>> dr4.get_data(model__time=50.,data_variable='item_size')
        [0.3, 0.4, 0.8, 0.4]
        >>> dr4.get_data(item_id=1, data_variable='grid_element')
        ['node']
        """
        if data_variable not in self.variable_names:
            raise KeyError("the variable '{}' is not in the "
                           "Datarecord".format(data_variable))
        if model__time is None:
            if item_id is None:
                return self[data_variable].values.tolist()
            else:
                return self.isel(
                        item_id=item_id)[data_variable].values.tolist()
        else:
            time_index=int(self.time_coordinates.index(model__time))
            if item_id is None:
                return self.isel(
                        time=time_index)[data_variable].values.tolist()
            else:
                return self.isel(time=time_index,
                                 item_id= item_id)[
                                 data_variable].values.tolist()

    def set_data(self,
                 model__time=None,
                 item_id=None,
                 data_variable=None,
                 new_value=np.nan):
        """Set the value of a variable in the record at a model time
        and/or for an item to a new value.

        Parameters
        ----------
        model__time : integer or float
            The time coordinate of the record to set.
        item_id : integer
            The item id of the record to set.
        variable : string
            The label of the variable to set.

        Returns
        -------
        Datarecord with updated data.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.data_record import DataRecord
        >>> grid = RasterModelGrid((3,3))

        Example of a Datarecord with dimensions time and item_id:
        >>> my_items4 = {'grid_element' : 'node',
        ...              'element_id': np.array([[1],[3],[3],[7]])}
        Note that both arrays have 2 dimensions as they vary along dimensions
        'time' and 'item_id'.
        >>> my_data4 = {'item_size': (['item_id', 'time'], np.array(
        ...                                 [[0.3],[0.4],[0.8],[0.4]]))}
        >>> dr4=DataRecord(grid,
        ...                time=[50.],
        ...                items=my_items4,
        ...                data_vars=my_data4)
        >>> dr4.to_dataframe()
                     grid_element  element_id  item_size
        item_id time
        0       50.0         node           1        0.3
        1       50.0         node           3        0.4
        2       50.0         node           3        0.8
        3       50.0         node           7        0.4
        >>> dr4.set_data(50.,2,'item_size', 0.5)
                     grid_element  element_id  item_size
        item_id time
        0       50.0         node           1        0.3
        1       50.0         node           3        0.4
        2       50.0         node           3        0.5
        3       50.0         node           7        0.4
        """


        # TO ADD TO THIS METHOD: if the modified data is in 'grid_element' or
        # 'element_id', check that the new_value is consistent with grid

        if data_variable not in self.variable_names:
            raise KeyError("the variable '{}' is not in the "
                           "Datarecord".format(data_variable))

        # If record to be changed is 'grid_element' or 'element_id',
        # check that provided grid_element is valid and that new
        # grid_element+element_id combination exist on the grid and
        # have valid format:
        if data_variable in ('grid_element', 'element_id'):
            if data_variable == 'grid_element':
                assoc_grid_element = new_value
                assoc_element_id = self.get_data(model__time,
                                                 item_id,
                                                 'element_id')
            if data_variable == 'element_id':
                assoc_element_id = new_value
                assoc_grid_element = self.get_data(model__time,
                                                   item_id,
                                                   'grid_element')
            _ = (self._check_grid_element_and_id(assoc_grid_element,
                                                 assoc_element_id,
                                                 flag=1))
            if assoc_element_id >= self._grid[assoc_grid_element].size:
                raise ValueError(('The location ' + assoc_grid_element + ' ' +
                                  str(assoc_element_id) +
                                  ' does not exist on this grid.'))
            if assoc_element_id < 0:
                raise ValueError(('You have passed an element id below zero. '
                                  ' This is not permitted.'))
            if type(assoc_element_id) != int:
                raise ValueError(('You have passed a non-integer element_id'
                                  ' to DataRecord, this is not permitted.'))

        if model__time is None:
            self[data_variable].loc[dict(item_id=item_id)] = new_value
        else:
            if item_id is None:
                self[data_variable].loc[dict(time=model__time)] = new_value
            else:
                self[data_variable].loc[
                        dict(time=model__time, item_id=item_id)] = new_value

    @property
    def variable_names(self):
        """Return the name(s) of the data variable(s) in the Datarecord.
        """
        return self.data_vars

    @property
    def number_of_timesteps(self):
        """Return the number of time steps in the Datarecord.
        """
        try:
            return self.number_of_times
        except AttributeError:
            print('This Datarecord does not record time')

    @property
    def time_coordinates(self):
        """Return a list of the time coordinates in the Datarecord.
        """
        try:
            return self.time.values.tolist()
        except AttributeError:
            print('This Datarecord does not record time')

    @property
    def earliest_time(self):
        """Return the earliest time coordinate in the Datarecord.
        """
        try:
            return min(self.time.values)
        except AttributeError:
            print('This Datarecord does not record time')

    @property
    def latest_time(self):
        """Return the latest time coordinate in the Datarecord.
        """
        try:
            return max(self.time.values)
        except AttributeError:
            print('This Datarecord does not record time')

    @property
    def prior_time(self):
        """Return the penultimate time coordinate in the Datarecord.
        """

        try:
            if self.number_of_timesteps < 2:
                return np.nan
            else:
                return sorted(self.time_coordinates)[-2]
        except AttributeError:
            print('This Datarecord does not record time')

