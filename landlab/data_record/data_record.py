
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from six import string_types

import xarray as xr
from xarray import Dataset

from landlab.field import GroupError

_LOCATIONS = {'node': 'number_of_nodes',
              'patch': 'number_of_patches',
              'link': 'number_of_links',
              'corner': 'number_of_corners',
              'face': 'number_of_faces',
              'cell': 'number_of_cells'}

_FILL_VALUE = np.nan

class DataRecord(Dataset):
    """

    Datastructure to hold variables relating to the grid or to generic items
    that live on grid elements.

    This is a base class to contain the majority of core ItemCollection
    functionality. It inherits from the xarray Dataset.

    Data variables can vary along one or both of the following dimensions:
        - time (model time)
        - item_id: variables can characterize a set of items (each identified
            by an individual id) that reside on the grid.

    Some examples:
        - the variable 'mean_elevation' characterizes the grid and
            varies with time,
        - the variable 'clast__rock_type' characterizes a set of items (clasts)
            and varies with item_id,
        - the variable 'clast__size' can vary with both time and item_id

    If an item or set of items is defined, each item must be defined by the
    grid element and the element id at which they each reside: e.g.,
    grid_element = 'node', element_id = 9. In this case, 'grid_element'
    and 'element_id' are default data variables (in addition to any
    user-specified variables).
    For each item, the element_id must be less than the number of this item's
    grid_element that exist on the grid: e.g, if you have a grid with 100
    links, no item can live at link 100 or link -3.


    """
    _name = "DataRecord"

    _cite_as = """"""

    _input_var_names = ("",)

    _output_var_names = ("",)

    _var_units = {
        "": "",
    }

    _var_mapping = {
        "": "node",
    }

    _var_doc = {
        "": "",
    }

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
            - if length=1: corresponds to first recorded timestep (later
                timesteps to be added manually by user)
            - if length=3, required format is [first_timestep,
                                            last_step,
                                            timestep_size]

        items : dict (optional)
            If None: no item is created
            Structure is:
                {'grid_element' : grid_element,
                 'element_id' : element_id}
            If items are to be created, items should be a 2 * number-of-items
            long array with [grid_element, element_id]:
                - grid_element is str or number-of-items long array containing
                strings of the grid element(s) on which the items live. Valid
                locations depend on the grid type. If provided as a string
                it is assumed that all items live on the same type of grid
                element.
                - element_id is an array of integers identifying the grid
                element ID on which each item resides.
            An example argument would be: [[['node'], ['node'], ['link']]
                                           [[1],      [5],      [1]     ]]

        data_vars : dict (optional)
            Dictionary of the data variables to be recorded.
            Structure is:
                {'variable_name_1' : (['dimensions'], variable_data_1),
                 'variable_name_2' : ...                                }
                - 'variable_name' is a string of the variable name
                - ['dimensions'] is the dimension(s) over which the variable
                exists, can be ['time'], ['item_id'] or ['time', 'item_id'].
                - variable_data is an array containing the data (can be empty),
                its size must match that of the variable dimension(s).

        attrs : dict (optional)
            Dictionary of global attributes on the Dataset (metadata).
            Suggested: {'time_units' : 'y'}

        compat: str (optional)
            String indicating how to compare variables of the same name
            for potential conflicts:
                - ‘broadcast_equals’: all values must be equal when variables
                are broadcast against each other to ensure common dimensions.
                - ‘equals’: all values and dimensions must be the same.
                - ‘identical’: all values, dimensions and attributes must be
                the same.
                Default value is 'broadcast_equals'.


        Examples
        --------
        >>> import numpy as np
        ...

        """

        # save a reference to the grid
        self._grid = grid

        # depending on the grid type, permitted locations (for items) vary:
        self.permitted_locations = self._grid.groups

        # TIME
        if time is not None:
            # create coordinates for the 'time' dimension:
            if len(time) == 1:
                self.times = np.array(time)
                self.number_of_times = 1
            elif isinstance(time, (list, np.ndarray)):
                self.times = np.array(range(time[0], time[1]+time[2], time[2])) #so that last timestep is not excluded
                self.number_of_times = len(self.times)
            else:
                raise TypeError(('time must be array')) #resolve imput type


        # ITEMS
        if items is not None:

            if isinstance(items, dict):
                try:
                    self.grid_element, self.element_id = items['grid_element'], items['element_id']
                except KeyError:
                    print('You must provide an ''items'' dictionnary,'
                          '(see documentation for required format)')
                self.number_of_items = len(self.element_id)

                # check that grid_element and element_id exist on the grid and
                #have valid format:
                self.grid_element = self._check_grid_element_and_id(
                        self.grid_element, self.element_id)

                # check that element IDs do not exceed number of elements
                # on this grid
                self._check_element_id_values(self.grid_element,
                                              self.element_id)

                self.item_ids = np.array(range(self.number_of_items))

                # create initial variables dictionary and coordinates:
                if time is not None:
                    data_vars_dict = {'grid_element' : (['item_id', 'time'], self.grid_element), #grid_element_init
                             'element_id' : (['item_id', 'time'], self.element_id)} #element_id_init
                    coords = {'time' : self.times,
                              'item_id' : self.item_ids}

                else: # no time
                    data_vars_dict = {'grid_element' : (['item_id'], self.grid_element),
                             'element_id' : (['item_id'], self.element_id)}
                    coords = {'item_id' : self.item_ids}

            else:
                raise TypeError(('You must provide an ''items'' dictionary '
                                 '(see documentation for required format)'))
        else:
            # no items, initial variables dictionary is empty:
            data_vars_dict = {}
            if time is not None:
                coords = {'time' : self.times}
            else:
                coords = {}


        # VARIABLES
        if data_vars is not None:

            # TO DO: CHECK FORMAT OF THE PASSED DATA_VARS DICTIONARY

            # create complete variables dictionary (=data_vars_dict+data_vars):
            data_vars_dict.update(data_vars)


# no need to check that length of data variables corresponds to dimension time
# and item_id, this is done internally when creating the Dataset, raises
# ValueError with problematic variables if conflict : TO TEST

        # ATTRIBUTES
        if attrs is not None:
            try:
                attrs.keys()
            except AttributeError:
                raise TypeError(('Attributes (attrs) passed to DataRecord'
                                'must be a dictionary'))

        # initialize the xarray Dataset:
        #self.Dataset = Dataset(data_vars_dict, coords, attrs, compat)
        super().__init__(data_vars=data_vars_dict,
                                         coords=coords,
                                         attrs=attrs,
                                         compat=compat)

# Modified from ItemCollection:
    def _check_grid_element_and_id(self, grid_element, element_id):
        """Check that grid_element and element_id are the right size."""
        # make sure that grid element is of a permitted type and the correct size.
        if isinstance(grid_element, string_types):
            #all items are on same type of grid_element
            if grid_element in self.permitted_locations:
                pass
            else:
                raise ValueError(('Location index provided: ' + grid_element +
                                  ' is not a permitted location for this grid '
                                  'type.'))

            #create list of grid_element for all items:
            ge_name = grid_element
            grid_element = np.empty((self.number_of_items, ), dtype=object)
            grid_element.fill(ge_name)

        else:
            for loc in grid_element:
                if loc[0] in self.permitted_locations:
                    pass
                else:
                    raise ValueError((
                            'Grid element provided: ' + loc[0] + ' is not'
                            ' a permitted location for this grid type.'))

        if len(grid_element) != self.number_of_items:
            raise ValueError(('grid_element passed to DataRecord must be '
                              ' the same length as the number of items or 1.'))

        if len(element_id) != self.number_of_items:
            raise ValueError(('element_id passed to DataRecord must be '
                              ' the same length as the number of items.'))

        return grid_element

# Modified from ItemCollection:
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



    def add_record(self, model__time = None, item_id=None, new_record=None, **kwargs):
        """ Add a time-related record to an existing DataRecord.
        New record can be a variable
        model__time must be list or array of size 1
        New record must be a dictionnary:
            {'variable_name_1' : (['dimensions'], variable_data_1)}
        Concatenate along one dimension at a time (e.g., new record of a given
        variable for all items, or new item with all time records(?), new
        non-item related variable record for a given timestep)
        """

        if model__time is not None:
            if isinstance(model__time, (list, np.array)) == False:
                raise TypeError(('You have passed a model__time that is not'
                             'permitted, must be list or array.'))
            else:
                if item_id is not None:
                    coords_to_add = {'time' : np.array(model__time),
                                     'item_id' : item_id}
                    # check that new_record has grid_element and element_id
                    # and that they are valid?

                else: # no item
                    coords_to_add = {'time' : np.array(model__time)}

        else: # no time
            if item_id is not None:
                coords_to_add = {'item_id' : item_id}
                # check that new_record has grid_element and element_id
                # and that they are valid?


        # Dataset of new record:
        ds_to_add = Dataset(data_vars=new_record,
                            coords=coords_to_add,
                            **kwargs)
        # Merge new record and original dataset:
        #new_ds = xr.concat((self, ds_to_add), compat='equals')
        self.merge(ds_to_add, inplace='True', compat='no_conflicts')


    def add_items(self,
                  model__time = None,
                  new_item=None,
                  new_item_spec=None,
                  **kwargs):

        """ Add new items to the current Datarecord.

        new_item : dict
            If None: no item is created
            Structure is:
                {'grid_element' : grid_element,
                 'element_id' : element_id}
            If items are to be created, items should be a 2 * number-of-items
            long array with [grid_element, element_id]:
                - grid_element is str or number-of-items long array containing
                strings of the grid element(s) on which the items live. Valid
                locations depend on the grid type. If provided as a string
                it is assumed that all items live on the same type of grid
                element.
                - element_id is an array of integers identifying the grid
                element ID on which each item resides.
            An example argument would be: [[['node'], ['node'], ['link']]
                                           [[1],      [5],      [1]     ]]
        new_item_spec : dict (optional)
            Dictionary containing any data variables (other than
            'grid_element' and 'element_id') relating to the new item(s) to be
            added.
        """

        if isinstance(new_item, dict):
                try:
                    self.grid_element, self.element_id = (
                            new_item['grid_element'], new_item['element_id'])
                except KeyError:
                    print('You must provide an ''items'' dictionnary,'
                          '(see documentation for required format)')
        else:
            print('You must provide an ''items'' dictionnary,'
                          '(see documentation for required format)')

        number_of_new_items = len(new_item['element_id'])
        # first id of new item = last item in existing datarecord+1:
        new_first_item_id = self['item_id'][-1] +1
        new_item_ids = range(new_first_item_id,
                             new_first_item_id+number_of_new_items)


        if model__time is not None:
            if isinstance(model__time, (list, np.array)) == False:
                raise TypeError(('You have passed a model__time that is not'
                             'permitted, must be list or array.'))
            else:
                if new_item is not None:
                    coords_to_add = {'time' : np.array(model__time),
                                     'item_id' : new_item_ids}
                    # check that grid_element and element_id exist on the grid
                    # and have valid format:
                    self.grid_element = self._check_grid_element_and_id(
                            self.grid_element, self.element_id)

                    # check that element IDs do not exceed number of elements
                    # on this grid
                    self._check_element_id_values(self.grid_element,
                                                  self.element_id)

        else: # no time
            if new_item_ids is not None:
                coords_to_add = {'item_id' : new_item_ids}
                # check that new_record has grid_element and element_id
                # and that they are valid?

        # other variables:
        if new_item_spec is not None:
            data_vars_dict.update(data_vars)

        # Dataset of new record:
        ds_to_add = Dataset(data_vars=new_record,
                            coords=coords_to_add,
                            **kwargs)
        # Merge new record and original dataset:
        #new_ds = xr.concat((self, ds_to_add), compat='equals')
        self.merge(ds_to_add, inplace='True', compat='no_conflicts')

    @property
    def variable_names(self):
        return self.Dataset.data_vars
