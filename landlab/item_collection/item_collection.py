#! /usr/bin/env python
"""

"""
import numpy as np
from six import string_types

from pandas import DataFrame

from landlab.field import GroupError
from landlab import BAD_INDEX_VALUE

_LOCATIONS = {'node': 'number_of_nodes',
              'patch': 'number_of_patches',
              'link': 'number_of_links',
              'corner': 'number_of_corners',
              'face': 'number_of_faces',
              'cell': 'number_of_cells'}


class ItemCollection(object):
    """
    Datastructure to hold generic items that live on grid elements.

    This is a base class to contain the majority of core ItemCollection
    functionality. It inherits from the Pandas Dataframe.

    It requires that all data live on a grid element (e.g. node, link), with
    an element id that is less than the number of those elements. For example
    if you have a grid with only 100 links, no item can live at link 100 or
    link -3.

    """

    def __init__(self, grid, data=None, grid_element=None, element_id=None):
        """
        Parameters
        ----------
        grid : ModelGrid
        data : dictionary
            A group of number-of-items long arrays. All arrays must be the same
            length. The dictionary keys will become the column names of the
            ItemCollection
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
        >>> particle_volumes = [1, 2, 3, 4, 5, 6, 7]
        >>> particle_ages = [10, 11, 12, 13, 14, 15, 16]  
        >>> data = {'particle_ages': particle_ages,
        ...         'particle_volumes': particle_volumes}
        >>> ic = ItemCollection(grid, 
        ...                     data = data, 
        ...                     grid_element ='node', 
        ...                     element_id = element_id)
        
        
        Grid element can be any of the valid grid elements for the grid type
        provided.
        
        >>> grid_element = ['node', 'link', 'node', 'link', 'node', 'link', 'node']
        >>> ic = ItemCollection(grid, 
        ...                     data = data, 
        ...                     grid_element = grid_element, 
        ...                     element_id = element_id)
        
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
        

        # add grid element and element ID to data frame
        data['grid_element'] = grid_element
        data['element_id'] = element_id

        # initialized the PD dataframe now that we've done checks
        self.DataFrame = DataFrame(data)

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
                    bad_inds = selected_elements[less_than_zero] == BAD_INDEX_VALUE
                    if sum(bad_inds) != sum(less_than_zero):
                        raise ValueError(('An item residing at ' + at + ' has '
                                          'an element id below zero that is '
                                          'not BAD_INDEX_VALUE. This is not '
                                          'permitted.'))
                        
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
        """Add a new variable to the ItemCollection."""
        if isinstance(variable, string_types):
            if np.array(values).size == self.DataFrame.shape[0]:
                self.DataFrame[variable] = values
            else:
                raise ValueError(('Values passed to add_variable must have '
                                  'the same size as the current ItemCollection.'))
        else:
            raise ValueError(('Variable name passed to add_variable must be of '
                              'type string.'))
            
    def add_items(self, data, grid_element, element_id):
        """Add new items to the ItemCollection"""
        
        self._check_sizes(data)
        
        grid_element = self._check_grid_element_and_id(grid_element, element_id)

        # add grid element and element ID to data frame
        data['grid_element'] = grid_element
        data['element_id'] = element_id

        new_data = DataFrame(data)
        
        self.DataFrame = self.DataFrame.append(new_data)
        
        # check that element IDs do not exceed number of elements on this grid
        self._check_element_id_values()

    def calc_aggregate_sum(self, var, at='node'):
        """Sum of variable at grid elements.

        Parameters
        ----------
        var : str
            Column name of variable to sum
        at : str, optional
            Name of grid element at which to sum. Default is "node".

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
        >>> element_id = [0, 0, 1, 1, 2, 3, 5]
        >>> particle_volumes = [1, 2, 3, 4, 5, 6, 7]
        >>> particle_ages = [10, 11, 12, 13, 14, 15, 16]  
        >>> grid_element = 'node'
        >>> data = {'particle_ages': particle_ages,
        ...         'particle_volumes':particle_volumes}
        >>> ic = ItemCollection(grid, 
        ...                     data = data, 
        ...                     grid_element ='node', 
        ...                     element_id = element_id)
        >>> s = ic.calc_aggregate_sum('particle_ages')
        >>> print(s)
        [ 21.  25.  14.  15.  nan  16.  nan  nan  nan]
        >>> len(s) == grid.number_of_nodes
        True
        
        """
        # select those items located on the correct type of element,
        # group by element_id and sum.
        vals = self.DataFrame.loc[self.DataFrame['grid_element'] == at].groupby('element_id').sum()
        
        # create a nan array that we will fill with the results of the sum
        # this should be the size of the number of elements, even if there are
        # no items living at some grid elements.
        out = np.nan * np.ones(self._grid[at].size)
        
        # put the values of the specified variable in to the correct location
        # of the out array.
        out[vals.index] = vals[var]

        return out
    
    def calc_aggregate_mean(self, var, at='node', nanmean=False):
        """Mean of variable at grid elements.

        Parameters
        ----------
        var : str
            Column name of variable to sum
        at : str, optional
            Name of grid element at which to sum. Default is "node".

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
        >>> element_id = [0, 0, 0, 0, 2, 3, 5]
        >>> particle_volumes = [1, 2, 3, 4, 5, 6, 7]
        >>> particle_ages = [10, 11, 12, 13, 14, 15, 16]  
        >>> grid_element = ['node', 'link', 'node', 'link', 'node', 'link', 'node']
        >>> data = {'particle_ages': particle_ages,
        ...         'particle_volumes':particle_volumes}
        >>> ic = ItemCollection(grid, 
        ...                     data = data, 
        ...                     grid_element = grid_element, 
        ...                     element_id = element_id)
        >>> s = ic.calc_aggregate_mean('particle_ages', at = 'link')
        >>> print(s)
        [ 12.  nan  nan  15.  nan  nan  nan  nan  nan  nan  nan  nan]
        >>> len(s) == grid.number_of_links
        True
        
        If ingoring nan's in item data is necessary, there is a nanmean option.
        
        
        
        """
        # select those items located on the correct type of element,
        # group by element_id and sum.
        if nanmean:
            pass
        else:
            vals = self.DataFrame.loc[self.DataFrame['grid_element'] == at].groupby('element_id').mean()
        
        # create a nan array that we will fill with the results of the sum
        # this should be the size of the number of elements, even if there are
        # no items living at some grid elements.
        out = np.nan * np.ones(self._grid[at].size)
        
        # put the values of the specified variable in to the correct location
        # of the out array.
        out[vals.index] = vals[var]

        return out

    def calc_aggregate_min(self, var, at='node'):
        """ """
        pass

    def calc_aggregate_max(self, var, at='node'):
        """ """
        pass

    def calc_aggregate_median(self, var, at='node'):
        """ """
        pass

    def calc_aggregate_var(self, var, at='node'):
        """ """
        pass
    
    def calc_aggregate_std(self, var, at='node'):
        """ """
        pass