from numpy import array, nan
from pandas import DataFrame


class RecordCollection(object):
    """Data structure that stores data of model time steps.

    Data is stored in the *DataFrame* property that is a Pandas DataFrame. Each
    row is a record of data at a model time.

    The first column, *model__time* is entered when data is added or modified
    with RecordCollection methods. Only 1 record can exist for a model time,
    although a record can be modified and removed after it is created.

    The data variables of records are the columns of a RecordCollection
    DataFrame. A variable may or may not have a value for every record. A value
    of NaN indicates the variable has no value at the model time of the record.

    RecordCollection facilitates some functionality of Pandas DataFrames.
    Additional Pandas functionality can be used by calling the methods of
    *DataFrame* directly.
    """

    def __init__(self, maintain_insert_order=False):
        """
        Parameters
        ----------
        maintain_insert_order : boolean, optional
            When 'True', model time entries are stored in the order that they
            were added/modified. When 'False' (default), time entries are
            sorted by time ascendingly. Sorting by time in the 'False' case
            occurs when entires are added and modified.
        """
        self.DataFrame = DataFrame(columns=['model__time'], dtype='object')

        self.maintain_insert_order = maintain_insert_order

    def _reset_index(self):
        if not self.maintain_insert_order:
            self.DataFrame = self.DataFrame.sort_values('model__time')
        self.DataFrame.reset_index(inplace=True, drop=True)

    def insert_time(self, model__time, data=None, clobber=False):
        """Insert a record.

        Only model times that do not exist in the RecordCollection may be
        entered using this method.

        Parameters
        ----------
        model__time : integer or float
            The model time of the record to insert.
        data : dictionary
            The data that will be inserted in the record at *model__time*. The
            dictionary keys will be the column labels. The dictionary values
            will be the values of the corresponding keys/columns.
        clobber : boolean

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()

        Insert a record for model time 0 with no data.
        >>> rc.insert_time(0)

        Insert a record for time 1 with variable label, 'test'.
        >>> rc.insert_time(100000, {'max_of_elevation': 1000})

        Insert a record for time 2 with variable label, 'test'.
        >>> rc.insert_time(200000, {'max_of_elevation': 1100,
                                    'a_list_variable': [1, 2, 3]})

        >>> print(rc.DataFrame)
          model__time  max_of_elevation a_list_variable
        0           0               NaN             NaN
        1      100000            1000.0             NaN
        2      200000            1100.0       [1, 2, 3]
        """
        if model__time in self.DataFrame.model__time.tolist():
            if clobber:
                self.remove_time(model__time)
            else:
                raise ValueError('the model time, {} already exists in '
                                 'DataFrame'.format(model__time))

        # Prepare data dictionary.

        if data == None:
            data = {}

        data['model__time'] = model__time

        if model__time in self.model__times:
            # Overwrite

            new_df = DataFrame([data])

            # Ensure data with lists are stored as lists.
            for variable_label, value in data.items():

                d = new_df.loc[0, variable_label]

                list_con = isinstance(d, list) == isinstance(data[variable_label], list)

                if not list_con:
                    data[variable_label] = [value]
                    new_df = DataFrame(data)

            self.DataFrame = new_df[~new_df.duplicated('model__time',
                                                       keep='last')]

        else:
            self.DataFrame = self.DataFrame.append(data, ignore_index=True,
                                                   sort=True)

        self._reset_index()

    def modify_time(self, model__time, data=None):
        """Modify an existing record.

        Parameters
        ----------
        model__time : integer or float
            The model time of the record to modify.
        data : dictionary
            The data that will be inserted in the record at *model__time*. The
            dictionary keys will be the column labels. The dictionary values
            will be the values of the corresponding keys/columns. Columns no

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(100000, data={'max_of_elevation': 400})
        >>> print(rc.DataFrame)
          model__time  max_of_elevation
        0      100000             400.0

        >>> rc.modify_time(100000, data={'a_variable': 500})
        >>> print(rc.DataFrame)
          model__time  max_of_elevation  a_variable
        0      100000             400.0       500.0
        """
        if model__time not in self.DataFrame.model__time.tolist():
            raise ValueError('the model time, {} is not in '
                             'RecordCollection'.format(model__time))

        # Get the row of *time*.
        old_df = self.DataFrame[self.DataFrame.model__time == model__time]
        old_df.reset_index(inplace=True, drop=True)

        # Create a DataFrame of the data to place at *model__time*.
        new_df = DataFrame([data])

        # Get columns only in old DataFrame.
        old_cols = old_df.columns.difference(new_df.columns)

        # Get columns shared in old and new DataFrames.
        intersected_cols = new_df.columns.intersection(old_df.columns).tolist()

        # Get columns not in the old DataFrame.
        new_cols = new_df.columns.difference(old_df.columns).tolist()

        adjusted_cols = [intersected_cols, new_cols]

        for cols in adjusted_cols:
            if cols:
                # Merge old and new DataFrames.
                df = old_df[old_cols].merge(new_df[cols], left_index=True,
                           right_index=True, how='outer')
                df = df[df.model__time.notnull()]
                df = self.DataFrame.append(df, ignore_index=True, sort=False)

                # Remove the copy of the row at *model__time* with unmodified
                # data.
                self.DataFrame = df[~df.duplicated('model__time', keep='last')]

        self._reset_index()

    def remove_time(self, model__time):
        """Remove a record.

        Parameters
        ----------
        model__time : integer or float
            The model time of the record to remove.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.remove_time(1)
        >>> print(rc.DataFrame)
          model__time
        0           0
        1           2
        """
        time_mask = self.DataFrame['model__time'] == model__time
        index = self.DataFrame.index[time_mask]
        self.DataFrame.drop(index, inplace=True)

        self._reset_index()

    def get_value(self, model__time, variable):
        """Get the value of a variable at a model time in the record.

        Parameters
        ----------
        model__time : integer or float
            The model time of the record to get.
        variable : string
            The label of the variable to get.

        Returns
        -------
        object
            The value of *variable* at *model__time*. The type of the returned
            object is dependent on the type of the variable value.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0, data={'variable_a': 1.5})
        >>> rc.insert_time(1, data={'variable_a': 2.0, 'variable_b': [3, 2]})
        >>> rc.get_value(0, 'variable_a')
        1.5
        >>> rc.get_value(1, 'variable_b')
        [3, 2]
        >>> rc.get_value(0, 'variable_b')
        nan
        """
        if variable not in self.DataFrame.columns:
            raise KeyError("the variable '{}' is not in the "
                           "RecordCollection".format(variable))

        time_index = self.model__times.index(model__time)
        return self.DataFrame.loc[time_index, variable]

    def get_model_time_prior_to_time(self, input_time):
        """Get the model time most immediately prior to an input time.

        The *input_time* does not need to be a time in *model__time* of
        *DataFrame*. Whichever *model__time* value is both nearest and prior to
        *input_time* is returned.

        Parameters
        ----------
        input_time : integer or float
            A time subsequent to the model time to get.

        Returns
        -------
        integer or float
            The model time prior to *input_time*. The type of the returned time
            is dependent on the number type initially inserted into the
            DataFrame.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.get_model_time_prior_to_time(2)
        1
        """
        times = array(self.model__times)
        input_time_greater_than_model_time = times < input_time

        # Return the model time most immediately prior to *input_time* if there
        # is such a model time.
        if any(input_time_greater_than_model_time):
            return times[input_time_greater_than_model_time].max()
        else:
            return nan

    @property
    def number_of_time_steps(self):
        """Get the number of model time steps in the record.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.number_of_time_steps
        0
        >>> rc.insert_time(1)
        >>> rc.insert_time(10)
        >>> rc.number_of_time_steps
        2
        """
        return len(self.DataFrame)

    @property
    def model__times(self):
        """Get a list of the model times in the record.

        Examples
        --------
        Get model times in the default, time ascending order.

        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(10000)
        >>> rc.insert_time(1)
        >>> rc.insert_time(1.5)
        >>> rc.model__times
        [1, 1.5, 10000]

        Get times in the order they were inserted.

        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection(maintain_insert_order=True)
        >>> rc.insert_time(10000)
        >>> rc.insert_time(1)
        >>> rc.insert_time(1.5)
        >>> rc.model__times
        [10000, 1, 1.5]
        """
        return self.DataFrame.model__time.tolist()

    @property
    def model__earliest_time(self):
        """Get the earliest model time in the record.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.model__earliest_time
        0
        """
        return self.DataFrame.model__time.min()

    @property
    def model__latest_time(self):
        """Get the model latest time in the record.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.model__latest_time
        2
        """
        return self.DataFrame.model__time.max()

    @property
    def model__prior_time(self):
        """Get the model time prior to the latest model time in the record.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.model__prior_time
        1
        """
        if self.number_of_time_steps < 2:
            return nan
        else:
            return sorted(self.model__times)[-2]
