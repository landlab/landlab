from numpy import array, nan
from pandas import DataFrame


class RecordCollection(object):
    """Data structure that stores data of time steps.

    Data is stored in the *DataFrame* property that is a Pandas DataFrame. Each
    row is a record of data at a time.

    The first column, *time* is entered when data is added or modified with
    RecordCollection methods. Only 1 record can exist for a time. A record can
    be modified after its created.

    The data variables of records are the columns of a RecordCollection. An
    variable may or may not have a value for every record. A value of NaN
    indicates the variable has no value at the time of the record.

    This object facilitates some DataFrame functionality. Access the methods of
    *DataFrame* to use other Pandas functionality (e.g., filtering, sorting).
    """

    def __init__(self, maintain_insert_order=False):
        """
        Parameters
        ----------
        maintain_insert_order : boolean, optional
            When 'True', time entries are stored in the order that they were
            added/modified. When 'False' (default), time entries are sorted by
            time ascendingly. Sorting by time in the 'False' case occurs when
            entires are added and modified.
        """
        self.DataFrame = DataFrame(columns=['time'], dtype='object')

        self.maintain_insert_order = maintain_insert_order

    def _reset_index(self):
        if not self.maintain_insert_order:
            self.DataFrame = self.DataFrame.sort_values('time')
        self.DataFrame.reset_index(inplace=True, drop=True)

    def insert_time(self, time, data=None):
        """Insert a record.

        Only times that do not exist in the RecordCollection may be entered
        using this method.

        Parameters
        ----------
        time : integer or float
            The time of the record to insert.
        data : dictionary
            The data that will be inserted in the record at *time*. The
            dictionary keys will be the column labels. The dictionary values
            will be the values of the corresponding keys/columns.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()

        Insert a record for time 0 with no data.
        >>> rc.insert_time(0)

        Insert a record for time 1 with variable label, 'test'.
        >>> rc.insert_time(100000, {'max_of_elevation': 1000})

        Insert a record for time 2 with variable label, 'test'.
        >>> rc.insert_time(200000, {'max_of_elevation': 1100,
                                    'a_list_variable': [1, 2, 3]})

        >>> print(rc.DataFrame)
             time  max_of_elevation a_list_variable
        0       0               NaN             NaN
        1  100000            1000.0             NaN
        2  200000            1100.0       [1, 2, 3]
        """
        if time in self.DataFrame.time.tolist():
            raise ValueError('*time* already exists in DataFrame.')

        if data == None:
            data = {}

        data['time'] = time

        if time in self.times:
            # Overwrite
            new_df = DataFrame([data])

            # Ensure data with lists are stored as lists.
            for k,v in data.items():

                d = new_df.loc[0, k]

                list_con = isinstance(d, list) == isinstance(data[k], list)

                if not list_con:
                    data[k] = [v]
                    new_df = DataFrame(data)

            self.DataFrame = new_df[~new_df.duplicated('time', keep='last')]

        else:
            self.DataFrame = self.DataFrame.append(data, ignore_index=True)

        self._reset_index()

    def modify_time(self, time, data=None):
        """Modify an existing record.

        Parameters
        ----------
        time : integer or float
            The time of the record to modify.
        data : dictionary
            The data that will be inserted in the record at *time*. The
            dictionary keys will be the column labels. The dictionary values
            will be the values of the corresponding keys/columns. Columns no

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(100000, data={'max_of_elevation': 400})
        >>> print(rc.DataFrame)
             time  max_of_elevation
        0  100000             400.0

        >>> rc.modify_time(100000, data={'a': 500})
        >>> print(rc.DataFrame)
             time  max_of_elevation      a
        0  100000             400.0  500.0
        """
        if time not in self.DataFrame.time.tolist():
            raise ValueError('*time* not in RecordCollection.')

        # Get the row of *time*.
        old_df = self.DataFrame[self.DataFrame.time == time]
        old_df.reset_index(inplace=True, drop=True)

        # Create a DataFrame of the data to place at *time*.
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
                df = df[df.time.notnull()]
                df = self.DataFrame.append(df, sort=False)

                # Remove the copy of the row at time with unmodified data.
                self.DataFrame = df[~df.duplicated('time', keep='last')]

        self._reset_index()

    def remove_time(self, time):
        """Remove a record.

        Parameters
        ----------
        time : integer or float
            The time of the record to remove.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.remove_time(1)
        >>> print(rc.DataFrame)
          time
        0    0
        1    2
        """
        index = self.DataFrame.index[self.DataFrame['time'] == time]
        self.DataFrame.drop(index, inplace=True)

        self._reset_index()

    def get_value(self, time, variable):
        """Get the value of a variable at a time in the record.

        Parameters
        ----------
        time : integer or float
            The time of the record to get.
        variable : string
            The label of the variable to get.

        Returns
        -------
        object
            The value of *variable* at *time*. The type of the returned object
            is dependent on the type of the variable value.

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

        time_index = self.times.index(time)
        return self.DataFrame.loc[time_index, variable]

    def get_time_prior_to_time(self, time):
        times = array(self.times)
        return times[times < time].max()

    @property
    def number_of_time_steps(self):
        """Get the number of time steps in the record.

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
    def times(self):
        """Get a list of the times in the record.

        Examples
        --------
        Get times in the default, time ascending order.

        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(10000)
        >>> rc.insert_time(1)
        >>> rc.insert_time(1.5)
        >>> rc.times
        [1, 1.5, 10000]

        Get times in the order they were inserted.

        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection(maintain_insert_order=True)
        >>> rc.insert_time(1000000)
        >>> rc.insert_time(1)
        >>> rc.insert_time(1.5)
        >>> rc.times
        [10000, 1, 1.5]
        """
        return self.DataFrame.time.tolist()

    @property
    def time__earliest(self):
        """Get the earliest time in the record.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.time__earliest
        0
        """
        return self.DataFrame.time.min()

    @property
    def time__latest(self):
        """Get the latest time in the record.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.time__latest
        2
        """
        return self.DataFrame.time.max()

    @property
    def time__prior(self):
        """Get the time prior to the latest time in the record.

        Examples
        --------
        >>> from landlab.components.species_evolution import RecordCollection
        >>> rc = RecordCollection()
        >>> rc.insert_time(0)
        >>> rc.insert_time(1)
        >>> rc.insert_time(2)
        >>> rc.time__prior
        1
        """
        if self.number_of_time_steps < 2:
            return nan
        else:
            return sorted(self.times)[-2]
