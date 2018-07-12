from numpy import array, where
from pandas import DataFrame, merge, concat


class Record(object):
    """Data structure to store attributes over time.

    Inherits from Pandas DataFrame.

    Timestep is the index and 'attributes' are the values. The attributes are
    in turn dictionaries where the key is the of the attribute and values are
    the values of that key.
    """

    def __init__(self, allow_time_duplicates=False, sort_time=True,
                 reset_index=True):
        """
        """
        self.DataFrame = DataFrame(columns=['time'], dtype='object')

        self.allow_time_duplicates = allow_time_duplicates

    def insert_time(self, time, data=None):
        """Append an entry to the record.

        Times previously entered will be overwritten.

        Parametersr,
        ----------
        time : integer or float

        data : dictionary
            A group of number-of-items long arrays....
        """
        if not isinstance(time, (int, float)):
            raise TypeError('*time* must be an integer or a float value.')

        if data == None:
            data = {}

        data['time'] = time


        if time in self.DataFrame.time:
            # Overwrite
            new_df = DataFrame([data])

            # Ensure data with lists are stored as lists.
            for k,v in data.items():

                d = new_df.loc[0, k]

                list_con = isinstance(d, list) == isinstance(data[k], list)

                if not list_con:
                    data[k] = [v]
                    new_df = DataFrame(data)

            new_df = self.DataFrame.append(new_df, sort=False)

            if self.allow_time_duplicates:
                self.DataFrame = new_df
            else:
                self.DataFrame = new_df[~new_df.duplicated('time',
                                                           keep='last')]

        else:
            self.DataFrame = self.DataFrame.append(data, ignore_index=True)

        self.DataFrame = self.DataFrame.sort_values('time')
        self.DataFrame.reset_index(inplace=True, drop=True)

    def modify_time(self, time, data=None):
        if time not in self.DataFrame.time.tolist():
            raise ValueError('*time* not in DataFrame.')

        data['time'] = time

        old_df = self.DataFrame
        new_df = DataFrame([data])

        i = where(old_df.time == time)[0]

        old_columns = old_df.columns.tolist()
        new_columns = new_df.columns.tolist()

        all_columns = old_df.columns.union(new_df.columns).tolist()

        new_data = {}

        print(merge(old_df, new_df, on='time', right_index=True))

        for col in all_columns:
            if col in new_columns:
#                old_df.at[i, col] = new_df.loc[0, col]
                new_data[col] = new_df.loc[0, col]
            else:
#                old_df.at[i, col] = old_df.loc[i, col]
                new_data[col] = old_df.loc[i, col]

        new_df = DataFrame([new_data])
        new_df = self.DataFrame.append(new_df, sort=False)
        self.DataFrame = new_df[~new_df.duplicated('time', keep='last')]

        self.DataFrame = self.DataFrame.sort_values('time')
        self.DataFrame.reset_index(inplace=True, drop=True)

#        print(new_data)

#        df = old_df.loc[:, ['time', 'a', 'v']]

#        print(df)
#
#        df.update(new_df)

#        self.DataFrame = df.reset_index(drop=True)

#        cols_to_use = old_df.columns.difference(new_df.columns)
#        print(cols_to_use)
#        new_df = merge(new_df, old_df, on='time', how='outer')

#        new_df = concat([old_df, new_df], sort=False, axis=0).sort_index().reset_index(drop=True)

#        print(new_df)
#        self.DataFrame = self.DataFrame.append(new_df, sort=False)
#        self.DataFrame.loc[i] = new_df

#        print(new_df)
#        new_df = self.DataFrame[i].append(DataFrame([data]), sort=False)
#        self.DataFrame.iloc[i] = self.DataFrame.iloc[i].combine_first(new_df)
#        print(new_df)
#        df = new_df.combine(old_df, lambda s1, s2: s1 if s1 and s2 else s2)
#        self.DataFrame = merge(old_df, new_df, how='right', on='time')
#        self.DataFrame = old_df.update(new_df)


    def get_time_prior_to_time(self, time):
        times = array(self.times)
        return times[times < time].max()

    @property
    def times(self):
        """Get a list of the times in the record."""
        return self.DataFrame.time.tolist()

    @property
    def time__earliest(self):
        """Get the earliest time in the record."""
        return self.DataFrame.time.min()

    @property
    def time__latest(self):
        """Get the latest time in the record."""
        return self.DataFrame.time.max()

    @property
    def time__prior(self):
        """Get the time prior to the latest time in the record."""
        if self.DataFrame.count == 1:
            return None
        else:
            return self.times[-2]
