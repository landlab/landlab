from numpy import array
from pandas import DataFrame


class Record(DataFrame):
    """Data structure to store attributes over time.

    Inherits from Pandas DataFrame.

    Timestep is the index and 'attributes' are the values. The attributes are
    in turn dictionaries where the key is the of the attribute and values are
    the values of that key.
    """

    def __init__(self, *args, **kw):
        if 'columns' in kw:
            kw['columns'].insert(0, 'time')
        else:
            kw['columns'] = ['time']
        super(Record, self).__init__(*args, **kw)

    def append_entry(self, time, dictionary=None):
        if dictionary == None:
            dictionary = {}

        dictionary['time'] = time

        # Handle unincluded columns.
        unset_cols = list(set(self.columns.tolist()) - set(dictionary.keys()))
        for c in unset_cols:
            dictionary[c] = None

        # Handle new columns.
        unset_cols = list(set(dictionary.keys()) - set(self.columns.tolist()))
        for c in unset_cols:
            self.loc[:len(self), c] = None

        i = len(self)

        for k, v in dictionary.items():
            self.loc[i, k] = v

    def get_time_prior_to_time(self, time):
        times = array(self.times)
        return times[times < time].max()

    @property
    def times(self):
        return sorted(self.time.tolist())

    @property
    def time__earliest(self):
        return self.time.min()

    @property
    def time__latest(self):
        return self.time.max()

    @property
    def time__prior(self):
        if self.count == 1:
            return None
        else:
            return self.times[-2]
