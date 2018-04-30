import pandas as pd


class Record(pd.DataFrame):
    """Pandas DataFrame to store attributes over time.

    Timestep is the index and 'attributes' are the values. The attributes are
    in turn dictionaries where the key is the of the attribute and values are
    the values of that key.
    """

    def __init__(self, *args, **kw):
        super(Record, self).__init__(*args, **kw)
        self.index.name = 'step'

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
