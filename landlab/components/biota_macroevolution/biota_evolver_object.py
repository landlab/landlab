"""BiotaEvolver base objects.
"""

import pandas as pd

class Record(pd.DataFrame):
    """ Data frame-like object to store attributes over time.

    Time is the index and 'attributes' are the values. The attributes are in turn
    dictionaries where the key is the of the attribute and values are the
    values of that key.
    """

    def __init__(self, *args, **kw):
        super(Record, self).__init__(*args, **kw)
        self.index.name = 'time'

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


class BiotaEvolverObject(object):
    """
    """

    def __init__(self):
        """Initialize a BiotaEvolverObject.
        """
        self._identifier = None
        self.record = Record()

    @property
    def identifier(self):
        return self._identifier
