"""BiotaEvolver base objects.
"""

class Record(dict):
    """ Dictionary-like data structure to store attributes over time.

    Time is the key and 'attributes' are the values. The attributes are in turn
    dictionaries where the key is the of the attribute and values are the
    values of that key.
    """

    def at_time(self, time):
        record_at_time = {}
        for t, data in self.items():
            if t <= time:
                record_at_time[t] = data

        return record_at_time

    def values_of_attribute(self, attribute):
        values = []
        for t in self.times:
            if attribute in self[t].keys():
                values.append(self[t][attribute])
        return values

    @property
    def count(self):
        return len(self.keys())

    @property
    def times(self):
        return sorted(self.keys())

    @property
    def time__first(self):
        return self.times[0]

    @property
    def time__last(self):
        return self.times[-1]

    @property
    def time__prior(self):
        if self.count == 1:
            return None
        else:
            return self.times[-2]

    @property
    def number_of_species(self):
        return len(self.values_of_attribute('species'))


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
