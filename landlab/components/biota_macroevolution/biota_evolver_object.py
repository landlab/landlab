"""BiotaEvolver base object.
"""


class BiotaEvolverObject(object):
    """
    """

    def __init__(self):
        """Initialize a BiotaEvolverObject.
        """
        self.record = {}

    def get_record_key_list(self, key):
        key_list = []
        for t in self.record.keys():
            key_list.append(self.record[t][key])
        return key_list

    @property
    def time__list(self):
        return sorted(self.record.keys())

    @property
    def time__last(self):
        return self.time__list[-1]

    @property
    def time__prior(self):
        if len(self.time__list) == 1:
            return None
        else:
            return self.time__list[-2]

    @property
    def number_of_records(self):
        return len(self.time__list)
