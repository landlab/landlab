#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 12:07:40 2019

@author: njlyons
"""
from collections import defaultdict, OrderedDict

import numpy as np
from pandas import DataFrame

class Record(object):

    def __init__(self, initial_time=0):
        self._dict = OrderedDict([('time', [initial_time])])

    @property
    def times(self):
        return self._dict['time']

    @property
    def variables(self):
        return list(self._dict.keys())

    @property
    def prior_time(self):
        if len(self.times) < 2:
            return self.times[0]
        else:
            return sorted(self.times)[-2]

    @property
    def latest_time(self):
        return np.nanmax(self.times)

    @property
    def record_count(self):
        return len(self._dict['time'])

    @property
    def dataframe(self):
        return DataFrame(self._dict)

    def iterate_time(self, dt):
        self._dict['time'].append(self.latest_time + dt)

    def insert_add_on(self, add_on):
        for key, value in add_on.items():
            if key not in self._dict.keys():
                self._dict[key] = [np.nan] * (self.record_count - 1)

            self._dict[key].append(value)

    def get_empty_add_on(self):
        return defaultdict(int)

    def change_value(self, time, var_name, new_value):
        idx = np.where(np.array(self.times) == time)[0][0]

        if len(self._dict[var_name]) - 1 < idx:
            self._dict[var_name].append(new_value)
        else:
            self._dict[var_name][idx] = new_value

