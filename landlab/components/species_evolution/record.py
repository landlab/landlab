#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Structure to store data over time for SpeciesEvolver."""
from collections import OrderedDict

import numpy as np
from pandas import DataFrame


class Record(object):
    """Structure to store data over time for SpeciesEvolver.

    This object is intended to be used internally by the SpeciesEvolver
    component and its objects.

    The record attribute, ``_dict`` is a dictionary that stores the data of the
    record. The dictionary key, 'time' is created when the record is
    initialized. All other dictionary keys are variables of the record.

    The values of variables correspond to the time values elementwise. For
    example, if the time values are [0, 1, 2] and the values of var1 are
    [5, 6, 7] then the values of var 1 at time 0, 1, and 2 are 5, 6, and 7,
    respectively.

    A new record entry is created by advancing the record time with the method,
    ``advance_time``. Variables can be added to the record with the method,
    ``set_value``. Variable values not set at a time have a value of `nan`.
    """

    def __init__(self, initial_time=0):
        """Instantiate Record.

        Parameters
        ----------
        initial_time : float or int
            The initial time in the record.
        """
        self._dict = OrderedDict([("time", [initial_time])])

    @property
    def times(self):
        """The times stored in the record."""
        return self._dict["time"]

    @property
    def prior_time(self):
        """The penultimate time in the record.

        nan is returned if the record stores only one time step.
        """
        if len(self.times) < 2:
            return np.nan
        else:
            return sorted(self.times)[-2]

    @property
    def earliest_time(self):
        """The earliest time in the record."""
        return min(self.times)

    @property
    def latest_time(self):
        """The latest time in the record."""
        return max(self.times)

    @property
    def count_of_time_steps(self):
        """The count of record time steps."""
        return len(self._dict["time"])

    @property
    def variables(self):
        """The variables in the record."""
        variables = list(self._dict.keys())
        variables.remove("time")
        return variables

    @property
    def data_frame(self):
        """A Pandas DataFrame of the record."""
        return DataFrame(self._dict)

    def __len__(self):
        """The count of record time steps."""
        return self.count_of_time_steps

    def advance_time(self, dt):
        """Advance the time in the record by a time step duration.

        This method creates a new record entry. The new entry includes the time
        and all variables have a value of nan for this time. The time of the
        entry is the sum of the latest time of the record and the time step
        duration, ``dt``.

        Parameters
        ----------
        dt : float or int
            The time step duration to advance time in the record.
        """
        self._dict["time"].append(self.latest_time + dt)

        for var in self.variables:
            self._dict[var].append(np.nan)

    def get_value(self, var_name, time=np.nan):
        """Get the value of a variable.

        Parameters
        ----------
        var_name : string
            The name of the variable to get.
        time : float or int, optional
            The time of the variable to get. The latest time in the record is
            the default.

        Returns
        -------
        T
            The value of `var_name` that has the type, T. `nan` is returned if
            `var_name` does not exist in the record.
        """
        if var_name not in self.variables:
            return np.nan

        time = self._check_time(time)

        idx = self._get_time_index(time)
        return self._dict[var_name][idx]

    def set_value(self, var_name, value, time=np.nan):
        """Set the value of a variable.

        Parameters
        ----------
        var_name : string
            The name of the variable to set. A new variable is created if
            `var_name` does not exist in the record.
        value : T
            The value of `var_name`. The type, T should be an appropriate type
            for `var_name`.
        time : float or int, optional
            The time in the record to change a variable. The latest time in the
            record is the default.
        """
        time = self._check_time(time)

        if var_name not in self.variables:
            self._dict[var_name] = [np.nan] * (self.count_of_time_steps)

        idx = self._get_time_index(time)
        self._dict[var_name][idx] = value

    def increment_value(self, var_name, increase, time=np.nan):
        """Increment the value of a variable.

        Parameters
        ----------
        var_name : string
            The name of the variable to set. A new variable is created if
            `var_name` does not exist in the record.
        increase : T
            The value to increase `var_name`. The type, T should be an
            appropriate type for `var_name`.
        time : float or int, optional
            The time in the record to change a variable. The latest time in the
            record is the default.
        """
        time = self._check_time(time)
        idx = self._get_time_index(time)

        if var_name in self.variables:
            vals = self._dict[var_name]

            if np.isnan(vals[idx]):
                vals[idx] = increase
            else:
                vals[idx] += increase
        else:
            self._dict[var_name] = [np.nan] * self.count_of_time_steps
            self._dict[var_name][idx] = increase

    def _get_time_index(self, time):
        return np.where(np.array(self.times) == time)[0][0]

    def _check_time(self, time):
        if np.isnan(time):
            time = self.latest_time
        elif time not in self.times:
            raise ValueError("the time, {} not in record".format(time))

        return time
