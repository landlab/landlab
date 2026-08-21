#! /usr/bin/env python

# # Base class for a grid-based Landlab model
#
# This code defines LandlabModel, a Python class that is designed to make it easier
# to create a standalone model code using Landlab. The model developer writes
# a class that inherits from LandlabModel and adds the functionality needed to
# implement their model. LandlabModel provides code to handle formatted user input,
# in the form of either a Python dictionary or the name of a yaml-format input
# file (given as a string). The LandlabModel __init__() method will combine the
# user inputs with a set of default parameter values defined in the derived
# class header (for parameters whose value has not been specified by the user).
# For model execution, the user simply needs to override the built-in update()
# method. LandlabModel calls this via a built-in run() method (which runs the model
# from start to finish) and a built-in update_until() method (which calls
# update() until the either the run is complete or it is time to pause and
# generate output).
#
# *(Greg Tucker, University of Colorado Boulder)*
#

import sys

import numpy as np


def verify_input_file_and_load_params(input_file: str) -> dict:
    """
    If the named file exists, read it and return a dict of parameters.
    Otherwise, raise FileNotFoundError.

    Note: this is meant to be a temporary fn because currently
    load_params does not handle missing files gracefully. Once it does,
    then this could simply be replaced by load_params.

    Parameters
    ----------
    input_file : str
        name of input file containing parameters in yaml format

    Returns
    -------
    dict containing parameter names and values
    """

    from landlab import load_params

    try:
        f = open(input_file)
        f.close()
        params = load_params(input_file)
        return params
    except FileNotFoundError:
        raise


def merge_user_and_default_params(user_params: dict, default_params: dict) -> None:
    """
    Merge default parameters into the user-parameter dictionary, adding
    defaults where user values are absent.

    Parameters
    ----------
    user_params : dict
        dict containing names and values of user-defined parameters
    default_params : dict
        dict containing all parameter names and their default values

    Examples
    --------
    >>> u = {"a": 1, "d": {"da": 4}, "e": 5, "grid": {"RasterModelGrid": []}}
    >>> d = {"a": 2, "b": 3, "d": {"db": 6}, "grid": {"HexModelGrid": []}}
    >>> merge_user_and_default_params(u, d)
    >>> u["a"]
    1
    >>> u["b"]
    3
    >>> u["d"]
    {'da': 4, 'db': 6}
    >>> u["grid"]
    {'RasterModelGrid': []}
    """
    for k in default_params.keys():
        if k not in user_params.keys():
            user_params[k] = default_params[k]
        elif (
            isinstance(user_params[k], dict)
            and isinstance(default_params[k], dict)
            and k != "grid"
        ):
            merge_user_and_default_params(user_params[k], default_params[k])


def get_or_create_node_field(grid, name: str, dtype: str = "float64") -> np.ndarray:
    """
    Get handle to a grid field if it exists, otherwise create it.

    Parameters
    ----------
    grid : Landlab grid object
        grid to operate on
    name : str
        name of field to get or create
    dtype : str (default "float64")
        data type for field if created

    Returns
    -------
    ndarray : the field as an array

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> g = RasterModelGrid((3, 3))
    >>> f = g.add_field("field1", np.arange(9), at="node")
    >>> get_or_create_node_field(g, "field1")
    array([0, 1, 2, 3, 4, 5, 6, 7, 8])
    >>> get_or_create_node_field(g, "field2")
    array([0., 0., 0., 0., 0., 0., 0., 0., 0.])
    >>> "field2" in g.at_node.keys()
    True
    """
    try:
        return grid.at_node[name]
    except KeyError:
        return grid.add_zeros(name, at="node", dtype=dtype, clobber=True)


def read_arrays_from_files(params):
    """
    Given a parameter dictionary params, identify any items that are dictionaries
    with the key "_filepath" and replace the item with the contents of the
    specified file path.

    Parameters
    ----------
    params : dict
        dict to analyze

    Returns:
        dict : the modified dict

    Examples
    --------
    >>> import numpy as np
    >>> np.save("test1", 0.1 * np.arange(3))
    >>> np.save("test2", np.arange(4).reshape((2, 2)) / 2)
    >>> np.save("test3", 1.5 * np.arange(3).reshape((3, 1)))
    >>> p = {
    ...     "a": 123,
    ...     "b": {"c": 456, "d": {"_filepath": "test1.npy"}},
    ...     "e": {"_filepath": "test2.npy"},
    ...     "f": {"_filepath": "test3.npy"},
    ... }
    >>> p = read_arrays_from_files(p)
    >>> p["b"]["d"]
    array([0. , 0.1, 0.2])
    >>> p["e"]
    array([[0. , 0.5],
           [1. , 1.5]])
    >>> p["f"]
    array([[0. ],
           [1.5],
           [3. ]])
    """
    for item in params:
        if isinstance(params[item], dict):
            if "_filepath" in params[item]:
                params[item] = np.load(params[item]["_filepath"])
            else:
                params[item] = read_arrays_from_files(params[item])
    return params


def _get_pause_time_list_and_next(time_info, clock_dict, no_first_pause=False):
    """
    Given a float or iterable as ``time_info``, return a list of times to pause the
    simulation to perform an action.

    Return a list of times at which to pause the simulation, including an item at
    the end that is after the termination of the run so that there is always an item
    to be popped.

    Parameters
    ----------
    time_info : float or list of float
        Interval for pausing (if float) or list of individual times
    clock_dict : dict
        Contains ``start`` and ``stop`` as float items, with ``stop`` > ``start``
    no_first_pause : bool
        Flag indicating whether to include the start time as pause time (default False)

    Returns
    -------
    list : list of simulation times at which to pause for a given action
    float : the next time at which to pause

    Examples
    --------
    >>> cldict = {"start": 0.0, "step": 1.0, "stop": 4.0}
    >>> _get_pause_time_list_and_next(1.0, cldict)
        ([1.0, 2.0, 3.0, 4.0], 0.0)
    >>> _get_pause_time_list_and_next(1.0, cldict, no_first_pause=True)
    ([2.0, 3.0, 4.0], 1.0)
    >>> _get_pause_time_list_and_next([0.0, 0.5, 2.0, 4.0], cldict)
    ([0.5, 2.0, 4.0], 0.0)
    """
    if isinstance(time_info, float) or isinstance(time_info, int):
        start = clock_dict["start"]
        if no_first_pause:
            start += time_info
        pause_times = list(np.arange(start, clock_dict["stop"] + time_info, time_info))
    elif isinstance(time_info, list):
        pause_times = time_info.copy()
    else:
        print("time_info must be of type float or list, not", type(time_info))
        raise (TypeError)
    next_pause = pause_times.pop(0)
    return pause_times, next_pause


class LandlabModel:
    """
    Base class for a generic Landlab grid-based model.

    Examples
    --------
    >>> from landlab.utils import LandlabModel
    >>> class MyModel(LandlabModel):
    ...     def __init__(self, params):
    ...         super().__init__(params=params)
    ...
    >>> p = {"grid": {"source": "create"}}
    >>> p["grid"]["create_grid"] = {
    ...     "RasterModelGrid": {"shape": (4, 5), "xy_spacing": 2.0}
    ... }
    >>> model = MyModel(params=p)
    >>> model.grid.shape
    (4, 5)
    """

    # Default parameters, to be overridden in derived classes
    DEFAULT_PARAMS = {
        "grid": {
            "source": "create",
            "create_grid": {
                "RasterModelGrid": [
                    {"shape": (5, 5)},
                    {"xy_spacing": 1.0},
                ],
            },
        },
        "clock": {"start": 0.0, "stop": 2.0, "step": 1.0},
        "output": {
            "plot_times": 10.0,  # float or list
            "save_times": 10.0,  # float or list
            "report_times": 1.0,  # float or list
            "save_path": "model_output",
            "clobber": True,
            "fields": None,
            "plot_to_file": True,
        },
    }

    def __init__(self, params: dict = {}, input_file: str = "") -> dict:
        """
        Initialize the model.

        Parameters
        ----------
        params : dict
            Dictionary containing names and values of model parameters
        input_file : str
            Name of yaml-format file containing names and values

        Notes
        -----
        User should pass either params or input_file. If input_file is
        passed, params will be ignored.
        """
        if len(input_file) > 0:
            params = verify_input_file_and_load_params(input_file)
        self.params = params
        merge_user_and_default_params(self.params, self.DEFAULT_PARAMS)
        read_arrays_from_files(self.params)
        self.setup_grid(self.params["grid"])
        self.setup_for_output(self.params)
        self.setup_run_control(self.params["clock"])

    def setup_grid(self, grid_params: dict) -> None:
        """Load or create the grid.

        Parameters
        ----------
        grid_params : dict
            Dictionary containing parameters related grid setup.

        Notes
        -----
        Must include an item "source" for which the valid values are
        "create" (create a new grid), "file" (read a grid from file), or
        "grid_object" (indicating that a grid object is included
        directly in the parameter dictionary).

        If "create", then there must be an item "create_grid" that
        contains a dict in which the key is the name of the grid type
        ("RasterModelGrid", "HexModelGrid") and the value is a dict
        containing the names and values for the grid object's
        parameters (such as a tuple for "shape", etc.)

        If "file", then there must be an item "grid_file_name" that
        contains the file name as a string.

        If "grid_object", then there must be an item called "grid_object"
        containing the grid object!

        Examples
        --------
        >>> p = {"grid": {"source": "create"}}
        >>> p["grid"]["create_grid"] = {
        ...     "RasterModelGrid": {"shape": (4, 5), "xy_spacing": 2.0}
        ... }
        >>> sim = LandlabModel(params=p)
        >>> sim.grid.shape
        (4, 5)
        >>> from landlab.io.native_landlab import save_grid
        >>> save_grid(sim.grid, "test.grid", clobber=True)
        >>> p = {"grid": {"source": "file", "grid_file_name": "test.grid"}}
        >>> sim = LandlabModel(params=p)
        >>> sim.grid.shape
        (4, 5)
        >>> from landlab import RasterModelGrid
        >>> p = {"grid": {"source": "grid_object"}}
        >>> p["grid"]["grid_object"] = RasterModelGrid((3, 3))
        >>> sim = LandlabModel(params=p)
        >>> sim.grid.shape
        (3, 3)
        >>> from numpy.testing import assert_raises
        >>> p["grid"]["grid_object"] = "spam"
        >>> assert_raises(ValueError, LandlabModel, p)
        grid_object must be a Landlab grid.
        """
        from landlab import ModelGrid
        from landlab import create_grid
        from landlab.io.native_landlab import load_grid

        if grid_params["source"] == "create":
            # print("Create here...")
            # print(grid_params["create_grid"])
            self.grid = create_grid(grid_params, section="create_grid")
        elif grid_params["source"] == "file":
            self.grid = load_grid(grid_params["grid_file_name"])
        elif grid_params["source"] == "grid_object":
            if isinstance(grid_params["grid_object"], ModelGrid):
                self.grid = grid_params["grid_object"]
            else:
                print("grid_object must be a Landlab grid.")
                raise ValueError

    def setup_for_output(self, params: dict) -> None:
        """
        Setup variables for control of plotting and saving.

        Parameters
        ----------
        params : dict
            Parameter dictionary. Must include a key ``output`` with a dictionary
        containing values for ``plot_times``, ``save_times``, and ``report_times``.
        Each of these should be either a ``float`` or a ``list``. If a list, the value
        is interpreted as a list of model times for plotting, saving, or reporting.
        If a single float, the value is interpreted as the (regular) time
        interval (in model time) for plotting, saving, or reporting.
            Should also contain a key ``clock`` as a dictionary that has values
        for ``start`` and ``stop``.

        Notes
        -----
        The "format" parameter can be "grid" (native Landlab grid format), "netcdf",
        or "vtk". The default is "grid".
        """
        op_params = params["output"]
        clock_params = params["clock"]

        self.plot_times, self.next_plot = _get_pause_time_list_and_next(
            op_params["plot_times"], clock_params
        )
        self.save_times, self.next_save = _get_pause_time_list_and_next(
            op_params["save_times"], clock_params, no_first_pause=True
        )
        self.report_times, self.next_report = _get_pause_time_list_and_next(
            op_params["report_times"], clock_params
        )

        self.ndigits_for_save_files = int(np.ceil(np.log10(len(self.save_times) + 1)))
        self.save_num = 0  # current save file frame number
        self.save_path = op_params["save_path"]
        if op_params["plot_to_file"]:
            self.ndigits_for_plot_files = int(
                np.ceil(np.log10(len(self.plot_times) + 1))
            )
            self.plot_num = 0  # current plot image frame number
        self.display_params = params

        if "format" in op_params:
            save_fmt = op_params["format"]
        else:
            save_fmt = "grid"

        if save_fmt == "grid":
            self.save_state = self.save_state_grid_format
        elif save_fmt == "vtk":
            self.save_state = self.save_state_vtk_format
        elif save_fmt == "netcdf":
            self.save_state = self.save_state_netcdf_format
        else:
            print("Unrecognized save format '" + save_fmt + "'.")
            print("Valid formats are: grid, vtk, netcdf")
            raise ValueError

    def setup_run_control(self, clock_params: dict) -> None:
        """
        Initialize variables related to control of run timing.

        Parameters
        ----------
        clock_params : dict
            Dictionary with items "start", "step", and "stop"
        """
        self.run_duration = clock_params["stop"] - clock_params["start"]
        self.dt = clock_params["step"]
        self.current_time = clock_params["start"]

    def report(self, current_time: float) -> None:
        """Issue a text update on status."""
        print(self.__class__.__name__, "time =", current_time)

    def plot(self, current_time: float = 0.0) -> None:
        """Virtual function for plotting; to be overridden."""
        print("Base class placeholder for plot() at time", current_time)

    def save_state_grid_format(
        self, save_path: str, save_num: int, ndigits: int
    ) -> None:
        """
        Save the grid and its fields in native Landlab format.

        Override this function to add to or modify what gets saved.
        """
        from landlab.io.native_landlab import save_grid

        save_grid(
            self.grid, save_path + str(save_num).zfill(ndigits) + ".grid", clobber=True
        )

    def save_state_vtk_format(
        self, save_path: str, save_num: int, ndigits: int
    ) -> None:
        """
        Save grid fields in legacy VTK format.

        Override this function to add to or modify what gets saved.
        """
        from landlab.io.legacy_vtk import write_legacy_vtk

        write_legacy_vtk(
            save_path + str(save_num).zfill(ndigits) + ".vtk", self.grid, clobber=True
        )

    def save_state_netcdf_format(
        self, save_path: str, save_num: int, ndigits: int
    ) -> None:
        """
        Save grid fields in NetCDF format.

        Override this function to add to or modify what gets saved.
        """
        from landlab.io.netcdf import write_netcdf

        write_netcdf(save_path + str(save_num).zfill(ndigits) + ".nc", self.grid)

    def update(self, dt: float) -> None:
        """
        Advance the model by one time step of duration dt.

        The derived class should override this function.
        """
        self.current_time += dt

    def update_until(self, update_to_time: float, dt: float) -> None:
        """Iterate up to given time, using time-step duration dt."""
        remaining_time = update_to_time - self.current_time
        while remaining_time > 0.0:
            dt = min(dt, remaining_time)
            self.update(dt)
            remaining_time -= dt

    def run(self, run_duration: float | None = None, dt: float | None = None) -> None:
        """Run the model for given duration, or self.run_duration if none
        given.

        Includes file output of images and model state at user-specified
        intervals.
        """
        if run_duration is None:
            run_duration = self.run_duration
        if dt is None:
            dt = self.dt

        stop_time = run_duration + self.current_time
        while self.current_time < stop_time:
            next_pause = min(self.next_plot, self.next_save)
            next_pause = min(next_pause, self.next_report)
            next_pause = min(next_pause, stop_time)
            self.update_until(next_pause, dt)
            if self.current_time >= self.next_report:
                self.report(self.current_time)
                self.next_report = self.report_times.pop(0)
            if self.current_time >= self.next_plot:
                self.plot(self.current_time)
                self.next_plot = self.plot_times.pop(0)
            if self.current_time >= self.next_save:
                self.save_num += 1
                self.save_state(
                    self.save_path, self.save_num, self.ndigits_for_save_files
                )
                self.next_save = self.save_times.pop(0)


if __name__ == "__main__":
    """Launch a run.

    Optional command-line argument is the name of a yaml-format text file with
    parameters. File should include sections for "grid_setup", "process",
    "run_control", and "output". Each of these should have the format shown in
    the defaults defined above in the class header.
    """
    if len(sys.argv) > 1:
        from landlab import load_params

        params = load_params(sys.argv[1])
        sim = LandlabModel(params)
    else:
        sim = LandlabModel()  # use default params
    sim.run()
