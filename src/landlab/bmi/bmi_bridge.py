"""
========================================================================================
Wrap landlab component with the Basic Modeling Interface (:mod:`landlab.bmi.bmi_bridge`)
========================================================================================

.. sectionauthor:: Eric Hutton

Function reference
------------------

The `wrap_as_bmi` function wraps a landlab component class so that it
exposes a Basic Modelling Interface.

"""

import inspect

import numpy as np
from bmipy import Bmi

from ..core import load_params
from ..core.model_component import Component
from ..framework.decorators import snake_case
from ..grid import HexModelGrid
from ..grid import RasterModelGrid
from ..grid.create import create_grid

BMI_LOCATION = {
    "node": "node",
    "link": "edge",
    "patch": "face",
    "corner": "node",
    "face": "edge",
    "cell": "face",
    "grid": "none",
}

BMI_GRID = {
    "node": 0,
    "link": 0,
    "patch": 0,
    "corner": 1,
    "face": 1,
    "cell": 1,
    "grid": 2,
}


class TimeStepper:
    """Step through time.

    Parameters
    ----------
    start : float, optional
        Clock start time.
    stop : float, optional
        Stop time.
    step : float, optional
        Time step.

    Examples
    --------
    >>> from landlab.bmi import TimeStepper
    >>> time_stepper = TimeStepper()
    >>> time_stepper.start
    0.0
    >>> time_stepper.stop is None
    True
    >>> time_stepper.step
    1.0
    >>> time_stepper.time
    0.0
    >>> for _ in range(10):
    ...     time_stepper.advance()
    ...
    >>> time_stepper.time
    10.0
    >>> time_stepper = TimeStepper(1.0, 13.0, 2.0)
    >>> [time for time in time_stepper]
    [1.0, 3.0, 5.0, 7.0, 9.0, 11.0]
    """

    def __init__(self, start=0.0, stop=None, step=1.0, units="s"):
        self._start = start
        self._stop = stop
        self._step = step
        self._units = units

        self._time = start

    def __iter__(self):
        if self.stop is None:
            while 1:
                yield self._time
                self._time += self._step
        else:
            while self._time < self._stop:
                yield self._time
                self._time += self._step
        return

    @property
    def time(self):
        """Current time."""
        return self._time

    @property
    def start(self):
        """Start time."""
        return self._start

    @property
    def stop(self):
        """Stop time."""
        return self._stop

    @property
    def step(self):
        """Time Step."""
        return self._step

    @step.setter
    def step(self, new_val):
        """Change the time step."""
        self._step = new_val

    @property
    def units(self):
        """Time units."""
        return self._units

    def advance(self):
        """Advance the time stepper by one time step."""
        self._time += self.step
        if self._stop is not None and self._time > self._stop:
            raise StopIteration()


def wrap_as_bmi(cls):
    """Wrap a landlab class so it exposes a BMI.

    Give a landlab component a Basic Model Interface (BMI). Since landlab
    components have an interface that is already in the style of BMI,
    this function adds just a light wrapping to landlab components. There
    are a number of differences that may cause some confusion to
    landlab users.

    1.  Because BMI doesn't have a concept of a dual grid, it only
        defines *nodes* (points), *edges* (vectors), and *faces*
        (areas). The dual-graph of landlab is considered as two
        separate grids by BMI.

    2. It is important to note that BMI has only three grid elements
       (*node*, *edge*, and *face*) while landlab has 6. The names
       used by landlab and BMI are also different.

       Thus, a BMI-wrapped landlab component will always have two
       grids with grid identifiers 0, and 1. Grid 0 will contain
       the landlab *nodes*, *links*, and *patches* while grid 1 will
       contain *corners*, *faces*, and *cells*.landlab and BMI
       refer to grid elements by different names. The mapping from
       landlab to BMI nomenclature is the following:

        Grid 0:
        *  *node*: *node*
        *  *link*: *edge*
        *  *patch*: *face*

        Grid 1:
        *  *corner*: *node*
        *  *face*: *edge*
        *  *cell*: *face*

    3.  In BMI, the *initialize* method requires an input file that is
        used to create and setup the model for time-stepping. landlab
        components generally do not have anything like this; instead
        this task is usually done programmatically. Thus, the
        input file that is used by the BMI *initialize* method is
        a standard landlab input file as used by the landlab *create_grid*
        function.

    Parameters
    ----------
    cls : class
        A landlab class that inherits from `Component`.

    Returns
    -------
    class
        A wrapped class that exposes a BMI.

    Examples
    --------
    >>> from landlab.bmi import wrap_as_bmi
    >>> from landlab.components.flexure import Flexure

    >>> BmiFlexure = wrap_as_bmi(Flexure)
    >>> flexure = BmiFlexure()
    >>> sorted(flexure.get_input_var_names())
    ['boundary_condition_flag', 'lithosphere__overlying_pressure_increment']
    >>> flexure.get_var_units("lithosphere__overlying_pressure_increment")
    'Pa'

    >>> config = '''
    ... flexure:
    ...     eet: 10.e+3
    ...     method: flexure
    ... clock:
    ...     start: 0.
    ...     stop: 10.
    ...     step: 2.
    ... grid:
    ...     RasterModelGrid:
    ...     - [20, 40]
    ...     - xy_spacing: [2000., 1000.]
    ...     - fields:
    ...        node:
    ...          lithosphere__overlying_pressure_increment:
    ...            constant:
    ...              - value: 0.0
    ... '''
    >>> flexure.initialize(config)
    >>> sorted(flexure.get_output_var_names())
    ['boundary_condition_flag', 'lithosphere_surface__elevation_increment']
    >>> flexure.get_var_grid("lithosphere_surface__elevation_increment")
    0
    >>> flexure.get_grid_shape(0, np.empty(flexure.get_grid_rank(0), dtype=int))
    array([20, 40])
    >>> dz = np.empty(flexure.get_grid_size(0))
    >>> _ = flexure.get_value("lithosphere_surface__elevation_increment", dz)

    >>> np.all(dz == 0.0)
    True
    >>> flexure.get_current_time()
    0.0

    >>> sorted(flexure.get_input_var_names())
    ['boundary_condition_flag', 'lithosphere__overlying_pressure_increment']
    >>> load = np.zeros((20, 40), dtype=float)
    >>> load[0, 0] = 1.0
    >>> flexure.set_value("lithosphere__overlying_pressure_increment", load)
    >>> flexure.update()
    >>> flexure.get_current_time()
    2.0
    >>> _ = flexure.get_value("lithosphere_surface__elevation_increment", dz)
    >>> np.all(dz == 0.0)
    False
    """
    if not issubclass(cls, Component):
        raise TypeError("class must inherit from Component")

    class BmiWrapper(Bmi):
        __doc__ = """
        Basic Modeling Interface for the {name} component.
        """.format(
            name=cls.__name__
        ).strip()

        _cls = cls

        def __init__(self):
            self._base = None
            self._clock = None
            super().__init__()

            self._input_var_names = tuple(
                set(self._cls.input_var_names) | {"boundary_condition_flag"}
            )
            self._output_var_names = tuple(
                set(self._cls.output_var_names) | {"boundary_condition_flag"}
            )
            self._info = self._cls._info.copy()

            self._info["boundary_condition_flag"] = {
                "mapping": "node",
                "units": "",
                "dtype": int,
                "intent": None,
                "doc": "boundary condition flag of grid nodes",
            }

        def get_component_name(self):
            """Name of the component."""
            return self._cls.name

        def get_input_var_names(self):
            """Names of the input exchange items."""
            return self._input_var_names

        def get_input_item_count(self):
            return len(self._input_var_names)

        def get_output_var_names(self):
            """Names of the output exchange items."""
            return self._output_var_names

        def get_output_item_count(self):
            return len(self._output_var_names)

        def get_current_time(self):
            """Current component time."""
            return self._clock.time

        def get_end_time(self):
            """Stop time for the component."""
            return self._clock.stop

        def get_start_time(self):
            """Start time of the component."""
            return self._clock.start

        def get_time_step(self):
            """Component time step."""
            return self._clock.step

        def get_time_units(self):
            """Time units used by the component."""
            return self._clock.units

        def initialize(self, config_file):
            """Initialize the component from a file.

            BMI-wrapped Landlab components use input files in YAML format.
            Component-specific parameters are listed at the top level,
            followed by grid and then time information. An example input
            file looks like::

                flexure:
                    eet: 15.e+3
                clock:
                    start: 0
                    stop: 100.
                    step: 2.
                grid:
                    type: raster
                    shape: [20, 40]
                    spacing: [1000., 2000.]

            In this case, a `RasterModelGrid` is created (with the given shape
            and spacing) and passed to the underlying landlab component. The
            `eet=15000.` is also given to the component but as a keyword
            parameter. The BMI clock is initialized with the given parameters.

            Parameters
            ----------
            config_file : str or file_like
                YAML-formatted input file for the component.
            """
            grid = create_grid(config_file, section="grid")

            if not grid:
                raise ValueError(f"no grid in config file ({config_file})")
            elif isinstance(grid, list):
                raise ValueError(f"multiple grids in config file ({config_file})")

            params = load_params(config_file)
            params.pop("grid")
            clock_params = params.pop("clock")
            self._clock = TimeStepper(**clock_params)

            self._base = self._cls(grid, **params.pop(snake_case(cls.__name__), {}))
            self._base.grid.at_node["boundary_condition_flag"] = (
                self._base.grid.status_at_node
            )

        def update(self):
            """Update the component one time step."""
            if hasattr(self._base, "update"):
                self._base.update()
            elif hasattr(self._base, "run_one_step"):
                args = []
                for name, arg in inspect.signature(
                    self._base.run_one_step
                ).parameters.items():
                    if arg.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD:
                        args.append(name)

                if len(args) == 0 or "dt" not in args:
                    self._base.run_one_step()
                else:
                    self._base.run_one_step(self._clock.step)

            self._clock.advance()

        def update_frac(self, frac):
            """Update the component a fraction of a time step."""
            time_step = self.get_time_step()
            self._clock.step = time_step * frac
            self.update()
            self._clock.step = time_step

        def update_until(self, then):
            """Update the component until a given time."""
            n_steps = (then - self.get_current_time()) / self.get_time_step()
            for _ in range(int(n_steps)):
                self.update()
            self.update_frac(n_steps - int(n_steps))

        def finalize(self):
            """Clean-up the component."""
            pass

        def get_var_grid(self, name):
            """Get the grid id for a variable."""
            at = self._info[name]["mapping"]
            return BMI_GRID[at]

        def get_var_itemsize(self, name):
            """Get the size of elements of a variable."""
            at = self._info[name]["mapping"]
            return self._base.grid[at][name].itemsize

        def get_var_nbytes(self, name):
            """Get the total number of bytes used by a variable."""
            at = self._info[name]["mapping"]
            return self._base.grid[at][name].nbytes

        def get_var_type(self, name):
            """Get the data type for a variable."""
            at = self._info[name]["mapping"]
            return str(self._base.grid[at][name].dtype)

        def get_var_units(self, name):
            """Get the unit used by a variable."""
            return self._info[name]["units"]

        def get_value_ref(self, name):
            """Get a reference to a variable's data."""
            at = self._info[name]["mapping"]
            return self._base.grid[at][name]

        def get_value(self, name, dest):
            """Get a copy of a variable's data."""
            at = self._info[name]["mapping"]
            dest[:] = self._base.grid[at][name]
            return dest

        def set_value(self, name, values):
            """Set the values of a variable."""
            if name in self.get_input_var_names():
                if name == "boundary_condition_flag":
                    self._base.grid.status_at_node = values
                else:
                    at = self._info[name]["mapping"]
                    self._base.grid[at][name][:] = values.flat
            else:
                raise KeyError(f"{name} is not an input item")

        def get_grid_origin(self, grid, origin):
            """Get the origin for a structured grid."""
            if grid == 0:
                origin[:] = (self._base.grid.node_y[0], self._base.grid.node_x[0])
            elif grid == 1:
                origin[:] = (
                    self._base.grid.node_y[0] + self._base.grid.dy * 0.5,
                    self._base.grid.node_x[0] + self._base.grid.dx * 0.5,
                )
            return origin

        def get_grid_rank(self, grid):
            """Get the number of dimensions of a grid."""
            if grid in (0, 1):
                return 2
            else:
                return 0

        def get_grid_shape(self, grid, shape):
            """Get the shape of a structured grid."""
            if grid == 0:
                shape[:] = (
                    self._base.grid.number_of_node_rows,
                    self._base.grid.number_of_node_columns,
                )
            elif grid == 1:
                shape[:] = (
                    self._base.grid.number_of_node_rows - 1,
                    self._base.grid.number_of_node_columns - 1,
                )
            return shape

        def get_grid_spacing(self, grid, spacing):
            """Get the row and column spacing of a structured grid."""
            spacing[:] = (self._base.grid.dy, self._base.grid.dx)
            return spacing

        def get_grid_type(self, grid):
            """Get the type of grid."""
            if grid == 2:
                return "scalar"
            elif isinstance(self._base.grid, RasterModelGrid):
                return "uniform_rectilinear"
            else:
                return "unstructured"

        def get_grid_edge_count(self, grid):
            if grid == 0:
                return self._base.grid.number_of_links
            elif grid == 1:
                return self._base.grid.number_of_faces

        def get_grid_edge_nodes(self, grid, edge_nodes):
            if grid == 0:
                return self._base.grid.nodes_at_link.reshape((-1,))
            elif grid == 1:
                return self._base.grid.corners_at_face.reshape((-1,))

        def get_grid_face_count(self, grid):
            if grid == 0:
                return self._base.grid.number_of_patches
            elif grid == 1:
                return self._base.grid.number_of_cells

        def get_grid_face_nodes(self, grid, face_nodes):
            if grid == 0:
                return self._base.grid.nodes_at_patch
            elif grid == 1:
                return self._base.grid.corners_at_cell

        def get_grid_face_edges(self, grid):
            if grid == 0:
                return self._base.grid.links_at_patch
            elif grid == 1:
                return self._base.grid.faces_at_cell

        def get_grid_node_count(self, grid):
            if grid == 0:
                return self._base.grid.number_of_nodes
            elif grid == 1:
                return self._base.grid.number_of_corners

        def get_grid_nodes_per_face(self, grid, nodes_per_face):
            if grid == 0:
                return np.full(self._base.grid.number_of_nodes, 3, dtype=int)
            elif grid == 1 and isinstance(self._base.grid, HexModelGrid):
                return np.full(self._base.grid.number_of_faces, 6, dtype=int)

        def get_grid_size(self, grid):
            if grid == 0:
                return self._base.grid.number_of_nodes
            elif grid == 1:
                return self._base.grid.number_of_corners

        def get_grid_x(self, grid, x):
            if grid == 0:
                return self._base.grid.x_of_node
            elif grid == 1:
                return self._base.grid.x_of_corner

        def get_grid_y(self, grid, y):
            if grid == 0:
                return self._base.grid.y_of_node
            elif grid == 1:
                return self._base.grid.y_of_corner

        def get_grid_z(self, grid, z):
            raise NotImplementedError("get_grid_z")
            # Only should be implemented for presently non-existant 3D grids.

        def get_value_at_indices(self, name, dest, inds):
            at = self._info[name]["mapping"]
            dest[:] = self._base.grid[at][name][inds]
            return dest

        def get_value_ptr(self, name):
            at = self._info[name]["mapping"]
            return self._base.grid[at][name]

        def get_var_location(self, name):
            return BMI_LOCATION[self._info[name]["mapping"]]

        def set_value_at_indices(self, name, inds, src):
            at = self._info[name]["mapping"]
            self._base.grid[at][name][inds] = src

    BmiWrapper.__name__ = f"{cls.__name__}BMI"
    return BmiWrapper
