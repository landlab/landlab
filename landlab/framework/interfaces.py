"""
The Basic Modeling Interface.
"""


class Error(Exception):
    """Base class for BMI exceptions"""

    pass


class FatalError(Exception):
    """
    Raise this exception if an unrecoverable error was found
    """

    pass


class BadVarNameError(Error):
    """Exception to indicate a bad input/output variable name"""

    def __init__(self, name):
        super(BadVarNameError, self).__init__()
        self.name = name

    def __str__(self):
        return self.name


class MissingModelAttributeError(Error):
    """
    Raise this exception if a component is missing a required attribute.
    """

    def __init__(self, attrib):
        super(MissingModelAttributeError, self).__init__()
        self.attrib = attrib

    def __str__(self):
        return self.attrib


class TimeBoundsError(Error):
    """
    Raise this exception if a component updates beyond its time horizon
    """

    pass


class BmiGridType(int):
    """
    Base type to indicate the type of a BMI model's grid.

    :code: Grid type code as an int
    :name: Name of the grid type as a string
    """

    def __new__(cls, code, name):
        obj = super(BmiGridType, cls).__new__(cls, code)
        obj.name = name
        return obj

    def __str__(self):
        return self.name

    def __repr__(self):
        return "BmiGridType(%d, %s)" % (self, self.name)


GRID_TYPE_UNKNOWN = BmiGridType(-1, "Unknown")
GRID_TYPE_NONE = BmiGridType(0, "No grid")
GRID_TYPE_UNIFORM = BmiGridType(1, "Uniform rectilinear")
GRID_TYPE_RECTILINEAR = BmiGridType(2, "Rectilinear")
GRID_TYPE_STRUCTURED = BmiGridType(3, "Structured")
GRID_TYPE_UNSTRUCTURED = BmiGridType(4, "Unstructured")


class BmiBase(object):
    """
    Definition of the Basic Modeling Interface
    """

    def initialize(self, file_name):
        """
        Initialize model.

        :file_name: String of configuration file
        """
        pass

    def update(self, **kwds):
        """
        Update model by one time step.
        """
        pass

    def finalize(self):
        """
        Clean-up model
        """
        pass

    def get_input_var_names(self):
        """
        Get names of input variables to the model as standard names.

        :returns: A list of input standard names as strings
        """
        pass

    def get_output_var_names(self):
        """
        Get names of output variables to the model as standard names.

        :returns: A list of output standard names as strings
        """
        pass

    def get_var_type(self, var_name):
        """
        Get type of an exchange item.
        """
        pass

    def get_var_units(self, var_name):
        """
        Get units of an exchange item.
        """
        pass

    def get_var_rank(self, var_name):
        """
        Rank of exchange item.
        """
        pass

    def get_time_step(self):
        """
        Model time step.
        """
        pass

    def get_start_time(self):
        """
        Model start time.
        """
        pass

    def get_current_time(self):
        """
        Current time of model.
        """
        pass

    def get_end_time(self):
        """
        Model stop time.
        """
        pass


class BmiExtendedBase(object):
    """
    An extension interface for a BMI.
    """

    def update_until(self, time):
        """
        Update model until some time.

        :time: Update duration
        """
        pass

    def run_model(self):
        """
        Initialize, run, and finalize a model.
        """
        pass


class BmiUnstructured(object):
    """
    BMI for a model that uses an unstructured grid.
    """

    def get_x(self, name):
        """
        Get x-coordinates of grid nodes.
        """
        pass

    def get_y(self, name):
        """
        Get y-coordinates of grid nodes.
        """
        pass

    def get_connectivity(self, name):
        """
        Get cell connectivity.
        """
        pass

    def get_offset(self, name):
        """
        Get cell offset.
        """
        pass


class BmiStructured(object):
    """
    BMI for a model that uses a structured grid.
    """

    def get_grid_shape(self, name):
        """
        Get shape of grid for variable, name.

        :name: Standard name
        """
        pass

    def get_x(self, name):
        """
        Get x-coordinates of grid nodes.
        """
        pass

    def get_y(self, name):
        """
        Get y-coordinates of grid nodes.
        """
        pass


class BmiRectilinear(object):
    """
    BMI for a model that uses a rectilinear grid.
    """

    def get_grid_shape(self, name):
        """
        Get shape of grid for variable, name.

        :name: Standard name
        """
        pass

    def get_columns(self, name):
        """
        Get coordinates of grid columns.
        """
        pass

    def get_rows(self, name):
        """
        Get coordinates of grid rows.
        """
        pass


class BmiUniformRectilinear(object):
    """
    BMI for a model that exposes a uniform rectilinear grid.
    """

    def get_grid_shape(self, name):
        """
        Get shape of grid for variable, name.

        :name: Standard name
        """
        pass

    def get_grid_spacing(self, name):
        """
        Get spacing of grid for variable, name.

        :name: Standard name
        """
        pass

    def get_grid_origin(self, name):
        """
        Get origin of grid for variable, name.

        :name: Standard name
        """
        pass


class BmiNoGrid(object):
    """
    BMI for a model that does not have a grid.
    """
