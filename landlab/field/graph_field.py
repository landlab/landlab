import numpy as np
import xarray as xr

from .grouped import GroupError
from .scalar_data_fields import FieldError


def reshape_for_storage(array, field_size=None):
    """Reshape an array to be stored as a field.

    For reshaping rules, see ``shape_for_storage``.

    Parameters
    ----------
    array : ndarray
        The array to be stored.
    field_size : int, optional
        The size of the field.

    Returns
    -------
    ndarray
        The (possibly) reshaped array.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.field.graph_field import reshape_for_storage

    The shape will be such that the first dimension in the field size.

    >>> data = np.arange(6)
    >>> reshape_for_storage(data, 3)
    array([[0, 1],
           [2, 3],
           [4, 5]])
    >>> reshape_for_storage(data, 2)
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> reshape_for_storage(data, 6)
    array([0, 1, 2, 3, 4, 5])

    If the array is already the correct shape, just return that array.

    >>> data = np.arange(6).reshape((2, 3))
    >>> reshape_for_storage(data, 2) is data
    True
    """
    shape = shape_for_storage(array, field_size)
    if shape == array.shape or array.ndim == 0:
        return array
    else:
        return array.reshape(shape)


def shape_for_storage(array, field_size=None):
    """The shape an array will be stored as.

    Parameters
    ----------
    array : ndarray
        The array to be stored.
    field_size : int, optional
        The size of the field.

    Returns
    -------
    tuple of int
        The shape the array will be stored as.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.field.graph_field import shape_for_storage

    The shape will be such that the first dimension in the field size.

    >>> data = np.arange(6)
    >>> shape_for_storage(data, 3) == (3, 2)
    True
    >>> shape_for_storage(data, 2) == (2, 3)
    True
    >>> shape_for_storage(data, 6) == (6, )
    True

    If a field size is not given, the array will be stored as a
    flattened array.

    >>> shape_for_storage(data) == (6, )
    True
    >>> data = np.arange(6).reshape((3, 2))
    >>> shape_for_storage(data) == (6, )
    True

    For field sizes of 1, the array is always flattened.

    >>> shape_for_storage(data, 1) == (1, 6)
    True

    For scalar arrays, the field size must be 1.

    >>> data = np.array(1.)
    >>> shape_for_storage(data) == (1, )
    True
    >>> shape_for_storage(data, field_size=1) == (1, )
    True

    If the array cannot be shaped into a storage shape, a `ValueError`
    is raised.

    >>> data = np.array(1.)
    >>> shape_for_storage(data, field_size=4) # DOCTEST: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    ValueError: unable to reshape array to field size
    >>> data = np.arange(6.)
    >>> shape_for_storage(data, field_size=4) # DOCTEST: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    ValueError: unable to reshape array to field size
    """
    if field_size is None:
        field_size = array.size

    if array.size % field_size != 0:
        raise ValueError(
            "unable to reshape array to field size ({0} != {1})".format(
                array.size, field_size
            )
        )

    if field_size == array.size:
        shape = (array.size,)
    else:
        shape = (field_size, array.size // field_size)

    return shape


class FieldDataset(dict):

    """Wrap an xarray.Dataset as a landlab field.

    This is a light wrapping of xarray.Dataset. The main differences
    are that a `FieldDataset` can be created with a size but not
    allocate any memory for data arrays until an array is actually
    needed. The setitem method is also overriden so that when arrays
    are added they are stored reshaped in the landlab style. That
    is, shaped as `(n_elements, values_per_element)`.

    Examples
    --------
    >>> from landlab.field.graph_field import FieldDataset

    >>> ds = FieldDataset("node")
    >>> ds.size is None
    True
    >>> ds.set_value("air_temperature", [1.0, 1.0, 1.0, 1.0])
    >>> ds["air_temperature"]
    array([ 1.,  1.,  1.,  1.])
    >>> ds.size
    4
    >>> ds.set_value("air_temperature", [1.0, 1.0, 1.0])
    Traceback (most recent call last):
    ValueError: unable to reshape array to field size (3 != 4)

    >>> ds = FieldDataset("node", fixed_size=False)
    >>> ds.size is None
    True
    >>> ds.set_value("air_temperature", [1.0, 1.0, 1.0, 1.0])
    >>> ds.size
    4
    >>> ds.set_value("air_temperature", [0.0, 0.0])
    >>> ds.size
    2
    >>> ds["air_temperature"]
    array([ 0.,  0.])
    """

    def __init__(self, name, size=None, fixed_size=True):
        self._name = name
        self._size = None
        self._fixed_size = bool(fixed_size)
        self._ds = xr.Dataset()
        self._units = {}

        self.size = size

        super(FieldDataset, self).__init__()

    @property
    def size(self):
        """Size of the field dataset as number of elements.

        Examples
        --------
        >>> from landlab.field.graph_field import FieldDataset

        >>> ds = FieldDataset("grid", size=1)
        >>> ds.set_value("air_temperature", [1.0, 1.0, 1.0])
        >>> ds.set_value("ground_temperature", [0.0, 0.0])
        >>> ds["ground_temperature"]
        array([[ 0.,  0.]])

        >>> ds = FieldDataset("grid", size=1)
        >>> ds.set_value("air_temperature", 0.1)
        >>> ds["air_temperature"]
        array(0.1)
        """
        return self._size

    @size.setter
    def size(self, size):
        if self._size != size:
            if self._size is not None and self.fixed_size:
                raise ValueError(
                    "size has already been set ({size}) and fixed_size is True".format(
                        size=self._size
                    )
                )
            elif not isinstance(size, int) or size < 0:
                raise ValueError(
                    "size must be a positive integer or None ({size})".format(size=size)
                )
            self._size = size

    @property
    def fixed_size(self):
        """Flag that indicates if arrays added to the dataset must be of a
        fixed size.

        Examples
        --------
        >>> from landlab.field.graph_field import FieldDataset

        >>> ds = FieldDataset("node", fixed_size=False)
        >>> ds.set_value("air_temperature", [1.0, 1.0, 1.0, 1.0])
        >>> ds.set_value("air_temperature", [0.0, 0.0])
        >>> ds["air_temperature"]
        array([ 0.,  0.])

        >>> ds.fixed_size = True
        >>> ds.size
        2
        >>> ds.set_value("air_temperature", [1.0, 1.0, 1.0])
        Traceback (most recent call last):
        ValueError: unable to reshape array to field size (4 != 2)
        """
        return self._fixed_size

    @fixed_size.setter
    def fixed_size(self, fixed_size):
        self._fixed_size = bool(fixed_size)
        if self._fixed_size:
            self.size = self._ds.dims[self._name]

    @property
    def units(self):
        return self._units

    @property
    def dataset(self):
        return self._ds

    def keys(self):
        return list(self._ds.variables)

    def set_value(self, name, value_array, attrs=None):
        attrs = attrs or {}
        attrs.setdefault("units", "?")

        value_array = np.asarray(value_array)

        if name in self._ds and self._ds[name].values is value_array:
            self._ds[name].values.shape = shape_for_storage(value_array, self.size)
            return

        if self.fixed_size:
            value_array = reshape_for_storage(value_array, self.size)
        else:
            value_array = reshape_for_storage(value_array, None)

        if value_array.ndim > 0:
            self.size = value_array.shape[0]
        else:
            self.size = value_array.size

        if value_array.ndim == 0:
            dims = ()
        elif value_array.ndim == 1:
            dims = (self._name,)
        else:
            dims = (self._name, name + "_per_" + self._name)

        if name in self._ds:
            self._ds = self._ds.drop(name)

        self._ds.update({name: xr.DataArray(value_array, dims=dims, attrs=attrs)})
        self._units[name] = attrs["units"]

    def pop(self, name):
        array = self._ds[name].values
        self._ds = self._ds.drop(name)
        return array

    def __getitem__(self, name):
        if isinstance(name, str):
            try:
                return self._ds[name].values
            except KeyError:
                raise FieldError(name)
        else:
            raise TypeError("field name not a string")

    def __setitem__(self, name, value_array):
        self.set_value(name, value_array)

    def __contains__(self, name):
        return name in self._ds

    def __str__(self):
        return str(self._ds)

    def __repr__(self):
        return "FieldDataset({name}, size={size}, fixed_size={fixed_size})".format(
            name=repr(self._name),
            size=repr(self.size),
            fixed_size=repr(self.fixed_size),
        )

    def __len__(self):
        return len(self._ds.variables)

    def __iter__(self):
        return iter(self._ds.variables)


class GraphFields(object):

    """Collection of grouped data-fields.

    The `GraphFields` class holds a set of data fields that are
    separated into *groups*. A typical use for this class would be to define
    the groups as being locations on a grid where the values are defined.
    For instance, the groups could be *node*, *cell*, *link*, and *face*.

    Attributes
    ----------
    groups

    See Also
    --------
    landlab.field.ScalarDataFields : Data fields within a *group* are
        stored as :class:`landlab.field.ScalarDataFields`.
    landlab.field.ModelDataFields : Equivalent data structure for
        old-style fields.

    Examples
    --------

    Create two groups of data fields defined at *node* and *cell*. Each set can
    have a differenct number of values.

    >>> from landlab.field import GraphFields
    >>> fields = GraphFields()
    >>> fields.new_field_location('node', 12)
    >>> fields.new_field_location('cell', 2)

    Create some new value arrays for each of the data fields.

    >>> fields.ones('node')
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
    >>> fields.zeros('cell')
    array([ 0.,  0.])

    Create new value arrays and add them to the data fields. Because the data
    fields are in different groups (node and cell), they can have the same
    name.

    >>> fields.add_ones('topographic__elevation', at='node')
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])
    >>> fields.at_node['topographic__elevation']
    array([ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

    >>> fields.add_ones('topographic__elevation', at='cell')
    array([ 1.,  1.])
    >>> fields.at_cell['topographic__elevation']
    array([ 1.,  1.])

    Each group acts as a `dict` so, for instance, to get the variables names
    in a group use the `keys` method,

    >>> list(fields.at_cell.keys())
    ['topographic__elevation']

    If the size of the new field location is ``None``, the field will be
    unsized. This means that fields added to this location can be of any
    size.

    >>> fields = GraphFields()
    >>> fields.new_field_location('grid', None)
    >>> fields.at_grid['g'] = 9.81
    >>> fields.at_grid['g']
    array(9.81)
    >>> fields.at_grid['w'] = (3., 4.)
    >>> fields.at_grid['w']
    array([ 3.,  4.])

    The dimensions of groups can also be specified when the object is
    instantiated. In this case, group sizes are specified as a dictionary
    with keys being group names and values group sizes.

    >>> fields = GraphFields({'node': 6, 'grid': None})
    >>> fields.at_grid['g'] = 9.81
    >>> fields.at_node['x'] = [0, 1, 2, 3, 4, 5]
    >>> fields.at_grid['g']
    array(9.81)
    >>> fields.at_node['x']
    array([0, 1, 2, 3, 4, 5])
    """

    def __init__(self, *args, **kwds):
        try:
            dims = args[0]
        except IndexError:
            dims = {}

        self._groups = set()
        for loc in dims:
            self.new_field_location(loc, dims[loc])

        self.default_group = kwds.get("default_group", None)

    def __getitem__(self, name):
        try:
            return getattr(self, "at_" + name)
        except AttributeError:
            raise GroupError(name)

    @property
    def default_group(self):
        return self._default_group

    @default_group.setter
    def default_group(self, loc):
        if self.has_group(loc) or loc is None:
            self._default_group = loc
        else:
            raise ValueError("{loc} is not a valid group name".format(loc=loc))

    def new_field_location(self, loc, size=None):
        """Add a new quantity to a field.

        Create an empty group into which new fields can be added. The new group
        is created but no memory allocated yet. The dictionary of the new group
        can be through a new *at_* attribute of the class instance.

        Parameters
        ----------
        group: str
            Name of the new group to add to the field.
        size: int, optional
            Number of elements in the new quantity. If not provided, the
            size is set to be the size of the first field added to the goup.

        Raises
        ------
        ValueError
            If the field already contains the group.

        Examples
        --------
        Create a collection of fields and add two groups, *node* and *cell*,
        to it.

        >>> from landlab.field import GraphFields
        >>> fields = GraphFields()
        >>> fields.new_field_location("node", 12)
        >>> fields.new_field_location("cell", 2)

        The group names in the collection are retrieved with the *groups*
        attribute as a `set`.

        >>> names = list(fields.groups)
        >>> names.sort()
        >>> names
        ['cell', 'node']

        Access the new (empty) groups with the *at_* attributes.

        >>> fields.at_cell
        FieldDataset('cell', size=2, fixed_size=True)
        >>> fields.at_node
        FieldDataset('node', size=12, fixed_size=True)

        >>> fields.new_field_location("core_node")
        >>> fields.at_core_node.size is None
        True
        >>> fields.at_core_node["air__temperature"] = [0, 1]
        >>> fields.at_core_node.size
        2

        LLCATS: FIELDCR
        """
        dataset_name = "at_" + loc
        if loc not in self._groups:
            setattr(
                self, dataset_name, FieldDataset(loc, size, fixed_size=size is not None)
            )
            self._groups.add(loc)
        else:
            raise ValueError("{loc} location already exists".format(loc=loc))

    @property
    def groups(self):
        """List of group names.

        Returns
        -------
        set
            Set of quantity names.
        """
        return self._groups

    def has_group(self, name):
        """Check if a group exists.

        Parameters
        ----------
        group: str
            Name of the group.

        Returns
        -------
        boolean
            True if the field contains *group*, otherwise False.

        Examples
        --------
        Check if the field has the groups named *node* or *cell*.

        >>> from landlab.field import GraphFields
        >>> fields = GraphFields()
        >>> fields.new_field_location('node', 12)
        >>> fields.has_group('node')
        True
        >>> fields.has_group('cell')
        False

        LLCATS: FIELDINF
        """
        return name in self._groups

    def has_field(self, group, field):
        """Check if a field is in a group.

        Parameters
        ----------
        group: str
            Name of the group.
        field: str
            Name of the field.

        Returns
        -------
        boolean
            ``True`` if the group contains the field, otherwise ``False``.

        Examples
        --------
        Check if the field named ``topographic__elevation`` is contained
        in a group.

        >>> from landlab.field import GraphFields
        >>> fields = GraphFields()
        >>> fields.new_field_location('node', 12)
        >>> _ = fields.add_ones("topographic__elevation", at="node")
        >>> fields.has_field('node', 'topographic__elevation')
        True
        >>> fields.has_field('cell', 'topographic__elevation')
        False

        LLCATS: FIELDINF
        """
        try:
            return field in self[group]
        except KeyError:
            return False

    def keys(self, group):
        """List of field names in a group.

        Returns a list of the field names as a list of strings.

        Parameters
        ----------
        group : str
            Group name.

        Returns
        -------
        list
            List of field names.

        Examples
        --------
        >>> from landlab.field import GraphFields
        >>> fields = GraphFields()
        >>> fields.new_field_location('node', 4)
        >>> list(fields.keys('node'))
        []
        >>> _ = fields.add_empty('topographic__elevation', at='node')
        >>> list(fields.keys('node'))
        ['topographic__elevation']

        LLCATS: FIELDINF
        """
        return self[group].keys()

    def size(self, group):
        """Size of the arrays stored in a group.

        Parameters
        ----------
        group : str
            Group name.

        Returns
        -------
        int
            Array size.

        Examples
        --------
        >>> from landlab.field import GraphFields
        >>> fields = GraphFields()
        >>> fields.new_field_location('node', 4)
        >>> fields.size('node')
        4

        LLCATS: GINF FIELDINF
        """
        return self[group].size

    def field_values(self, group, field):
        """Get values of a field.

        Given a *group* and a *field*, return a reference to the associated
        data array.

        Parameters
        ----------
        group: str
            Name of the group.
        field: str
            Name of the field withing *group*.

        Returns
        -------
        array
            The values of the field.

        Raises
        ------
        GroupError
            If *group* does not exits
        FieldError
            If *field* does not exits

        Examples
        --------
        Create a group of fields called *node*.

        >>> from landlab.field import GraphFields
        >>> fields = GraphFields()
        >>> fields.new_field_location('node', 4)

        Add a field, initialized to ones, called *topographic__elevation*
        to the *node* group. The *field_values* method returns a reference
        to the field's data.

        >>> _ = fields.add_ones("topographic__elevation", at="node")
        >>> fields.field_values('node', 'topographic__elevation')
        array([ 1.,  1.,  1.,  1.])

        Raise FieldError if *field* does not exist in *group*.

        >>> fields.field_values('node', 'planet_surface__temperature')
        ...     # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        FieldError: planet_surface__temperature

        If *group* does not exists, Raise GroupError.

        >>> fields.field_values('cell', 'topographic__elevation')
        ...     # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        GroupError: cell

        LLCATS: FIELDIO
        """
        try:
            fields = self[group]
        except KeyError:
            raise GroupError(group)
        try:
            return fields[field]
        except KeyError:
            raise FieldError(group)

    def return_array_or_field_values(self, group, field):
        """Return field given a field name, or array of with the correct shape.

        Given a *group* and a *field*, return a reference to the associated
        data array. *field* is either a string that is a field in the group
        or an array of the correct size.

        This function is meant to serve like the ``use_field_name_or_array``
        decorator for bound functions.

        Parameters
        ----------
        group: str
            Name of the group.
        field: str or array
            Name of the field withing *group*.

        Returns
        -------
        array
            The values of the field.

        Raises
        ------
        GroupError
            If *group* does not exits
        FieldError
            If *field* does not exits

        Examples
        --------
        Create a group of fields called *node*.

        >>> import numpy as np
        >>> from landlab.field import ModelDataFields
        >>> fields = ModelDataFields()
        >>> fields.new_field_location('node', 4)

        Add a field, initialized to ones, called *topographic__elevation*
        to the *node* group. The *field_values* method returns a reference
        to the field's data.

        >>> _ = fields.add_ones("topographic__elevation", at="node")
        >>> fields.field_values('node', 'topographic__elevation')
        array([ 1.,  1.,  1.,  1.])

        Alternatively, if the second argument is an array, its size is
        checked and returned if correct.

        >>> vals = np.array([4., 5., 7., 3.])
        >>> fields.return_array_or_field_values('node', vals)
        array([ 4.,  5.,  7.,  3.])

        Raise FieldError if *field* does not exist in *group*.

        >>> fields.return_array_or_field_values('node', 'surface__temperature')
        ...     # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        FieldError: surface__temperature

        If *group* does not exists, Raise GroupError.

        >>> fields.return_array_or_field_values('cell', 'topographic__elevation')
        ...     # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        GroupError: cell

        And if the array of values provided is incorrect, raise a ValueError.

        >>> vals = np.array([3., 2., 1.])
        >>> fields.return_array_or_field_values('node', vals)
        ...     # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ValueError: Array has incorrect size.

        LLCATS: FIELDIO
        """
        if isinstance(field, str):
            vals = self.field_values(group, field)
        else:
            vals = np.asarray(field)
            if vals.size != self[group].size:
                msg = "Array has incorrect size."
                raise ValueError(msg)
        return vals

    def field_units(self, group, field):
        """Get units for a field.

        Returns the unit string associated with the data array in *group* and
        *field*.

        Parameters
        ----------
        group: str
            Name of the group.
        field: str
            Name of the field withing *group*.

        Returns
        -------
        str
            The units of the field.

        Raises
        ------
        KeyError
            If either *field* or *group* does not exist.

        LLCATS: FIELDINF
        """
        return self[group]._ds[field].attrs["units"]

    def empty(self, *args, **kwds):
        """Uninitialized array whose size is that of the field.

        Return a new array of the data field size, without initializing
        entries. Keyword arguments are the same as that for the equivalent
        numpy function.

        Parameters
        ----------
        group : str
            Name of the group.

        See Also
        --------
        numpy.empty : See for a description of optional keywords.
        landlab.field.GraphFields.ones : Equivalent method that
            initializes the data to 1.
        landlab.field.GraphFields.zeros : Equivalent method that
            initializes the data to 0.

        Examples
        --------
        >>> from landlab.field import GraphFields
        >>> field = GraphFields()
        >>> field.new_field_location('node', 4)
        >>> field.empty('node') # doctest: +SKIP
        array([  2.31584178e+077,  -2.68156175e+154,   9.88131292e-324,
        ... 2.78134232e-309]) # Uninitialized memory

        Note that a new field is *not* added to the collection of fields.

        >>> list(field.keys('node'))
        []

        LLCATS: FIELDCR
        """
        if len(args) == 0:
            group = kwds.pop("at", kwds.pop("centering", "node"))
        else:
            group = args[0]

        if group == "grid":
            raise ValueError(
                "ones is not supported for at='grid', if you "
                "want to create a field at the grid, use\n"
                "grid.at_grid['value_name']=value\n"
                "instead.\nAlternatively, if you want ones"
                "of the shape stored at_grid, use np.array(1)."
            )

        size = getattr(self, "at_{group}".format(group=group)).size
        if size is None:
            raise ValueError("group is not yet sized.")

        return np.empty(size, **kwds)

    def ones(self, *args, **kwds):
        """Array, initialized to 1, whose size is that of the field.

        Return a new array of the data field size, filled with ones. Keyword
        arguments are the same as that for the equivalent numpy function.

        Parameters
        ----------
        group : str
            Name of the group.

        See Also
        --------
        numpy.ones : See for a description of optional keywords.
        landlab.field.GraphFields.empty : Equivalent method that
            does not initialize the new array.
        landlab.field.GraphFields.zeros : Equivalent method that
            initializes the data to 0.

        Examples
        --------
        >>> from landlab.field import GraphFields
        >>> field = GraphFields()
        >>> field.new_field_location('node', 4)
        >>> field.ones('node')
        array([ 1.,  1.,  1.,  1.])
        >>> field.ones('node', dtype=int)
        array([1, 1, 1, 1])

        Note that a new field is *not* added to the collection of fields.

        >>> list(field.keys('node'))
        []

        LLCATS: FIELDCR
        """
        allocated = self.empty(*args, **kwds)
        allocated.fill(1)
        return allocated

    def zeros(self, *args, **kwds):
        """Array, initialized to 0, whose size is that of the field.

        Parameters
        ----------
        group : str
            Name of the group.

        Return a new array of the data field size, filled with zeros. Keyword
        arguments are the same as that for the equivalent numpy function.

        This method is not valid for the group *grid*.

        See Also
        --------
        numpy.zeros : See for a description of optional keywords.
        landlab.field.GraphFields.empty : Equivalent method that does not
            initialize the new array.
        landlab.field.GraphFields.ones : Equivalent
            method that initializes the data to 1.

        Examples
        --------
        >>> from landlab.field import GraphFields
        >>> field = GraphFields()
        >>> field.new_field_location('node', 4)
        >>> field.zeros('node')
        array([ 0.,  0.,  0.,  0.])

        Note that a new field is *not* added to the collection of fields.

        >>> list(field.keys('node'))
        []

        LLCATS: FIELDCR
        """
        allocated = self.empty(*args, **kwds)
        allocated.fill(0)
        return allocated

    def add_field(self, *args, **kwds):
        """add_field(name, value_array, at='node', units='-', copy=False,
        clobber=False)

        Add an array of values to the field.

        Add an array of data values to a collection of fields and associate it
        with the key, *name*. Use the *copy* keyword to, optionally, add a
        copy of the provided array.

        In the case of adding to the collection *grid*, the added field is a
        numpy scalar rather than a numpy array.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        value_array : numpy.array
            Array of values to add to the field.
        at : str, optional
            Grid location to store values. If not given, values are
            assumed to be on `node`.
        units : str, optional
            Optionally specify the units of the field.
        copy : boolean, optional
            If True, add a *copy* of the array to the field. Otherwise save add
            a reference to the array.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        numpy.array
            The data array added to the field. Depending on the *copy*
            keyword, this could be a copy of *value_array* or *value_array*
            itself.

        Raises
        ------
        ValueError :
            If *value_array* has a size different from the field.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.field import GraphFields
        >>> field = GraphFields()
        >>> field.new_field_location('node', 4)
        >>> values = np.ones(4, dtype=int)
        >>> field.add_field("topographic__elevation", values, at="node")
        array([1, 1, 1, 1])

        A new field is added to the collection of fields. The saved value
        array is the same as the one initially created.

        >>> field.at_node['topographic__elevation'] is values
        True

        If you want to save a copy of the array, use the *copy* keyword. In
        addition, adding values to an existing field will remove the reference
        to the previously saved array. The *clobber=False* keyword changes this
        behavior to raise an exception in such a case.

        >>> field.add_field(
        ...     "topographic__elevation", values, at="node", copy=True, clobber=True
        ... )
        array([1, 1, 1, 1])
        >>> field.at_node['topographic__elevation'] is values
        False
        >>> field.add_field(
        ...     "topographic__elevation", values, at="node", clobber=False
        ... ) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        FieldError: topographic__elevation

        LLCATS: FIELDCR
        """
        if len(args) == 3:
            at, name, value_array = args
        elif len(args) == 2:
            at, name, value_array = (kwds.pop("at", None), args[0], args[1])
        else:
            raise ValueError("number of arguments must be 2 or 3")

        units = kwds.get("units", "?")
        copy = kwds.get("copy", False)
        clobber = kwds.get("clobber", False)
        value_array = np.asarray(value_array)

        at = at or self.default_group
        if at is None:
            raise ValueError("no group specified")

        attrs = {"long_name": name}
        attrs["units"] = units

        if copy:
            value_array = value_array.copy()

        ds = getattr(self, "at_" + at)

        if not clobber and name in ds:
            raise FieldError("{name}@{at}".format(name=name, at=at))

        dims = (at,)
        if value_array.ndim > 1:
            dims += (name + "_per_" + at,)
            value_array = value_array.reshape((value_array.shape[0], -1))

        ds.set_value(name, value_array, attrs=attrs)
        # ds[name] = value_array
        return ds[name]

    def delete_field(self, loc, name):
        """Erases an existing field.

        Parameters
        ----------
        group : str
            Name of the group.
        name: str
            Name of the field.

        Raises
        ------
        KeyError
            If the named field does not exist.

        LLCATS: FIELDCR
        """
        try:
            ds = getattr(self, "at_" + loc)
        except AttributeError:
            raise KeyError(loc)
        ds._ds = ds._ds.drop(name)

    def add_empty(self, *args, **kwds):
        """add_empty(name, at='node', units='-', clobber=False)

        Create and add an uninitialized array of values to the field.

        Create a new array of the data field size, without initializing
        entries, and add it to the field as *name*. The *units* keyword gives
        the units of the new fields as a string. Remaining keyword arguments
        are the same as that for the equivalent numpy function.

        This method is not valid for the group *grid*.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        at : str, optional
            Grid location to store values. If not given, values are
            assumed to be on `node`.
        units : str, optional
            Optionally specify the units of the field.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        array :
            A reference to the newly-created array.

        See Also
        --------
        numpy.empty : See for a description of optional keywords.
        landlab.field.GraphFields.empty : Equivalent method that
            does not initialize the new array.
        landlab.field.GraphFields.zeros : Equivalent method that
            initializes the data to 0.

        LLCATS: FIELDCR
        """
        if len(args) == 2:
            loc, name = args
        elif len(args) == 1:
            loc, name = kwds.pop("at"), args[0]
        else:
            raise ValueError("number of arguments must be 1 or 2")
        units = kwds.pop("units", "?")
        copy = kwds.pop("copy", False)
        clobber = kwds.pop("clobber", False)
        return self.add_field(
            name,
            self.empty(at=loc, **kwds),
            at=loc,
            units=units,
            copy=copy,
            clobber=clobber,
        )

    def add_ones(self, *args, **kwds):
        """add_ones(name, at='node', units='-', clobber=False)

        Create and add an array of values, initialized to 1, to the field.

        Create a new array of the data field size, filled with ones, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        This method is not valid for the group *grid*.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        at : str, optional
            Grid location to store values. If not given, values are
            assumed to be on `node`.
        units : str, optional
            Optionally specify the units of the field.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        array :
            A reference to the newly-created array.

        See Also
        --------
        numpy.ones : See for a description of optional keywords.
        andlab.field.GraphFields.add_empty : Equivalent method that
            does not initialize the new array.
        andlab.field.GraphFields.add_zeros : Equivalent method that
            initializes the data to 0.

        Examples
        --------
        Add a new, named field to a collection of fields.

        >>> from landlab.field import GraphFields
        >>> field = GraphFields()
        >>> field.new_field_location('node', 4)
        >>> field.add_ones("topographic__elevation", at="node")
        array([ 1.,  1.,  1.,  1.])
        >>> list(field.keys('node'))
        ['topographic__elevation']
        >>> field['node']['topographic__elevation']
        array([ 1.,  1.,  1.,  1.])
        >>> field.at_node['topographic__elevation']
        array([ 1.,  1.,  1.,  1.])

        LLCATS: FIELDCR
        """
        data = self.add_empty(*args, **kwds)
        data.fill(1)
        return data

    def add_zeros(self, *args, **kwds):
        """add_zeros(name, at='node', units='-', clobber=False)

        Create and add an array of values, initialized to 0, to the field.

        Create a new array of the data field size, filled with zeros, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        at : str, optional
            Grid location to store values. If not given, values are
            assumed to be on `node`.
        units : str, optional
            Optionally specify the units of the field.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        array :
            A reference to the newly-created array.

        See also
        --------
        numpy.zeros : See for a description of optional keywords.
        landlab.field.GraphFields.add_empty : Equivalent method that
            does not initialize the new array.
        landlab.field.GraphFields.add_ones : Equivalent method that
            initializes the data to 1.

        LLCATS: FIELDCR
        """
        data = self.add_empty(*args, **kwds)
        data.fill(0)
        return data

    def add_full(self, *args, **kwds):
        """Create and add an array of values, fill with `fill_value`.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        fill_value : scalar
            Fill value.
        at : str, optional
            Grid location to store values. If not given, values are
            assumed to be on `node`.
        units : str, optional
            Optionally specify the units of the field.
        copy : boolean, optional
            If True, add a *copy* of the array to the field. Otherwise save add
            a reference to the array.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        array :
            A reference to the newly-created array.

        LLCATS: FIELDCR
        """
        if len(args) == 3:
            at, name, fill_value = args
        elif len(args) == 2:
            at = kwds.pop("at", "node")
            name, fill_value = args
        else:
            raise ValueError("number of arguments must be 2 or 3")

        data = self.add_empty(name, at=at, **kwds)
        data.fill(fill_value)
        return data
