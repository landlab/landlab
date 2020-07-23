import os

import numpy as np


def _deposit_or_erode(layers, n_layers, dz):
    """Update the array that contains layers with deposition or erosion.

    This function operates on the entire array that contain the layers (active,
    and allocated but not yet active). The number of active layers includes the
    layer that is currently being added. Thus the row with the index
    ``n_layers - 1`` is the layer that is currently being added as an active
    layer.

    Note that in EventLayers, layers represent an event, and not necessarily
    material.

    This means that if only erosion occurs, the array elements in the row with
    the index ``n_layers - 1`` will be set to zero and thickness will be
    removed from lower layers. Note that lower layers have smaller row indicies
    as the bottom of the layer stack has row index zero.

    If deposition occurs, the array elements in the row with index
    ``n_layers - 1`` will be set to the value of dz.

    Parameters
    ----------
    layers : ndarray of shape `(n_layers, n_nodes)`
        Array of layer thicknesses. This array is the datastructure that
        contains all allocated layers, active or inactive.
    n_layers : int
        Number of active layers.
    dz : ndarray of shape `(n_nodes, )`
        Thickness of the new layer. Negative thicknesses mean
        erode the top-most layers.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.layers.eventlayers import _deposit_or_erode

    First we create a numpy array allocated to contain layers. We fill it with
    -1. These -1.'s do not represent negative layer thickness. In EventLayers
    this array is created with np.empty, but that creates different numbers
    every time and doesn't work for testing.

    >>> allocated_layers_array = np.full((4, 3), 0.)

    Next we add a layer with spatially variable thickness. We specify that the
    number of active layers (including the one being added) is 1.

    >>> dz = np.array([1., 2., 3.])
    >>> _deposit_or_erode(allocated_layers_array, 1, dz)
    >>> allocated_layers_array
    array([[ 1.,  2.,  3.],
           [ 0.,  0.,  0.],
           [ 0.,  0.,  0.],
           [ 0.,  0.,  0.]])

    As you can see, this changes the value of the first row in the array. The
    remainder of the array represents space in the datatastructure that has
    been allocated to contain layers, but does not yet contain active layers.

    Next we add a layer of thickness 1. To do this, we now need to specify that
    the number of active layers is 2.

    >>> dz = np.array([1., 1., 1.])
    >>> _deposit_or_erode(allocated_layers_array, 2, dz)
    >>> allocated_layers_array
    array([[ 1.,  2.,  3.],
           [ 1.,  1.,  1.],
           [ 0.,  0.,  0.],
           [ 0.,  0.,  0.]])

    Finally, we do some erosion. We specify that the number of active layers is
    3 and give a spatially variable field of erosion and deposition.

    >>> _deposit_or_erode(allocated_layers_array, 3, [1., -1., -2.])
    >>> allocated_layers_array
    array([[ 1.,  2.,  2.],
           [ 1.,  0.,  0.],
           [ 1.,  0.,  0.],
           [ 0.,  0.,  0.]])

    >>> _deposit_or_erode(allocated_layers_array, 3, [1., -1., -2.])
    >>> allocated_layers_array
    array([[ 1.,  1.,  0.],
           [ 1.,  0.,  0.],
           [ 2.,  0.,  0.],
           [ 0.,  0.,  0.]])
    """
    from .ext.eventlayers import deposit_or_erode

    layers = layers.reshape((layers.shape[0], -1))
    try:
        dz = dz.reshape((layers.shape[1],))
    except (AttributeError, ValueError):
        dz = np.broadcast_to(dz, (layers.shape[1],))
    finally:
        dz = np.asfarray(dz)

    deposit_or_erode(layers, n_layers, dz)


def _get_surface_index(layers, n_layers, surface_index):
    """Get index within each stack of the layer at the topographic surface.

    Parameters
    ----------
    layers : ndarray of shape `(n_layers, n_nodes)`
        Array of layer thicknesses. This array is the datastructure that
        contains all allocated layers, active or inactive.
    n_layers : int
        Number of active layers.
    surface_index : ndarray of shape `(n_nodes, )`

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.layers.eventlayers import (_deposit_or_erode,
    ...                                         _get_surface_index)

    >>> layers = np.full((5, 3), 1.)
    >>> dz = np.array([-1., -2., -3.])

    Note here, if you are very confused by the use of ``_deposit_or_erode``
    we recommend you read the docstring associated with that function.

    >>> _deposit_or_erode(layers, 5, dz)
    >>> layers
    array([[ 1.,  1.,  1.],
           [ 1.,  1.,  1.],
           [ 1.,  1.,  0.],
           [ 1.,  0.,  0.],
           [ 0.,  0.,  0.]])

    >>> surface_index = np.empty(3, dtype=int)
    >>> _get_surface_index(layers, 5, surface_index)
    >>> surface_index
    array([3, 2, 1])
    """
    from .ext.eventlayers import get_surface_index

    layers = layers.reshape((layers.shape[0], -1))

    get_surface_index(layers, n_layers, surface_index)


def _reduce_matrix(array, step, reducer):
    """Combine rows of a 2D matrix.

    Parameters
    ----------
    array : ndarray, shape (m, n)
        Matrix to reduce.
    step : int
        Number of rows in each block to reduce.
    reducer : ufunc
        Function to use for combining rows.

    Returns
    -------
    ndarray
        Matrix with rows combined.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.layers.eventlayers import _reduce_matrix
    >>> array = np.arange(12).reshape((4, 3))
    >>> array
    array([[ 0,  1,  2],
           [ 3,  4,  5],
           [ 6,  7,  8],
           [ 9, 10, 11]])
    >>> _reduce_matrix(array, 4, np.sum)
    array([[18, 22, 26]])
    >>> _reduce_matrix(array, 2, np.sum)
    array([[ 3,  5,  7],
           [15, 17, 19]])
    """
    return reducer(array.reshape((-1, step) + array.shape[1:]), axis=1).reshape(
        (-1,) + array.shape[1:]
    )


class _BlockSlice:
    """Slices that divide a matrix into equally sized blocks."""

    def __init__(self, *args):
        """_BlockSlice([start], stop, [step])"""
        if len(args) > 3:
            raise TypeError(
                "_BlockSlice expected at most 3 arguments, got {0}".format(len(args))
            )

        self._args = tuple(args)

        self._start = 0
        self._stop = None
        self._step = None

        if len(args) == 1:
            self._stop = args[0]
        elif len(args) == 2:
            self._start, self._stop = args
        elif len(args) == 3:
            self._start, self._stop, self._step = args
            self._step = min(self._stop - self._start, self._step)

        if self._stop is not None and self._stop < self._start:
            raise ValueError(
                "stop ({0}) must be greater than start ({1})".format(
                    self._stop, self._start
                )
            )

    def __repr__(self):
        return "_BlockSlice({0})".format(", ".join([repr(arg) for arg in self._args]))

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def step(self):
        return self._step

    def indices(self, n_rows):
        """Row indices to blocks within a matrix.

        Parameters
        ----------
        n_rows : int
            The number of rows in the matrix.

        Returns
        -------
        (start, stop, step)
            Tuple of (int* that gives the row of the first block, row of the
            last block, and the number of rows in each block.

        Examples
        --------
        >>> from landlab.layers.eventlayers import _BlockSlice

        The default is one single block that encomapses all the rows.

        >>> _BlockSlice().indices(4)
        (0, 4, 4)

        >>> _BlockSlice(3).indices(4)
        (0, 3, 3)

        >>> _BlockSlice(1, 3).indices(4)
        (1, 3, 2)

        >>> _BlockSlice(1, 7, 2).indices(8)
        (1, 7, 2)
        """
        start, stop, step = self.start, self.stop, self.step
        if stop is None:
            stop = n_rows

        start, stop, _ = slice(start, stop).indices(n_rows)

        if step is None:
            step = stop - start

        if step != 0 and (stop - start) % step != 0:
            stop = (stop - start) // step * step + start

        return start, stop, step


def _valid_keywords_or_raise(kwds, required=(), optional=()):
    """Check for valid keyword arguments.

    Parameters
    ----------
    kwds : iterable of str
        Keywords to check for validity.
    required : iterable of str, optional
        Keywords that are required.
    optional : iterable of str, optional
        Keywords that are optional.

    Examples
    --------
    >>> from landlab.layers.eventlayers import _valid_keywords_or_raise
    >>> _valid_keywords_or_raise(["foo"], optional=["foo", "bar"])
    >>> _valid_keywords_or_raise(["foo"], required=["foo", "bar"])
    Traceback (most recent call last):
        ...
    TypeError: missing keyword arguments ('bar')
    >>> _valid_keywords_or_raise(["baz"], optional=["foo", "bar"])
    Traceback (most recent call last):
        ...
    TypeError: invalid keyword arguments ('baz' not in {'bar', 'foo'})
    """
    keys = set(kwds)
    required = set(required)
    optional = required | set(optional)

    unknown = keys - optional
    if unknown:
        raise TypeError(
            "invalid keyword arguments ({0} not in {{{1}}})".format(
                ", ".join(sorted([repr(name) for name in unknown])),
                ", ".join(sorted([repr(name) for name in optional])),
            )
        )

    missing = required - keys
    if missing:
        raise TypeError(
            "missing keyword arguments ({0})".format(
                ", ".join(sorted([repr(name) for name in missing]))
            )
        )


def resize_array(array, newsize, exact=False):
    """Increase the size of an array, leaving room to grow.

    Parameters
    ----------
    array : ndarray
        The array to resize.
    newsize : int
        Size of the zero-th dimension of the resized array.
    exact : bool, optional
        Should the new array have the exact size provided or
        at least that size.

    Returns
    -------
    ndarray
        Copy of the input array resized.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.layers.eventlayers import resize_array

    >>> x = np.arange(6)
    >>> bigger_x = resize_array(x, 10)
    >>> bigger_x.size
    17
    >>> np.all(x[:6] == bigger_x[:6])
    True
    >>> x is bigger_x
    False

    >>> x = np.arange(6).reshape((2, 3))
    >>> bigger_x = resize_array(x, 4)
    >>> bigger_x.shape == (10, 3)
    True

    >>> bigger_x = resize_array(x, 4, exact=True)
    >>> bigger_x.shape == (4, 3)
    True
    """
    newsize = int(newsize)
    allocated = array.shape[0]

    if newsize <= allocated:
        return array

    if exact:
        new_allocated = newsize
    else:
        new_allocated = (newsize >> 3) + 6 + newsize

    larger_array = np.empty((new_allocated,) + array.shape[1:], dtype=array.dtype)
    larger_array[:allocated] = array

    return larger_array


def _allocate_layers_for(array, number_of_layers, number_of_stacks):
    """Allocate a layer matrix.

    Parameters
    ----------
    array : number or ndarray
        Array of layer properties to track.
    number_of_layers : int
        Number of layers to allocate.
    number_of_stacks : int
        Number of stacks to allocate.

    Returns
    -------
    ndarray of size `(number_of_layers, number_of_stacks, values_per_stack)`
        Newly allocated matrix for storing layer properties.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.layers.eventlayers import _allocate_layers_for

    >>> layers = _allocate_layers_for(3, 2, 4)
    >>> layers.shape == (2, 4)
    True
    >>> layers.dtype.kind == 'i'
    True

    >>> layers = _allocate_layers_for(np.zeros(4), 2, 4)
    >>> layers.shape == (2, 4)
    True
    >>> layers.dtype.kind == 'f'
    True

    >>> layers = _allocate_layers_for(np.zeros(2), 2, 4)
    >>> layers.shape == (2, 4, 2)
    True
    >>> layers.dtype.kind == 'f'
    True
    """
    array = np.asarray(array)

    if array.ndim > 0 and len(array) != number_of_stacks:
        values_per_stack = array.shape
    else:
        values_per_stack = array.shape[1:]

    return np.empty(
        (number_of_layers, number_of_stacks) + values_per_stack, dtype=array.dtype
    )


class EventLayersMixIn:

    """MixIn that adds a EventLayers attribute to a ModelGrid."""

    @property
    def event_layers(self):
        """EventLayers for each cell."""
        try:
            self._event_layers
        except AttributeError:
            self._event_layers = EventLayers(self.number_of_cells)
        finally:
            return self._event_layers

    @property
    def at_layer(self):
        """EventLayers for each cell."""
        return self.event_layers


class EventLayers:

    """Track EventLayers where each event is its own layer.

    EventLayers are meant to represent a layered object in which each layer
    represents a event. Thus they are likely the most appropriate tool to use
    if the user is interested in chronostratigraphy. If erosion occurs, a new
    layer with zero thickness is created. Thus, EventLayers may not be the most
    memory efficent layers datastructure.

    EventLayers exists in contrast to the MaterialLayers object which does not
    make a new layer if only erosion occurs and if the attributes of the new
    layer are equivalent to the attributes of the material at the surface of the
    layer stack.

    Attributes
    ----------
    allocated
    dz
    thickness
    number_of_layers
    number_of_stacks
    surface_index
    z

    Methods
    -------
    add
    get_surface_values

    Parameters
    ----------
    number_of_stacks : int
        Number of layer stacks to track.

    Examples
    --------
    >>> from landlab.layers.eventlayers import EventLayers

    Create an empty layer stack with 5 stacks.

    >>> layers = EventLayers(5)
    >>> layers.number_of_stacks
    5
    >>> layers.number_of_layers
    0

    Add a layer with a uniform thickness.

    >>> layers.add(1.5)
    >>> layers.number_of_layers
    1
    >>> layers.dz
    array([[ 1.5,  1.5,  1.5,  1.5,  1.5]])

    Add a second layer with uneven thickness.

    >>> layers.add([1., 2., .5, 5., 0.])
    >>> layers.dz
    array([[ 1.5,  1.5,  1.5,  1.5,  1.5],
           [ 1. ,  2. ,  0.5,  5. ,  0. ]])

    Adding a layer with negative thickness will remove
    existing layers for the top of the stack. Note that
    this will create a new layer with thickness zero
    that represents this 'event'. If instead your
    application would prefer that no new row is added to
    the layers datastructure, you may want to consider
    the MaterialLayers object.

    >>> layers.add(-1)
    >>> layers.dz
    array([[ 1.5,  1.5,  1. ,  1.5,  0.5],
           [ 0. ,  1. ,  0. ,  4. ,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ]])

    Get the index value of the layer within each stack
    at the topographic surface.

    >>> layers.surface_index
    array([0, 1, 0, 1, 0])
    """

    def __init__(self, number_of_stacks, allocated=0):
        self._number_of_layers = 0
        self._number_of_stacks = number_of_stacks
        self._surface_index = np.zeros(number_of_stacks, dtype=int)
        self._attrs = dict()

        dims = (self.number_of_layers, self.number_of_stacks)
        self._attrs["_dz"] = np.empty(dims, dtype=float)
        self._resize(allocated, exact=True)

    def __getitem__(self, name):
        return self._attrs[name][: self.number_of_layers]

    def __setitem__(self, name, values):
        dims = (self.allocated, self.number_of_stacks)
        values = np.asarray(values)
        if values.ndim == 1:
            values = np.expand_dims(values, 1)
        values = np.broadcast_to(values, (self.number_of_layers, self.number_of_stacks))
        self._attrs[name] = _allocate_layers_for(values.flatten()[0], *dims)
        self._attrs[name][: self.number_of_layers] = values

    def __iter__(self):
        return (name for name in self._attrs if not name.startswith("_"))

    def __str__(self):
        lines = [
            "number_of_layers: {number_of_layers}",
            "number_of_stacks: {number_of_stacks}",
            "tracking: {attrs}",
        ]
        return os.linesep.join(lines).format(
            number_of_layers=self.number_of_layers,
            number_of_stacks=self.number_of_stacks,
            attrs=", ".join(self.tracking) or "null",
        )

    def __repr__(self):
        return self.__class__.__name__ + "({number_of_stacks})".format(
            number_of_stacks=self.number_of_stacks
        )

    @property
    def tracking(self):
        """Layer properties being tracked.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers
        >>> layers = EventLayers(3)
        >>> layers.tracking
        []
        >>> layers.add(1., age=1.)
        >>> layers.tracking
        ['age']
        """
        return [name for name in self._attrs if not name.startswith("_")]

    def _setup_layers(self, **kwds):
        dims = (self.allocated, self.number_of_stacks)
        for name, array in kwds.items():
            self._attrs[name] = _allocate_layers_for(array, *dims)

    @property
    def number_of_stacks(self):
        """Number of stacks."""
        return self._number_of_stacks

    @property
    def thickness(self):
        """Total thickness of the columns.

        The sum of all layer thicknesses for each stack as an array
        of shape `(number_of_stacks, )`.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        Initially there are no layers so the total thickness is 0.

        >>> layers = EventLayers(3)
        >>> layers.thickness
        array([ 0.,  0.,  0.])

        After adding some layers, the stacks have varying thicknesses.

        >>> layers.add(15.)
        >>> layers.add([1., -1., 2.])
        >>> layers.thickness
        array([ 16.,  14.,  17.])
        """
        return np.sum(self.dz, axis=0)

    @property
    def z(self):
        """Thickness to top of each layer.

        Thickness from the bottom of each stack to the top of each layer
        as an array of shape `(number_of_layers, number_of_stacks)`.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        Initially there are no layers so the elevation to the top
        is 0.

        >>> layers = EventLayers(3)
        >>> layers.z.shape == (0, 3)
        True

        After adding some layers, elevations are to the top of each layer.

        >>> layers.add(15.)
        >>> layers.add([1., -1., 2.])
        >>> layers.dz
        array([[ 15.,  14.,  15.],
               [  1.,   0.,   2.]])
        >>> layers.z
        array([[ 15.,  14.,  15.],
               [ 16.,  14.,  17.]])
        """
        return np.cumsum(self.dz, axis=0)

    @property
    def dz(self):
        """Thickness of each layer.

        The thickness of each layer at each stack as an array of shape
        `(number_of_layers, number_of_stacks)`.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        Initially there are no layers so there are not thicknesses.

        >>> layers = EventLayers(3)
        >>> layers.dz.shape == (0, 3)
        True

        Now add two layers, the first of uniform thickness and the
        second non-uniform and with some erosion.

        >>> layers.add(15.)
        >>> layers.add([1., -1., 2.])
        >>> layers.dz
        array([[ 15.,  14.,  15.],
               [  1.,   0.,   2.]])
        """
        return self._attrs["_dz"][: self.number_of_layers]

    @property
    def number_of_layers(self):
        """Total number of layers.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        >>> layers = EventLayers(3)
        >>> layers.number_of_layers
        0

        >>> layers.add(15.)
        >>> layers.add([1., -1., 2.])
        >>> layers.number_of_layers
        2
        """
        return self._number_of_layers

    @property
    def allocated(self):
        """Total number of allocated layers.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        >>> layers = EventLayers(3)
        >>> layers.number_of_layers
        0
        >>> layers.allocated == 0
        True

        >>> layers.add(15.)
        >>> layers.number_of_layers
        1
        >>> layers.allocated == 7
        True
        >>> for _ in range(layers.allocated): layers.add(0.)
        >>> layers.number_of_layers
        8
        >>> layers.allocated == 15
        True

        If you know how many layers you will ultimately have, you
        can allocated enough memory for them when you create your
        layer stacks.

        >>> layers = EventLayers(3, allocated=15)
        >>> layers.number_of_layers
        0
        >>> layers.allocated == 15
        True

        >>> layers.add(15.)
        >>> layers.number_of_layers
        1
        >>> layers.allocated == 15
        True
        >>> for _ in range(layers.allocated): layers.add(0.)
        >>> layers.number_of_layers
        16
        >>> layers.allocated == 24
        True
        """
        return self._attrs["_dz"].shape[0]

    def add(self, dz, **kwds):
        """Add a layer to the stacks.

        Parameters
        ----------
        dz : float or array_like
            Thickness to add to each stack.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        Create an empty layer stack with 3 stacks.

        >>> layers = EventLayers(3)
        >>> layers.number_of_layers
        0

        To add a layer of uniform thickness to every stack.

        >>> layers.add(1.5)
        >>> layers.number_of_layers
        1
        >>> layers.dz
        array([[ 1.5,  1.5,  1.5]])

        Add a second layer with uneven thickness.

        >>> layers.add([1., 2., .5])
        >>> layers.dz
        array([[ 1.5,  1.5,  1.5],
               [ 1. ,  2. ,  0.5]])

        Adding a layer with negative thickness will remove
        existing layers for the top of the stack.

        >>> layers.add(-1)
        >>> layers.dz
        array([[ 1.5,  1.5,  1. ],
               [ 0. ,  1. ,  0. ],
               [ 0. ,  0. ,  0. ]])

        Use keywords to track properties of each layer. For instance,
        here we create a new stack and add a layer with a particular
        *age*. You can access the layer properties as if the object
        were a dictionary.

        >>> layers = EventLayers(3)
        >>> layers.add(1., age=3.)
        >>> layers.dz
        array([[ 1.,  1.,  1.]])
        >>> layers['age']
        array([[ 3.,  3.,  3.]])
        >>> layers.add(2., age=6.)
        >>> layers['age']
        array([[ 3.,  3.,  3.],
               [ 6.,  6.,  6.]])

        Attributes for each layer will exist even if the the layer is
        associated with erosion.

        >>> layers.add([-2, -1, 1], age=8.)
        >>> layers.dz
        array([[ 1.,  1.,  1.],
               [ 0.,  1.,  2.],
               [ 0.,  0.,  1.]])
        >>> layers['age']
        array([[ 3.,  3.,  3.],
               [ 6.,  6.,  6.],
               [ 8.,  8.,  8.]])

        To get the values at the surface of the layer stack:

        >>> layers.get_surface_values('age')
        array([ 3.,  6.,  8.])
        """
        if self.number_of_layers == 0:
            self._setup_layers(**kwds)

        self._add_empty_layer()

        _deposit_or_erode(self._attrs["_dz"], self.number_of_layers, dz)
        _get_surface_index(
            self._attrs["_dz"], self.number_of_layers, self._surface_index
        )

        for name in kwds:
            try:
                self[name][-1] = kwds[name]
            except KeyError:
                raise ValueError(
                    "EventLayers: {0} is not being tracked. Error in adding.".format(
                        name
                    )
                )

    def reduce(self, *args, **kwds):
        """reduce([start], stop, [step])
        Combine layers.

        Reduce adjacent layers into a single layer.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        Create an empty layer stack with 3 stacks.

        >>> layers = EventLayers(3)
        >>> layers.number_of_layers
        0

        To add a layer of uniform thickness to every stack.

        >>> layers.add(1.5)
        >>> layers.number_of_layers
        1
        >>> layers.dz
        array([[ 1.5,  1.5,  1.5]])

        Add a second layer with uneven thickness.

        >>> layers.add([1., 2., .5])
        >>> layers.dz
        array([[ 1.5,  1.5,  1.5],
               [ 1. ,  2. ,  0.5]])

        Combine all of the layers into a single layer.

        >>> layers.reduce()
        >>> layers.dz
        array([[ 2.5,  3.5,  2. ]])

        Add two additional layers to the top. The bottom-most layer is row
        0, and the two new layers are rows 1 and 2.

        >>> layers.add([1., 2., .5])
        >>> layers.add([1., 2., .5])
        >>> layers.dz
        array([[ 2.5,  3.5,  2. ],
               [ 1. ,  2. ,  0.5],
               [ 1. ,  2. ,  0.5]])

        Combine the two new layers (layers 1 and 2) into a single layer.

        >>> layers.reduce(1, 3)
        >>> layers.dz
        array([[ 2.5,  3.5,  2. ],
               [ 2. ,  4. ,  1. ]])

        >>> layers.add([1., 2., .5])
        >>> layers.add([1., 2., .5])
        >>> layers.dz
        array([[ 2.5,  3.5,  2. ],
               [ 2. ,  4. ,  1. ],
               [ 1. ,  2. ,  0.5],
               [ 1. ,  2. ,  0.5]])

        Combine the middle two layers.

        >>> layers.reduce(1, 3)
        >>> layers.dz
        array([[ 2.5,  3.5,  2. ],
               [ 3. ,  6. ,  1.5],
               [ 1. ,  2. ,  0.5]])
        >>> layers.add([1., 1., 1.])
        >>> layers.dz
        array([[ 2.5,  3.5,  2. ],
               [ 3. ,  6. ,  1.5],
               [ 1. ,  2. ,  0.5],
               [ 1. ,  1. ,  1. ]])

        Combine every two layers (layers 0 and 1 and combined, and layers
        1 and 2 are combined).

        >>> layers.reduce(0, 4, 2)
        >>> layers.dz
        array([[ 5.5,  9.5,  3.5],
               [ 2. ,  3. ,  1.5]])

        When layers are combined, thicknesses are summed but layer attributes
        can be combined in other ways (e.g. max, or mean)

        >>> layers = EventLayers(3)
        >>> layers.add([1, 1, 1], age=0.)
        >>> layers.add([1, 2, 5], age=1.)
        >>> layers.add([2, 2, 2], age=2.)
        >>> layers.reduce(age=np.max)
        >>> layers["age"]
        array([[ 2.,  2.,  2.]])

        >>> layers.add([2, 2, 2], age=3.)
        >>> layers.add([2, 2, 2], age=4.)
        >>> layers.reduce(1, 3, age=np.mean)
        >>> layers["age"]
        array([[ 2. ,  2. ,  2. ],
               [ 3.5,  3.5,  3.5]])
        """
        _valid_keywords_or_raise(kwds, required=self.tracking, optional=self._attrs)

        start, stop, step = _BlockSlice(*args).indices(self._number_of_layers)

        if step <= 1:
            return

        n_blocks = (stop - start) // step
        n_removed = n_blocks * (step - 1)
        for name, array in self._attrs.items():
            middle = _reduce_matrix(array[start:stop, :], step, kwds.get(name, np.sum))
            top = array[stop : self._number_of_layers, :]

            array[start : start + n_blocks, :] = middle
            array[start + n_blocks : start + n_blocks + len(top)] = top

        self._number_of_layers -= n_removed
        self._surface_index[:] -= n_removed

    @property
    def surface_index(self):
        """Index to the top non-empty layer.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayers

        Create an empty layer stack with 5 stacks.

        >>> layers = EventLayers(3)
        >>> layers.surface_index
        array([0, 0, 0])

        Add a layer with a uniform thickness.

        >>> for _ in range(5): layers.add(1.0)
        >>> layers.surface_index
        array([4, 4, 4])

        Add a layer with varying thickness. Negative thickness
        removes thickness from underlying layers, zero thickness adds a
        layer but doesn't change the surface index.

        >>> layers.add([-1.0, 0.0, 1.0])
        >>> layers.surface_index
        array([3, 4, 5])
        """
        return self._surface_index

    def get_surface_values(self, name):
        """Values of a field on the surface layer."""
        return self._attrs[name][self.surface_index, np.arange(self._number_of_stacks)]

    def _add_empty_layer(self):
        """Add a new empty layer to the stacks."""
        if self.number_of_layers >= self.allocated:
            self._resize(self.allocated + 1)

        self._number_of_layers += 1
        self._attrs["_dz"][self.number_of_layers - 1, :] = 0.0
        for name in self._attrs:
            self._attrs[name][self.number_of_layers - 1] = 0.0

    def _resize(self, newsize, exact=False):
        """Allocate more memory for the layers."""
        for name in self._attrs:
            self._attrs[name] = resize_array(self._attrs[name], newsize, exact=exact)
