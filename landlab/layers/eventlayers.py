from __future__ import print_function

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


class EventLayersMixIn(object):

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


class EventLayers(object):

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
                msg = "MaterialLayers: {0} is not being tracked. Error in adding.".format(
                    name
                )
                raise ValueError(msg)

    @property
    def surface_index(self):
        return self._surface_index

    def get_surface_values(self, name):
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
