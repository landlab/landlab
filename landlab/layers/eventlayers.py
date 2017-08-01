from __future__ import print_function

import os
import sys

import numpy as np


def deposit_or_erode(layers, n_layers, dz):
    from .ext.eventlayers import deposit_or_erode as _deposit_or_erode

    layers = layers.reshape((layers.shape[0], -1))
    try:
        dz = dz.reshape((layers.shape[1], ))
    except AttributeError:
        dz = np.broadcast_to(dz, (layers.shape[1], ))
    finally:
        dz = np.asfarray(dz)

    _deposit_or_erode(layers, n_layers, dz)


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
    >>> bigger_x.shape
    (10, 3)

    >>> bigger_x = resize_array(x, 4, exact=True)
    >>> bigger_x.shape
    (4, 3)
    """
    newsize = int(newsize)
    allocated = array.shape[0]

    if newsize <= allocated:
        return array

    if exact:
        new_allocated = newsize
    else:
        new_allocated = (newsize >> 3) + 6 + newsize

    larger_array = np.empty((new_allocated, ) + array.shape[1:],
                            dtype=array.dtype)
    larger_array[:allocated] = array

    return larger_array


def allocate_layers_for(array, nlayers, nstacks):
    array = np.asarray(array)

    if array.ndim > 0 and len(array) != nstacks:
        values_per_stack = array.shape
    else:
        values_per_stack = array.shape[1:]
        
    return np.empty((nlayers, nstacks) + values_per_stack , dtype=array.dtype)
    

class EventLayersMixIn(object):

    """MixIn that adds a layers attribute to a ModelGrid."""

    @property
    def layers(self):
        """Layers for each cell."""
        try:
            self._layers
        except AttributeError:
            self._layers = EventLayerStack(self.number_of_cells)
        finally:
            return self._layers


class EventLayerStack(object):

    """Track layers where each event is its own layer.
    
    Parameters
    ----------
    nstacks : int
        Number of layer stacks to track.

    Examples
    --------
    >>> from landlab.layers.eventlayers import EventLayerStack

    Create an empty layer stack with 5 stacks.

    >>> layers = EventLayerStack(5)
    >>> layers.nstacks
    5
    >>> layers.nlayers
    0

    Add a layer with a uniform thickness.

    >>> layers.add(1.5)
    >>> layers.nlayers
    1
    >>> layers.dz
    array([[ 1.5,  1.5,  1.5,  1.5,  1.5]])

    Add a second layer with uneven thickness.

    >>> layers.add([1., 2., .5, 5., 0.])
    >>> layers.dz
    array([[ 1.5,  1.5,  1.5,  1.5,  1.5],
           [ 1. ,  2. ,  0.5,  5. ,  0. ]])

    Adding a layer with negative thickness will remove
    existing layers for the top of the stack.

    >>> layers.add(-1)
    >>> layers.dz
    array([[ 1.5,  1.5,  1. ,  1.5,  0.5],
           [ 0. ,  1. ,  0. ,  4. ,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ]])
    """

    def __init__(self, nstacks, allocated=0):
        self._nlayers = 0
        self._nstacks = nstacks

        self._attrs = dict()

        dims = (self.nlayers, self.nstacks)
        self._attrs['_dz'] = np.empty(dims , dtype=float)
        self._resize(allocated, exact=True)

    def __getitem__(self, name):
        return self._attrs[name][:self.nlayers]

    def __str__(self):
        lines = [
            "nlayers: {nlayers}",
            "nstacks: {nstacks}",
            "tracking: {attrs}",
        ]
        return os.linesep.join(lines).format(
            nlayers=self.nlayers,
            nstacks=self.nstacks,
            attrs=', '.join(self.tracking) or 'null')

    def __repr__(self):
        return 'EventLayerStack({nstacks})'.format(nstacks=self.nstacks)

    @property
    def tracking(self):
        """Layer properties being tracked.
        
        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayerStack
        >>> layers = EventLayerStack(3)
        >>> layers.tracking
        []
        >>> layers.add(1., age=1.)
        >>> layers.tracking
        ['age']
        """
        return [name for name in self._attrs if not name.startswith('_')]

    def _setup_layers(self, **kwds):
        dims = (self.allocated, self.nstacks)
        for name, array in kwds.items():
            self._attrs[name] = allocate_layers_for(array, *dims)

    @property
    def nstacks(self):
        """Number of stacks."""
        return self._nstacks

    @property
    def thickness(self):
        """Total thickness of the columns.
        
        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayerStack

        Initially there are no layers so the total thickness is 0.

        >>> layers = EventLayerStack(3)
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
        
        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayerStack

        Initially there are no layers so the elevation to the top
        is 0.

        >>> layers = EventLayerStack(3)
        >>> layers.z
        array([], shape=(0, 3), dtype=float64)

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
        
        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayerStack

        Initially there are no layers so there are not thicknesses.

        >>> layers = EventLayerStack(3)
        >>> layers.dz
        array([], shape=(0, 3), dtype=float64)

        Now add two layers, the first of uniform thickness and the
        second non-uniform and with some erosion.

        >>> layers.add(15.)
        >>> layers.add([1., -1., 2.])
        >>> layers.dz
        array([[ 15.,  14.,  15.],
               [  1.,   0.,   2.]])
        """
        return self._attrs['_dz'][:self.nlayers]

    @property
    def nlayers(self):
        """Total number of layers.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayerStack

        >>> layers = EventLayerStack(3)
        >>> layers.nlayers
        0

        >>> layers.add(15.)
        >>> layers.add([1., -1., 2.])
        >>> layers.nlayers
        2
        """
        return self._nlayers

    @property
    def allocated(self):
        """Total number of allocated layers.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayerStack

        >>> layers = EventLayerStack(3)
        >>> layers.nlayers, layers.allocated
        (0, 0)

        >>> layers.add(15.)
        >>> layers.nlayers, layers.allocated
        (1, 7)
        >>> for _ in range(layers.allocated): layers.add(0.)
        >>> layers.nlayers, layers.allocated
        (8, 15)

        If you know how many layers you will ultimately have, you
        can allocated enough memory for them when you create your
        layer stacks.

        >>> layers = EventLayerStack(3, allocated=15)
        >>> layers.nlayers, layers.allocated
        (0, 15)

        >>> layers.add(15.)
        >>> layers.nlayers, layers.allocated
        (1, 15)
        >>> for _ in range(layers.allocated): layers.add(0.)
        >>> layers.nlayers, layers.allocated
        (16, 24)
        """
        return self._attrs['_dz'].shape[0]

    def add(self, dz, **kwds):
        """Add a layer to the stacks.

        Parameters
        ----------
        dz : float or array_like
            Thickness to add to each stack.

        Examples
        --------
        >>> from landlab.layers.eventlayers import EventLayerStack

        Create an empty layer stack with 3 stacks.

        >>> layers = EventLayerStack(3)
        >>> layers.nlayers
        0

        To add a layer of uniform thickness to every stack.

        >>> layers.add(1.5)
        >>> layers.nlayers
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

        >>> layers = EventLayerStack(3)
        >>> layers.add(1., age=3.)
        >>> layers.dz
        array([[ 1.,  1.,  1.]])
        >>> layers['age']
        array([[ 3.,  3.,  3.]])
        >>> layers.add(2., age=6.)
        >>> layers['age']
        array([[ 3.,  3.,  3.],
               [ 6.,  6.,  6.]])
        """
        if self.nlayers == 0:
            self._setup_layers(**kwds)

        self._add_empty_layer()

        deposit_or_erode(self._attrs['_dz'], self.nlayers, dz)

        for name in kwds:
            try:
                self[name][-1] = kwds[name]
            except KeyError:
                print('{0} is not being tracked. Ignoring'.format(name),
                      file=sys.stderr)

    def _add_empty_layer(self):
        """Add a new empty layer to the stacks."""
        if self.nlayers >= self.allocated:
            self._resize(self.allocated + 1)

        self._nlayers += 1
        for name in self._attrs:
            self._attrs[name][self.nlayers - 1] = 0.

    def _resize(self, newsize, exact=False):
        """Allocate more memory for the layers."""
        for name in self._attrs:
            self._attrs[name] = resize_array(self._attrs[name], newsize,
                                             exact=exact)
