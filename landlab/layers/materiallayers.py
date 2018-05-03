from __future__ import print_function

import os
import sys

import numpy as np

from landlab.layers.eventlayers import (EventLayers,
                                        _deposit_or_erode,
                                        _get_surface_index,
                                        resize_array,
                                        _allocate_layers_for)



class MaterialLayersMixIn(object):

    """MixIn that adds a layers attribute to a ModelGrid."""

    @property
    def layers(self):
        """Layers for each cell."""
        try:
            self._layers
        except AttributeError:
            self._layers = MaterialLayers(self.number_of_cells)
        finally:
            return self._layers


class MaterialLayers(EventLayers):

    """Track layers where each layer has some material in it.

    Parameters
    ----------
    number_of_stacks : int
        Number of layer stacks to track.

    Examples
    --------
    >>> from landlab.layers.materiallayers import MaterialLayers

    Create an empty layer stack with 5 stacks.

    >>> layers = MaterialLayers(5)
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
    existing layers for the top of the stack.

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
        self._surface_index = np.empty(number_of_stacks, dtype=int)
        self._attrs = dict()

        dims = (self.number_of_layers, self.number_of_stacks)
        self._attrs['_dz'] = np.empty(dims , dtype=float)
        self._resize(allocated, exact=True)

    def __getitem__(self, name):
        return self._attrs[name][:self.number_of_layers]

    def __setitem__(self, name, values):
        dims = (self.allocated, self.number_of_stacks)
        self._attrs[name] = _allocate_layers_for(values.flatten()[0], *dims)
        self._attrs[name][:self.number_of_layers] = values

    def __str__(self):
        lines = [
            "number_of_layers: {number_of_layers}",
            "number_of_stacks: {number_of_stacks}",
            "tracking: {attrs}",
        ]
        return os.linesep.join(lines).format(
            number_of_layers=self.number_of_layers,
            number_of_stacks=self.number_of_stacks,
            attrs=', '.join(self.tracking) or 'null')

    def __repr__(self):
        return 'EventLayers({number_of_stacks})'.format(
            number_of_stacks=self.number_of_stacks)

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
        return [name for name in self._attrs if not name.startswith('_')]

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
        return self._attrs['_dz'][:self.number_of_layers]

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
        return self._attrs['_dz'].shape[0]

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

        >>> layers.surface_values('age')
        array([ 3.,  6.,  8.])
        """
        if self.number_of_layers == 0:
            self._setup_layers(**kwds)

        self._add_empty_layer()

        _deposit_or_erode(self._attrs['_dz'], self.number_of_layers, dz)

        for name in kwds:
            try:
                self[name][-1] = kwds[name]
            except KeyError:
                print('{0} is not being tracked. Ignoring'.format(name),
                      file=sys.stderr)

    @property
    def surface_index(self):
        _get_surface_index(self._attrs['_dz'], self.number_of_layers, self._surface_index)
        return self._surface_index

    def surface_values(self, name):
        _get_surface_index(self._attrs['_dz'], self.number_of_layers, self._surface_index)

        return self._attrs[name][self._surface_index,
                                 np.arange(self._number_of_stacks)]

    def _add_empty_layer(self):
        """Add a new empty layer to the stacks."""
        if self.number_of_layers >= self.allocated:
            self._resize(self.allocated + 1)

        self._number_of_layers += 1
        for name in self._attrs:
            self._attrs[name][self.number_of_layers - 1] = 0.

    def _resize(self, newsize, exact=False):
        """Allocate more memory for the layers."""
        for name in self._attrs:
            self._attrs[name] = resize_array(self._attrs[name], newsize,
                                             exact=exact)
