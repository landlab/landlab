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


def resize_array(array, newsize):
    newsize = int(newsize)
    allocated = array.shape[0]

    if newsize < allocated:
        return array

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

    @property
    def layers(self):
        try:
            self._layers
        except AttributeError:
            self._layers = EventLayerStack(self.number_of_cells)
        finally:
            return self._layers


class EventLayerStack(object):

    def __init__(self, nstacks):
        self._nlayers = 0
        self._nstacks = nstacks

        self._attrs = dict()

        dims = (self.nlayers, self.nstacks)
        self._attrs['_dz'] = np.empty(dims , dtype=float)

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
        return [name for name in self._attrs if not name.startswith('_')]

    def setup_layers(self, **kwds):
        dims = (self.allocated, self.nstacks)
        for name, array in kwds.items():
            self._attrs[name] = allocate_layers_for(array, *dims)

    @property
    def nstacks(self):
        return self._nstacks

    @property
    def thickness(self):
        """Total thickness of the columns."""
        return np.sum(self.dz, axis=0)

    @property
    def depth(self):
        """Burial depth to bottom of each layer."""
        return self.thickness - self.z

    @property
    def z(self):
        """Thickness to top of each layer."""
        return np.cumsum(self.dz, axis=0)

    @property
    def dz(self):
        """Thickness to bottom of each layer."""
        return self._attrs['_dz'][:self.nlayers]

    @property
    def nlayers(self):
        return self._nlayers

    @property
    def allocated(self):
        return self._attrs['_dz'].shape[0]

    def add(self, dz, **kwds):
        if self.nlayers == 0:
            self.setup_layers(**kwds)

        self.add_empty_layer()

        deposit_or_erode(self._attrs['_dz'], self.nlayers, dz)

        for name in kwds:
            try:
                self[name][-1] = kwds[name]
            except KeyError:
                print('{0} is not being tracked. Ignoring'.format(name),
                      file=sys.stderr)

    def add_empty_layer(self):
        if self.nlayers >= self.allocated:
            self.resize(self.allocated + 1)

        self._nlayers += 1
        for name in self._attrs:
            self._attrs[name][self.nlayers - 1] = 0.

    def resize(self, newsize):
        for name in self._attrs:
            self._attrs[name] = resize_array(self._attrs[name], newsize)
