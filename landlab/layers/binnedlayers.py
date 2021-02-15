import numpy as np

from .eventlayers import EventLayers


class BinnedLayers(EventLayers):
    def __init__(self, number_of_stacks, allocated=0, bin_size=1.0):
        """Layers of a maximum bin size.

        Parameters
        ----------
        number_of_stacks : int
            Number of layer stacks to track.
        bin_size : float or array of float
            The maximum bin size for each layer stack.

        """
        EventLayers.__init__(self, number_of_stacks, allocated=allocated)

        self._bin_size = np.broadcast_to(bin_size, number_of_stacks)

    def add(self, dz, **kwds):
        """Add a layer to the stacks.

        Parameters
        ----------
        dz : float or array_like
            Thickness to add to each stack.
        """
        available = self._bin_size - top_layer_thickness

        to_average = np.clip(dz, None, available)
        new_layer = dz - to_average

        EventLayers.add(to_average, **kwds)
        self.reduce(self.number_of_layers - 1, **kwds)
        EventLayers.add(new_layer, **kwds)
