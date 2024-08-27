import numpy as np

from landlab.layers.eventlayers import EventLayers
from landlab.layers.eventlayers import _deposit_or_erode
from landlab.layers.eventlayers import _get_surface_index


class MaterialLayersMixIn:
    """MixIn that adds a MaterialLayers attribute to a ModelGrid."""

    @property
    def material_layers(self):
        """MaterialLayers for each cell."""
        try:
            self._material_layers
        except AttributeError:
            self._material_layers = MaterialLayers(self.number_of_cells)
        return self._material_layers


class MaterialLayers(EventLayers):
    """Track MaterialLayers where each layer has some material in it.

    MaterialLayers are meant to represent a layered object in which each layer
    has some material in it. If erosion occurs, no new layer is created. These
    layers stand in contrast to the EventLayers for which each event is
    represented by a layer.

    MaterialLayers is likely a more memory efficent data structure than
    EventLayers as it does not record erosion as an array of zeros.

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
    array([[1.5,  1.5,  1.5,  1.5,  1.5]])

    MaterialLayers will combine layers if they have the same attributes.
    Adding a second layer with uneven thickness. Will increment the
    first layers thickness. This stands in contrast with EventLayers
    which will track each addition as a separate entry in the layers
    datastructure.

    >>> layers.add([1.0, 2.0, 3.0, 5.0, 0.0])
    >>> layers.dz
    array([[2.5,  3.5,  4.5,  6.5,  1.5]])

    Adding a layer with negative thickness will remove
    material from the layers. Unlike EventLayers, it will not add a
    layer of zeros that represent an event with no deposition.

    >>> layers.add(-1)
    >>> layers.dz
    array([[1.5,  2.5,  3.5,  5.5,  0.5]])

    Get the index value of the layer within each stack
    at the topographic surface.

    >>> layers.surface_index
    array([0, 0, 0, 0, 0])

    See the example in the ``add`` method to learn how MaterialLayers
    behaves if material properties are also tracked.
    """

    def add(self, dz, **kwds):
        """Add a layer to the  MaterialLayers stacks.

        Parameters
        ----------
        dz : float or array_like
            Thickness to add to each stack.

        Examples
        --------
        >>> from landlab.layers.materiallayers import MaterialLayers

        Create an empty layer stack with 3 stacks.

        >>> layers = MaterialLayers(3)
        >>> layers.number_of_layers
        0

        To add a layer of uniform thickness to every stack.

        >>> layers.add(1.5)
        >>> layers.number_of_layers
        1
        >>> layers.dz
        array([[1.5,  1.5,  1.5]])

        Add a second layer with uneven thickness.

        >>> layers.add([1.0, 2.0, 0.5])
        >>> layers.dz
        array([[2.5,  3.5,  2. ]])

        Because the attributes of this layer and the previous layer
        are the same (e.g. they don't exist), MaterialLayer will combine
        them. This is the primary difference between MaterialLayers and
        EventLayers.

        Adding a layer with negative thickness will remove material from
        the top of the stack.

        >>> layers.add(-1)
        >>> layers.dz
        array([[1.5,  2.5,  1. ]])
        >>> layers.number_of_layers
        1

        Use keywords to track properties of each layer. For instance,
        here we create a new stack and add a layer with a particular
        *type* and a particular *size*. You can access the layer properties as
        if the object were a dictionary.

        >>> layers = MaterialLayers(3)
        >>> layers.add(1.0, type=3.0, size="sand")
        >>> layers.dz
        array([[1.,  1.,  1.]])
        >>> layers["type"]
        array([[3.,  3.,  3.]])

        As you can see, there is no rule that says you can't use a string as
        the value of an attribute.

        Adding a layer with the same attributes as the entire surface of the
        MaterialLayers will result in the layers being combined.

        >>> layers.add(1.0, type=3.0, size="sand")
        >>> layers.add([2, -1, 0], type=3.0, size="sand")
        >>> layers.dz
        array([[4.,  1.,  2.]])

        Adding material with different attributes results in the creation of
        a new layer.

        >>> layers.add(2.0, type=6.0, size="sand")
        >>> layers.dz
        array([[4.,  1.,  2.],
               [2.,  2.,  2.]])
        >>> layers["type"]
        array([[3.,  3.,  3.],
               [6.,  6.,  6.]])
        >>> np.all(
        ...     layers["size"] == [["sand", "sand", "sand"], ["sand", "sand", "sand"]]
        ... )
        True

        Attributes for each layer will exist even if part the the layer is
        associated with erosion.

        >>> layers.add([-2, -1, 1], type=8.0, size="gravel")
        >>> layers.dz
        array([[4.,  1.,  2.],
               [0.,  1.,  2.],
               [0.,  0.,  1.]])
        >>> layers["type"]
        array([[3.,  3.,  3.],
               [6.,  6.,  6.],
               [8.,  8.,  8.]])

        To get the values at the surface of the layer stack:

        >>> layers.get_surface_values("type")
        array([3.,  6.,  8.])

        Removing enough material such that an entire layer's
        thickness is no longer present, results in that layer
        no longer being tracked. This is another difference
        between MaterialLayers and EventLayers.

        >>> layers.add([0.0, 0.0, -1.0])
        >>> layers.dz
        array([[4.,  1.,  2.],
               [0.,  1.,  2.]])
        >>> layers["type"]
        array([[3.,  3.,  3.],
               [6.,  6.,  6.]])
        >>> np.all(
        ...     layers["size"] == [["sand", "sand", "sand"], ["sand", "sand", "sand"]]
        ... )
        True
        >>> layers.number_of_layers
        2

        If attributes (like age and size in this example) are tracked, a layer
        will be combined with the surface layer only if all attributes are the
        same across the entire layer. Right now, the surface values vary.

        >>> layers.get_surface_values("type")
        array([3.,  6.,  6.])
        >>> np.all(layers.get_surface_values("size") == ["sand", "sand", "sand"])
        True

        Since the surface has different types, adding material will create a
        new layer.

        >>> layers.add(3.0, type=6.0, size="sand")
        >>> layers.dz
        array([[4.,  1.,  2.],
               [0.,  1.,  2.],
               [3.,  3.,  3.]])
        >>> layers["type"]
        array([[3.,  3.,  3.],
               [6.,  6.,  6.],
               [6.,  6.,  6.]])
        >>> layers.number_of_layers
        3

        But now, the entire surface has the qualities of type = 6. and size =
        'sand', so layers will be combined. This even works if the thickness of
        the new layer includes both erosion and deposition.

        >>> layers.add([-3.5, 0.0, 2.0], type=6.0, size="sand")
        >>> layers.dz
        array([[3.5,  1. ,  2. ],
               [0. ,  1. ,  2. ],
               [0. ,  3. ,  5. ]])
        >>> layers["type"]
        array([[3.,  3.,  3.],
               [6.,  6.,  6.],
               [6.,  6.,  6.]])
        >>> layers.number_of_layers
        3
        """
        dz = np.asarray(dz)

        if self.number_of_layers == 0:
            self._setup_layers(**kwds)

        compatible = self.number_of_layers > 0 and self.is_compatible(dz, **kwds)

        if not compatible:
            self._add_empty_layer()

        _deposit_or_erode(self._attrs["_dz"], self.number_of_layers, dz)
        _get_surface_index(
            self._attrs["_dz"], self.number_of_layers, self._surface_index
        )

        self._remove_empty_layers()

        if not compatible:
            for name in kwds:
                self[name][-1] = kwds[name]

    def _remove_empty_layers(self):
        number_of_filled_layers = self.surface_index.max() + 1
        if number_of_filled_layers < self.number_of_layers:
            self._number_of_layers = number_of_filled_layers

    def is_compatible(self, dz, **kwds):
        """Check if a new layer is compatible with the existing top layer.

        Parameters
        ----------
        dz : float or array_like
            Thickness to add to each stack.

        Returns
        -------
        bool
            ``True`` if the new layer is compatible, otherwise ``False``.
        """
        where_deposition = np.where(dz > 0.0)[0]
        if len(where_deposition) > 0:
            for name in kwds:
                try:
                    is_compatible = self[name][self.surface_index] == kwds[name]
                except KeyError as exc:
                    raise ValueError(
                        f"{name!r} is not being tracked. Error in adding."
                    ) from exc

                if not np.all(is_compatible[where_deposition]):
                    return False
        return True
