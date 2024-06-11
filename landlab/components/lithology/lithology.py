#!/usr/bin/env python3
"""Create a Lithology object with different properties."""

import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

from landlab import Component
from landlab.layers import EventLayers
from landlab.layers import MaterialLayers
from landlab.utils.return_array import return_array_at_node


class Lithology(Component):
    """Create a Lithology object.

    A Lithology is a three dimentional representation of material operated on
    by landlab components. Material can be removed through erosion or added to
    through deposition. Rock types can have multiple attributes (e.g. age,
    erodability or other parameter values, etc).

    If the tracked properties are model grid fields, they will be updated to
    the surface values of the Lithology. If the properties are not grid fields
    then at-node grid fields will be created with their names. Lithology and
    its derived versions will make a at-node grid field called `rock_type__id`
    to store the rock type id.

    Lithology was designed to be used on its own and to be inherited from and
    improved. Currently one other Lithology variant exists: LithoLayers
    which makes it easy to specify parallel layers of rock with generic layer
    geometries.

    It is constructed by specifying a series of thicknesses and a series of
    rock type IDs. Thicknesses and IDs are both specified in order of closest
    to the surface to furthest from the surface. Thicknesses can either be a
    single value (corresponding to a layer of uniform thickness) or a number-of
    -nodes length array (corresponding to a non-uniform layer).

    Additionally, an attribute dictionary specifies the properties of each
    rock type. This dictionary is expected to have the form of:

    .. code-block:: python

        attrs = {"K_sp": {1: 0.001, 2: 0.0001}, "D": {1: 0.01, 2: 0.001}}

    Where ``'K_sp'`` and ``'D'`` are properties to track, and ``1`` and ``2``
    are rock type IDs. The rock type IDs can be any type that is valid as a
    python dictionary key.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Barnhart, K., Hutton, E., Gasparini, N., Tucker, G. (2018). Lithology: A
    Landlab submodule for spatially variable rock properties. Journal of Open
    Source Software  3(30), 979 - 2. https://dx.doi.org/10.21105/joss.00979

    **Additional References**

    None Listed

    """

    _name = "Lithology"

    _unit_agnostic = True

    _cite_as = """
    @article{barnhart2018lithology,
        title = "Lithology: A Landlab submodule for spatially variable rock properties",
        journal = "Journal of Open Source Software",
        volume = "",
        pages = "",
        year = "2018",
        doi = "10.21105/joss.00979",
        author = {Katherine R. Barnhart and Eric Hutton and Nicole M. Gasparini
                  and Gregory E. Tucker},
    }"""

    _info = {}

    def __init__(
        self,
        grid,
        thicknesses,
        ids,
        attrs,
        layer_type="MaterialLayers",
        dz_advection=0,
        rock_id=None,
    ):
        """Create a new instance of Lithology.

        Parameters
        ----------
        grid : Landlab ModelGrid
        thicknesses : ndarray of shape `(n_layers, )` or `(n_layers, n_nodes)`
            Values of layer thicknesses from surface to depth. Layers do not
            have to have constant thickness. Layer thickness can be zero,
            though the entirety of Lithology must have non-zero thickness.
        ids : ndarray of shape `(n_layers, )` or `(n_layers, n_nodes)`
            Values of rock type IDs corresponding to each layer specified in
            **thicknesses**. A single layer may have multiple rock types if
            specified by the user.
        attrs : dict
            Rock type property dictionary. See class docstring for example of
            required format.
        layer_type : str, optional
            Type of Landlab layers object used to store the layers. If
            MaterialLayers (default) is specified, then erosion removes material
            and does not create a layer of thickness zero. If EventLayers is
            used, then erosion removes material and creates layers of thickness
            zero. Thus, EventLayers may be appropriate if the user is interested
            in chronostratigraphy.
        dz_advection : float, `(n_nodes, )` shape array, or at-node field array optional
            Change in rock elevation due to advection by some external process.
            This can be changed using the property setter. Dimensions are in
            length, not length per time.
        rock_id : value or `(n_nodes, )` shape array, optional
            Rock type id for new material if deposited.
            This can be changed using the property setter.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")

        Create a Lithology with uniform thicknesses that alternates between
        layers of type 1 and type 2 rock.

        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)

        After creating a Lithology, the model grid will have an at-node grid
        field set to the surface values of 'K_sp'.

        >>> mg.at_node["K_sp"]
        array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])

        The surface values are also properties of the Lithology.

        >>> lith["K_sp"]
        array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])

        We can access information about the Lithology like the total thickness
        or layer thicknesses.

        >>> lith.thickness
        array([8., 8., 8., 8., 8., 8., 8., 8., 8.])
        >>> lith.dz
        array([[1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [4., 4., 4., 4., 4., 4., 4., 4., 4.],
               [2., 2., 2., 2., 2., 2., 2., 2., 2.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.]])

        This might look confusing -- that the layers are in reverse order, but
        it is OK. The last layers in the Lithology are those that are closest
        to the surface.

        The layers don't all have to have the same thickness as in the prior
        example. If the layers have non-uniform thickness, then they must be
        specified in an array of shape `(n_layer, n_nodes)`. In this case, the
        layer IDs must be specified in either an array of `(n_layer)` or
        `(n_layer, n_nodes)`.

        Here we make a layer that gets thicker as a function of the x value of
        the model grid.

        >>> layer_pattern = (0.5 * mg.x_of_node) + 1.0
        >>> thicknesses = [1 * layer_pattern, 2 * layer_pattern, 4 * layer_pattern]
        >>> ids = [1, 2, 1]
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.thickness
        array([ 7. , 10.5, 14. ,  7. , 10.5, 14. ,  7. , 10.5, 14. ])
        >>> lith.dz
        array([[4. , 6. , 8. , 4. , 6. , 8. , 4. , 6. , 8. ],
               [2. , 3. , 4. , 2. , 3. , 4. , 2. , 3. , 4. ],
               [1. , 1.5, 2. , 1. , 1.5, 2. , 1. , 1.5, 2. ]])
        """
        super().__init__(grid)

        try:
            self._last_elevation = self._grid["node"]["topographic__elevation"][
                :
            ].copy()
        except KeyError as exc:
            raise ValueError(
                "Lithology requires that topographic__elevation already "
                "exists as an at-node field."
            ) from exc

        # save inital information about thicknesses, layers, attributes, and ids.
        self._init_thicknesses = np.asarray(thicknesses)
        self._attrs = attrs
        self._number_of_init_layers = self._init_thicknesses.shape[0]
        self._properties = list(attrs.keys())
        self._rock_id_name = "rock_type__id"
        # assert that thicknesses and ids are correct and consistent shapes

        # if thickness is a 2d array.
        if self._init_thicknesses.ndim == 2:
            # assert that the 2nd dimension is the same as the number of nodes.
            if self._init_thicknesses.shape[1] != self._grid.number_of_nodes:
                raise ValueError(
                    "Thicknesses provided to Lithology are ",
                    "inconsistent with the ModelGrid.",
                )

            # if IDs is a 2d array assert that it is the same size as thicknesses
            if np.asarray(ids).ndim == 2:
                if self._init_thicknesses.shape != np.asarray(ids).shape:
                    raise ValueError(
                        "Thicknesses and IDs provided to Lithology are ",
                        "inconsistent with each other.",
                    )
                # if tests pass set value of IDs.
                self._layer_ids = np.asarray(ids)

            # if IDS is a 1d array
            elif np.asarray(ids).ndim == 1:
                if np.asarray(ids).size != self._number_of_init_layers:
                    raise ValueError(
                        "Number of IDs provided to Lithology is ",
                        "inconsistent with number of layers provided in "
                        "thicknesses.",
                    )
                # if tests pass, broadcast ids to correct shape.
                self._layer_ids = np.broadcast_to(
                    np.atleast_2d(np.asarray(ids)).T, self._init_thicknesses.shape
                )

            else:
                raise ValueError(
                    "IDs must be of shape `(n_layers, )` or `(n_layers, n_nodes)`. "
                    "Passed array has more than 2 dimensions."
                )

        elif self._init_thicknesses.ndim == 1:
            if self._init_thicknesses.shape != np.asarray(ids).shape:
                raise ValueError(
                    "Thicknesses and IDs provided to Lithology are ",
                    "inconsistent with each other.",
                )
            self._layer_ids = np.asarray(ids)
        else:
            raise ValueError(
                "Thicknesses must be of shape `(n_layers, )` or "
                "`(n_layers, n_nodes)`. Passed array has more than 2 dimensions."
            )

        # assert that attrs are pointing to fields (or create them)
        for at in self._properties:
            if at not in grid.at_node:
                self._grid.add_empty(at, at="node")

        # add a field for the rock type id
        if self._rock_id_name not in self._grid.at_node:
            self._grid.add_empty(self._rock_id_name, at="node")

        # verify that all IDs have attributes.
        self._check_property_dictionary()

        # create a EventLayers instance
        if layer_type == "EventLayers":
            self._layers = EventLayers(
                grid.number_of_nodes, self._number_of_init_layers
            )
        elif layer_type == "MaterialLayers":
            self._layers = MaterialLayers(
                grid.number_of_nodes, self._number_of_init_layers
            )
        else:
            raise ValueError("Lithology passed an invalid option for " "layer type.")

        # From bottom to top, add layers to the Lithology with attributes.
        for i in range(self._number_of_init_layers - 1, -1, -1):
            try:
                self.add_layer(self._init_thicknesses[i, :], self._layer_ids[i, :])
            except IndexError:
                self.add_layer(self._init_thicknesses[i], self._layer_ids[i])

        self.dz_advection = dz_advection
        self.rock_id = rock_id

    def __getitem__(self, name):
        return self._get_surface_values(name)

    @property
    def dz_advection(self):
        """Rate of vertical advection.

        Parameters
        ----------
        dz_advection : float, `(n_nodes, )` shape array, or at-node field array optional
            Change in rock elevation due to advection by some external process.
            This can be changed using the property setter. Dimensions are in
            length, not length per time.

        Returns
        -------
        current rate of vertical advection
        """
        return return_array_at_node(self._grid, self._dz_advection)

    @dz_advection.setter
    def dz_advection(self, dz_advection):
        return_array_at_node(self._grid, dz_advection)  # verify that this will work.
        self._dz_advection = dz_advection

    @property
    def rock_id(self):
        """Rock type for deposition.

        Parameters
        ----------
        rock_id : value or `(n_nodes, )` shape array, optional
            Rock type id for new material if deposited.
            This can be changed using the property setter.

        Returns
        -------
        current type of rock being deposited (if deposition occurs)
        """
        if self._rock_id is None:
            return None
        else:
            return return_array_at_node(self._grid, self._rock_id)

    @rock_id.setter
    def rock_id(self, rock_id):
        return_array_at_node(self._grid, rock_id)  # verify that this will work.
        # verify that all rock types are valid
        self._rock_id = rock_id

    @property
    def ids(self):
        """Rock type IDs used by Lithology."""
        return list(self._ids)

    @property
    def tracked_properties(self):
        """Properties tracked by Lithology.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.tracked_properties
        ['K_sp']
        """
        self._properties.sort()
        return self._properties

    @property
    def properties(self):
        """Properties dictionary used by Lithology.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.properties
        {'K_sp': {1: 0.001, 2: 0.0001}}
        """
        return self._attrs

    @property
    def thickness(self):
        """Total thickness of the Lithology at each node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.thickness
        array([8., 8., 8., 8., 8., 8., 8., 8., 8.])
        """
        return self._layers.thickness

    @property
    def dz(self):
        """Thickness of each layer in the Lithology at each node.

        The thickness of each layer in the Lithology as an array of shape
        `(number_of_layers, number_of_nodes)`.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.dz
        array([[1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [4., 4., 4., 4., 4., 4., 4., 4., 4.],
               [2., 2., 2., 2., 2., 2., 2., 2., 2.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.]])
        """
        return self._layers.dz

    @property
    def z_bottom(self):
        """Thickness from the surface to the bottom of each layer in Lithology.

        Thickness from the topographic surface to the bottom of each layer as
        an array of shape `(number_of_layers, number_of_nodes)`.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.z_bottom
        array([[8., 8., 8., 8., 8., 8., 8., 8., 8.],
               [7., 7., 7., 7., 7., 7., 7., 7., 7.],
               [3., 3., 3., 3., 3., 3., 3., 3., 3.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.]])
        """
        thick = np.broadcast_to(self._layers.thickness, self._layers.z.shape)
        return thick - self._layers.z + self._layers.dz

    @property
    def z_top(self):
        """Thickness from the surface to the top of each layer in Lithology.

        Thickness from the topographic surface to the top of each layer as
        an array of shape `(number_of_layers, number_of_nodes)`.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.z_top
        array([[7., 7., 7., 7., 7., 7., 7., 7., 7.],
               [3., 3., 3., 3., 3., 3., 3., 3., 3.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.]])
        """
        thick = np.broadcast_to(self._layers.thickness, self._layers.z.shape)
        return thick - self._layers.z

    def _check_property_dictionary(self):
        """Check compatibility of Lithology and property dictionary."""
        ids = []
        for at in self._properties:
            ids.extend(self._attrs[at].keys())
        self._ids = frozenset(np.unique(ids))

        for at in self._properties:
            for i in self._ids:
                if i not in self._attrs[at]:
                    raise ValueError(
                        f"A rock type with ID value {i} was specified in Lithology. "
                        f"No value for this ID was provided in property {at}."
                    )

    def _update_surface_values(self):
        """Update Lithology surface values."""
        # Update surface values for each attribute.
        self._grid["node"][self._rock_id_name][:] = self._surface_rock_type
        for at in self._properties:
            self._grid["node"][at][:] = self[at]

    def add_layer(self, thickness, rock_id=None):
        """Add a new layer to Lithology.

        Parameters
        ----------
        thickness : float or `(n_nodes,)` array
            Positive values deposit material on to Lithology while negative
            values erode Lithology.
        rock_id : single value or `n_nodes` long itterable, optional if only erosion occurs
            Rock type ID for new deposits. Can be single value or an number-
            of-nodes array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]

        We can instantiate Lithology with rock type properties we know we will
        use in the future.

        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001, 3: 0.01}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)

        Add a layer of thickness 3 and rock type 3.

        >>> lith.add_layer(3, rock_id=3)

        The value of `K_sp` at node is now updated to the value of rock type 3

        >>> mg.at_node["K_sp"]
        array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

        A negative value will erode. We can also pass a `(n_nodes,) long array
        to erode unevenly. If all parts of the layer erode, then no `rock_id`
        needs to be passed.

        >>> erosion_amount = [-2.0, -2.0, -2.0, -4.0, -4.0, -4.0, -6.0, -6.0, -6.0]
        >>> lith.add_layer(erosion_amount)
        >>> mg.at_node["K_sp"]
        array([0.01  , 0.01  , 0.01  , 0.0001, 0.0001, 0.0001, 0.001 ,
               0.001 , 0.001 ])

        Now different layers are exposed at the surface and the value of `K_sp`
        is spatially variable.
        """
        thickness = np.array(thickness)

        # verify that Lithology will still have thickness after change
        if np.any((self._layers.thickness + thickness) <= 0):
            raise ValueError(
                "add_layer will result in Lithology having a thickness of "
                "zero at at least one node."
            )

        # verify that rock type added exists.
        try:
            all_ids_present = self._ids.issuperset(rock_id)
            new_ids = rock_id
        except TypeError:
            all_ids_present = self._ids.issuperset([rock_id])
            new_ids = [rock_id]

        if not all_ids_present:
            missing_ids = set(new_ids).difference(self._ids)

            if np.any(thickness > 0):
                raise ValueError(
                    "Lithology add_layer was given a rock type id that does "
                    "not yet exist and will need to deposit. Use a valid "
                    f"rock type or add_rock_type. {missing_ids}"
                )

        # add_rock_type
        if rock_id is not None:
            # add layer
            attributes = {self._rock_id_name: rock_id}
            self._layers.add(thickness, **attributes)
        else:
            self._layers.add(thickness)

        # update surface rock type
        self._surface_rock_type = self._layers.get_surface_values(self._rock_id_name)

        # update surface values
        self._update_surface_values()

    def add_property(self, attrs):
        """Add new property to Lithology.

        Parameters
        ----------
        attrs : dict
            Rock attribute dictionary for the new property(s).

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.add_property({"D": {1: 0.03, 2: 0.004}})
        >>> lith.tracked_properties
        ['D', 'K_sp']
        >>> mg.at_node["D"]
        array([0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03])
        """
        for at in attrs:
            if at in self._properties:
                raise ValueError(
                    "add_property is trying to add an existing "
                    f"attribute, this is not permitted. {at}"
                )

            new_rids = attrs[at].keys()
            for rid in new_rids:
                if rid not in self._ids:
                    raise ValueError(
                        f"add_property has an attribute({at}) for rock type {rid!s} "
                        "that no other rock type has. This is not permitted."
                    )

            for rid in self._ids:
                if rid not in new_rids:
                    raise ValueError(
                        "add_property needs a value for id {rid!s} and attribute {at}."
                    )

        for at in attrs:
            if at not in self._grid.at_node:
                self._grid.add_empty(at, at="node")
            self._attrs[at] = attrs[at]
            self._properties.append(at)

        # update surface values
        self._update_surface_values()

    def add_rock_type(self, attrs):
        """Add rock type to Lithology.

        Parameters
        ----------
        attrs : dict
            Rock attribute dictionary for the new rock type(s).

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.add_rock_type({"K_sp": {4: 0.03, 6: 0.004}})
        >>> lith.ids
        [1, 2, 4, 6]
        >>> lith.properties
        {'K_sp': {1: 0.001, 2: 0.0001, 4: 0.03, 6: 0.004}}
        """
        # Check that the new rock type has all existing attributes
        for at in self._properties:
            if at not in attrs:
                raise ValueError(f"The new rock type is missing attribute {at!s}.")
        # And no new attributes
        for at in attrs:
            if at not in self._properties:
                raise ValueError(
                    "The new rock type has an attribute (e{at!s}) "
                    "that no other rock type has. This is not permitted."
                )

        new_ids = []
        for at in attrs:
            att_dict = attrs[at]
            rids = att_dict.keys()
            for rid in rids:
                if rid in self._layer_ids:
                    raise ValueError(
                        "Rock type ID {rid!s} for attribute {at!s} "
                        "has already been added. This is not allowed"
                    )
                else:
                    new_ids.append(rid)
                    self._attrs[at][rid] = att_dict[rid]
        self._ids = self._ids.union(new_ids)

        # update surface values
        self._update_surface_values()

    def update_rock_properties(self, at, rock_id, value):
        """Update rock type attribute.

        Parameters
        ----------
        at : str
            Attribute name
        rock_id : value
            Rock type ID
        value : value
            New value for rock type attribute

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)

        >>> mg.at_node["K_sp"]
        array([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])

        >>> lith.update_rock_properties("K_sp", 1, 0.03)

        >>> mg.at_node["K_sp"]
        array([0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03])
        """
        if at not in self._properties:
            raise ValueError(
                f"Lithology cannot update the value of {at!s} as "
                "this attribute does not exist."
            )

        if not self._ids.issuperset([rock_id]):
            raise ValueError(
                f"Lithology cannot update the value of rock type {rock_id!s} "
                f"for attribute {at!s} as this rock type is not yet defined."
            )

        # set the value in the attribute dictionary
        self._attrs[at][rock_id] = value

        # update surface values
        self._update_surface_values()

    def _get_surface_values(self, at):
        """Get surface values for attribute."""
        return np.array(list(map(self._attrs[at].get, self._surface_rock_type)))

    def rock_cube_to_xarray(self, depths):
        """Construct a 3D rock cube of rock type ID as an xarray dataset.

        Create an xarray dataset in (x, y, z) that shows the rock type with
        depth relative to the current topographic surface.

        Here the z dimension is depth relative to the current topographic
        surface, NOT depth relative to an absolute datum.

        Note also that when this method is called, it will construct the current
        values of lithology with depth, NOT the initial values.

        Parameters
        ----------
        depths : array

        Returns
        -------
        ds : xarray dataset
        """
        depths = np.asarray(depths)
        rock_type = self._layers[self._rock_id_name]
        rock_cube = np.empty((depths.size, self._grid.shape[0], self._grid.shape[1]))

        # at each node point, interpolate between ztop/bottomo correct.y
        for sid in range(self._layers.number_of_stacks):
            coord = np.unravel_index(sid, (self._grid.shape[0], self._grid.shape[1]))
            real_layers = self.dz[:, sid] > 0
            f = interp1d(
                np.flipud(self.z_top[real_layers, sid]),
                np.flipud(rock_type[real_layers, sid]),
                kind="previous",
            )
            vals = f(depths)
            rock_cube[:, coord[0], coord[1]] = vals

        ds = xr.Dataset(
            data_vars={
                "rock_type__id": (
                    ("z", "y", "x"),
                    rock_cube,
                    {"units": "-", "long_name": "Rock Type ID Code"},
                )
            },
            coords={
                "x": (
                    ("x"),
                    self._grid.x_of_node.reshape(self._grid.shape)[0, :],
                    {"units": "meters"},
                ),
                "y": (
                    ("y"),
                    self._grid.y_of_node.reshape(self._grid.shape)[:, 1],
                    {"units": "meters"},
                ),
                "z": (
                    ("z"),
                    depths,
                    {"units": "meters", "long_name": "Depth Below Topographic Surface"},
                ),
            },
        )

        return ds

    def run_one_step(self):
        """Update Lithology.

        The ``run_one_step`` method calculates elevation change of the
        Lithology surface (taking into account any advection due to external
        processes) and then either deposits or erodes based on elevation
        change.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import Lithology
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_ones("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs)
        >>> lith.dz
        array([[1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [4., 4., 4., 4., 4., 4., 4., 4., 4.],
               [2., 2., 2., 2., 2., 2., 2., 2., 2.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.]])
        >>> lith.thickness
        array([8., 8., 8., 8., 8., 8., 8., 8., 8.])

        If we erode the surface, and then update Lithology, the thickness will
        change.

        >>> z -= 0.5
        >>> lith.run_one_step()
        >>> lith.thickness
        array([7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5])

        The default of Lithology is to use MaterialLayers from the Landlab
        layers submodule. This means that when we erode, we will remove a layer
        from the layers datastructure if it has no material anywhere.

        >>> lith.dz
        array([[1. , 1. , 1. , 1. , 1. , 1. , 1. , 1. , 1. ],
               [4. , 4. , 4. , 4. , 4. , 4. , 4. , 4. , 4. ],
               [2. , 2. , 2. , 2. , 2. , 2. , 2. , 2. , 2. ],
               [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]])

        We can see the value of the rock type at the surface.

        >>> mg.at_node["rock_type__id"]
        array([1., 1., 1., 1., 1., 1., 1., 1., 1.])

        If you deposit, a valid rock_id must be provided. If the rock type
        is the same as the current surface value everywhere, then the layers
        will be combined. This rock_id can be provided as part of the init of
        Lithology or by setting a property (as shown below).

        >>> z += 1.5
        >>> lith.rock_id = 1
        >>> lith.run_one_step()
        >>> lith.thickness
        array([9., 9., 9., 9., 9., 9., 9., 9., 9.])

        >>> lith.dz
        array([[1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [4., 4., 4., 4., 4., 4., 4., 4., 4.],
               [2., 2., 2., 2., 2., 2., 2., 2., 2.],
               [2., 2., 2., 2., 2., 2., 2., 2., 2.]])

        This contrasts with the behavior of Lithology if we use EventLayers.
        Next we repeat this example with EventLayers. Note that no matter which
        method you use, the values of the model grid fields will be the same.
        These two methods differ only in the details of the data structure they
        use to store the layer information.

        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_ones("topographic__elevation", at="node")
        >>> thicknesses = [1, 2, 4, 1]
        >>> ids = [1, 2, 1, 2]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = Lithology(mg, thicknesses, ids, attrs, layer_type="EventLayers")
        >>> lith.dz
        array([[1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [4., 4., 4., 4., 4., 4., 4., 4., 4.],
               [2., 2., 2., 2., 2., 2., 2., 2., 2.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.]])
        >>> lith.thickness
        array([8., 8., 8., 8., 8., 8., 8., 8., 8.])

        If we erode the surface, and then update Lithology, the thickness
        will change. However, with EventLayers, the ``lith.dz`` structure
        will be different. It will have a layer with thickness zero that
        represents the event of erosion.

        >>> z -= 0.5
        >>> lith.run_one_step()
        >>> lith.thickness
        array([7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 7.5])
        >>> lith.dz
        array([[1. , 1. , 1. , 1. , 1. , 1. , 1. , 1. , 1. ],
               [4. , 4. , 4. , 4. , 4. , 4. , 4. , 4. , 4. ],
               [2. , 2. , 2. , 2. , 2. , 2. , 2. , 2. , 2. ],
               [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
               [0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ]])

        We can see the value of the rock type at the surface. As expected,
        it is just the same as if we used MaterialLayers.

        >>> mg.at_node["rock_type__id"]
        array([1., 1., 1., 1., 1., 1., 1., 1., 1.])

        If you deposit, a valid rock_id must be provided. Unlike
        MaterialLayers, these two layers will not be combined, even if they
        have the same properties.

        >>> z += 1.5
        >>> lith.rock_id = 1
        >>> lith.run_one_step()
        >>> lith.thickness
        array([9., 9., 9., 9., 9., 9., 9., 9., 9.])

        >>> lith.dz
        array([[1. , 1. , 1. , 1. , 1. , 1. , 1. , 1. , 1. ],
               [4. , 4. , 4. , 4. , 4. , 4. , 4. , 4. , 4. ],
               [2. , 2. , 2. , 2. , 2. , 2. , 2. , 2. , 2. ],
               [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
               [0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. ],
               [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]])
        """
        # calculate amount of erosion
        elevation_change = self._grid["node"]["topographic__elevation"] - (
            self._last_elevation + self.dz_advection
        )

        # add layer
        self.add_layer(elevation_change, rock_id=self.rock_id)

        # update the last elevation.
        self._last_elevation = self._grid["node"]["topographic__elevation"][:].copy()
