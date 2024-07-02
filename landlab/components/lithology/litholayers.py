#!/usr/bin/env python3
"""Create a LithoLayers component with different properties."""

import numpy as np

from landlab.components.lithology.lithology import Lithology


class LithoLayers(Lithology):
    """Create LithoLayers component.

    A LithoLayers is a three dimentional representation of material operated on
    by landlab components. Material can be removed through erosion or added to
    through deposition. Rock types can have multiple attributes (e.g. age,
    erodability or other parameter values, etc).

    If the tracked properties are model grid fields, they will be updated to
    the surface values of the Lithology. If the properties are not grid fields
    then at-node grid fields will be created with their names.

    It is constructed by specifying a series of depths below the surface, an
    anchor point, a series of rock type ids, and the functional form of a
    surface. Depths and IDs are both specified in order of closest
    to the surface to furthest from the surface.

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

    _name = "LithoLayers"

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
        z0s,
        ids,
        attrs,
        x0=0,
        y0=0,
        function=lambda x, y: 0 * x + 0 * y,
        layer_type="EventLayers",
        dz_advection=0,
        rock_id=None,
    ):
        """Create a new instance of a LithoLayers.

        Parameters
        ----------
        grid : Landlab ModelGrid
        z0s : ndarray of shape `(n_layers, )`
            Values of layer depth from surface at horizontal location (x0, y0).
        ids : ndarray of shape `(n_layers, )`
            Values of rock type IDs corresponding to each layer specified in
            **z0s**.
        attrs : dict
            Rock type property dictionary. See class docstring for example of
            required format.
        x0 : float, optional
            x value of anchor point for all layers.
        y0 : float, optional
            y value of anchor point for all layers.
        function : function, optional
            Functional form of layers as a function of two variables, x and y.
            Default value is `lambda x, y: 0*x + 0*y` for flatlying layers.
        layer_type : str, optional
            Type of Landlab layers object used to store the layers. If
            MaterialLayers (default) is specified, then erosion removes material
            and does not create a layer of thickness zero. If EventLayers is
            used, then erosion removes material and creates layers of thickness
            zero. Thus, EventLayers may be appropriate if the user is interested
            in chronostratigraphy.
        dz_advection : float, `(n_nodes, )` shape array, or at-node field array optional
            Change in rock elevation due to advection by some external process.
            This can be changed using the property setter.
        rock_id : value or `(n_nodes, )` shape array, optional
            Rock type id for new material if deposited.
            This can be changed using the property setter.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import LithoLayers
        >>> mg = RasterModelGrid((3, 3))
        >>> z = mg.add_zeros("node", "topographic__elevation")

        Create a LithoLayers with flatlying layers that altrnate between
        layers of type 1 and type 2 rock.

        >>> z0s = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
        >>> ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
        >>> attrs = {"K_sp": {1: 0.001, 2: 0.0001}}
        >>> lith = LithoLayers(mg, z0s, ids, attrs)
        >>> lith.dz
        array([[1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.]])

        Now create a set of layers that dip. Our anchor point will be the
        default value of (x0, y0) = (0, 0)

        >>> lith = LithoLayers(mg, z0s, ids, attrs, function=lambda x, y: x + y)
        >>> lith.dz
        array([[1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [1., 1., 1., 1., 1., 1., 1., 1., 1.],
               [0., 1., 1., 1., 1., 1., 1., 1., 1.],
               [0., 0., 1., 0., 1., 1., 1., 1., 1.],
               [0., 0., 0., 0., 0., 1., 0., 1., 1.],
               [0., 0., 0., 0., 0., 0., 0., 0., 1.],
               [0., 0., 0., 0., 0., 0., 0., 0., 0.]])

        We can get the surface values, and as we'd expect, they alternate as
        the dipping layers are exposed at the surface.

        >>> lith["K_sp"]
        array([0.0001, 0.001 , 0.0001, 0.001 , 0.0001, 0.001 , 0.0001, 0.001 , 0.0001])
        """

        function_args = function.__code__.co_varnames
        if len(function_args) != 2:
            raise ValueError(
                "LithoLayers: function must take exactly two arguments, x and y."
            )

        if np.asarray(z0s).size != np.asarray(ids).size:
            raise ValueError(
                "LithoLayers: Size of layer depths and layer IDs must be the same"
            )

        if np.any(np.diff(z0s) < 0):
            raise ValueError("LithoLayers: Bad layer depth order passed.")

        z_surf = function(grid.x_of_node - x0, grid.y_of_node - y0)

        if hasattr(z_surf, "shape"):
            if z_surf.shape != grid.x_of_node.shape:
                raise ValueError(
                    "LithoLayers: function must return an array of shape (n_nodes,)"
                )
        else:
            raise ValueError(
                "LithoLayers: function must return an array of shape (n_nodes,)"
            )

        layer_thicknesses = []
        layer_ids = []

        num_layers = np.asarray(z0s).size

        last_layer_elev = np.zeros(grid.number_of_nodes)

        # create layers (here listed from the top to the bottom.)
        for i in range(num_layers):
            layer_depth = z_surf + z0s[i]
            layer_depth[layer_depth < 0] = 0

            layer_thickness = layer_depth.copy() - last_layer_elev.copy()

            last_layer_elev = layer_depth.copy()

            layer_thicknesses.append(layer_thickness)
            layer_ids.append(ids[i] * np.ones(z_surf.size))

        super().__init__(
            grid,
            layer_thicknesses,
            layer_ids,
            attrs,
            layer_type=layer_type,
            dz_advection=dz_advection,
            rock_id=rock_id,
        )
