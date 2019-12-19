import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
from six.moves import range

from landlab import Component, MissingKeyError, ModelParameterDictionary
from landlab.utils.decorators import use_file_name_or_kwds

# Things to add: 1. Explicit stability check.
# 2. Implicit handling of scenarios where kappa*dt exceeds critical step -
#    subdivide dt automatically.


class PerronNLDiffuse(Component):

    """Nonlinear diffusion, following Perron (2011).

    This module uses Taylor Perron's implicit (2011) method to solve the
    nonlinear hillslope diffusion equation across a rectangular, regular grid
    for a single timestep. Note it works with the mass flux implicitly, and
    thus does not actually calculate it. Grid must be at least 5x5.

    Boundary condition handling assumes each edge uses the same BC for each of
    its nodes.
    This component cannot yet handle looped boundary conditions, but all others
    should be fine.

    This component has KNOWN STABILITY ISSUES which will be resolved in a
    future release; use at your own risk.

    The primary method of this class is :func:`run_one_step`.

    Examples
    --------
    >>> from landlab.components import PerronNLDiffuse
    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> mg = RasterModelGrid((5, 5))
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> nl = PerronNLDiffuse(mg, nonlinear_diffusivity=1.)
    >>> dt = 100.
    >>> nt = 20
    >>> uplift_rate = 0.001
    >>> for i in range(nt):
    ...     z[mg.core_nodes] += uplift_rate*dt
    ...     nl.run_one_step(dt)
    >>> z_target = np.array(
    ...     [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
    ...       0.        ,  0.00778637,  0.0075553 ,  0.00778637,  0.        ,
    ...       0.        ,  0.0075553 ,  0.0078053 ,  0.0075553 ,  0.        ,
    ...       0.        ,  0.00778637,  0.0075553 ,  0.00778637,  0.        ,
    ...       0.        ,  0.        ,  0.        ,  0.        ,  0.        ])
    >>> np.allclose(z, z_target)
    True
    """

    _name = "PerronNLDiffuse"

    _input_var_names = ("topographic__elevation",)

    _output_var_names = ("topographic__elevation",)

    _var_units = {"topographic__elevation": "m"}

    _var_mapping = {"topographic__elevation": "node"}

    _var_doc = {
        "topographic__elevation": (
            "Land surface topographic elevation; can "
            + "be overwritten in initialization"
        )
    }

    @use_file_name_or_kwds
    def __init__(
        self,
        grid,
        nonlinear_diffusivity=None,
        S_crit=33.0 * np.pi / 180.0,
        rock_density=2700.0,
        sed_density=2700.0,
        **kwds
    ):
        """
        Parameters
        ----------
        grid : RasterModelGrid
            A Landlab raster grid
        nonlinear_diffusivity : float, array or field name
            The nonlinear diffusivity
        S_crit : float (radians)
            The critical hillslope angle
        rock_density : float (kg*m**-3)
            The density of intact rock
        sed_density : float (kg*m**-3)
            The density of the mobile (sediment) layer
        """
        # disable internal_uplift option:
        internal_uplift = None
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code
        self.values_to_diffuse = "topographic__elevation"
        if nonlinear_diffusivity is not None:
            if nonlinear_diffusivity is not str:
                self._kappa = nonlinear_diffusivity
            else:
                self._kappa = self.grid.at_node[nonlinear_diffusivity]
        else:
            try:
                self._kappa = kwds.pop("kappa", None)
            except KeyError:
                raise KeyError(
                    "nonlinear_diffusivity must be provided to "
                    + "the PerronNLDiffuse component"
                )
        if internal_uplift is None:
            self.internal_uplifts = False
            self._uplift = 0.0
        else:
            self.internal_uplifts = True
            self._uplift = float(internal_uplift)
            # self._uplift = self.grid.zeros('node', dtype=float)
            # self._uplift[self.grid.core_nodes] = internal_uplift
        self._rock_density = rock_density
        self._sed_density = sed_density
        self._S_crit = S_crit

        # for component back compatibility (undocumented):
        # ###
        self.timestep_in = kwds.pop("dt", None)
        if "values_to_diffuse" in kwds.keys():
            self.values_to_diffuse = kwds.pop("values_to_diffuse")
            for mytups in (self._input_var_names, self._output_var_names):
                myset = set(mytups)
                myset.remove("topographic__elevation")
                myset.add(self.values_to_diffuse)
                mytups = tuple(myset)
            for mydicts in (self._var_units, self._var_mapping, self._var_doc):
                mydicts[self.values_to_diffuse] = mydicts.pop("topographic__elevation")

        self._delta_x = grid.dx
        self._delta_y = grid.dy
        self._one_over_delta_x = 1.0 / self._delta_x
        self._one_over_delta_y = 1.0 / self._delta_y
        self._one_over_delta_x_sqd = self._one_over_delta_x ** 2.0
        self._one_over_delta_y_sqd = self._one_over_delta_y ** 2.0
        self._b = 1.0 / self._S_crit ** 2.0

        ncols = grid.number_of_node_columns
        self.ncols = ncols
        nrows = grid.number_of_node_rows
        self.nrows = nrows
        nnodes = grid.number_of_nodes
        self.nnodes = nnodes
        ninteriornodes = grid.number_of_interior_nodes
        ncorenodes = ninteriornodes - 2 * (ncols + nrows - 6)
        self.ninteriornodes = ninteriornodes
        self.interior_grid_width = ncols - 2
        self.core_cell_width = ncols - 4

        self._interior_corners = np.array(
            [ncols + 1, 2 * ncols - 2, nnodes - 2 * ncols + 1, nnodes - ncols - 2]
        )
        _left_list = np.array(range(2 * ncols + 1, nnodes - 2 * ncols, ncols))
        # ^these are still real IDs
        _right_list = np.array(range(3 * ncols - 2, nnodes - 2 * ncols, ncols))
        _bottom_list = np.array(range(ncols + 2, 2 * ncols - 2))
        _top_list = np.array(range(nnodes - 2 * ncols + 2, nnodes - ncols - 2))
        self._left_list = _left_list
        self._right_list = _right_list
        self._bottom_list = _bottom_list
        self._top_list = _top_list

        self._core_nodes = self._coreIDtoreal(np.arange(ncorenodes, dtype=int))
        self.corenodesbyintIDs = self._realIDtointerior(self._core_nodes)
        self.ncorenodes = len(self._core_nodes)

        self.corner_interior_IDs = self._realIDtointerior(self._interior_corners)
        # ^i.e., interior corners as interior IDs
        self.bottom_interior_IDs = self._realIDtointerior(np.array(_bottom_list))
        self.top_interior_IDs = self._realIDtointerior(np.array(_top_list))
        self.left_interior_IDs = self._realIDtointerior(np.array(_left_list))
        self.right_interior_IDs = self._realIDtointerior(np.array(_right_list))

        # build an ID map to let us easily map the variables of the core nodes
        # onto the operating matrix:
        # This array is ninteriornodes long, but the IDs it contains are
        # REAL IDs
        operating_matrix_ID_map = np.empty((ninteriornodes, 9))
        self.interior_IDs_as_real = self._interiorIDtoreal(np.arange(ninteriornodes))
        for j in range(ninteriornodes):
            i = self.interior_IDs_as_real[j]
            operating_matrix_ID_map[j, :] = np.array(
                [
                    (i - ncols - 1),
                    (i - ncols),
                    (i - ncols + 1),
                    (i - 1),
                    i,
                    (i + 1),
                    (i + ncols - 1),
                    (i + ncols),
                    (i + ncols + 1),
                ]
            )
        self.operating_matrix_ID_map = operating_matrix_ID_map
        self.operating_matrix_core_int_IDs = self._realIDtointerior(
            operating_matrix_ID_map[self.corenodesbyintIDs, :]
        )
        # ^shape(ncorenodes,9)
        # see below for corner and edge maps

        # Build masks for the edges and corners to be applied to the operating
        # matrix map.
        # Antimasks are the boundary nodes, masks are "normal"
        self.topleft_mask = [1, 2, 4, 5]
        topleft_antimask = [0, 3, 6, 7, 8]
        self.topright_mask = [0, 1, 3, 4]
        topright_antimask = [2, 5, 6, 7, 8]
        self.bottomleft_mask = [4, 5, 7, 8]
        bottomleft_antimask = [0, 1, 2, 3, 6]
        self.bottomright_mask = [3, 4, 6, 7]
        bottomright_antimask = [0, 1, 2, 5, 8]
        self.corners_masks = np.vstack(
            (
                self.bottomleft_mask,
                self.bottomright_mask,
                self.topleft_mask,
                self.topright_mask,
            )
        )
        # ^(each_corner,mask_for_each_corner)
        self.corners_antimasks = np.vstack(
            (
                bottomleft_antimask,
                bottomright_antimask,
                topleft_antimask,
                topright_antimask,
            )
        )
        # ^so shape becomes (4,5)
        self.left_mask = [1, 2, 4, 5, 7, 8]
        self.left_antimask = [0, 3, 6]
        self.top_mask = [0, 1, 2, 3, 4, 5]
        self.top_antimask = [6, 7, 8]
        self.right_mask = [0, 1, 3, 4, 6, 7]
        self.right_antimask = [2, 5, 8]
        self.bottom_mask = [3, 4, 5, 6, 7, 8]
        self.bottom_antimask = [0, 1, 2]
        self.antimask_corner_position = [0, 2, 2, 4]
        # ^this is the position w/i the corner antimasks that the true corner
        # actually occupies

        self.modulator_mask = np.array(
            [-ncols - 1, -ncols, -ncols + 1, -1, 0, 1, ncols - 1, ncols, ncols + 1]
        )

        self.updated_boundary_conditions()

    def updated_boundary_conditions(self):
        """Call if grid BCs are updated after component instantiation.
        """
        grid = self.grid
        nrows = self.nrows
        ncols = self.ncols
        # ^Set up terms for BC handling (still feels very clumsy)
        bottom_edge = grid.nodes_at_bottom_edge[1:-1]
        top_edge = grid.nodes_at_top_edge[1:-1]
        left_edge = grid.nodes_at_left_edge[1:-1]
        right_edge = grid.nodes_at_right_edge[1:-1]
        self.bottom_flag = 1
        self.top_flag = 1
        self.left_flag = 1
        self.right_flag = 1
        # self.corner_flags = [1,1,1,1] #In ID order, so BL,BR,TL,TR
        if np.all(grid.status_at_node[bottom_edge] == 4):
            # ^This should be all of them, or none of them
            self.bottom_flag = 4
        elif np.all(grid.status_at_node[bottom_edge] == 3):
            self.bottom_flag = 3
        elif np.all(grid.status_at_node[bottom_edge] == 2):
            self.bottom_flag = 2
        elif np.all(grid.status_at_node[bottom_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
            # Note this could get fraught if we need to open a cell to let
            # water flow out...
        if np.all(grid.status_at_node[top_edge] == 4):
            self.top_flag = 4
        elif np.all(grid.status_at_node[top_edge] == 3):
            self.top_flag = 3
        elif np.all(grid.status_at_node[top_edge] == 2):
            self.top_flag = 2
        elif np.all(grid.status_at_node[top_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
        if np.all(grid.status_at_node[left_edge] == 4):
            self.left_flag = 4
        elif np.all(grid.status_at_node[left_edge] == 3):
            self.left_flag = 3
        elif np.all(grid.status_at_node[left_edge] == 2):
            self.left_flag = 2
        elif np.all(grid.status_at_node[left_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
        if np.all(grid.status_at_node[right_edge] == 4):
            self.right_flag = 4
        elif np.all(grid.status_at_node[right_edge] == 3):
            self.right_flag = 3
        elif np.all(grid.status_at_node[right_edge] == 2):
            self.right_flag = 2
        elif np.all(grid.status_at_node[right_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )

        self.fixed_grad_BCs_present = (
            self.bottom_flag == 2
            or self.top_flag == 2
            or self.left_flag == 2
            or self.right_flag == 2
        )
        self.looped_BCs_present = (
            self.bottom_flag == 3
            or self.top_flag == 3
            or self.left_flag == 3
            or self.right_flag == 3
        )
        if self.fixed_grad_BCs_present:
            if self.values_to_diffuse != grid.fixed_gradient_of:
                raise ValueError(
                    "Boundary conditions set in the grid don't "
                    "apply to the data the diffuser is trying to "
                    "work with"
                )

        if np.any(grid.status_at_node == 2):
            self.fixed_grad_offset_map = np.empty(nrows * ncols, dtype=float)
            self.fixed_grad_anchor_map = np.empty_like(self.fixed_grad_offset_map)
            self.fixed_grad_offset_map[
                grid.fixed_gradient_node_properties["boundary_node_IDs"]
            ] = grid.fixed_gradient_node_properties["values_to_add"]

        self.corner_flags = grid.status_at_node[[0, ncols - 1, -ncols, -1]]

        op_mat_just_corners = self.operating_matrix_ID_map[self.corner_interior_IDs, :]
        op_mat_cnr0 = op_mat_just_corners[0, self.bottomleft_mask]
        op_mat_cnr1 = op_mat_just_corners[1, self.bottomright_mask]
        op_mat_cnr2 = op_mat_just_corners[2, self.topleft_mask]
        op_mat_cnr3 = op_mat_just_corners[3, self.topright_mask]
        op_mat_just_active_cnrs = np.vstack(
            (op_mat_cnr0, op_mat_cnr1, op_mat_cnr2, op_mat_cnr3)
        )
        self.operating_matrix_corner_int_IDs = self._realIDtointerior(
            op_mat_just_active_cnrs
        )
        # ^(4corners,4nodesactivepercorner)
        self.operating_matrix_bottom_int_IDs = self._realIDtointerior(
            self.operating_matrix_ID_map[self.bottom_interior_IDs, :][
                :, self.bottom_mask
            ]
        )
        # ^(nbottomnodes,6activenodeseach)
        self.operating_matrix_top_int_IDs = self._realIDtointerior(
            self.operating_matrix_ID_map[self.top_interior_IDs, :][:, self.top_mask]
        )
        self.operating_matrix_left_int_IDs = self._realIDtointerior(
            self.operating_matrix_ID_map[self.left_interior_IDs, :][:, self.left_mask]
        )
        self.operating_matrix_right_int_IDs = self._realIDtointerior(
            self.operating_matrix_ID_map[self.right_interior_IDs, :][:, self.right_mask]
        )

    def _initialize(self, grid, input_stream):
        inputs = ModelParameterDictionary(input_stream)
        self.inputs = inputs
        self.grid = grid

        self.internal_uplifts = False

        if self.internal_uplifts:
            try:
                self._uplift = inputs.read_float("uplift")
            except MissingKeyError:
                self._uplift = inputs.read_float("uplift_rate")
        else:
            self._uplift = 0.0
        self._rock_density = inputs.read_float("rock_density")
        self._sed_density = inputs.read_float("sed_density")
        self._kappa = inputs.read_float("kappa")  # ==_a
        self._S_crit = inputs.read_float("S_crit")
        try:
            self.values_to_diffuse = inputs.read_str("values_to_diffuse")
        except MissingKeyError:
            self.values_to_diffuse = "topographic__elevation"
        try:
            self.timestep_in = inputs.read_float("dt")
        except MissingKeyError:
            raise NameError(
                """No fixed timestep supplied, it must be set
                       dynamically somewhere else. Be sure to call
                       input_timestep(timestep_in) as part of your run
                       loop."""
            )

        self._delta_x = grid.dx
        self._delta_y = grid.dy
        self._one_over_delta_x = 1.0 / self._delta_x
        self._one_over_delta_y = 1.0 / self._delta_y
        self._one_over_delta_x_sqd = self._one_over_delta_x ** 2.0
        self._one_over_delta_y_sqd = self._one_over_delta_y ** 2.0
        self._b = 1.0 / self._S_crit ** 2.0

        ncols = grid.number_of_node_columns
        self.ncols = ncols
        nrows = grid.number_of_node_rows
        self.nrows = nrows
        nnodes = grid.number_of_nodes
        self.nnodes = nnodes
        ninteriornodes = grid.number_of_interior_nodes
        ncorenodes = ninteriornodes - 2 * (ncols + nrows - 6)
        self.ninteriornodes = ninteriornodes
        self.interior_grid_width = ncols - 2
        self.core_cell_width = ncols - 4

        self._interior_corners = np.array(
            [ncols + 1, 2 * ncols - 2, nnodes - 2 * ncols + 1, nnodes - ncols - 2]
        )
        _left_list = np.array(range(2 * ncols + 1, nnodes - 2 * ncols, ncols))
        # ^these are still real IDs
        _right_list = np.array(range(3 * ncols - 2, nnodes - 2 * ncols, ncols))
        _bottom_list = np.array(range(ncols + 2, 2 * ncols - 2))
        _top_list = np.array(range(nnodes - 2 * ncols + 2, nnodes - ncols - 2))
        self._left_list = _left_list
        self._right_list = _right_list
        self._bottom_list = _bottom_list
        self._top_list = _top_list

        self._core_nodes = self._coreIDtoreal(np.arange(ncorenodes, dtype=int))
        self.corenodesbyintIDs = self._realIDtointerior(self._core_nodes)
        self.ncorenodes = len(self._core_nodes)

        self.corner_interior_IDs = self._realIDtointerior(self._interior_corners)
        # ^i.e., interior corners as interior IDs
        self.bottom_interior_IDs = self._realIDtointerior(np.array(_bottom_list))
        self.top_interior_IDs = self._realIDtointerior(np.array(_top_list))
        self.left_interior_IDs = self._realIDtointerior(np.array(_left_list))
        self.right_interior_IDs = self._realIDtointerior(np.array(_right_list))

        # build an ID map to let us easily map the variables of the core nodes
        # onto the operating matrix:
        # This array is ninteriornodes long, but the IDs it contains are
        # REAL IDs
        operating_matrix_ID_map = np.empty((ninteriornodes, 9))
        self.interior_IDs_as_real = self._interiorIDtoreal(np.arange(ninteriornodes))
        for j in range(ninteriornodes):
            i = self.interior_IDs_as_real[j]
            operating_matrix_ID_map[j, :] = np.array(
                [
                    (i - ncols - 1),
                    (i - ncols),
                    (i - ncols + 1),
                    (i - 1),
                    i,
                    (i + 1),
                    (i + ncols - 1),
                    (i + ncols),
                    (i + ncols + 1),
                ]
            )
        self.operating_matrix_ID_map = operating_matrix_ID_map
        self.operating_matrix_core_int_IDs = self._realIDtointerior(
            operating_matrix_ID_map[self.corenodesbyintIDs, :]
        )
        # ^shape(ncorenodes,9)
        # see below for corner and edge maps

        # Build masks for the edges and corners to be applied to the operating
        # matrix map.
        # Antimasks are the boundary nodes, masks are "normal"
        topleft_mask = [1, 2, 4, 5]
        topleft_antimask = [0, 3, 6, 7, 8]
        topright_mask = [0, 1, 3, 4]
        topright_antimask = [2, 5, 6, 7, 8]
        bottomleft_mask = [4, 5, 7, 8]
        bottomleft_antimask = [0, 1, 2, 3, 6]
        bottomright_mask = [3, 4, 6, 7]
        bottomright_antimask = [0, 1, 2, 5, 8]
        self.corners_masks = np.vstack(
            (bottomleft_mask, bottomright_mask, topleft_mask, topright_mask)
        )
        # ^(each_corner,mask_for_each_corner)
        self.corners_antimasks = np.vstack(
            (
                bottomleft_antimask,
                bottomright_antimask,
                topleft_antimask,
                topright_antimask,
            )
        )
        # ^so shape becomes (4,5)
        self.left_mask = [1, 2, 4, 5, 7, 8]
        self.left_antimask = [0, 3, 6]
        self.top_mask = [0, 1, 2, 3, 4, 5]
        self.top_antimask = [6, 7, 8]
        self.right_mask = [0, 1, 3, 4, 6, 7]
        self.right_antimask = [2, 5, 8]
        self.bottom_mask = [3, 4, 5, 6, 7, 8]
        self.bottom_antimask = [0, 1, 2]
        self.antimask_corner_position = [0, 2, 2, 4]
        # ^this is the position w/i the corner antimasks that the true corner
        # actually occupies

        self.modulator_mask = np.array(
            [-ncols - 1, -ncols, -ncols + 1, -1, 0, 1, ncols - 1, ncols, ncols + 1]
        )

        # ^Set up terms for BC handling (still feels very clumsy)
        bottom_edge = grid.nodes_at_bottom_edge[1:-1]
        top_edge = grid.nodes_at_top_edge[1:-1]
        left_edge = grid.nodes_at_left_edge[1:-1]
        right_edge = grid.nodes_at_right_edge[1:-1]
        self.bottom_flag = 1
        self.top_flag = 1
        self.left_flag = 1
        self.right_flag = 1
        # self.corner_flags = [1,1,1,1] #In ID order, so BL,BR,TL,TR
        if np.all(grid.status_at_node[bottom_edge] == 4):
            # ^This should be all of them, or none of them
            self.bottom_flag = 4
        elif np.all(grid.status_at_node[bottom_edge] == 3):
            self.bottom_flag = 3
        elif np.all(grid.status_at_node[bottom_edge] == 2):
            self.bottom_flag = 2
        elif np.all(grid.status_at_node[bottom_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
            # Note this could get fraught if we need to open a cell to let
            # water flow out...
        if np.all(grid.status_at_node[top_edge] == 4):
            self.top_flag = 4
        elif np.all(grid.status_at_node[top_edge] == 3):
            self.top_flag = 3
        elif np.all(grid.status_at_node[top_edge] == 2):
            self.top_flag = 2
        elif np.all(grid.status_at_node[top_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
        if np.all(grid.status_at_node[left_edge] == 4):
            self.left_flag = 4
        elif np.all(grid.status_at_node[left_edge] == 3):
            self.left_flag = 3
        elif np.all(grid.status_at_node[left_edge] == 2):
            self.left_flag = 2
        elif np.all(grid.status_at_node[left_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
        if np.all(grid.status_at_node[right_edge] == 4):
            self.right_flag = 4
        elif np.all(grid.status_at_node[right_edge] == 3):
            self.right_flag = 3
        elif np.all(grid.status_at_node[right_edge] == 2):
            self.right_flag = 2
        elif np.all(grid.status_at_node[right_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )

        self.fixed_grad_BCs_present = (
            self.bottom_flag == 2
            or self.top_flag == 2
            or self.left_flag == 2
            or self.right_flag == 2
        )
        self.looped_BCs_present = (
            self.bottom_flag == 3
            or self.top_flag == 3
            or self.left_flag == 3
            or self.right_flag == 3
        )
        if self.fixed_grad_BCs_present:
            if self.values_to_diffuse != grid.fixed_gradient_of:
                raise ValueError(
                    "Boundary conditions set in the grid don't "
                    "apply to the data the diffuser is trying to "
                    "work with"
                )

        if np.any(grid.status_at_node == 2):
            self.fixed_grad_offset_map = np.empty(nrows * ncols, dtype=float)
            self.fixed_grad_anchor_map = np.empty_like(self.fixed_grad_offset_map)
            self.fixed_grad_offset_map[
                grid.fixed_gradient_node_properties["boundary_node_IDs"]
            ] = grid.fixed_gradient_node_properties["values_to_add"]

        self.corner_flags = grid.status_at_node[[0, ncols - 1, -ncols, -1]]

        op_mat_just_corners = operating_matrix_ID_map[self.corner_interior_IDs, :]
        op_mat_cnr0 = op_mat_just_corners[0, bottomleft_mask]
        op_mat_cnr1 = op_mat_just_corners[1, bottomright_mask]
        op_mat_cnr2 = op_mat_just_corners[2, topleft_mask]
        op_mat_cnr3 = op_mat_just_corners[3, topright_mask]
        op_mat_just_active_cnrs = np.vstack(
            (op_mat_cnr0, op_mat_cnr1, op_mat_cnr2, op_mat_cnr3)
        )
        self.operating_matrix_corner_int_IDs = self._realIDtointerior(
            op_mat_just_active_cnrs
        )
        # ^(4corners,4nodesactivepercorner)
        self.operating_matrix_bottom_int_IDs = self._realIDtointerior(
            operating_matrix_ID_map[self.bottom_interior_IDs, :][:, self.bottom_mask]
        )
        # ^(nbottomnodes,6activenodeseach)
        self.operating_matrix_top_int_IDs = self._realIDtointerior(
            operating_matrix_ID_map[self.top_interior_IDs, :][:, self.top_mask]
        )
        self.operating_matrix_left_int_IDs = self._realIDtointerior(
            operating_matrix_ID_map[self.left_interior_IDs, :][:, self.left_mask]
        )
        self.operating_matrix_right_int_IDs = self._realIDtointerior(
            operating_matrix_ID_map[self.right_interior_IDs, :][:, self.right_mask]
        )

    def input_timestep(self, timestep_in):
        """
        Allows the user to set a dynamic (evolving) timestep manually as part
        of a run loop.
        """
        self.timestep_in = timestep_in

    def _gear_timestep(self, timestep_in, new_grid):
        """
        This method allows the gearing between the model run step and the
        component (shorter) step.
        The method becomes unstable if S>Scrit, so we test to prevent this.
        We implicitly assume the initial condition does not contain
        slopes > Scrit. If the method persistently explodes, this may be the
        problem.
        """
        extended_elevs = np.empty(self.grid.number_of_nodes + 1, dtype=float)
        extended_elevs[-1] = np.nan
        node_neighbors = self.grid.active_adjacent_nodes_at_node
        extended_elevs[:-1] = new_grid["node"][self.values_to_diffuse]
        max_offset = np.nanmax(
            np.fabs(
                extended_elevs[:-1][node_neighbors]
                - extended_elevs[:-1].reshape((self.grid.number_of_nodes, 1))
            )
        )
        if max_offset > np.tan(self._S_crit) * min(self.grid.dx, self.grid.dy):
            # ^using S not tan(S) adds a buffer - but not appropriate
            self.internal_repeats = (
                int(
                    max_offset
                    // (np.tan(self._S_crit) * min(self.grid.dx, self.grid.dy))
                )
                + 1
            )
            # now we rig it so the actual timestep is an integer divisor
            # of T_in:
            self._delta_t = timestep_in / self.internal_repeats
            self.uplift_per_step = (
                new_grid["node"][self.values_to_diffuse]
                - self.grid["node"][self.values_to_diffuse]
            ) / self.internal_repeats
            if self.internal_repeats > 10000:
                raise ValueError(
                    """Uplift rate is too high; solution is not
                                 stable!!"""
                )
        else:
            self.internal_repeats = 1
            self._delta_t = timestep_in
            self.uplift_per_step = (
                new_grid["node"][self.values_to_diffuse]
                - self.grid["node"][self.values_to_diffuse]
            )
        return self._delta_t

    def _set_variables(self, grid):
        """
        This function sets the variables needed for update().
        Now vectorized, shouold run faster.
        At the moment, this method can only handle fixed value BCs.
        """
        n_interior_nodes = grid.number_of_interior_nodes

        # Initialize the local builder lists
        _mat_RHS = np.zeros(n_interior_nodes)

        try:
            elev = grid["node"][self.values_to_diffuse]
        except KeyError:
            raise NameError("elevations not found in grid!")
        try:
            _delta_t = self._delta_t
        except AttributeError:
            raise NameError(
                """Timestep not set! Call _gear_timestep(tstep)
                            after initializing the component, but before
                            running it."""
            )
        _one_over_delta_x = self._one_over_delta_x
        _one_over_delta_x_sqd = self._one_over_delta_x_sqd
        _one_over_delta_y = self._one_over_delta_y
        _one_over_delta_y_sqd = self._one_over_delta_y_sqd
        _kappa = self._kappa
        _b = self._b
        _S_crit = self._S_crit
        _core_nodes = self._core_nodes
        corenodesbyintIDs = self.corenodesbyintIDs
        operating_matrix_core_int_IDs = self.operating_matrix_core_int_IDs
        operating_matrix_corner_int_IDs = self.operating_matrix_corner_int_IDs
        _interior_corners = self._interior_corners
        corners_antimasks = self.corners_antimasks
        corner_interior_IDs = self.corner_interior_IDs
        modulator_mask = self.modulator_mask
        corner_flags = self.corner_flags
        bottom_interior_IDs = self.bottom_interior_IDs
        top_interior_IDs = self.top_interior_IDs
        left_interior_IDs = self.left_interior_IDs
        right_interior_IDs = self.right_interior_IDs
        bottom_antimask = self.bottom_antimask
        _bottom_list = self._bottom_list
        top_antimask = self.top_antimask
        _top_list = self._top_list
        left_antimask = self.left_antimask
        _left_list = self._left_list
        right_antimask = self.right_antimask
        _right_list = self._right_list

        # Need to modify the "effective" values of the edge nodes if any of
        # the edges are inactive:
        if self.bottom_flag == 4:
            bottom_edge, inside_bottom_edge = grid.nodes[(0, 1), :]
            elev[bottom_edge] = elev[inside_bottom_edge]
            # corners are special cases, and assumed linked to the bottom and
            # top edge BCs...
            elev[bottom_edge[0]] = elev[inside_bottom_edge[1]]
            elev[bottom_edge[-1]] = elev[inside_bottom_edge[-2]]
        if self.top_flag == 4:
            top_edge, inside_top_edge = grid.nodes[(-1, -2), :]
            elev[top_edge] = elev[inside_top_edge]
            # corners are special cases, and assumed linked to the bottom and
            # top edge BCs...
            elev[top_edge[0]] = elev[inside_top_edge[1]]
            elev[top_edge[-1]] = elev[inside_top_edge[-2]]
        if self.left_flag == 4:
            left_edge = grid.nodes[1:-1, 0]
            inside_left_edge = grid.nodes[1:-1, 1]
            elev[left_edge] = elev[inside_left_edge]
        if self.right_flag == 4:
            right_edge = grid.nodes[1:-1, -1]
            inside_right_edge = grid.nodes[1:-1, -2]
            elev[right_edge] = elev[inside_right_edge]

        # replacing loop:
        cell_neighbors = grid.active_adjacent_nodes_at_node
        # ^E,N,W,S
        cell_diagonals = grid.diagonal_adjacent_nodes_at_node  # NE,NW,SW,SE
        # ^this should be dealt with by active_neighbors... (skips bad nodes)

        _z_x = (
            (elev[cell_neighbors[:, 0]] - elev[cell_neighbors[:, 2]])
            * 0.5
            * _one_over_delta_x
        )
        _z_y = (
            (elev[cell_neighbors[:, 1]] - elev[cell_neighbors[:, 3]])
            * 0.5
            * _one_over_delta_y
        )
        _z_xx = (
            elev[cell_neighbors[:, 0]] - 2.0 * elev + elev[cell_neighbors[:, 2]]
        ) * _one_over_delta_x_sqd
        _z_yy = (
            elev[cell_neighbors[:, 1]] - 2.0 * elev + elev[cell_neighbors[:, 3]]
        ) * _one_over_delta_y_sqd
        _z_xy = (
            (
                elev[cell_diagonals[:, 0]]
                - elev[cell_diagonals[:, 1]]
                - elev[cell_diagonals[:, 3]]
                + elev[cell_diagonals[:, 2]]
            )
            * 0.25
            * _one_over_delta_x
            * _one_over_delta_y
        )
        _d = 1.0 / (1.0 - _b * (_z_x * _z_x + _z_y * _z_y))

        _abd_sqd = _kappa * _b * _d * _d
        _F_ij = -2.0 * _kappa * _d * (
            _one_over_delta_x_sqd + _one_over_delta_y_sqd
        ) - 4.0 * _abd_sqd * (
            _z_x * _z_x * _one_over_delta_x_sqd + _z_y * _z_y * _one_over_delta_y_sqd
        )
        _F_ijminus1 = (
            _kappa * _d * _one_over_delta_x_sqd
            - _abd_sqd * _z_x * (_z_xx + _z_yy) * _one_over_delta_x
            - 4.0
            * _abd_sqd
            * _b
            * _d
            * (_z_x * _z_x * _z_xx + _z_y * _z_y * _z_yy + 2.0 * _z_x * _z_y * _z_xy)
            * _z_x
            * _one_over_delta_x
            - 2.0
            * _abd_sqd
            * (
                _z_x * _z_xx * _one_over_delta_x
                - _z_x * _z_x * _one_over_delta_x_sqd
                + _z_y * _z_xy * _one_over_delta_x
            )
        )
        _F_ijplus1 = (
            _kappa * _d * _one_over_delta_x_sqd
            + _abd_sqd * _z_x * (_z_xx + _z_yy) * _one_over_delta_x
            + 4.0
            * _abd_sqd
            * _b
            * _d
            * (_z_x * _z_x * _z_xx + _z_y * _z_y * _z_yy + 2.0 * _z_x * _z_y * _z_xy)
            * _z_x
            * _one_over_delta_x
            + 2.0
            * _abd_sqd
            * (
                _z_x * _z_xx * _one_over_delta_x
                + _z_x * _z_x * _one_over_delta_x_sqd
                + _z_y * _z_xy * _one_over_delta_x
            )
        )
        _F_iminus1j = (
            _kappa * _d * _one_over_delta_y_sqd
            - _abd_sqd * _z_y * (_z_xx + _z_yy) * _one_over_delta_y
            - 4.0
            * _abd_sqd
            * _b
            * _d
            * (_z_x * _z_x * _z_xx + _z_y * _z_y * _z_yy + 2.0 * _z_x * _z_y * _z_xy)
            * _z_y
            * _one_over_delta_y
            - 2.0
            * _abd_sqd
            * (
                _z_y * _z_yy * _one_over_delta_y
                - _z_y * _z_y * _one_over_delta_y_sqd
                + _z_x * _z_xy * _one_over_delta_y
            )
        )
        _F_iplus1j = (
            _kappa * _d * _one_over_delta_y_sqd
            + _abd_sqd * _z_y * (_z_xx + _z_yy) * _one_over_delta_y
            + 4.0
            * _abd_sqd
            * _b
            * _d
            * (_z_x * _z_x * _z_xx + _z_y * _z_y * _z_yy + 2.0 * _z_x * _z_y * _z_xy)
            * _z_y
            * _one_over_delta_y
            + 2.0
            * _abd_sqd
            * (
                _z_y * _z_yy * _one_over_delta_y
                + _z_y * _z_y * _one_over_delta_y_sqd
                + _z_x * _z_xy * _one_over_delta_y
            )
        )
        _F_iplus1jplus1 = _abd_sqd * _z_x * _z_y * _one_over_delta_x * _one_over_delta_y
        _F_iminus1jminus1 = _F_iplus1jplus1
        _F_iplus1jminus1 = -_F_iplus1jplus1
        _F_iminus1jplus1 = _F_iplus1jminus1

        _equ_RHS_calc_frag = (
            _F_ij * elev
            + _F_ijminus1 * elev[cell_neighbors[:, 2]]
            + _F_ijplus1 * elev[cell_neighbors[:, 0]]
            + _F_iminus1j * elev[cell_neighbors[:, 3]]
            + _F_iplus1j * elev[cell_neighbors[:, 1]]
            + _F_iminus1jminus1 * elev[cell_diagonals[:, 2]]
            + _F_iplus1jplus1 * elev[cell_diagonals[:, 0]]
            + _F_iplus1jminus1 * elev[cell_diagonals[:, 1]]
            + _F_iminus1jplus1 * elev[cell_diagonals[:, 3]]
        )

        # NB- all _z_... and _F_... variables are nnodes long, and thus use
        # real IDs (tho calcs will be flawed for Bnodes)

        # RHS of equ 6 (see para [20])
        _func_on_z = self._rock_density / self._sed_density * self._uplift + _kappa * (
            (_z_xx + _z_yy) / (1.0 - (_z_x * _z_x + _z_y * _z_y) / _S_crit * _S_crit)
            + 2.0
            * (_z_x * _z_x * _z_xx + _z_y * _z_y * _z_yy + 2.0 * _z_x * _z_y * _z_xy)
            / (
                _S_crit
                * _S_crit
                * (1.0 - (_z_x * _z_x + _z_y * _z_y) / _S_crit * _S_crit) ** 2.0
            )
        )

        # Remember, the RHS is getting wiped each loop as part of
        # self._set_variables()
        # _mat_RHS is ninteriornodes long, but were only working on a
        # ncorenodes long subset here
        _mat_RHS[corenodesbyintIDs] += elev[_core_nodes] + _delta_t * (
            _func_on_z[_core_nodes] - _equ_RHS_calc_frag[_core_nodes]
        )
        low_row = (
            np.vstack((_F_iminus1jminus1, _F_iminus1j, _F_iminus1jplus1)) * -_delta_t
        )
        mid_row = np.vstack(
            (-_delta_t * _F_ijminus1, 1.0 - _delta_t * _F_ij, -_delta_t * _F_ijplus1)
        )
        top_row = np.vstack((_F_iplus1jminus1, _F_iplus1j, _F_iplus1jplus1)) * -_delta_t
        nine_node_map = np.vstack((low_row, mid_row, top_row)).T
        # ^Note shape is (nnodes,9); it's realID indexed
        core_op_mat_row = np.repeat(corenodesbyintIDs, 9)
        core_op_mat_col = operating_matrix_core_int_IDs.astype(int).flatten()
        core_op_mat_data = nine_node_map[_core_nodes, :].flatten()

        # Now the interior corners; BL,BR,TL,TR
        _mat_RHS[corner_interior_IDs] += elev[_interior_corners] + _delta_t * (
            _func_on_z[_interior_corners] - _equ_RHS_calc_frag[_interior_corners]
        )
        corners_op_mat_row = np.repeat(self.corner_interior_IDs, 4)
        corners_op_mat_col = operating_matrix_corner_int_IDs.astype(int).flatten()
        corners_op_mat_data = nine_node_map[_interior_corners, :][
            (np.arange(4).reshape((4, 1)), self.corners_masks)
        ].flatten()
        # ^1st index gives (4,9), 2nd reduces to (4,4), then flattened
        for i in range(4):  # loop over each corner, as so few
            # Note that this ONLY ADDS THE VALUES FOR THE TRUE GRID CORNERS.
            # The sides get done in the edge tests, below.
            if corner_flags[i] == 1:
                true_corner = self.antimask_corner_position[i]
                _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                    nine_node_map[_interior_corners[i], :][
                        corners_antimasks[i, true_corner]
                    ]
                    * elev[
                        _interior_corners[i]
                        + modulator_mask[corners_antimasks[i, true_corner]]
                    ]
                )
            elif corner_flags[i] == 4 or corner_flags[i] == 3:
                # ^inactive boundary cell
                # Actually the easiest case! Equivalent to fixed gradient,
                # but the gradient is zero, so material only goes in the linked
                # cell. And because it's a true corner, that linked cell
                # doesn't appear in the interior matrix at all!
                pass
            elif corner_flags[i] == 2:
                true_corner = self.antimask_corner_position[i]
                _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                    nine_node_map[_interior_corners[i], :][
                        corners_antimasks[i, true_corner]
                    ]
                    * self.fixed_gradient_offset_map[
                        _interior_corners[i]
                        + modulator_mask[corners_antimasks[i, true_corner]]
                    ]
                )
            else:
                raise NameError(
                    """Sorry! This module cannot yet handle fixed
                    gradient or looped BCs..."""
                )
            # Todo: handle these BCs properly, once the grid itself works with
            # them.
            # Can follow old routines; see self.set_bc_cell() commented out
            # method below.

        # Now the edges
        _mat_RHS[bottom_interior_IDs] += elev[_bottom_list] + _delta_t * (
            _func_on_z[_bottom_list] - _equ_RHS_calc_frag[_bottom_list]
        )
        _mat_RHS[top_interior_IDs] += elev[_top_list] + _delta_t * (
            _func_on_z[_top_list] - _equ_RHS_calc_frag[_top_list]
        )
        _mat_RHS[left_interior_IDs] += elev[_left_list] + _delta_t * (
            _func_on_z[_left_list] - _equ_RHS_calc_frag[_left_list]
        )
        _mat_RHS[right_interior_IDs] += elev[_right_list] + _delta_t * (
            _func_on_z[_right_list] - _equ_RHS_calc_frag[_right_list]
        )
        bottom_op_mat_row = np.repeat(bottom_interior_IDs, 6)
        top_op_mat_row = np.repeat(top_interior_IDs, 6)
        left_op_mat_row = np.repeat(left_interior_IDs, 6)
        right_op_mat_row = np.repeat(right_interior_IDs, 6)
        bottom_op_mat_col = self.operating_matrix_bottom_int_IDs.astype(int).flatten()
        top_op_mat_col = self.operating_matrix_top_int_IDs.astype(int).flatten()
        left_op_mat_col = self.operating_matrix_left_int_IDs.astype(int).flatten()
        right_op_mat_col = self.operating_matrix_right_int_IDs.astype(int).flatten()
        bottom_op_mat_data = nine_node_map[_bottom_list, :][
            :, self.bottom_mask
        ].flatten()
        top_op_mat_data = nine_node_map[_top_list, :][:, self.top_mask].flatten()
        left_op_mat_data = nine_node_map[_left_list, :][:, self.left_mask].flatten()
        right_op_mat_data = nine_node_map[_right_list, :][:, self.right_mask].flatten()

        if self.bottom_flag == 1:
            # goes to RHS only
            _mat_RHS[bottom_interior_IDs] -= _delta_t * np.sum(
                nine_node_map[_bottom_list, :][:, bottom_antimask]
                * elev[
                    _bottom_list.reshape((len(_bottom_list), 1))
                    + (modulator_mask[bottom_antimask]).reshape((1, 3))
                ],
                axis=1,
            )
            # ^note the broadcasting to (nedge,3) in final fancy index
            # ...& the corners
            edges = [(1, 2), (0, 1), (0, 0), (0, 0)]
            for i in [0, 1]:
                edge_list = edges[i]
                _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                    nine_node_map[_interior_corners[i], :][
                        corners_antimasks[i, edge_list]
                    ]
                    * elev[
                        _interior_corners[i]
                        + modulator_mask[corners_antimasks[i, edge_list]]
                    ]
                )
            # make dummy array objects for the x,y coords in coo creation of
            # _operating_matrix
            bottom_op_mat_row_add = np.empty(0)
            bottom_op_mat_col_add = np.empty(0)
            bottom_op_mat_data_add = np.empty(0)
        elif self.bottom_flag == 4 or self.bottom_flag == 2:
            # ^i.e., fixed zero gradient (4) or more general case...
            bottom_op_mat_row_add = np.empty((bottom_interior_IDs.size * 3 + 6))
            bottom_op_mat_col_add = np.empty((bottom_interior_IDs.size * 3 + 6))
            bottom_op_mat_data_add = np.empty((bottom_interior_IDs.size * 3 + 6))
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            bottom_op_mat_row_add[: (bottom_interior_IDs.size * 3)] = np.repeat(
                bottom_interior_IDs, 3
            )
            bottom_op_mat_col_add[
                : (bottom_interior_IDs.size * 3)
            ] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.bottom_interior_IDs, :][
                    :, self.bottom_mask[0:3]
                ]
            ).flatten()
            bottom_op_mat_data_add[: (bottom_interior_IDs.size * 3)] = (
                _delta_t
                * (nine_node_map[_bottom_list, :][:, bottom_antimask]).flatten()
            )
            # ...& the corners
            this_corner_coords = np.array([0, 1])
            # order is bottom 2 lower left, bottom 2 lower right, lower left
            # true corner, lower right true corner.
            bottom_op_mat_row_add[-6:-2] = np.repeat(
                corner_interior_IDs[this_corner_coords], 2
            )
            bottom_op_mat_col_add[-6:-2] = self.operating_matrix_corner_int_IDs[
                this_corner_coords.reshape((2, 1)), this_corner_coords
            ].flatten()
            bottom_op_mat_row_add[-2:] = corner_interior_IDs[this_corner_coords]
            bottom_op_mat_col_add[-2:] = self.operating_matrix_corner_int_IDs[
                (this_corner_coords[0], this_corner_coords[0]),
                (this_corner_coords[1], this_corner_coords[1]),
            ].flatten()
            bottom_op_mat_data_add[-6:-4] = (
                _delta_t
                * nine_node_map[_interior_corners[0], :][
                    corners_antimasks[0, [1, 2]]
                ].flatten()
            )
            bottom_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[1], :][
                    corners_antimasks[1, [0, 1]]
                ].flatten()
            )
            bottom_op_mat_data_add[-2] = (
                _delta_t
                * nine_node_map[_interior_corners[0], :][corners_antimasks[0, 0]]
            )
            bottom_op_mat_data_add[-1] = (
                _delta_t
                * nine_node_map[_interior_corners[1], :][corners_antimasks[1, 2]]
            )
            if self.bottom_flag == 2:
                # Read the offsets from the map we made in the __init__,
                # use them as constant terms, incorporated into RHS
                _mat_RHS[bottom_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_bottom_list, :][:, bottom_antimask]
                    * self.fixed_gradient_offset_map[
                        _bottom_list.reshape((len(_bottom_list), 1))
                        + (modulator_mask[bottom_antimask]).reshape((1, 3))
                    ],
                    axis=1,
                )
                # ^note the broadcasting to (nedge,3) in final fancy index
                # ...& the corners
                edges = [(1, 2), (0, 1), (0, 0), (0, 0)]
                for i in [0, 1]:
                    edge_list = edges[i]
                    _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                        nine_node_map[_interior_corners[i], :][
                            corners_antimasks[i, edge_list]
                        ]
                        * self.fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self.bottom_flag == 3:
            # This will handle both top and bottom BCs...
            bottom_op_mat_row_add = np.empty((bottom_interior_IDs.size * 3 + 6))
            bottom_op_mat_col_add = np.empty((bottom_interior_IDs.size * 3 + 6))
            bottom_op_mat_data_add = np.empty((bottom_interior_IDs.size * 3 + 6))
            bottom_op_mat_row_add[: (bottom_interior_IDs.size * 3)] = np.repeat(
                bottom_interior_IDs, 3
            )
            # ^...put the values in the same places in the operating matrix...
            bottom_op_mat_col_add[
                : (bottom_interior_IDs.size * 3)
            ] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.top_interior_IDs, :][
                    :, self.top_mask[3:6]
                ]
            ).flatten()
            bottom_op_mat_data_add[: (bottom_interior_IDs.size * 3)] = (
                _delta_t
                * (nine_node_map[_bottom_list, :][:, bottom_antimask]).flatten()
            )
            # ^...but the values refer to the TOP of the grid
            top_op_mat_row_add = np.empty((top_interior_IDs.size * 3 + 6))
            top_op_mat_col_add = np.empty((top_interior_IDs.size * 3 + 6))
            top_op_mat_data_add = np.empty((top_interior_IDs.size * 3 + 6))
            top_op_mat_row_add[: (top_interior_IDs.size * 3)] = np.repeat(
                top_interior_IDs, 3
            )
            top_op_mat_col_add[: (top_interior_IDs.size * 3)] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.bottom_interior_IDs, :][
                    :, self.bottom_mask[0:3]
                ]
            ).flatten()
            top_op_mat_data_add[: (top_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_top_list, :][:, top_antimask]).flatten()
            )
            # & the corners
            bottom_corner_coords = np.array([0, 1])
            top_corner_coords = np.array([2, 3])
            bottom_op_mat_row_add[-6:-2] = np.repeat(
                corner_interior_IDs[bottom_corner_coords], 2
            )
            bottom_op_mat_col_add[-6:-2] = self.operating_matrix_corner_int_IDs[
                top_corner_coords.reshape((2, 1)), top_corner_coords
            ].flatten()
            bottom_op_mat_row_add[-2:] = corner_interior_IDs[bottom_corner_coords]
            bottom_op_mat_col_add[-2:] = self.operating_matrix_corner_int_IDs[
                (top_corner_coords[0], top_corner_coords[0]),
                (top_corner_coords[1], top_corner_coords[1]),
            ].flatten()
            bottom_op_mat_data_add[-6:-4] = (
                _delta_t
                * nine_node_map[_interior_corners[0], :][
                    corners_antimasks[0, [1, 2]]
                ].flatten()
            )
            bottom_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[1], :][
                    corners_antimasks[1, [0, 1]]
                ].flatten()
            )
            bottom_op_mat_data_add[-2] = (
                _delta_t
                * nine_node_map[_interior_corners[0], :][corners_antimasks[0, 0]]
            )
            bottom_op_mat_data_add[-1] = (
                _delta_t
                * nine_node_map[_interior_corners[1], :][corners_antimasks[1, 2]]
            )
            top_op_mat_row_add[-6:-2] = np.repeat(
                corner_interior_IDs[top_corner_coords], 2
            )
            top_op_mat_col_add[-6:-2] = self.operating_matrix_corner_int_IDs[
                bottom_corner_coords.reshape((2, 1)), bottom_corner_coords
            ].flatten()
            top_op_mat_row_add[-2:] = corner_interior_IDs[top_corner_coords]
            top_op_mat_col_add[-2:] = self.operating_matrix_corner_int_IDs[
                (bottom_corner_coords[0], bottom_corner_coords[0]),
                (bottom_corner_coords[1], bottom_corner_coords[1]),
            ].flatten()
            top_op_mat_data_add[-6:-4] = (
                _delta_t
                * nine_node_map[_interior_corners[2], :][
                    corners_antimasks[2, [3, 4]]
                ].flatten()
            )
            top_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[3], :][
                    corners_antimasks[3, [2, 3]]
                ].flatten()
            )
            top_op_mat_data_add[-2] = (
                _delta_t
                * nine_node_map[_interior_corners[2], :][corners_antimasks[2, 2]]
            )
            top_op_mat_data_add[-1] = (
                _delta_t
                * nine_node_map[_interior_corners[3], :][corners_antimasks[3, 4]]
            )
        else:
            raise NameError(
                """Something is very wrong with your boundary
                            conditions...!"""
            )

        if self.top_flag == 1:
            # goes to RHS only
            _mat_RHS[top_interior_IDs] -= _delta_t * np.sum(
                nine_node_map[_top_list, :][:, top_antimask]
                * elev[
                    _top_list.reshape((len(_top_list), 1))
                    + (modulator_mask[top_antimask]).reshape((1, 3))
                ],
                axis=1,
            )
            # ...& the corners
            edges = [(0, 0), (0, 0), (3, 4), (2, 3)]
            for i in [2, 3]:
                edge_list = edges[i]
                _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                    nine_node_map[_interior_corners[i], :][
                        corners_antimasks[i, edge_list]
                    ]
                    * elev[
                        _interior_corners[i]
                        + modulator_mask[corners_antimasks[i, edge_list]]
                    ]
                )
            top_op_mat_row_add = np.empty(0)
            top_op_mat_col_add = np.empty(0)
            top_op_mat_data_add = np.empty(0)
        elif self.top_flag == 4 or self.top_flag == 2:
            top_op_mat_row_add = np.empty((top_interior_IDs.size * 3 + 6))
            top_op_mat_col_add = np.empty((top_interior_IDs.size * 3 + 6))
            top_op_mat_data_add = np.empty((top_interior_IDs.size * 3 + 6))
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            top_op_mat_row_add[: (top_interior_IDs.size * 3)] = np.repeat(
                top_interior_IDs, 3
            )
            top_op_mat_col_add[: (top_interior_IDs.size * 3)] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.top_interior_IDs, :][
                    :, self.top_mask[3:6]
                ]
            ).flatten()
            top_op_mat_data_add[: (top_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_top_list, :][:, top_antimask]).flatten()
            )
            # ...& the corners
            this_corner_coords = np.array([2, 3])
            top_op_mat_row_add[-6:-2] = np.repeat(
                corner_interior_IDs[this_corner_coords], 2
            )
            top_op_mat_col_add[-6:-2] = self.operating_matrix_corner_int_IDs[
                this_corner_coords.reshape((2, 1)), this_corner_coords
            ].flatten()
            top_op_mat_row_add[-2:] = corner_interior_IDs[this_corner_coords]
            top_op_mat_col_add[-2:] = self.operating_matrix_corner_int_IDs[
                (this_corner_coords[0], this_corner_coords[0]),
                (this_corner_coords[1], this_corner_coords[1]),
            ].flatten()
            top_op_mat_data_add[-6:-4] = (
                _delta_t
                * nine_node_map[_interior_corners[2], :][
                    corners_antimasks[2, [3, 4]]
                ].flatten()
            )
            top_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[3], :][
                    corners_antimasks[3, [2, 3]]
                ].flatten()
            )
            top_op_mat_data_add[-2] = (
                _delta_t
                * nine_node_map[_interior_corners[2], :][corners_antimasks[2, 2]]
            )
            top_op_mat_data_add[-1] = (
                _delta_t
                * nine_node_map[_interior_corners[3], :][corners_antimasks[3, 4]]
            )
            if self.top_flag == 2:
                _mat_RHS[top_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_top_list, :][:, top_antimask]
                    * self.fixed_gradient_offset_map[
                        _top_list.reshape((len(_top_list), 1))
                        + (modulator_mask[top_antimask]).reshape((1, 3))
                    ],
                    axis=1,
                )
                # ...& the corners
                edges = [(0, 0), (0, 0), (3, 4), (2, 3)]
                for i in [2, 3]:
                    edge_list = edges[i]
                    _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                        nine_node_map[_interior_corners[i], :][
                            corners_antimasks[i, edge_list]
                        ]
                        * self.fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self.top_flag == 3:
            pass  # dealt with above
        else:
            raise NameError(
                """Something is very wrong with your boundary
                            conditions...!"""
            )

        if self.left_flag == 1:
            # goes to RHS only
            _mat_RHS[left_interior_IDs] -= _delta_t * np.sum(
                nine_node_map[_left_list, :][:, left_antimask]
                * elev[
                    _left_list.reshape((len(_left_list), 1))
                    + (modulator_mask[left_antimask]).reshape((1, 3))
                ],
                axis=1,
            )
            # ...& the corners
            edges = [(3, 4), (0, 0), (0, 1), (0, 0)]
            for i in [0, 2]:
                edge_list = edges[i]
                _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                    nine_node_map[_interior_corners[i], :][
                        corners_antimasks[i, edge_list]
                    ]
                    * elev[
                        _interior_corners[i]
                        + modulator_mask[corners_antimasks[i, edge_list]]
                    ]
                )
            left_op_mat_row_add = np.empty(0)
            left_op_mat_col_add = np.empty(0)
            left_op_mat_data_add = np.empty(0)
        elif self.left_flag == 4 or self.left_flag == 2:
            left_op_mat_row_add = np.empty((left_interior_IDs.size * 3 + 4))
            left_op_mat_col_add = np.empty((left_interior_IDs.size * 3 + 4))
            left_op_mat_data_add = np.empty((left_interior_IDs.size * 3 + 4))
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            left_op_mat_row_add[: (left_interior_IDs.size * 3)] = np.repeat(
                left_interior_IDs, 3
            )
            left_op_mat_col_add[
                : (left_interior_IDs.size * 3)
            ] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.left_interior_IDs, :][
                    :, self.left_mask[::2]
                ]
            ).flatten()
            left_op_mat_data_add[: (left_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_left_list, :][:, left_antimask]).flatten()
            )
            # ...& the corners
            this_corner_coords = np.array([0, 2])
            left_op_mat_row_add[-4:] = np.repeat(
                corner_interior_IDs[this_corner_coords], 2
            )
            left_op_mat_col_add[-4:] = self.operating_matrix_corner_int_IDs[
                this_corner_coords.reshape((2, 1)), this_corner_coords
            ].flatten()
            left_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[0], :][
                    corners_antimasks[0, [3, 4]]
                ].flatten()
            )
            left_op_mat_data_add[-2:] = (
                _delta_t
                * nine_node_map[_interior_corners[2], :][
                    corners_antimasks[2, [0, 1]]
                ].flatten()
            )
            if self.left_flag == 2:
                _mat_RHS[left_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_left_list, :][:, left_antimask]
                    * self.fixed_gradient_offset_map[
                        _left_list.reshape((len(_left_list), 1))
                        + (modulator_mask[left_antimask]).reshape((1, 3))
                    ],
                    axis=1,
                )
                # ...& the corners
                edges = [(3, 4), (0, 0), (0, 1), (0, 0)]
                for i in [0, 2]:
                    edge_list = edges[i]
                    _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                        nine_node_map[_interior_corners[i], :][
                            corners_antimasks[i, edge_list]
                        ]
                        * self.fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self.left_flag == 3:
            left_op_mat_row_add = np.empty((left_interior_IDs.size * 3 + 4))
            left_op_mat_col_add = np.empty((left_interior_IDs.size * 3 + 4))
            left_op_mat_data_add = np.empty((left_interior_IDs.size * 3 + 4))
            left_op_mat_row_add[: (left_interior_IDs.size * 3)] = np.repeat(
                left_interior_IDs, 3
            )
            left_op_mat_col_add[
                : (left_interior_IDs.size * 3)
            ] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.right_interior_IDs, :][
                    :, self.right_mask[1::2]
                ]
            ).flatten()
            left_op_mat_data_add[: (left_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_left_list, :][:, left_antimask]).flatten()
            )
            right_op_mat_row_add = np.empty((right_interior_IDs.size * 3 + 4))
            right_op_mat_col_add = np.empty((right_interior_IDs.size * 3 + 4))
            right_op_mat_data_add = np.empty((right_interior_IDs.size * 3 + 4))
            right_op_mat_row_add[: (right_interior_IDs.size * 3)] = np.repeat(
                right_interior_IDs, 3
            )
            right_op_mat_col_add[
                : (right_interior_IDs.size * 3)
            ] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.left_interior_IDs, :][
                    :, self.left_mask[::2]
                ]
            ).flatten()
            right_op_mat_data_add[: (right_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_right_list, :][:, right_antimask]).flatten()
            )
            # & the corners
            left_corner_coords = np.array([0, 2])
            right_corner_coords = np.array([1, 3])
            left_op_mat_row_add[-4:] = np.repeat(
                corner_interior_IDs[left_corner_coords], 2
            )
            left_op_mat_col_add[-4:] = self.operating_matrix_corner_int_IDs[
                right_corner_coords.reshape((2, 1)), right_corner_coords
            ].flatten()
            left_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[0], :][
                    corners_antimasks[0, [3, 4]]
                ].flatten()
            )
            left_op_mat_data_add[-2:] = (
                _delta_t
                * nine_node_map[_interior_corners[2], :][
                    corners_antimasks[2, [0, 1]]
                ].flatten()
            )
            right_op_mat_row_add[-4:] = np.repeat(
                corner_interior_IDs[right_corner_coords], 2
            )
            right_op_mat_col_add[-4:] = self.operating_matrix_corner_int_IDs[
                left_corner_coords.reshape((2, 1)), left_corner_coords
            ].flatten()
            right_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[1], :][
                    corners_antimasks[1, [3, 4]]
                ].flatten()
            )
            right_op_mat_data_add[-2:] = (
                _delta_t
                * nine_node_map[_interior_corners[3], :][
                    corners_antimasks[3, [0, 1]]
                ].flatten()
            )
        else:
            raise NameError(
                """Something is very wrong with your boundary
                            conditions...!"""
            )

        if self.right_flag == 1:
            # goes to RHS only
            _mat_RHS[right_interior_IDs] -= _delta_t * np.sum(
                nine_node_map[_right_list, :][:, right_antimask]
                * elev[
                    _right_list.reshape((len(_right_list), 1))
                    + (modulator_mask[right_antimask]).reshape((1, 3))
                ],
                axis=1,
            )
            # ...& the corners
            edges = [(0, 0), (3, 4), (0, 0), (0, 1)]
            for i in [1, 3]:
                edge_list = edges[i]
                _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                    nine_node_map[_interior_corners[i], :][
                        corners_antimasks[i, edge_list]
                    ]
                    * elev[
                        _interior_corners[i]
                        + modulator_mask[corners_antimasks[i, edge_list]]
                    ]
                )
            right_op_mat_row_add = np.empty(0)
            right_op_mat_col_add = np.empty(0)
            right_op_mat_data_add = np.empty(0)
        elif self.right_flag == 4 or self.right_flag == 2:
            right_op_mat_row_add = np.empty((right_interior_IDs.size * 3 + 4))
            right_op_mat_col_add = np.empty((right_interior_IDs.size * 3 + 4))
            right_op_mat_data_add = np.empty((right_interior_IDs.size * 3 + 4))
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            right_op_mat_row_add[: (right_interior_IDs.size * 3)] = np.repeat(
                right_interior_IDs, 3
            )
            right_op_mat_col_add[
                : (right_interior_IDs.size * 3)
            ] = self._realIDtointerior(
                self.operating_matrix_ID_map[self.right_interior_IDs, :][
                    :, self.right_mask[1::2]
                ]
            ).flatten()
            right_op_mat_data_add[: (right_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_right_list, :][:, right_antimask]).flatten()
            )
            # ...& the corners
            this_corner_coords = np.array([1, 3])
            right_op_mat_row_add[-4:] = np.repeat(
                corner_interior_IDs[this_corner_coords], 2
            )
            right_op_mat_col_add[-4:] = self.operating_matrix_corner_int_IDs[
                this_corner_coords.reshape((2, 1)), this_corner_coords
            ].flatten()
            right_op_mat_data_add[-4:-2] = (
                _delta_t
                * nine_node_map[_interior_corners[1], :][
                    corners_antimasks[1, [3, 4]]
                ].flatten()
            )
            right_op_mat_data_add[-2:] = (
                _delta_t
                * nine_node_map[_interior_corners[3], :][
                    corners_antimasks[3, [0, 1]]
                ].flatten()
            )
            if self.right_flag == 2:
                _mat_RHS[right_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_right_list, :][:, right_antimask]
                    * self.fixed_gradient_offset_map[
                        _right_list.reshape((len(_right_list), 1))
                        + (modulator_mask[right_antimask]).reshape((1, 3))
                    ],
                    axis=1,
                )
                # ...& the corners
                edges = [(0, 0), (3, 4), (0, 0), (0, 1)]
                for i in [1, 3]:
                    edge_list = edges[i]
                    _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                        nine_node_map[_interior_corners[i], :][
                            corners_antimasks[i, edge_list]
                        ]
                        * self.fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self.top_flag == 3:
            pass  # dealt with above
        else:
            raise NameError(
                """Something is very wrong with your boundary
                            conditions...!"""
            )

        # new approach using COO sparse matrix requires we build the matrix
        # only now...
        self._operating_matrix = sparse.coo_matrix(
            (
                np.concatenate(
                    (
                        core_op_mat_data,
                        corners_op_mat_data,
                        bottom_op_mat_data,
                        top_op_mat_data,
                        left_op_mat_data,
                        right_op_mat_data,
                        bottom_op_mat_data_add,
                        top_op_mat_data_add,
                        left_op_mat_data_add,
                        right_op_mat_data_add,
                    )
                ),
                (
                    np.concatenate(
                        (
                            core_op_mat_row,
                            corners_op_mat_row,
                            bottom_op_mat_row,
                            top_op_mat_row,
                            left_op_mat_row,
                            right_op_mat_row,
                            bottom_op_mat_row_add,
                            top_op_mat_row_add,
                            left_op_mat_row_add,
                            right_op_mat_row_add,
                        )
                    ),
                    np.concatenate(
                        (
                            core_op_mat_col,
                            corners_op_mat_col,
                            bottom_op_mat_col,
                            top_op_mat_col,
                            left_op_mat_col,
                            right_op_mat_col,
                            bottom_op_mat_col_add,
                            top_op_mat_col_add,
                            left_op_mat_col_add,
                            right_op_mat_col_add,
                        )
                    ),
                ),
            ),
            shape=(n_interior_nodes, n_interior_nodes),
        ).tocsr()
        self._mat_RHS = _mat_RHS

    # These methods translate ID numbers between arrays of differing sizes
    def _realIDtointerior(self, ID):
        ncols = self.ncols
        interior_ID = (ID // ncols - 1) * (ncols - 2) + (ID % ncols) - 1
        if np.any(interior_ID < 0) or np.any(interior_ID >= self.ninteriornodes):
            raise NameError(
                """One of the supplied nodes was outside the
                            interior grid!"""
            )
        else:
            return interior_ID.astype(int)

    def _interiorIDtoreal(self, ID):
        IGW = self.interior_grid_width
        real_ID = (ID // IGW + 1) * self.ncols + (ID % IGW) + 1
        assert np.all(real_ID < self.nnodes)
        return real_ID.astype(int)

    def _realIDtocore(self, ID):
        ncols = self.ncols
        core_ID = (ID // ncols - 2) * (ncols - 4) + (ID % ncols) - 2
        if np.any(core_ID < 0) or np.any(core_ID >= self.ncorenodes):
            raise NameError(
                """One of the supplied nodes was outside the
                            core grid!"""
            )
        else:
            return core_ID.astype(int)

    def _coreIDtoreal(self, ID):
        CCW = self.core_cell_width
        real_ID = (ID // CCW + 2) * self.ncols + (ID % CCW) + 2
        assert np.all(real_ID < self.nnodes)
        return real_ID.astype(int)

    def _interiorIDtocore(self, ID):
        IGW = self.interior_grid_width
        core_ID = (ID // IGW - 1) * (self.ncols - 4) + (ID % IGW) - 1
        if np.any(core_ID < 0) or np.any(core_ID >= self.ncorenodes):
            raise NameError(
                """One of the supplied nodes was outside the
                            core grid!"""
            )
        else:
            return core_ID.astype(int)

    def _coreIDtointerior(self, ID):
        CCW = self.core_cell_width
        interior_ID = (ID // CCW + 1) * (self.ncols - 2) + (ID % CCW) + 1
        assert np.all(interior_ID < self.ninteriornodes)
        return interior_ID.astype(int)

    def diffuse(self, grid_in, elapsed_time, num_uplift_implicit_comps=1):
        """
        This is the "old style" run method of the class, superceded by
        :func:`run_one_step`.
        Takes *grid_in*, the model grid, and *elapsed_time*, the
        total model time elapsed so far.

        *grid_in* must contain the field to diffuse, which defaults to
        'topographic__elevation'. This can be overridden with the
        values_to_diffuse property in the input file.

        See the class docstring for a list of the other properties necessary
        in the input file for this component to run.

        Note that the implicit nature of this component requires it to
        incorporate uplift into its execution in order to stay stable.
        If you only have one module that requires this, do not add uplift
        manually in your loop; this method will include uplift automatically.

        If more than one of your components has this requirement, set
        *num_uplift_implicit_comps* to the total number of components that
        do.
        """
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code
        if self.internal_uplifts:
            # this is adhoc to fix for the duration of Germany visit
            self._uplift = self.inputs.read_float("uplift_rate")
            self._delta_t = self.timestep_in
            self._set_variables(self.grid)
            _interior_elevs = linalg.spsolve(self._operating_matrix, self._mat_RHS)
            self.grid["node"][self.values_to_diffuse][
                self.interior_IDs_as_real
            ] = _interior_elevs
            grid_in = self.grid
        else:
            self._gear_timestep(self.timestep_in, grid_in)
            for i in range(self.internal_repeats):
                grid_in["node"][self.values_to_diffuse][:] = (
                    self.grid["node"][self.values_to_diffuse] + self.uplift_per_step
                )
                # Initialize the variables for the step:
                self._set_variables(grid_in)
                # Solve interior of grid:
                _interior_elevs = linalg.spsolve(self._operating_matrix, self._mat_RHS)
                # this fn solves Ax=B for x

                # Handle the BC cells; test common cases first for speed
                self.grid["node"][self.values_to_diffuse][
                    self.interior_IDs_as_real
                ] = _interior_elevs

                # if BC==1 or BC==4, don't need to take any action; in both
                # cases the values are unchanged.
                if self.fixed_grad_BCs_present:
                    self.grid["node"][self.values_to_diffuse][
                        grid_in.fixed_gradient_node_properties["boundary_node_IDs"]
                    ] = (
                        self.grid["node"][self.values_to_diffuse][
                            self.grid.fixed_gradient_node_properties["anchor_node_IDs"]
                        ]
                        + self.grid.fixed_gradient_node_properties["values_to_add"]
                    )
                if self.looped_BCs_present:
                    self.grid["node"][self.values_to_diffuse][
                        self.grid.looped_node_properties["boundary_node_IDs"]
                    ] = self.grid["node"][self.values_to_diffuse][
                        self.grid.looped_node_properties["linked_node_IDs"]
                    ]

        return self.grid

    def run_one_step(self, dt):
        """Run the diffuser for one timestep, dt.

        This is the primary method of the class.

        Parameters
        ----------
        dt : float (time)
            The imposed timestep.
        """
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code
        if self.internal_uplifts:
            self._delta_t = self.timestep_in
            self._set_variables(self.grid)
            _interior_elevs = linalg.spsolve(self._operating_matrix, self._mat_RHS)
            self.grid["node"][self.values_to_diffuse][
                self.interior_IDs_as_real
            ] = _interior_elevs
        else:
            self._gear_timestep(dt, self.grid)
            for i in range(self.internal_repeats):
                # Initialize the variables for the step:
                self._set_variables(self.grid)
                # Solve interior of grid:
                _interior_elevs = linalg.spsolve(self._operating_matrix, self._mat_RHS)
                # this fn solves Ax=B for x

                # Handle the BC cells; test common cases first for speed
                self.grid["node"][self.values_to_diffuse][
                    self.interior_IDs_as_real
                ] = _interior_elevs

        # if BC==1 or BC==4, don't need to take any action; in both
        # cases the values are unchanged.
        if self.fixed_grad_BCs_present:
            self.grid["node"][self.values_to_diffuse][
                self.grid.fixed_gradient_node_properties["boundary_node_IDs"]
            ] = (
                self.grid["node"][self.values_to_diffuse][
                    self.grid.fixed_gradient_node_properties["anchor_node_IDs"]
                ]
                + self.grid.fixed_gradient_node_properties["values_to_add"]
            )
        if self.looped_BCs_present:
            self.grid["node"][self.values_to_diffuse][
                self.grid.looped_node_properties["boundary_node_IDs"]
            ] = self.grid["node"][self.values_to_diffuse][
                self.grid.looped_node_properties["linked_node_IDs"]
            ]
