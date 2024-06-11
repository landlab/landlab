import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from landlab import Component

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
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> nl = PerronNLDiffuse(mg, nonlinear_diffusivity=1.0)
    >>> dt = 100.0
    >>> nt = 20
    >>> uplift_rate = 0.001
    >>> for i in range(nt):
    ...     z[mg.core_nodes] += uplift_rate * dt
    ...     nl.run_one_step(dt)
    ...
    >>> z_target = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 0.00778637, 0.0075553, 0.00778637, 0.0],
    ...     [0.0, 0.0075553, 0.0078053, 0.0075553, 0.0],
    ...     [0.0, 0.00778637, 0.0075553, 0.00778637, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> np.allclose(z.reshape(mg.shape), z_target)
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Perron, J. (2011). Numerical methods for nonlinear hillslope transport laws.
    Journal of Geophysical Research  116(F2), 23 - 13.
    https://dx.doi.org/10.1029/2010jf001801

    """

    _name = "PerronNLDiffuse"

    _unit_agnostic = True

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        }
    }

    def __init__(
        self,
        grid,
        nonlinear_diffusivity=0.01,
        S_crit=33.0 * np.pi / 180.0,
        rock_density=2700.0,
        sed_density=2700.0,
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
        super().__init__(grid)

        self._bc_set_code = self._grid.bc_set_code
        self._values_to_diffuse = "topographic__elevation"
        self._kappa = nonlinear_diffusivity
        self._rock_density = rock_density
        self._sed_density = sed_density
        self._S_crit = S_crit
        self._uplift = 0.0
        self._delta_x = grid.dx
        self._delta_y = grid.dy
        self._one_over_delta_x = 1.0 / self._delta_x
        self._one_over_delta_y = 1.0 / self._delta_y
        self._one_over_delta_x_sqd = self._one_over_delta_x**2.0
        self._one_over_delta_y_sqd = self._one_over_delta_y**2.0
        self._b = 1.0 / self._S_crit**2.0

        ncols = grid.number_of_node_columns
        self._ncols = ncols
        nrows = grid.number_of_node_rows
        self._nrows = nrows
        nnodes = grid.number_of_nodes
        self._nnodes = nnodes
        ninteriornodes = grid.number_of_interior_nodes
        ncorenodes = ninteriornodes - 2 * (ncols + nrows - 6)
        self._ninteriornodes = ninteriornodes
        self._interior_grid_width = ncols - 2
        self._core_cell_width = ncols - 4

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
        self._corenodesbyintIDs = self._realIDtointerior(self._core_nodes)
        self._ncorenodes = len(self._core_nodes)

        self._corner_interior_IDs = self._realIDtointerior(self._interior_corners)
        # ^i.e., interior corners as interior IDs
        self._bottom_interior_IDs = self._realIDtointerior(np.array(_bottom_list))
        self._top_interior_IDs = self._realIDtointerior(np.array(_top_list))
        self._left_interior_IDs = self._realIDtointerior(np.array(_left_list))
        self._right_interior_IDs = self._realIDtointerior(np.array(_right_list))

        # build an ID map to let us easily map the variables of the core nodes
        # onto the operating matrix:
        # This array is ninteriornodes long, but the IDs it contains are
        # REAL IDs
        operating_matrix_ID_map = np.empty((ninteriornodes, 9))
        self._interior_IDs_as_real = self._interiorIDtoreal(np.arange(ninteriornodes))
        for j in range(ninteriornodes):
            i = self._interior_IDs_as_real[j]
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
        self._operating_matrix_ID_map = operating_matrix_ID_map
        self._operating_matrix_core_int_IDs = self._realIDtointerior(
            operating_matrix_ID_map[self._corenodesbyintIDs, :]
        )
        # ^shape(ncorenodes,9)
        # see below for corner and edge maps

        # Build masks for the edges and corners to be applied to the operating
        # matrix map.
        # Antimasks are the boundary nodes, masks are "normal"
        self._topleft_mask = [1, 2, 4, 5]
        topleft_antimask = [0, 3, 6, 7, 8]
        self._topright_mask = [0, 1, 3, 4]
        topright_antimask = [2, 5, 6, 7, 8]
        self._bottomleft_mask = [4, 5, 7, 8]
        bottomleft_antimask = [0, 1, 2, 3, 6]
        self._bottomright_mask = [3, 4, 6, 7]
        bottomright_antimask = [0, 1, 2, 5, 8]
        self._corners_masks = np.vstack(
            (
                self._bottomleft_mask,
                self._bottomright_mask,
                self._topleft_mask,
                self._topright_mask,
            )
        )
        # ^(each_corner,mask_for_each_corner)
        self._corners_antimasks = np.vstack(
            (
                bottomleft_antimask,
                bottomright_antimask,
                topleft_antimask,
                topright_antimask,
            )
        )
        # ^so shape becomes (4,5)
        self._left_mask = [1, 2, 4, 5, 7, 8]
        self._left_antimask = [0, 3, 6]
        self._top_mask = [0, 1, 2, 3, 4, 5]
        self._top_antimask = [6, 7, 8]
        self._right_mask = [0, 1, 3, 4, 6, 7]
        self._right_antimask = [2, 5, 8]
        self._bottom_mask = [3, 4, 5, 6, 7, 8]
        self._bottom_antimask = [0, 1, 2]
        self._antimask_corner_position = [0, 2, 2, 4]
        # ^this is the position w/i the corner antimasks that the true corner
        # actually occupies

        self._modulator_mask = np.array(
            [-ncols - 1, -ncols, -ncols + 1, -1, 0, 1, ncols - 1, ncols, ncols + 1]
        )

        self.updated_boundary_conditions()

    def updated_boundary_conditions(self):
        """Call if grid BCs are updated after component instantiation."""
        grid = self._grid
        nrows = self._nrows
        ncols = self._ncols
        # ^Set up terms for BC handling (still feels very clumsy)
        bottom_edge = grid.nodes_at_bottom_edge[1:-1]
        top_edge = grid.nodes_at_top_edge[1:-1]
        left_edge = grid.nodes_at_left_edge[1:-1]
        right_edge = grid.nodes_at_right_edge[1:-1]
        self._bottom_flag = 1
        self._top_flag = 1
        self._left_flag = 1
        self._right_flag = 1
        # self._corner_flags = [1,1,1,1] #In ID order, so BL,BR,TL,TR
        if np.all(grid.status_at_node[bottom_edge] == 4):
            # ^This should be all of them, or none of them
            self._bottom_flag = 4
        elif np.all(grid.status_at_node[bottom_edge] == 3):
            self._bottom_flag = 3
        elif np.all(grid.status_at_node[bottom_edge] == 2):
            self._bottom_flag = 2
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
            self._top_flag = 4
        elif np.all(grid.status_at_node[top_edge] == 3):
            self._top_flag = 3
        elif np.all(grid.status_at_node[top_edge] == 2):
            self._top_flag = 2
        elif np.all(grid.status_at_node[top_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
        if np.all(grid.status_at_node[left_edge] == 4):
            self._left_flag = 4
        elif np.all(grid.status_at_node[left_edge] == 3):
            self._left_flag = 3
        elif np.all(grid.status_at_node[left_edge] == 2):
            self._left_flag = 2
        elif np.all(grid.status_at_node[left_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )
        if np.all(grid.status_at_node[right_edge] == 4):
            self._right_flag = 4
        elif np.all(grid.status_at_node[right_edge] == 3):
            self._right_flag = 3
        elif np.all(grid.status_at_node[right_edge] == 2):
            self._right_flag = 2
        elif np.all(grid.status_at_node[right_edge] == 1):
            pass
        else:
            raise NameError(
                "Different cells on the same grid edge have "
                "different boundary statuses"
            )

        self._fixed_grad_BCs_present = (
            self._bottom_flag == 2
            or self._top_flag == 2
            or self._left_flag == 2
            or self._right_flag == 2
        )
        self._looped_BCs_present = (
            self._bottom_flag == 3
            or self._top_flag == 3
            or self._left_flag == 3
            or self._right_flag == 3
        )
        if (
            self._fixed_grad_BCs_present
            and self._values_to_diffuse != grid.fixed_gradient_of
        ):
            raise ValueError(
                "Boundary conditions set in the grid don't "
                "apply to the data the diffuser is trying to "
                "work with"
            )

        if np.any(grid.status_at_node == 2):
            self._fixed_grad_offset_map = np.empty(nrows * ncols, dtype=float)
            self._fixed_grad_anchor_map = np.empty_like(self._fixed_grad_offset_map)
            self._fixed_grad_offset_map[
                grid.fixed_gradient_node_properties["boundary_node_IDs"]
            ] = grid.fixed_gradient_node_properties["values_to_add"]

        self._corner_flags = grid.status_at_node[[0, ncols - 1, -ncols, -1]]

        op_mat_just_corners = self._operating_matrix_ID_map[
            self._corner_interior_IDs, :
        ]
        op_mat_cnr0 = op_mat_just_corners[0, self._bottomleft_mask]
        op_mat_cnr1 = op_mat_just_corners[1, self._bottomright_mask]
        op_mat_cnr2 = op_mat_just_corners[2, self._topleft_mask]
        op_mat_cnr3 = op_mat_just_corners[3, self._topright_mask]
        op_mat_just_active_cnrs = np.vstack(
            (op_mat_cnr0, op_mat_cnr1, op_mat_cnr2, op_mat_cnr3)
        )
        self._operating_matrix_corner_int_IDs = self._realIDtointerior(
            op_mat_just_active_cnrs
        )
        # ^(4corners,4nodesactivepercorner)
        self._operating_matrix_bottom_int_IDs = self._realIDtointerior(
            self._operating_matrix_ID_map[self._bottom_interior_IDs, :][
                :, self._bottom_mask
            ]
        )
        # ^(nbottomnodes,6activenodeseach)
        self._operating_matrix_top_int_IDs = self._realIDtointerior(
            self._operating_matrix_ID_map[self._top_interior_IDs, :][:, self._top_mask]
        )
        self._operating_matrix_left_int_IDs = self._realIDtointerior(
            self._operating_matrix_ID_map[self._left_interior_IDs, :][
                :, self._left_mask
            ]
        )
        self._operating_matrix_right_int_IDs = self._realIDtointerior(
            self._operating_matrix_ID_map[self._right_interior_IDs, :][
                :, self._right_mask
            ]
        )

    def _gear_timestep(self, timestep_in, new_grid):
        """This method allows the gearing between the model run step and the
        component (shorter) step.

        The method becomes unstable if S>Scrit, so we test to prevent
        this. We implicitly assume the initial condition does not
        contain slopes > Scrit. If the method persistently explodes,
        this may be the problem.
        """
        extended_elevs = np.empty(self._grid.number_of_nodes + 1, dtype=float)
        extended_elevs[-1] = np.nan
        node_neighbors = self._grid.active_adjacent_nodes_at_node
        extended_elevs[:-1] = new_grid["node"][self._values_to_diffuse]
        max_offset = np.nanmax(
            np.fabs(
                extended_elevs[:-1][node_neighbors]
                - extended_elevs[:-1].reshape((self._grid.number_of_nodes, 1))
            )
        )
        if max_offset > np.tan(self._S_crit) * min(self._grid.dx, self._grid.dy):
            # ^using S not tan(S) adds a buffer - but not appropriate
            self._internal_repeats = (
                int(
                    max_offset
                    // (np.tan(self._S_crit) * min(self._grid.dx, self._grid.dy))
                )
                + 1
            )
            # now we rig it so the actual timestep is an integer divisor
            # of T_in:
            self._delta_t = timestep_in / self._internal_repeats
            self._uplift_per_step = (
                new_grid["node"][self._values_to_diffuse]
                - self._grid["node"][self._values_to_diffuse]
            ) / self._internal_repeats
            if self._internal_repeats > 10000:
                raise ValueError(
                    """Uplift rate is too high; solution is not
                                 stable!!"""
                )
        else:
            self._internal_repeats = 1
            self._delta_t = timestep_in
            self._uplift_per_step = (
                new_grid["node"][self._values_to_diffuse]
                - self._grid["node"][self._values_to_diffuse]
            )
        return self._delta_t

    def _set_variables(self, grid):
        """This function sets the variables needed for update().

        Now vectorized, shouold run faster. At the moment, this method
        can only handle fixed value BCs.
        """
        n_interior_nodes = grid.number_of_interior_nodes

        # Initialize the local builder lists
        _mat_RHS = np.zeros(n_interior_nodes)

        try:
            elev = grid["node"][self._values_to_diffuse]
        except KeyError as exc:
            raise NameError("elevations not found in grid!") from exc
        try:
            _delta_t = self._delta_t
        except AttributeError as exc:
            raise NameError(
                "Timestep not set! Call _gear_timestep(tstep) "
                "after initializing the component, but before running it."
            ) from exc
        _one_over_delta_x = self._one_over_delta_x
        _one_over_delta_x_sqd = self._one_over_delta_x_sqd
        _one_over_delta_y = self._one_over_delta_y
        _one_over_delta_y_sqd = self._one_over_delta_y_sqd
        _kappa = self._kappa
        _b = self._b
        _S_crit = self._S_crit
        _core_nodes = self._core_nodes
        corenodesbyintIDs = self._corenodesbyintIDs
        operating_matrix_core_int_IDs = self._operating_matrix_core_int_IDs
        operating_matrix_corner_int_IDs = self._operating_matrix_corner_int_IDs
        _interior_corners = self._interior_corners
        corners_antimasks = self._corners_antimasks
        corner_interior_IDs = self._corner_interior_IDs
        modulator_mask = self._modulator_mask
        corner_flags = self._corner_flags
        bottom_interior_IDs = self._bottom_interior_IDs
        top_interior_IDs = self._top_interior_IDs
        left_interior_IDs = self._left_interior_IDs
        right_interior_IDs = self._right_interior_IDs
        bottom_antimask = self._bottom_antimask
        _bottom_list = self._bottom_list
        top_antimask = self._top_antimask
        _top_list = self._top_list
        left_antimask = self._left_antimask
        _left_list = self._left_list
        right_antimask = self._right_antimask
        _right_list = self._right_list

        # Need to modify the "effective" values of the edge nodes if any of
        # the edges are inactive:
        if self._bottom_flag == 4:
            bottom_edge, inside_bottom_edge = grid.nodes[(0, 1), :]
            elev[bottom_edge] = elev[inside_bottom_edge]
            # corners are special cases, and assumed linked to the bottom and
            # top edge BCs...
            elev[bottom_edge[0]] = elev[inside_bottom_edge[1]]
            elev[bottom_edge[-1]] = elev[inside_bottom_edge[-2]]
        if self._top_flag == 4:
            top_edge, inside_top_edge = grid.nodes[(-1, -2), :]
            elev[top_edge] = elev[inside_top_edge]
            # corners are special cases, and assumed linked to the bottom and
            # top edge BCs...
            elev[top_edge[0]] = elev[inside_top_edge[1]]
            elev[top_edge[-1]] = elev[inside_top_edge[-2]]
        if self._left_flag == 4:
            left_edge = grid.nodes[1:-1, 0]
            inside_left_edge = grid.nodes[1:-1, 1]
            elev[left_edge] = elev[inside_left_edge]
        if self._right_flag == 4:
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
        corners_op_mat_row = np.repeat(self._corner_interior_IDs, 4)
        corners_op_mat_col = operating_matrix_corner_int_IDs.astype(int).flatten()
        corners_op_mat_data = nine_node_map[_interior_corners, :][
            (np.arange(4).reshape((4, 1)), self._corners_masks)
        ].flatten()
        # ^1st index gives (4,9), 2nd reduces to (4,4), then flattened
        for i in range(4):  # loop over each corner, as so few
            # Note that this ONLY ADDS THE VALUES FOR THE TRUE GRID CORNERS.
            # The sides get done in the edge tests, below.
            if corner_flags[i] == 1:
                true_corner = self._antimask_corner_position[i]
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
                true_corner = self._antimask_corner_position[i]
                _mat_RHS[corner_interior_IDs[i]] -= _delta_t * np.sum(
                    nine_node_map[_interior_corners[i], :][
                        corners_antimasks[i, true_corner]
                    ]
                    * self._fixed_gradient_offset_map[
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
        bottom_op_mat_col = self._operating_matrix_bottom_int_IDs.astype(int).flatten()
        top_op_mat_col = self._operating_matrix_top_int_IDs.astype(int).flatten()
        left_op_mat_col = self._operating_matrix_left_int_IDs.astype(int).flatten()
        right_op_mat_col = self._operating_matrix_right_int_IDs.astype(int).flatten()
        bottom_op_mat_data = nine_node_map[_bottom_list, :][
            :, self._bottom_mask
        ].flatten()
        top_op_mat_data = nine_node_map[_top_list, :][:, self._top_mask].flatten()
        left_op_mat_data = nine_node_map[_left_list, :][:, self._left_mask].flatten()
        right_op_mat_data = nine_node_map[_right_list, :][:, self._right_mask].flatten()

        if self._bottom_flag == 1:
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
        elif self._bottom_flag == 4 or self._bottom_flag == 2:
            # ^i.e., fixed zero gradient (4) or more general case...
            bottom_op_mat_row_add = np.empty(bottom_interior_IDs.size * 3 + 6)
            bottom_op_mat_col_add = np.empty(bottom_interior_IDs.size * 3 + 6)
            bottom_op_mat_data_add = np.empty(bottom_interior_IDs.size * 3 + 6)
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            bottom_op_mat_row_add[: (bottom_interior_IDs.size * 3)] = np.repeat(
                bottom_interior_IDs, 3
            )
            bottom_op_mat_col_add[: (bottom_interior_IDs.size * 3)] = (
                self._realIDtointerior(
                    self._operating_matrix_ID_map[self._bottom_interior_IDs, :][
                        :, self._bottom_mask[0:3]
                    ]
                ).flatten()
            )
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
            bottom_op_mat_col_add[-6:-2] = self._operating_matrix_corner_int_IDs[
                this_corner_coords.reshape((2, 1)), this_corner_coords
            ].flatten()
            bottom_op_mat_row_add[-2:] = corner_interior_IDs[this_corner_coords]
            bottom_op_mat_col_add[-2:] = self._operating_matrix_corner_int_IDs[
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
            if self._bottom_flag == 2:
                # Read the offsets from the map we made in the __init__,
                # use them as constant terms, incorporated into RHS
                _mat_RHS[bottom_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_bottom_list, :][:, bottom_antimask]
                    * self._fixed_gradient_offset_map[
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
                        * self._fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self._bottom_flag == 3:
            # This will handle both top and bottom BCs...
            bottom_op_mat_row_add = np.empty(bottom_interior_IDs.size * 3 + 6)
            bottom_op_mat_col_add = np.empty(bottom_interior_IDs.size * 3 + 6)
            bottom_op_mat_data_add = np.empty(bottom_interior_IDs.size * 3 + 6)
            bottom_op_mat_row_add[: (bottom_interior_IDs.size * 3)] = np.repeat(
                bottom_interior_IDs, 3
            )
            # ^...put the values in the same places in the operating matrix...
            bottom_op_mat_col_add[: (bottom_interior_IDs.size * 3)] = (
                self._realIDtointerior(
                    self._operating_matrix_ID_map[self._top_interior_IDs, :][
                        :, self._top_mask[3:6]
                    ]
                ).flatten()
            )
            bottom_op_mat_data_add[: (bottom_interior_IDs.size * 3)] = (
                _delta_t
                * (nine_node_map[_bottom_list, :][:, bottom_antimask]).flatten()
            )
            # ^...but the values refer to the TOP of the grid
            top_op_mat_row_add = np.empty(top_interior_IDs.size * 3 + 6)
            top_op_mat_col_add = np.empty(top_interior_IDs.size * 3 + 6)
            top_op_mat_data_add = np.empty(top_interior_IDs.size * 3 + 6)
            top_op_mat_row_add[: (top_interior_IDs.size * 3)] = np.repeat(
                top_interior_IDs, 3
            )
            top_op_mat_col_add[: (top_interior_IDs.size * 3)] = self._realIDtointerior(
                self._operating_matrix_ID_map[self._bottom_interior_IDs, :][
                    :, self._bottom_mask[0:3]
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
            bottom_op_mat_col_add[-6:-2] = self._operating_matrix_corner_int_IDs[
                top_corner_coords.reshape((2, 1)), top_corner_coords
            ].flatten()
            bottom_op_mat_row_add[-2:] = corner_interior_IDs[bottom_corner_coords]
            bottom_op_mat_col_add[-2:] = self._operating_matrix_corner_int_IDs[
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
            top_op_mat_col_add[-6:-2] = self._operating_matrix_corner_int_IDs[
                bottom_corner_coords.reshape((2, 1)), bottom_corner_coords
            ].flatten()
            top_op_mat_row_add[-2:] = corner_interior_IDs[top_corner_coords]
            top_op_mat_col_add[-2:] = self._operating_matrix_corner_int_IDs[
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

        if self._top_flag == 1:
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
        elif self._top_flag == 4 or self._top_flag == 2:
            top_op_mat_row_add = np.empty(top_interior_IDs.size * 3 + 6)
            top_op_mat_col_add = np.empty(top_interior_IDs.size * 3 + 6)
            top_op_mat_data_add = np.empty(top_interior_IDs.size * 3 + 6)
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            top_op_mat_row_add[: (top_interior_IDs.size * 3)] = np.repeat(
                top_interior_IDs, 3
            )
            top_op_mat_col_add[: (top_interior_IDs.size * 3)] = self._realIDtointerior(
                self._operating_matrix_ID_map[self._top_interior_IDs, :][
                    :, self._top_mask[3:6]
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
            top_op_mat_col_add[-6:-2] = self._operating_matrix_corner_int_IDs[
                this_corner_coords.reshape((2, 1)), this_corner_coords
            ].flatten()
            top_op_mat_row_add[-2:] = corner_interior_IDs[this_corner_coords]
            top_op_mat_col_add[-2:] = self._operating_matrix_corner_int_IDs[
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
            if self._top_flag == 2:
                _mat_RHS[top_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_top_list, :][:, top_antimask]
                    * self._fixed_gradient_offset_map[
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
                        * self._fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self._top_flag == 3:
            pass  # dealt with above
        else:
            raise NameError(
                """Something is very wrong with your boundary
                            conditions...!"""
            )

        if self._left_flag == 1:
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
        elif self._left_flag == 4 or self._left_flag == 2:
            left_op_mat_row_add = np.empty(left_interior_IDs.size * 3 + 4)
            left_op_mat_col_add = np.empty(left_interior_IDs.size * 3 + 4)
            left_op_mat_data_add = np.empty(left_interior_IDs.size * 3 + 4)
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            left_op_mat_row_add[: (left_interior_IDs.size * 3)] = np.repeat(
                left_interior_IDs, 3
            )
            left_op_mat_col_add[: (left_interior_IDs.size * 3)] = (
                self._realIDtointerior(
                    self._operating_matrix_ID_map[self._left_interior_IDs, :][
                        :, self._left_mask[::2]
                    ]
                ).flatten()
            )
            left_op_mat_data_add[: (left_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_left_list, :][:, left_antimask]).flatten()
            )
            # ...& the corners
            this_corner_coords = np.array([0, 2])
            left_op_mat_row_add[-4:] = np.repeat(
                corner_interior_IDs[this_corner_coords], 2
            )
            left_op_mat_col_add[-4:] = self._operating_matrix_corner_int_IDs[
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
            if self._left_flag == 2:
                _mat_RHS[left_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_left_list, :][:, left_antimask]
                    * self._fixed_gradient_offset_map[
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
                        * self._fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self._left_flag == 3:
            left_op_mat_row_add = np.empty(left_interior_IDs.size * 3 + 4)
            left_op_mat_col_add = np.empty(left_interior_IDs.size * 3 + 4)
            left_op_mat_data_add = np.empty(left_interior_IDs.size * 3 + 4)
            left_op_mat_row_add[: (left_interior_IDs.size * 3)] = np.repeat(
                left_interior_IDs, 3
            )
            left_op_mat_col_add[: (left_interior_IDs.size * 3)] = (
                self._realIDtointerior(
                    self._operating_matrix_ID_map[self._right_interior_IDs, :][
                        :, self._right_mask[1::2]
                    ]
                ).flatten()
            )
            left_op_mat_data_add[: (left_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_left_list, :][:, left_antimask]).flatten()
            )
            right_op_mat_row_add = np.empty(right_interior_IDs.size * 3 + 4)
            right_op_mat_col_add = np.empty(right_interior_IDs.size * 3 + 4)
            right_op_mat_data_add = np.empty(right_interior_IDs.size * 3 + 4)
            right_op_mat_row_add[: (right_interior_IDs.size * 3)] = np.repeat(
                right_interior_IDs, 3
            )
            right_op_mat_col_add[: (right_interior_IDs.size * 3)] = (
                self._realIDtointerior(
                    self._operating_matrix_ID_map[self._left_interior_IDs, :][
                        :, self._left_mask[::2]
                    ]
                ).flatten()
            )
            right_op_mat_data_add[: (right_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_right_list, :][:, right_antimask]).flatten()
            )
            # & the corners
            left_corner_coords = np.array([0, 2])
            right_corner_coords = np.array([1, 3])
            left_op_mat_row_add[-4:] = np.repeat(
                corner_interior_IDs[left_corner_coords], 2
            )
            left_op_mat_col_add[-4:] = self._operating_matrix_corner_int_IDs[
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
            right_op_mat_col_add[-4:] = self._operating_matrix_corner_int_IDs[
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

        if self._right_flag == 1:
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
        elif self._right_flag == 4 or self._right_flag == 2:
            right_op_mat_row_add = np.empty(right_interior_IDs.size * 3 + 4)
            right_op_mat_col_add = np.empty(right_interior_IDs.size * 3 + 4)
            right_op_mat_data_add = np.empty(right_interior_IDs.size * 3 + 4)
            # Equivalent to fixed gradient, but the gradient is zero, so
            # material only goes in the linked cell(i.e., each cell in the
            # op_mat edges points back to itself).
            right_op_mat_row_add[: (right_interior_IDs.size * 3)] = np.repeat(
                right_interior_IDs, 3
            )
            right_op_mat_col_add[: (right_interior_IDs.size * 3)] = (
                self._realIDtointerior(
                    self._operating_matrix_ID_map[self._right_interior_IDs, :][
                        :, self._right_mask[1::2]
                    ]
                ).flatten()
            )
            right_op_mat_data_add[: (right_interior_IDs.size * 3)] = (
                _delta_t * (nine_node_map[_right_list, :][:, right_antimask]).flatten()
            )
            # ...& the corners
            this_corner_coords = np.array([1, 3])
            right_op_mat_row_add[-4:] = np.repeat(
                corner_interior_IDs[this_corner_coords], 2
            )
            right_op_mat_col_add[-4:] = self._operating_matrix_corner_int_IDs[
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
            if self._right_flag == 2:
                _mat_RHS[right_interior_IDs] -= _delta_t * np.sum(
                    nine_node_map[_right_list, :][:, right_antimask]
                    * self._fixed_gradient_offset_map[
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
                        * self._fixed_gradient_offset_map[
                            _interior_corners[i]
                            + modulator_mask[corners_antimasks[i, edge_list]]
                        ]
                    )
        elif self._top_flag == 3:
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
        ncols = self._ncols
        interior_ID = (ID // ncols - 1) * (ncols - 2) + (ID % ncols) - 1
        if np.any(interior_ID < 0) or np.any(interior_ID >= self._ninteriornodes):
            raise NameError(
                """One of the supplied nodes was outside the
                            interior grid!"""
            )
        else:
            return interior_ID.astype(int)

    def _interiorIDtoreal(self, ID):
        IGW = self._interior_grid_width
        real_ID = (ID // IGW + 1) * self._ncols + (ID % IGW) + 1
        assert np.all(real_ID < self._nnodes)
        return real_ID.astype(int)

    def _realIDtocore(self, ID):
        ncols = self._ncols
        core_ID = (ID // ncols - 2) * (ncols - 4) + (ID % ncols) - 2
        if np.any(core_ID < 0) or np.any(core_ID >= self._ncorenodes):
            raise NameError(
                """One of the supplied nodes was outside the
                            core grid!"""
            )
        else:
            return core_ID.astype(int)

    def _coreIDtoreal(self, ID):
        CCW = self._core_cell_width
        real_ID = (ID // CCW + 2) * self._ncols + (ID % CCW) + 2
        assert np.all(real_ID < self._nnodes)
        return real_ID.astype(int)

    def _interiorIDtocore(self, ID):
        IGW = self._interior_grid_width
        core_ID = (ID // IGW - 1) * (self._ncols - 4) + (ID % IGW) - 1
        if np.any(core_ID < 0) or np.any(core_ID >= self._ncorenodes):
            raise NameError(
                """One of the supplied nodes was outside the
                            core grid!"""
            )
        else:
            return core_ID.astype(int)

    def _coreIDtointerior(self, ID):
        CCW = self._core_cell_width
        interior_ID = (ID // CCW + 1) * (self._ncols - 2) + (ID % CCW) + 1
        assert np.all(interior_ID < self._ninteriornodes)
        return interior_ID.astype(int)

    def run_one_step(self, dt):
        """Run the diffuser for one timestep, dt.

        This is the primary method of the class.

        Parameters
        ----------
        dt : float (time)
            The imposed timestep.
        """
        if self._bc_set_code != self._grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self._grid.bc_set_code
        else:
            self._gear_timestep(dt, self._grid)
            for _ in range(self._internal_repeats):
                # Initialize the variables for the step:
                self._set_variables(self._grid)
                # Solve interior of grid:
                _interior_elevs = linalg.spsolve(self._operating_matrix, self._mat_RHS)
                # this fn solves Ax=B for x

                # Handle the BC cells; test common cases first for speed
                self._grid["node"][self._values_to_diffuse][
                    self._interior_IDs_as_real
                ] = _interior_elevs

        # if BC==1 or BC==4, don't need to take any action; in both
        # cases the values are unchanged.
        if self._fixed_grad_BCs_present:
            self._grid["node"][self._values_to_diffuse][
                self._grid.fixed_gradient_node_properties["boundary_node_IDs"]
            ] = (
                self._grid["node"][self._values_to_diffuse][
                    self._grid.fixed_gradient_node_properties["anchor_node_IDs"]
                ]
                + self._grid.fixed_gradient_node_properties["values_to_add"]
            )
        if self._looped_BCs_present:
            self._grid["node"][self._values_to_diffuse][
                self._grid.looped_node_properties["boundary_node_IDs"]
            ] = self._grid["node"][self._values_to_diffuse][
                self._grid.looped_node_properties["linked_node_IDs"]
            ]
