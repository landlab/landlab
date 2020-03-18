#! /usr/env/python
"""Component that models 2D diffusion using an explicit finite-volume method.

Created July 2013 GT Last updated March 2016 DEJH with LL v1.0 component
style
"""


import numpy as np

from landlab import Component, FieldError, LinkStatus, NodeStatus, RasterModelGrid

_ALPHA = 0.15  # time-step stability factor
# ^0.25 not restrictive enough at meter scales w S~1 (possible cases)


class LinearDiffuser(Component):
    """This component implements linear diffusion of a Landlab field.

    Component assumes grid does not deform. If the boundary conditions on the
    grid change after component instantiation, be sure to also call
    :func:`updated_boundary_conditions` to ensure these are reflected in the
    component (especially if fixed_links are present).

    The method keyword allows control of the way the solver works. Options
    other than 'simple' will make the component run slower, but give second
    order solutions that incorporate information from more distal nodes (i.e.,
    the diagonals). Note that the option 'resolve_on_patches' can result in
    somewhat counterintuitive behaviour - this option tells the component to
    treat the diffusivity as a field **with directionality to it** (i.e., that
    the diffusivites are defined on links). Thus if all links have the same
    diffusivity value, with this flag active "effective" diffusivities
    at the nodes will be *higher* than this value (by a factor of root 2) as
    the diffusivity at each patch will be the mean vector sum of that at the
    bounding links.

    The primary method of this class is :func:`run_one_step`.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> mg = RasterModelGrid((9, 9))
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> z.reshape((9, 9))[4, 4] = 1.
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> ld = LinearDiffuser(mg, linear_diffusivity=1.)
    >>> for i in range(1):
    ...     ld.run_one_step(1.)
    >>> np.isclose(z[mg.core_nodes].sum(), 1.)
    True
    >>> mg2 = RasterModelGrid((5, 30))
    >>> z2 = mg2.add_zeros("topographic__elevation", at="node")
    >>> z2.reshape((5, 30))[2, 8] = 1.
    >>> z2.reshape((5, 30))[2, 22] = 1.
    >>> mg2.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> kd = mg2.node_x/mg2.node_x.mean()
    >>> ld2 = LinearDiffuser(mg2, linear_diffusivity=kd)
    >>> for i in range(10):
    ...     ld2.run_one_step(0.1)
    >>> z2[mg2.core_nodes].sum() == 2.
    True
    >>> z2.reshape((5, 30))[2, 8] > z2.reshape((5, 30))[2, 22]
    True

    An example using links:

    >>> mg1 = RasterModelGrid((10, 10), xy_spacing=100.)
    >>> mg2 = RasterModelGrid((10, 10), xy_spacing=100.)
    >>> z1 = mg1.add_zeros("topographic__elevation", at="node")
    >>> z2 = mg2.add_zeros("topographic__elevation", at="node")
    >>> dt = 1.
    >>> nt = 10
    >>> kappa_links = mg2.add_ones("surface_water__discharge", at="link")
    >>> kappa_links *= 10000.
    >>> dfn1 = LinearDiffuser(mg1, linear_diffusivity=10000.)
    >>> dfn2 = LinearDiffuser(mg2, linear_diffusivity='surface_water__discharge')
    >>> for i in range(nt):
    ...     z1[mg1.core_nodes] += 1.
    ...     z2[mg2.core_nodes] += 1.
    ...     dfn1.run_one_step(dt)
    ...     dfn2.run_one_step(dt)
    >>> np.allclose(z1, z2)
    True
    >>> z2.fill(0.)
    >>> dfn2 = LinearDiffuser(mg2, linear_diffusivity='surface_water__discharge',
    ...                       method='resolve_on_patches')
    >>> for i in range(nt):
    ...     z2[mg2.core_nodes] += 1.
    ...     dfn2.run_one_step(dt)
    >>> np.all(z2[mg2.core_nodes] < z1[mg2.core_nodes])
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Culling, W. (1963). Soil Creep and the Development of Hillside Slopes.
    The Journal of Geology  71(2), 127-161. https://dx.doi.org/10.1086/626891

    """

    _name = "LinearDiffuser"

    _unit_agnostic = True

    _info = {
        "hillslope_sediment__unit_volume_flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2/s",
            "mapping": "link",
            "doc": "Volume flux per unit width along links",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "link",
            "doc": "Gradient of the ground surface",
        },
    }

    def __init__(self, grid, linear_diffusivity=0.01, method="simple", deposit=True):
        """
        Parameters
        ----------
        grid : ModelGrid
            A grid.
        linear_diffusivity : float, array, or field name (m**2/time)
            The diffusivity. If an array or field name, these must be the
            diffusivities on either nodes or links - the component will
            distinguish which based on array length. Values on nodes will be
            mapped to links using an upwind scheme in the simple case.
        method : {'simple', 'resolve_on_patches', 'on_diagonals'}
            The method used to represent the fluxes. 'simple' solves a finite
            difference method with a simple staggered grid scheme onto the links.
            'resolve_on_patches' solves the scheme by mapping both slopes and
            diffusivities onto the patches and solving there before resolving
            values back to the nodes (and at the moment requires a raster grid).
            Note that this scheme enforces directionality in the diffusion field;
            it's no longer just a scalar field. Thus diffusivities must be defined
            *on links* when this option is chosen.
            'on_diagonals' permits Raster diagonals to carry the diffusional
            fluxes. These latter techniques are more computationally expensive,
            but can suppress cardinal direction artifacts if diffusion is
            performed on a raster. 'on_diagonals' pretends that the "faces" of a
            cell with 8 links are represented by a stretched regular octagon set
            within the true cell.
        deposit : {True, False}
            Whether diffusive material can be deposited. True means that diffusive
            material will be deposited if the divergence of sediment flux is
            negative. False means that even when the divergence of sediment flux is
            negative, no material is deposited. (No deposition ever.) The False
            case is a bit of a band-aid to account for cases when fluvial incision
            likely removes any material that would be deposited. If one couples
            fluvial detachment-limited incision with linear diffusion, the channels
            will not reach the predicted analytical solution unless deposit is set
            to False.
        """
        super().__init__(grid)

        self._bc_set_code = self._grid.bc_set_code
        assert method in ("simple", "resolve_on_patches", "on_diagonals")
        if method == "resolve_on_patches":
            assert isinstance(self._grid, RasterModelGrid)
            self._use_patches = True
        else:
            self._use_patches = False
        if method == "on_diagonals" and isinstance(self._grid, RasterModelGrid):
            self._use_diags = True
        else:
            self._use_diags = False
        self._current_time = 0.0
        self._run_before = False
        self._kd_on_links = False

        if isinstance(linear_diffusivity, str):
            try:
                self._kd = self._grid.at_link[linear_diffusivity]
                self._kd_on_links = True
            except KeyError:
                self._kd = self._grid.at_node[linear_diffusivity]
        else:
            self._kd = linear_diffusivity

            if isinstance(self._kd, (float, int)):
                self._kd = float(self._kd)
            else:
                if self._kd.size == self._grid.number_of_links:
                    self._kd_on_links = True
                else:
                    assert self._kd.size == self._grid.number_of_nodes

        if self._kd_on_links is True:
            assert isinstance(self._grid, RasterModelGrid)

        # if we're using patches, it is VITAL that diffusivity is defined on
        # links. The whole point of this functionality is that we honour
        # *directionality* in the diffusivities.
        if self._use_patches:
            assert self._kd_on_links
        # set _deposit flag to tell code whether or not diffusion can deposit.

        self._deposit = deposit

        self._values_to_diffuse = "topographic__elevation"

        # Set internal time step
        # ..todo:
        #   implement mechanism to compute time-steps dynamically if grid is
        #   adaptive/changing
        # as of modern componentization (Spring '16), this can take arrays
        # and irregular grids
        # CFL condition precalc:
        CFL_prefactor = (
            _ALPHA * self._grid.length_of_link[: self._grid.number_of_links] ** 2.0
        )
        # ^ link_length can include diags, if not careful...
        self._CFL_actives_prefactor = CFL_prefactor[self._grid.active_links]
        # ^note we can do this as topology shouldn't be changing

        # Get a list of interior cells
        self._interior_cells = self._grid.node_at_core_cell

        self._z = self._grid.at_node[self._values_to_diffuse]
        self._dqsds = self._grid.zeros("node", dtype=float)
        if not self._use_diags:
            g = self._grid.zeros(at="link")
            qs = self._grid.zeros(at="link")
            try:
                self._g = self._grid.add_field(
                    "topographic__gradient", g, at="link", clobber=False
                )
                # ^note this will object if this exists already
            except FieldError:  # keep a ref
                self._g = self._grid.at_link["topographic__gradient"]
            try:
                self._qs = self._grid.add_field(
                    "hillslope_sediment__unit_volume_flux", qs, at="link", clobber=False
                )
            except FieldError:
                self._qs = self._grid.at_link["hillslope_sediment__unit_volume_flux"]
            # note all these terms are deliberately loose, as we won't always
            # be dealing with topo
        else:
            g = np.zeros(self._grid.number_of_d8, dtype=float)
            qs = np.zeros(self._grid.number_of_d8, dtype=float)
            self._g = g
            self._qs = qs
            # now we have to choose what the face width of a diagonal is...
            # Adopt a regular octagon config if it's a square raster, and
            # stretch this geometry as needed.
            # Conceptually, this means we're passing mass across diamond-
            # shaped holes centered at the corners.
            # Note that this WON'T affect the inferred cell size - that's
            # still derived from the rectangle.
            self._d8width_face_at_link = np.empty(self._grid.number_of_d8)
            # note there will be null entries here
            # by our defs, every active link must have a face.
            # calc the length of a diag "face":
            rt2 = np.sqrt(2.0)
            horizontal_face = self._grid.dx / (1.0 + rt2)
            vertical_face = self._grid.dy / (1.0 + rt2)
            diag_face = np.sqrt(0.5 * (horizontal_face ** 2 + vertical_face ** 2))

            # NOTE: Do these need to be flattened?
            # self._hoz = self.grid.horizontal_links.flatten()
            # self._vert = self.grid.vertical_links.flatten()
            self._hoz = self.grid.horizontal_links
            self._vert = self.grid.vertical_links

            self._d8width_face_at_link[self._hoz] = vertical_face
            self._d8width_face_at_link[self._vert] = horizontal_face
            # ^ this operation pastes in faces where there are none, but
            # we'll never use them
            self._d8width_face_at_link[self._grid.number_of_links :] = diag_face

        self._vertlinkcomp = np.sin(self._grid.angle_of_link)
        self._hozlinkcomp = np.cos(self._grid.angle_of_link)

        if self._use_patches or self._kd_on_links:
            mg = self._grid
            try:
                self._hoz = self.grid.horizontal_links
                self._vert = self.grid.vertical_links
            except AttributeError:
                pass
            self._x_link_patches = mg.patches_at_link[self._hoz]
            x_link_patch_pres = mg.patches_present_at_link[self._hoz]
            self._x_link_patch_mask = np.logical_not(x_link_patch_pres)
            self._y_link_patches = mg.patches_at_link[self._vert]
            y_link_patch_pres = mg.patches_present_at_link[self._vert]
            self._y_link_patch_mask = np.logical_not(y_link_patch_pres)
            self._hoz_link_neighbors = np.empty((self._hoz.size, 4), dtype=int)
            self._vert_link_neighbors = np.empty((self._vert.size, 4), dtype=int)

        # do some pre-work to make fixed grad BC updating faster in the loop:
        self.updated_boundary_conditions()

    @property
    def fixed_grad_nodes(self):
        """Fixed gradient nodes."""
        return self._fixed_grad_nodes

    @property
    def fixed_grad_anchors(self):
        """Fixed gradient anchors."""
        return self._fixed_grad_anchors

    @property
    def fixed_grad_offsets(self):
        """Fixed gradient offsets."""
        return self._fixed_grad_offsets

    def updated_boundary_conditions(self):
        """Call if grid BCs are updated after component instantiation.

        Sets `fixed_grad_nodes`, `fixed_grad_anchors`, & `fixed_grad_offsets`,
        such that::

            value[fixed_grad_nodes] = value[fixed_grad_anchors] + offset

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> import numpy as np
        >>> mg = RasterModelGrid((4, 5))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> z[mg.core_nodes] = 1.
        >>> ld = LinearDiffuser(mg, linear_diffusivity=1.)
        >>> ld.fixed_grad_nodes.size == 0
        True
        >>> ld.fixed_grad_anchors.size == 0
        True
        >>> ld.fixed_grad_offsets.size == 0
        True
        >>> mg.at_link['topographic__slope'] = mg.calc_grad_at_link(
        ...     'topographic__elevation')
        >>> mg.status_at_node[mg.perimeter_nodes] = mg.BC_NODE_IS_FIXED_GRADIENT
        >>> ld.updated_boundary_conditions()
        >>> ld.fixed_grad_nodes
        array([ 1,  2,  3,  5,  9, 10, 14, 16, 17, 18])
        >>> ld.fixed_grad_anchors
        array([ 6,  7,  8,  6,  8, 11, 13, 11, 12, 13])
        >>> ld.fixed_grad_offsets
        array([-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.])
        >>> np.allclose(z[ld.fixed_grad_nodes],
        ...             z[ld.fixed_grad_anchors] + ld.fixed_grad_offsets)
        True
        """
        fixed_grad_nodes = np.where(
            self._grid.status_at_node == NodeStatus.FIXED_GRADIENT
        )[0]
        heads = self._grid.node_at_link_head[self._grid.fixed_links]
        tails = self._grid.node_at_link_tail[self._grid.fixed_links]
        head_is_fixed = np.in1d(heads, fixed_grad_nodes)
        self._fixed_grad_nodes = np.where(head_is_fixed, heads, tails)
        self._fixed_grad_anchors = np.where(head_is_fixed, tails, heads)
        vals = self._grid.at_node[self._values_to_diffuse]
        self._fixed_grad_offsets = (
            vals[self._fixed_grad_nodes] - vals[self._fixed_grad_anchors]
        )
        if self._use_diags:
            self._g.fill(0.0)

        if self._kd_on_links or self._use_patches:
            mg = self._grid
            x_link_patch_pres = mg.patches_present_at_link[self._hoz]
            self._x_link_patch_mask = np.logical_not(x_link_patch_pres)
            y_link_patch_pres = mg.patches_present_at_link[self._vert]
            self._y_link_patch_mask = np.logical_not(y_link_patch_pres)
            self._hoz_link_neighbors[:, :2] = mg.links_at_node[
                mg.node_at_link_head[self._hoz], 1:4:2
            ]
            self._hoz_link_neighbors[:, 2:] = mg.links_at_node[
                mg.node_at_link_tail[self._hoz], 1:4:2
            ]
            self._vert_link_neighbors[:, :2] = mg.links_at_node[
                mg.node_at_link_head[self._vert], 0:3:2
            ]
            self._vert_link_neighbors[:, 2:] = mg.links_at_node[
                mg.node_at_link_tail[self._vert], 0:3:2
            ]
            self._vert_link_badlinks = np.logical_or(
                mg.status_at_link[self._vert_link_neighbors] == LinkStatus.INACTIVE,
                self._vert_link_neighbors == -1,
            )
            self._hoz_link_badlinks = np.logical_or(
                mg.status_at_link[self._hoz_link_neighbors] == LinkStatus.INACTIVE,
                self._hoz_link_neighbors == -1,
            )

    def run_one_step(self, dt):
        """Run the diffuser for one timestep, dt.

        If the imposed timestep dt is longer than the Courant-Friedrichs-Lewy
        condition for the diffusion, this timestep will be internally divided
        as the component runs, as needed.

        Parameters
        ----------
        dt : float (time)
            The imposed timestep.
        """
        mg = self._grid
        z = self._grid.at_node[self._values_to_diffuse]

        if not self._run_before:
            self.updated_boundary_conditions()  # just in case
            self._run_before = True
        if self._bc_set_code != self._grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self._grid.bc_set_code

        core_nodes = self._grid.node_at_core_cell
        # do mapping of array kd here, in case it points at an updating
        # field:
        if isinstance(self._kd, np.ndarray):
            if not self._kd_on_links:
                kd_links = self._grid.map_max_of_link_nodes_to_link(self._kd)
                kd_activelinks = kd_links[self._grid.active_links]
                # re-derive CFL condition, as could change dynamically:
                dt_links = self._CFL_actives_prefactor / kd_activelinks
                self._dt = np.nanmin(dt_links)
            else:
                kd_links = self._kd
                kd_activelinks = self._kd[self._grid.active_links]
                dt_links = self._CFL_actives_prefactor / kd_activelinks
                self._dt_links = dt_links
                self._dt = np.nanmin(np.fabs(dt_links))
        else:
            kd_activelinks = self._kd
            # re-derive CFL condition, as could change dynamically:
            dt_links = self._CFL_actives_prefactor / kd_activelinks
            self._dt = np.nanmin(dt_links)

        if self._use_patches:
            # need this else diffusivities on inactive links deform off-angle
            # calculations
            kd_links = kd_links.copy()
            kd_links[self._grid.status_at_link == LinkStatus.INACTIVE] = 0.0

        # Take the smaller of delt or built-in time-step size self._dt
        self._tstep_ratio = dt / self._dt
        repeats = int(self._tstep_ratio // 1.0)
        extra_time = self._tstep_ratio - repeats

        # Can really get into trouble if no diffusivity happens but we run...
        if self._dt < np.inf:
            loops = repeats + 1
        else:
            loops = 0
        for i in range(loops):
            if not self._use_diags:
                grads = mg.calc_grad_at_link(z)
                self._g[mg.active_links] = grads[mg.active_links]
                if not self._use_patches:  # currently forbidden
                    # if diffusivity is an array, self._kd is already
                    # active_links-long
                    self._qs[mg.active_links] = (
                        -kd_activelinks * self._g[mg.active_links]
                    )
                    # Calculate the net deposition/erosion rate at each node
                    mg.calc_flux_div_at_node(self._qs, out=self._dqsds)
                else:  # project onto patches
                    slx = mg.zeros("link")
                    sly = mg.zeros("link")
                    slx[self._hoz] = self._g[self._hoz]
                    sly[self._vert] = self._g[self._vert]
                    patch_dx, patch_dy = mg.calc_grad_at_patch(z)
                    xvecs_vert = np.ma.array(
                        patch_dx[self._y_link_patches], mask=self._y_link_patch_mask
                    )
                    slx[self._vert] = xvecs_vert.mean()
                    yvecs_hoz = np.ma.array(
                        patch_dy[self._x_link_patches], mask=self._x_link_patch_mask
                    )
                    sly[self._hoz] = yvecs_hoz.mean()
                    # now map diffusivities (already on links, but we want
                    # more spatial averaging)
                    Kx = mg.zeros("link")
                    Ky = mg.zeros("link")
                    Kx[self._hoz] = kd_links[self._hoz]
                    Ky[self._vert] = kd_links[self._vert]
                    vert_link_crosslink_K = np.ma.array(
                        kd_links[self._vert_link_neighbors],
                        mask=self._vert_link_badlinks,
                    )
                    hoz_link_crosslink_K = np.ma.array(
                        kd_links[self._hoz_link_neighbors], mask=self._hoz_link_badlinks
                    )
                    Kx[self._vert] = vert_link_crosslink_K.mean(axis=1)
                    Ky[self._hoz] = hoz_link_crosslink_K.mean(axis=1)
                    Cslope = np.sqrt(slx ** 2 + sly ** 2)
                    v = np.sqrt(Kx ** 2 + Ky ** 2)
                    flux_links = v * Cslope
                    # NEW, to resolve issue with K being off angle to S:
                    # in fact, no. Doing this just makes this equivalent
                    # to the basic diffuser, but with a bunch more crap
                    # involved.
                    # flux_x = slx * Kx
                    # flux_y = sly * Ky
                    # flux_links = np.sqrt(flux_x*flux_x + flux_y*flux_y)
                    theta = np.arctan(np.fabs(sly) / (np.fabs(slx) + 1.0e-10))
                    flux_links[self._hoz] *= np.sign(slx[self._hoz]) * np.cos(
                        theta[self._hoz]
                    )
                    flux_links[self._vert] *= np.sign(sly[self._vert]) * np.sin(
                        theta[self._vert]
                    )
                    # zero out the inactive links
                    self._qs[mg.active_links] = -flux_links[mg.active_links]

                    self._grid.calc_flux_div_at_node(self._qs, out=self._dqsds)

            else:  # ..._use_diags
                # NB: this is dirty code. It uses the obsolete diagonal data
                # structures, and necessarily has to do a bunch of mapping
                # on the fly.

                # remap the kds onto the links, as necessary
                if isinstance(self._kd, np.ndarray):
                    d8link_kd = np.empty(self._grid.number_of_d8, dtype=float)
                    d8link_kd[self._grid.active_links] = kd_activelinks
                    d8link_kd[self._grid.active_diagonals] = np.amax(
                        self._kd[
                            self._grid.nodes_at_diagonal[self._grid.active_diagonals]
                        ],
                        axis=1,
                    ).flatten()
                else:
                    d8link_kd = self._kd
                self._g[self._grid.active_links] = self._grid.calc_grad_at_link(z)[
                    self._grid.active_links
                ]
                self._g[self._grid.active_diagonals] = (
                    z[self._grid._diag_activelink_tonode]
                    - z[self._grid._diag_activelink_fromnode]
                ) / self._grid.length_of_d8[self._grid.active_diagonals]
                self._qs[:] = -d8link_kd * self._g

                total_flux = self._qs * self._d8width_face_at_link  # nlinks
                totalflux_allnodes = (
                    total_flux[self._grid.links_at_node]
                    * self._grid.active_link_dirs_at_node
                ).sum(axis=1)
                totalflux_allnodes += (
                    total_flux[self._grid.d8s_at_node[:, 4:]]
                    * self._grid.active_diagonal_dirs_at_node
                ).sum(axis=1)
                self._dqsds[self._grid.node_at_cell] = (
                    -totalflux_allnodes[self._grid.node_at_cell]
                    / self._grid.area_of_cell
                )

            # Calculate the total rate of elevation change
            dzdt = -self._dqsds
            if not self._deposit:
                dzdt[np.where(dzdt > 0)] = 0.0
            # Update the elevations
            timestep = self._dt
            if i == (repeats):
                timestep *= extra_time
            else:
                pass
            self._grid.at_node[self._values_to_diffuse][core_nodes] += (
                dzdt[core_nodes] * timestep
            )

            # check the BCs, update if fixed gradient
            vals = self._grid.at_node[self._values_to_diffuse]
            vals[self._fixed_grad_nodes] = (
                vals[self._fixed_grad_anchors] + self._fixed_grad_offsets
            )

    @property
    def time_step(self):
        """Returns internal time-step size (as a property)."""
        return self._dt
