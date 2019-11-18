#! /usr/env/python
"""

Component that models 2D diffusion using an explicit finite-volume method.

Created July 2013 GT
Last updated March 2016 DEJH with LL v1.0 component style
"""

from __future__ import print_function

import numpy as np
from six.moves import range

from landlab import (
    FIXED_GRADIENT_BOUNDARY,
    INACTIVE_LINK,
    Component,
    FieldError,
    RasterModelGrid,
)
from landlab.utils.decorators import use_file_name_or_kwds

_ALPHA = 0.15  # time-step stability factor
# ^0.25 not restrictive enough at meter scales w S~1 (possible cases)


class LinearDiffuser(Component):
    """
    This component implements linear diffusion of a Landlab field.

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
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> z.reshape((9, 9))[4, 4] = 1.
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> ld = LinearDiffuser(mg, linear_diffusivity=1.)
    >>> for i in range(1):
    ...     ld.run_one_step(1.)
    >>> np.isclose(z[mg.core_nodes].sum(), 1.)
    True
    >>> mg2 = RasterModelGrid((5, 30))
    >>> z2 = mg2.add_zeros('node', 'topographic__elevation')
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
    >>> z1 = mg1.add_zeros('node', 'topographic__elevation')
    >>> z2 = mg2.add_zeros('node', 'topographic__elevation')
    >>> dt = 1.
    >>> nt = 10
    >>> kappa_links = mg2.add_ones('link', 'surface_water__discharge')
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
    """

    _name = "LinearDiffuser"

    _input_var_names = ("topographic__elevation",)

    _output_var_names = (
        "topographic__elevation",
        "topographic__gradient",
        "hillslope_sediment__unit_volume_flux",
    )

    _var_units = {
        "topographic__elevation": "m",
        "topographic__gradient": "-",
        "hillslope_sediment__unit_volume_flux": "m**2/s",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "topographic__gradient": "link",
        "hillslope_sediment__unit_volume_flux": "link",
    }

    _var_doc = {
        "topographic__elevation": (
            "Land surface topographic elevation; can "
            + "be overwritten in initialization"
        ),
        "topographic__gradient": "Gradient of surface, on links",
        "hillslope_sediment__unit_volume_flux": "Volume flux per unit width along links",
    }

    @use_file_name_or_kwds
    def __init__(
        self, grid, linear_diffusivity=None, method="simple", deposit=True, **kwds
    ):
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
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code
        assert method in ("simple", "resolve_on_patches", "on_diagonals")
        if method == "resolve_on_patches":
            assert isinstance(self.grid, RasterModelGrid)
            self._use_patches = True
        else:
            self._use_patches = False
        if method == "on_diagonals" and isinstance(self.grid, RasterModelGrid):
            self._use_diags = True
        else:
            self._use_diags = False
        self.current_time = 0.0
        self._run_before = False
        self._kd_on_links = False
        if linear_diffusivity is not None:
            if type(linear_diffusivity) is not str:
                self._kd = linear_diffusivity
                if type(self._kd) in (float, int):
                    if type(self._kd) is int:
                        self._kd = float(self._kd)
                else:
                    if self._kd.size == self.grid.number_of_links:
                        self._kd_on_links = True
                    else:
                        assert self._kd.size == self.grid.number_of_nodes
            else:
                try:
                    self._kd = self.grid.at_link[linear_diffusivity]
                    self._kd_on_links = True
                except KeyError:
                    self._kd = self.grid.at_node[linear_diffusivity]
        else:
            raise KeyError(
                "linear_diffusivity must be provided to the "
                + "LinearDiffuser component"
            )
        if self._kd_on_links is True:
            assert isinstance(self.grid, RasterModelGrid)

        # if we're using patches, it is VITAL that diffusivity is defined on
        # links. The whole point of this functionality is that we honour
        # *directionality* in the diffusivities.
        if self._use_patches:
            assert self._kd_on_links
        # set _deposit flag to tell code whether or not diffusion can deposit.

        self._deposit = deposit

        # for component back compatibility (undocumented):
        # note component can NO LONGER do internal uplift, at all.
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
        else:
            self.values_to_diffuse = "topographic__elevation"
        # Raise an error if somehow someone is using this weird functionality
        if self._grid is None:
            raise ValueError("You must now provide an existing grid!")
        # ###

        # Set internal time step
        # ..todo:
        #   implement mechanism to compute time-steps dynamically if grid is
        #   adaptive/changing
        # as of modern componentization (Spring '16), this can take arrays
        # and irregular grids
        # CFL condition precalc:
        CFL_prefactor = (
            _ALPHA * self.grid.length_of_link[: self.grid.number_of_links] ** 2.0
        )
        # ^ link_length can include diags, if not careful...
        self._CFL_actives_prefactor = CFL_prefactor[self.grid.active_links]
        # ^note we can do this as topology shouldn't be changing

        # Get a list of interior cells
        self.interior_cells = self.grid.node_at_core_cell

        self.z = self.grid.at_node[self.values_to_diffuse]
        self.dqsds = self.grid.zeros("node", dtype=float)
        if not self._use_diags:
            g = self.grid.zeros(at="link")
            qs = self.grid.zeros(at="link")
            try:
                self.g = self.grid.add_field(
                    "link", "topographic__gradient", g, noclobber=True
                )
                # ^note this will object if this exists already
            except FieldError:  # keep a ref
                self.g = self.grid.at_link["topographic__gradient"]
            try:
                self.qs = self.grid.add_field(
                    "link", "hillslope_sediment__unit_volume_flux", qs, noclobber=True
                )
            except FieldError:
                self.qs = self.grid.at_link["hillslope_sediment__unit_volume_flux"]
            # note all these terms are deliberately loose, as we won't always
            # be dealing with topo
        else:
            g = np.zeros(self.grid.number_of_d8, dtype=float)
            qs = np.zeros(self.grid.number_of_d8, dtype=float)
            self.g = g
            self.qs = qs
            # now we have to choose what the face width of a diagonal is...
            # Adopt a regular octagon config if it's a square raster, and
            # stretch this geometry as needed.
            # Conceptually, this means we're passing mass across diamond-
            # shaped holes centered at the corners.
            # Note that this WON'T affect the inferred cell size - that's
            # still derived from the rectangle.
            self._d8width_face_at_link = np.empty(self.grid.number_of_d8)
            # note there will be null entries here
            # by our defs, every active link must have a face.
            # calc the length of a diag "face":
            rt2 = np.sqrt(2.0)
            horizontal_face = self.grid.dx / (1.0 + rt2)
            vertical_face = self.grid.dy / (1.0 + rt2)
            diag_face = np.sqrt(0.5 * (horizontal_face ** 2 + vertical_face ** 2))
            self._hoz = self.grid.horizontal_links.flatten()
            self._vert = self.grid.vertical_links.flatten()
            self._d8width_face_at_link[self._hoz] = vertical_face
            self._d8width_face_at_link[self._vert] = horizontal_face
            # ^ this operation pastes in faces where there are none, but
            # we'll never use them
            self._d8width_face_at_link[self.grid.number_of_links :] = diag_face

        self._vertlinkcomp = np.sin(self.grid.angle_of_link)
        self._hozlinkcomp = np.cos(self.grid.angle_of_link)

        if self._use_patches or self._kd_on_links:
            mg = self.grid
            try:
                self._hoz = self.grid.horizontal_links.flatten()
                self._vert = self.grid.vertical_links.flatten()
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
        >>> z = mg.add_zeros('node', 'topographic__elevation')
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
        >>> mg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
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
            self.grid.status_at_node == FIXED_GRADIENT_BOUNDARY
        )[0]
        heads = self.grid.node_at_link_head[self.grid.fixed_links]
        tails = self.grid.node_at_link_tail[self.grid.fixed_links]
        head_is_fixed = np.in1d(heads, fixed_grad_nodes)
        self.fixed_grad_nodes = np.where(head_is_fixed, heads, tails)
        self.fixed_grad_anchors = np.where(head_is_fixed, tails, heads)
        vals = self.grid.at_node[self.values_to_diffuse]
        self.fixed_grad_offsets = (
            vals[self.fixed_grad_nodes] - vals[self.fixed_grad_anchors]
        )
        if self._use_diags:
            self.g.fill(0.0)

        if self._kd_on_links or self._use_patches:
            mg = self.grid
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
                mg.status_at_link[self._vert_link_neighbors] == INACTIVE_LINK,
                self._vert_link_neighbors == -1,
            )
            self._hoz_link_badlinks = np.logical_or(
                mg.status_at_link[self._hoz_link_neighbors] == INACTIVE_LINK,
                self._hoz_link_neighbors == -1,
            )

    def diffuse(self, dt, **kwds):
        """
        See :func:`run_one_step`.
        """
        if "internal_uplift" in kwds.keys():
            raise KeyError(
                "LinearDiffuser can no longer work with internal " + "uplift"
            )
        mg = self.grid
        z = self.grid.at_node[self.values_to_diffuse]

        if not self._run_before:
            self.updated_boundary_conditions()  # just in case
            self._run_before = True
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code

        core_nodes = self.grid.node_at_core_cell
        # do mapping of array kd here, in case it points at an updating
        # field:
        if type(self._kd) is np.ndarray:
            if not self._kd_on_links:
                kd_links = self.grid.map_max_of_link_nodes_to_link(self._kd)
                kd_activelinks = kd_links[self.grid.active_links]
                # re-derive CFL condition, as could change dynamically:
                dt_links = self._CFL_actives_prefactor / kd_activelinks
                self.dt = np.nanmin(dt_links)
            else:
                kd_links = self._kd
                kd_activelinks = self._kd[self.grid.active_links]
                dt_links = self._CFL_actives_prefactor / kd_activelinks
                self.dt_links = dt_links
                self.dt = np.nanmin(np.fabs(dt_links))
        else:
            kd_activelinks = self._kd
            # re-derive CFL condition, as could change dynamically:
            dt_links = self._CFL_actives_prefactor / kd_activelinks
            self.dt = np.nanmin(dt_links)

        if self._use_patches:
            # need this else diffusivities on inactive links deform off-angle
            # calculations
            kd_links = kd_links.copy()
            kd_links[self.grid.status_at_link == INACTIVE_LINK] = 0.0

        # Take the smaller of delt or built-in time-step size self.dt
        self.tstep_ratio = dt / self.dt
        repeats = int(self.tstep_ratio // 1.0)
        extra_time = self.tstep_ratio - repeats

        # Can really get into trouble if no diffusivity happens but we run...
        if self.dt < np.inf:
            loops = repeats + 1
        else:
            loops = 0
        for i in range(loops):
            if not self._use_diags:
                grads = mg.calc_grad_at_link(z)
                self.g[mg.active_links] = grads[mg.active_links]
                if not self._use_patches:  # currently forbidden
                    # if diffusivity is an array, self._kd is already
                    # active_links-long
                    self.qs[mg.active_links] = -kd_activelinks * self.g[mg.active_links]
                    # Calculate the net deposition/erosion rate at each node
                    mg.calc_flux_div_at_node(self.qs, out=self.dqsds)
                else:  # project onto patches
                    slx = mg.zeros("link")
                    sly = mg.zeros("link")
                    slx[self._hoz] = self.g[self._hoz]
                    sly[self._vert] = self.g[self._vert]
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
                    self.qs[mg.active_links] = -flux_links[mg.active_links]

                    self.grid.calc_flux_div_at_node(self.qs, out=self.dqsds)

            else:  # ..._use_diags
                # NB: this is dirty code. It uses the obsolete diagonal data
                # structures, and necessarily has to do a bunch of mapping
                # on the fly.

                # remap the kds onto the links, as necessary
                if isinstance(self._kd, np.ndarray):
                    d8link_kd = np.empty(self.grid.number_of_d8, dtype=float)
                    d8link_kd[self.grid.active_links] = kd_activelinks
                    d8link_kd[self.grid.active_diagonals] = np.amax(
                        self._kd[
                            self.grid.nodes_at_diagonal[self.grid.active_diagonals]
                        ],
                        axis=1,
                    ).flatten()
                else:
                    d8link_kd = self._kd
                self.g[self.grid.active_links] = self.grid.calc_grad_at_link(z)[
                    self.grid.active_links
                ]
                self.g[self.grid.active_diagonals] = (
                    z[self.grid._diag_activelink_tonode]
                    - z[self.grid._diag_activelink_fromnode]
                ) / self.grid.length_of_d8[self.grid.active_diagonals]
                self.qs[:] = -d8link_kd * self.g

                total_flux = self.qs * self._d8width_face_at_link  # nlinks
                totalflux_allnodes = (
                    total_flux[self.grid.links_at_node]
                    * self.grid.active_link_dirs_at_node
                ).sum(axis=1)
                totalflux_allnodes += (
                    total_flux[self.grid.d8s_at_node[:, 4:]]
                    * self.grid.active_diagonal_dirs_at_node
                ).sum(axis=1)
                self.dqsds[self.grid.node_at_cell] = (
                    -totalflux_allnodes[self.grid.node_at_cell] / self.grid.area_of_cell
                )

            # Calculate the total rate of elevation change
            dzdt = -self.dqsds
            if not self._deposit:
                dzdt[np.where(dzdt > 0)] = 0.0
            # Update the elevations
            timestep = self.dt
            if i == (repeats):
                timestep *= extra_time
            else:
                pass
            self.grid.at_node[self.values_to_diffuse][core_nodes] += (
                dzdt[core_nodes] * timestep
            )

            # check the BCs, update if fixed gradient
            vals = self.grid.at_node[self.values_to_diffuse]
            vals[self.fixed_grad_nodes] = (
                vals[self.fixed_grad_anchors] + self.fixed_grad_offsets
            )

        return self.grid

    def run_one_step(self, dt, **kwds):
        """Run the diffuser for one timestep, dt.

        If the imposed timestep dt is longer than the Courant-Friedrichs-Lewy
        condition for the diffusion, this timestep will be internally divided
        as the component runs, as needed.

        Parameters
        ----------
        dt : float (time)
            The imposed timestep.
        """
        self.diffuse(dt, **kwds)

    @property
    def time_step(self):
        """Returns internal time-step size (as a property).
        """
        return self.dt
