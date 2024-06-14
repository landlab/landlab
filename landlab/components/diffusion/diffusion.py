#! /usr/env/python
"""Component that models 2D diffusion using an explicit finite-volume method.

Created July 2013 GT Last updated March 2016 DEJH with LL v1.0 component
style
"""


import numpy as np

from landlab import Component
from landlab import LinkStatus
from landlab import NodeStatus
from landlab import RasterModelGrid

_ALPHA = 0.15  # time-step stability factor
# ^0.25 not restrictive enough at meter scales w S~1 (possible cases)


class LinearDiffuser(Component):
    """Linear diffusion of a Landlab field.

    Component assumes grid does not deform. If the boundary conditions on the
    grid change after component instantiation, be sure to also call
    :func:`updated_boundary_conditions` to ensure these are reflected in the
    component (especially if fixed_links are present).

    The ``method`` keyword allows control of the way the solver works.
    Note that the option 'resolve_on_patches' can result in
    somewhat counterintuitive behavior - this option tells the component to
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

    >>> grid = RasterModelGrid((9, 9))
    >>> z = grid.add_zeros("topographic__elevation", at="node")
    >>> z.reshape((9, 9))[4, 4] = 1.0
    >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> ld = LinearDiffuser(grid, linear_diffusivity=1.0)
    >>> ld.run_one_step(1.0)
    >>> np.isclose(z[grid.core_nodes].sum(), 1.0)
    True

    >>> grid = RasterModelGrid((5, 30))
    >>> z2 = grid.add_zeros("topographic__elevation", at="node")
    >>> z2.reshape((5, 30))[2, 8] = 1.0
    >>> z2.reshape((5, 30))[2, 22] = 1.0
    >>> grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> kd = grid.node_x / grid.node_x.mean()
    >>> ld2 = LinearDiffuser(grid, linear_diffusivity=kd)
    >>> for _ in range(10):
    ...     ld2.run_one_step(0.1)
    ...
    >>> z2[grid.core_nodes].sum() == 2.0
    True
    >>> z2.reshape((5, 30))[2, 8] > z2.reshape((5, 30))[2, 22]
    True

    An example using links:

    >>> grid1 = RasterModelGrid((10, 10), xy_spacing=100.0)
    >>> grid2 = RasterModelGrid((10, 10), xy_spacing=100.0)
    >>> z1 = grid1.add_zeros("topographic__elevation", at="node")
    >>> z2 = grid2.add_zeros("topographic__elevation", at="node")
    >>> dt = 1.0
    >>> nt = 10
    >>> grid2.at_link["surface_water__discharge"] = np.full(
    ...     grid2.number_of_links, 10000.0
    ... )
    >>> dfn1 = LinearDiffuser(grid1, linear_diffusivity=10000.0)
    >>> dfn2 = LinearDiffuser(grid2, linear_diffusivity="surface_water__discharge")
    >>> for i in range(nt):
    ...     z1[grid1.core_nodes] += 1.0
    ...     z2[grid2.core_nodes] += 1.0
    ...     dfn1.run_one_step(dt)
    ...     dfn2.run_one_step(dt)
    ...
    >>> np.allclose(z1, z2)
    True
    >>> z2.fill(0.0)
    >>> dfn2 = LinearDiffuser(
    ...     grid2,
    ...     linear_diffusivity="surface_water__discharge",
    ...     method="resolve_on_patches",
    ... )
    >>> for i in range(nt):
    ...     z2[grid2.core_nodes] += 1.0
    ...     dfn2.run_one_step(dt)
    ...
    >>> np.all(z2[grid2.core_nodes] < z1[grid2.core_nodes])
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
        method : {'simple', 'resolve_on_patches'}
            The method used to represent the fluxes. 'simple' solves a finite
            difference method with a simple staggered grid scheme onto the links.
            'resolve_on_patches' solves the scheme by mapping both slopes and
            diffusivities onto the patches and solving there before resolving
            values back to the nodes (and at the moment requires a raster grid).
            Note that this scheme enforces directionality in the diffusion field;
            it's no longer just a scalar field. Thus diffusivities must be defined
            *on links* when this option is chosen.
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
        method = self._validate_method(method)

        if method == "resolve_on_patches" and not isinstance(grid, RasterModelGrid):
            raise TypeError(
                "the resolve_on_patches method is only available for RasterModelGrid."
            )

        self._use_patches = method == "resolve_on_patches"

        self._current_time = 0.0
        self._run_before = False

        self._kd = self._validate_linear_diffusivity(grid, linear_diffusivity)

        if self._use_patches and np.ndim(self._kd) == 0:
            self._kd = np.broadcast_to(self._kd, grid.number_of_links)

        self._kd_on_links = np.size(self._kd) == grid.number_of_links

        if self._kd_on_links and not isinstance(grid, RasterModelGrid):
            raise TypeError(
                "linear_diffusivity defined at links is only available for"
                " RasterModelGrid."
            )

        # if we're using patches, it is VITAL that diffusivity is defined on
        # links. The whole point of this functionality is that we honour
        # *directionality* in the diffusivities.
        if self._use_patches and not self._kd_on_links:
            raise ValueError(
                "if using the resolve_on_patches method, linear_diffusivity"
                " must be defined at links."
            )

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

        self._CFL_actives_prefactor = CFL_prefactor[self._grid.active_links]
        # ^note we can do this as topology shouldn't be changing

        # Get a list of interior cells
        self._interior_cells = self._grid.node_at_core_cell

        self._z = self._grid.at_node[self._values_to_diffuse]
        self._dqsds = self._grid.zeros(at="node", dtype=float)

        for name in ("topographic__gradient", "hillslope_sediment__unit_volume_flux"):
            if name not in self._grid.at_link:
                self._grid.add_zeros(name, at="link")

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

    @staticmethod
    def _validate_method(method):
        valid_methods = {"simple", "resolve_on_patches"}

        if method not in valid_methods:
            raise ValueError(
                f"method {method} not understood"
                f" (must be one of {', '.join(sorted(valid_methods))})."
            )
        return method

    @staticmethod
    def _validate_linear_diffusivity(grid, linear_diffusivity):
        if isinstance(linear_diffusivity, str):
            if linear_diffusivity in grid.at_link:
                k = grid.at_link[linear_diffusivity]
            elif linear_diffusivity in grid.at_node:
                k = grid.at_node[linear_diffusivity]
            else:
                raise ValueError(
                    f"linear_diffusivity {linear_diffusivity!r}, it must be defined "
                    "at either nodes, or links."
                )
        elif np.ndim(linear_diffusivity) == 0:
            k = float(linear_diffusivity)
        else:
            k = np.asarray(linear_diffusivity)
            if k.size not in (grid.number_of_nodes, grid.number_of_links):
                raise ValueError(
                    "linear_diffusivity must be defined at either nodes, or links."
                )

        return k

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
        >>> z[mg.core_nodes] = 1.0
        >>> ld = LinearDiffuser(mg, linear_diffusivity=1.0)
        >>> ld.fixed_grad_nodes.size == 0
        True
        >>> ld.fixed_grad_anchors.size == 0
        True
        >>> ld.fixed_grad_offsets.size == 0
        True
        >>> mg.at_link["topographic__slope"] = mg.calc_grad_at_link(
        ...     "topographic__elevation"
        ... )
        >>> mg.status_at_node[mg.perimeter_nodes] = mg.BC_NODE_IS_FIXED_GRADIENT
        >>> ld.updated_boundary_conditions()
        >>> ld.fixed_grad_nodes
        array([ 1,  2,  3,  5,  9, 10, 14, 16, 17, 18])
        >>> ld.fixed_grad_anchors
        array([ 6,  7,  8,  6,  8, 11, 13, 11, 12, 13])
        >>> ld.fixed_grad_offsets
        array([-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.])
        >>> np.allclose(
        ...     z[ld.fixed_grad_nodes], z[ld.fixed_grad_anchors] + ld.fixed_grad_offsets
        ... )
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
        gradient = self._grid.at_link["topographic__gradient"]
        sediment_flux = self._grid.at_link["hillslope_sediment__unit_volume_flux"]

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
            grads = mg.calc_grad_at_link(z)
            gradient[mg.active_links] = grads[mg.active_links]
            if not self._use_patches:  # currently forbidden
                # if diffusivity is an array, self._kd is already
                # active_links-long
                sediment_flux[mg.active_links] = (
                    -kd_activelinks * gradient[mg.active_links]
                )
                # Calculate the net deposition/erosion rate at each node
                mg.calc_flux_div_at_node(sediment_flux, out=self._dqsds)
            else:  # project onto patches
                slx = mg.zeros("link")
                sly = mg.zeros("link")
                slx[self._hoz] = gradient[self._hoz]
                sly[self._vert] = gradient[self._vert]
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
                Cslope = np.sqrt(slx**2 + sly**2)
                v = np.sqrt(Kx**2 + Ky**2)
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
                sediment_flux[mg.active_links] = -flux_links[mg.active_links]

                self._grid.calc_flux_div_at_node(sediment_flux, out=self._dqsds)

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
