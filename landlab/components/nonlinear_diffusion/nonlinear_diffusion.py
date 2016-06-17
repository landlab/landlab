#! /usr/env/python
"""

Component that models 2D nonlinear diffusion using an
explicit finite-volume method.

Created June 2016 NMG
"""

from __future__ import print_function

import numpy as np
from six.moves import range
from inspect import getmro

from landlab import ModelParameterDictionary, Component, FieldError, \
    create_and_initialize_grid, FIXED_GRADIENT_BOUNDARY, FIXED_LINK, \
    RasterModelGrid, INACTIVE_LINK
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.utils.decorators import use_file_name_or_kwds

_ALPHA = 0.15   # time-step stability factor
# maximum value for for ratio of slope to critical slope
_max_slope_ratio = 0.999


class NonLinearDiffuser(Component):
    """
    This component implements nonlinear diffusion of a Landlab field.

    Component assumes grid does not deform. If the boundary conditions on the
    grid change after component instantiation, be sure to also call
    :func:`updated_boundary_conditions` to ensure these are reflected in the
    component (especially if fixed_links are present).

    The primary method of this class is :func:`run_one_step`.

    Construction::

        NonLinearDiffuser(grid, nonlinear_diffusivity=None,
                          critical_slope=None)

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    nonlinear_diffusivity : float, array, or field name (m**2/time)
        The diffusivity. If an array or field name, these must be the
        diffusivities on either nodes or links - the component will
        distinguish which based on array length. Values on nodes will be
        mapped to links using an upwind scheme in the simple case.
    critical_slope : float, array, or field name (dimensionless, -dz/dx)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> mg = RasterModelGrid((9, 9), 1.)
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> z.reshape((9, 9))[4, 4] = 1.
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> ld = LinearDiffuser(mg, linear_diffusivity=1.)
    >>> for i in range(1):
    ...     ld.run_one_step(1.)
    >>> np.isclose(z[mg.core_nodes].sum(), 1.)
    True
    >>> mg2 = RasterModelGrid((5, 30), 1.)
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
    >>> mg1 = RasterModelGrid((10, 10), 100.)
    >>> mg2 = RasterModelGrid((10, 10), 100.)
    >>> z1 = mg1.add_zeros('node', 'topographic__elevation')
    >>> z2 = mg2.add_zeros('node', 'topographic__elevation')
    >>> dt = 1.
    >>> nt = 10
    >>> kappa_links = mg2.add_ones('link', 'water__discharge')
    >>> kappa_links *= 10000.
    >>> dfn1 = LinearDiffuser(mg1, linear_diffusivity=10000.)
    >>> dfn2 = LinearDiffuser(mg2, linear_diffusivity='water__discharge')
    >>> for i in range(nt):
    ...     z1[mg1.core_nodes] += 1.
    ...     z2[mg2.core_nodes] += 1.
    ...     dfn1.run_one_step(dt)
    ...     dfn2.run_one_step(dt)
    >>> np.allclose(z1, z2)
    True
    >>> z2.fill(0.)
    >>> dfn2 = LinearDiffuser(mg2, linear_diffusivity='water__discharge',
    ...                       method='resolve_on_patches')
    >>> for i in range(nt):
    ...     z2[mg2.core_nodes] += 1.
    ...     dfn2.run_one_step(dt)
    >>> np.all(z2[mg2.core_nodes] < z1[mg2.core_nodes])
    True
    """

    _name = 'NonLinearDiffuser'

    _input_var_names = ('topographic__elevation',)

    _output_var_names = ('topographic__elevation',
                         'topographic__gradient',
                         'unit_flux',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'topographic__gradient': '-',
                  'unit_flux': 'm**3/s',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'topographic__gradient': 'link',
                    'unit_flux': 'link',
                    }

    _var_doc = {
        'topographic__elevation': ('Land surface topographic elevation; can ' +
                                   'be overwritten in initialization'),
        'topographic__gradient': 'Gradient of surface, on links',
        'unit_flux': 'Volume flux per unit width along links'
    }

    @use_file_name_or_kwds
    def __init__(self, grid, nonlinear_diffusivity=None,
                 critical_slope=None, **kwds):
        self._grid = grid
        self.current_time = 0.
        self._run_before = False
        self._kd_on_links = False
        if nonlinear_diffusivity is not None:
            if type(nonlinear_diffusivity) is not str:
                self._kd = nonlinear_diffusivity
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
                    self._kd = self.grid.at_link[nonlinear_diffusivity]
                    self._kd_on_links = True
                except KeyError:
                    self._kd = self.grid.at_node[nonlinear_diffusivity]
        else:
            raise KeyError("nonlinear_diffusivity must be provided to the " +
                           "NonLinearDiffuser component")
        if self._kd_on_links is True:
            assert isinstance(self.grid, RasterModelGrid)

        # NMG ADD SOME CHECKS ON CRITICAL SLOPE to make sure that something
        # appropriate is passed.

        self.timestep_in = kwds.pop('dt', None)
        if 'values_to_diffuse' in kwds.keys():
            self.values_to_diffuse = kwds.pop('values_to_diffuse')
            for mytups in (self._input_var_names, self._output_var_names):
                myset = set(mytups)
                myset.remove('topographic__elevation')
                myset.add(self.values_to_diffuse)
                mytups = tuple(myset)
            for mydicts in (self._var_units, self._var_mapping, self._var_doc):
                mydicts[self.values_to_diffuse] = mydicts.pop(
                    'topographic__elevation')
        else:
            self.values_to_diffuse = 'topographic__elevation'
        # Raise an error if somehow someone is using this weird functionality
        if self._grid is None:
            raise ValueError('You must now provide an existing grid!')

        # Set internal time step
        # CFL condition precalc:
        # Note that this will need a correction that incorporates the slope
        # but this prefactor should not change through time.
        CFL_prefactor = _ALPHA * self.grid.length_of_link[
            :self.grid.number_of_links] ** 2.
        # ^ link_length can include diags, if not careful...
        self._CFL_actives_prefactor = CFL_prefactor[self.grid.active_links]
        # ^note we can do this as topology shouldn't be changing

        # Get a list of interior cells
        self.interior_cells = self.grid.node_at_core_cell

        self.z = self.grid.at_node[self.values_to_diffuse]
        self.dqsds = self.grid.zeros('node', dtype=float)
        g = self.grid.zeros(at='link')
        qs = self.grid.zeros(at='link')
        # slope to critical slope ratio
        slope_ratio = self.grid.zeros(at='link')
        # denominator in nonlinear diffusion qs equation
        denom = self.grid.zeros(at='link')
        try:
            self.g = self.grid.add_field('link', 'topographic__gradient',
                                         g, noclobber=True)
            # ^note this will object if this exists already
        except FieldError:  # keep a ref
            self.g = self.grid.at_link['topographic__gradient']
        try:
            self.qs = self.grid.add_field('link', 'unit_flux', qs,
                                          noclobber=True)
        except FieldError:
            self.qs = self.grid.at_link['unit_flux']
            # note all these terms are deliberately loose, as we won't always
            # be dealing with topo

        self._angle_of_link = self.grid.angle_of_link()
        self._vertlinkcomp = np.sin(self._angle_of_link)
        self._hozlinkcomp = np.cos(self._angle_of_link)

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
        >>> mg = RasterModelGrid((4, 5), 1.)
        >>> z = mg.add_zeros('node', 'topographic__elevation')
        >>> z[mg.core_nodes] = 1.
        >>> nld = NonLinearDiffuser(mg, nonlinear_diffusivity=1.)
        >>> nld.fixed_grad_nodes.size == 0
        True
        >>> nld.fixed_grad_anchors.size == 0
        True
        >>> nld.fixed_grad_offsets.size == 0
        True
        >>> mg.at_link['topographic__slope'] = mg.calc_grad_at_link(
        ...     'topographic__elevation')
        >>> mg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
        >>> nld.updated_boundary_conditions()
        >>> nld.fixed_grad_nodes
        array([ 1,  2,  3,  5,  9, 10, 14, 16, 17, 18])
        >>> nld.fixed_grad_anchors
        array([ 6,  7,  8,  6,  8, 11, 13, 11, 12, 13])
        >>> nld.fixed_grad_offsets
        array([-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.])
        >>> np.allclose(z[ld.fixed_grad_nodes],
        ...             z[ld.fixed_grad_anchors] + nld.fixed_grad_offsets)
        True
        """
        fixed_grad_nodes = np.where(self.grid.status_at_node ==
                                    FIXED_GRADIENT_BOUNDARY)[0]
        heads = self.grid.node_at_link_head[self.grid.fixed_links]
        tails = self.grid.node_at_link_tail[self.grid.fixed_links]
        head_is_fixed = np.in1d(heads, fixed_grad_nodes)
        self.fixed_grad_nodes = np.where(head_is_fixed, heads, tails)
        self.fixed_grad_anchors = np.where(head_is_fixed, tails, heads)
        vals = self.grid.at_node[self.values_to_diffuse]
        self.fixed_grad_offsets = (vals[self.fixed_grad_nodes] -
                                   vals[self.fixed_grad_anchors])

        if self._kd_on_links:
            mg = self.grid
            x_link_patch_pres = mg.patches_present_at_link[self._hoz]
            self._x_link_patch_mask = np.logical_not(x_link_patch_pres)
            y_link_patch_pres = mg.patches_present_at_link[self._vert]
            self._y_link_patch_mask = np.logical_not(y_link_patch_pres)
            self._hoz_link_neighbors[:, :2] = mg.links_at_node[
                mg.node_at_link_head[self._hoz], 1:4:2]
            self._hoz_link_neighbors[:, 2:] = mg.links_at_node[
                mg.node_at_link_tail[self._hoz], 1:4:2]
            self._vert_link_neighbors[:, :2] = mg.links_at_node[
                mg.node_at_link_head[self._vert], 0:3:2]
            self._vert_link_neighbors[:, 2:] = mg.links_at_node[
                mg.node_at_link_tail[self._vert], 0:3:2]
            self._vert_link_badlinks = np.logical_or(
                mg.status_at_link[self._vert_link_neighbors] == INACTIVE_LINK,
                self._vert_link_neighbors == -1)
            self._hoz_link_badlinks = np.logical_or(
                mg.status_at_link[self._hoz_link_neighbors] == INACTIVE_LINK,
                self._hoz_link_neighbors == -1)

    def nonlinear_diffuse(self, dt, **kwds):
        """
        See :func:`run_one_step`.
        """
        mg = self.grid
        z = self.grid.at_node[self.values_to_diffuse]

        if not self._run_before:
            self.updated_boundary_conditions()  # just in case
            self._run_before = True

        core_nodes = self.grid.node_at_core_cell
        # do mapping of array kd here, in case it points at an updating
        # field:
        # All dt calculations need to be repeated in the loop
        # because of changing slopes.
        if type(self._kd) is np.ndarray:
            if not self._kd_on_links:
                kd_links = self.grid.map_max_of_link_nodes_to_link(self._kd)
                kd_activelinks = kd_links[self.grid.active_links]
                # gradients at all links
                grads = mg.calc_grad_at_link(z)
                # ratio of slope to critical slope
                slope_ratio[mg.active_links] = grads[mg.active_links]/
                                                critical_slope
                # limit the maximum slope ratio
                slope_ratio[np.where(slope_ratio > _max_slope_ratio)] =
                                                _max_slope_ratio
                # denominator in nonlinear qs equation
                denom[mg.active_links] = 1 - (slope_ratio ** 2)
                # re-derive CFL condition, as could change dynamically:
                dt_links = self._CFL_actives_prefactor * denom[mg.active_links]
                             ** 1.5 / kd_activelinks
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

        # Take the smaller of delt or built-in time-step size self.dt
        self.tstep_ratio = dt / self.dt
        repeats = int(self.tstep_ratio // 1.)
        extra_time = self.tstep_ratio - repeats

        # Can really get into trouble if no diffusivity happens but we run...
        if self.dt < np.inf:
            loops = repeats+1
        else:
            loops = 0
        for i in range(loops):
            # gradients at all links
            grads = mg.calc_grad_at_link(z)
            # gradients just at active links
            self.g[mg.active_links] = grads[mg.active_links]
            # denominator in non linear qs equation
            denom[mg.active_links] = 1 - ((self.g[mg.active_links]
                                           [:mg.number_of_active_links] /
                                           critical_slope) ** 2)
            # now calculating qs at links
            self.qs[mg.active_links] = (
                        -kd_activelinks * self.g[mg.active_links]/denom)
            # Calculate the net deposition/erosion rate at each node
            mg.calc_flux_div_at_node(self.qs, out=self.dqsds)

            # Calculate the total rate of elevation change
            dzdt = - self.dqsds
            # Update the elevations
            timestep = self.dt
            if i == (repeats):
                timestep *= extra_time
            else:
                pass
            self.grid.at_node[self.values_to_diffuse][core_nodes] += dzdt[
                core_nodes] * timestep

            # check the BCs, update if fixed gradient
            vals = self.grid.at_node[self.values_to_diffuse]
            vals[self.fixed_grad_nodes] = (vals[self.fixed_grad_anchors] +
                                           self.fixed_grad_offsets)

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
