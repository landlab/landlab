#! /usr/env/python
"""

Component that models 2D diffusion using an explicit finite-volume method.

Created July 2013 GT
Last updated March 2016 DEJH with LL v1.0 component style
"""

from __future__ import print_function

import numpy as np
from six.moves import range

from landlab import ModelParameterDictionary, Component, FieldError, \
    create_and_initialize_grid, FIXED_GRADIENT_BOUNDARY, FIXED_LINK
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.utils.decorators import use_file_name_or_kwds

_ALPHA = 0.15   # time-step stability factor
# ^0.25 not restrictive enough at meter scales w S~1 (possible cases)


class LinearDiffuser(Component):
    """
    This component implements linear diffusion of a Landlab field.

    Component assumes grid does not deform. If the boundary conditions on the
    grid change after component instantiation, be sure to also call
    :func:`updated_boundary_conditions` to ensure these are reflected in the
    component (especially if fixed_links are present).

    The primary method of this class is :func:`run_one_step`.

    Construction::

        LinearDiffuser(grid, linear_diffusivity=None)

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    linear_diffusivity : float, array, or field name (m**2/time)
        The diffusivity.

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
    """

    _name = 'LinearDiffuser'

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
    def __init__(self, grid, linear_diffusivity=None, **kwds):
        self._grid = grid
        self.current_time = 0.
        if linear_diffusivity is not None:
            if type(linear_diffusivity) is not str:
                self.kd = linear_diffusivity
            else:
                self.kd = self.grid.at_node[linear_diffusivity]
        else:
            raise KeyError("linear_diffusivity must be provided to the " +
                           "LinearDiffuser component")

        # for component back compatibility (undocumented):
        # note component can NO LONGER do internal uplift, at all.
        # ###
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
        # ###

        # Set internal time step
        # ..todo:
        #   implement mechanism to compute time-steps dynamically if grid is
        #   adaptive/changing
        # as of modern componentization (Spring '16), this can take arrays
        # and irregular grids
        if type(self.kd) not in (int, float):
            assert len(self.kd) == self.grid.number_of_nodes, self.kd
            kd_links = self.grid.map_max_of_link_nodes_to_link(self.kd)
        else:
            kd_links = float(self.kd)
        # assert CFL condition:
        CFL_prefactor = _ALPHA * self.grid.length_of_link ** 2.
        self._CFL_actives_prefactor = CFL_prefactor[self.grid.active_links]
        # ^note we can do this as topology shouldn't be changing
        dt_links = CFL_prefactor / kd_links
        self.dt = np.amin(dt_links[self.grid.active_links])

        # Get a list of interior cells
        self.interior_cells = self.grid.node_at_core_cell

        self.z = self.grid.at_node[self.values_to_diffuse]
        g = self.grid.zeros(at='link')
        qs = self.grid.zeros(at='link')
        try:
            self.g = self.grid.add_field('link', 'topographic__gradient', g,
                                         noclobber=True)
            # ^note this will object if this exists already
        except FieldError:
            self.g = self.grid.at_link['topographic__gradient'] # keep a ref
        try:
            self.qs = self.grid.add_field('link', 'unit_flux', qs,
                                          noclobber=True)
        except FieldError:
            self.qs = self.grid.at_link['unit_flux']
        # note all these terms are deliberately loose, as we won't always be
        # dealing with topo

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

    def diffuse(self, dt, **kwds):
        """
        See :func:`run_one_step`.
        """
        if 'internal_uplift' in kwds.keys():
            raise KeyError('LinearDiffuser can no longer work with internal ' +
                           'uplift')
        # Take the smaller of delt or built-in time-step size self.dt
        self.tstep_ratio = dt/self.dt
        repeats = int(self.tstep_ratio//1.)
        extra_time = self.tstep_ratio - repeats
        z = self.grid.at_node[self.values_to_diffuse]

        core_nodes = self.grid.node_at_core_cell
        # do mapping of array kd here, in case it points at an updating field:
        if type(self.kd) is np.ndarray:
            kd_activelinks = self.grid.map_max_of_link_nodes_to_link(
                self.kd)[self.grid.active_links]
            # re-derive CFL condition, as could change dynamically:
            dt_links = self._CFL_actives_prefactor / kd_activelinks
            self.dt = dt_links.min()
        else:
            kd_activelinks = self.kd

        for i in range(repeats+1):
            # Calculate the gradients and sediment fluxes
            self.g[self.grid.active_links] = \
                    self.grid.calc_grad_at_link(z)[self.grid.active_links]
            # if diffusivity is an array, self.kd is already active_links-long
            self.qs[self.grid.active_links] = (-kd_activelinks *
                                               self.g[self.grid.active_links])

            # Calculate the net deposition/erosion rate at each node
            self.dqsds = self.grid.calc_flux_div_at_node(self.qs)
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
