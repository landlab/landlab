#! /usr/env/python
"""

Component that models 2D diffusion using an explicit finite-volume method.

Created July 2013 GT
Last updated May 2015 DEJH

"""

##############DEJH is unsure if the uplift is correctly (not) incorporated here


from __future__ import print_function

import numpy as np
from six.moves import range

from landlab import ModelParameterDictionary, Component, FieldError, \
    create_and_initialize_grid, FIXED_GRADIENT_BOUNDARY, FIXED_LINK
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.utils.decorators import use_file_name_or_kwds

_ALPHA = 0.25   # time-step stability factor

#_VERSION = 'make_all_data'
#_VERSION = 'explicit'
_VERSION = 'pass_grid'


class LinearDiffuser(Component):
    """
    This component implements linear diffusion of a field in the supplied
    ModelGrid.

    This components requires the following parameters be set in the input file,
    *input_stream*, set in the component initialization:

        'linear_diffusivity', the diffusivity to use

    Optional inputs are:

        'uplift_rate', if you want this component to include the uplift
            internally
        'dt', the model timestep (assumed constant)
        'values_to_diffuse', a string giving the name of the grid field
        containing the data to diffuse.

    Supply *dt* to the diffuser through the diffuse() argument.
    This allows you to set a dynamic timestep for this class.
    If 'values_to_diffuse' is not provided, defaults to
    'topographic__elevation'.

    No particular units are necessary where they are not specified, as long as
    all units are internally consistent.

    The component takes *grid*, the ModelGrid object, and (optionally)
    *current_time* and *input_stream*. If *current_time* is not set, it
    defaults to 0.0. If *input_stream* is not set in instantiation of the
    class, :func:`initialize` with *input_stream* as in input must be called
    instead.
    *Input_stream* is the filename of (& optionally, path to) the parameter
    file.

    At the moment, this diffuser can only work with constant diffusivity.
    Spatially variable diffusivity hopefully coming soon.

    Component assumes grid does not deform, and BCs are in place before
    initialization.

    The primary method of this class is :func:`diffuse`.
    """

    _name = 'LinearDiffuser'

    _input_var_names = ('topographic__elevation',)

#############################UPDATE ME
    _output_var_names = ('topographic__elevation',
                         'surface_gradient',
                         'unit_flux')

    _var_units = {'topographic__elevation': 'm',
                  'surface_gradient': '-',
                  'unit_flux': 'm**3/s',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'surface_gradient': 'link',
                    'unit_flux': 'link',
                    }

    _var_doc = {
        'topographic__elevation': ('Land surface topographic elevation; can ' +
                                   'be overwritten in initialization'),
        'surface_gradient': 'Gradient of surface, on links',
        'unit_flux': 'Volume flux per unit width along links'
    }

    @use_file_name_or_kwds
    def __init__(self, grid, linear_diffusivity=None, **kwds):
        self._grid = grid
        self.current_time = 0.
        if linear_diffusivity is not None:
            self.kd = linear_diffusivity
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
            assert len(self.kd) is self.grid.number_of_nodes, self.kd
            kd_links = self.grid.map_max_of_link_nodes_to_link(self.kd)
        else:
            kd_links = float(self.kd)
        # assert CFL condition:
        dt_links = _ALPHA * self.grid.link_length ** 2. / kd_links
        self.dt = np.amin(dt_links[self.grid.active_links])

        # Get a list of interior cells
        self.interior_cells = self.grid.node_at_core_cell

        self.z = self.grid.at_node[self.values_to_diffuse]
        g = self.grid.zeros(centering='link')
        qs = self.grid.zeros(centering='link')
        try:
            self.g = self.grid.add_field('link', 'surface__gradient', g,
                                         noclobber=True)
            # ^note this will object if this exists already
        except FieldError:
            pass  # field exists, so no problem
        try:
            self.qs = self.grid.add_field('link', 'unit_flux', qs,
                                          noclobber=True)
        except FieldError:
            pass
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
        >>> ld.fixed_grad_nodes.size == 0 and ld.fixed_grad_anchors.size == 0
        ...     and ld.fixed_grad_offsets.size == 0
        True
        >>> mg.at_link['topographic__slope'] = mg.calculate_gradients_at_links(
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
        This is the primary method of the class. Call it to perform an iteration
        of the model. Takes *dt*, the current timestep.

        The modelgrid must contain the field to diffuse, which defaults to
        'topographic__elevation'. This can be overridden with the
        values_to_diffuse property in the input file.

        See the class docstring for a list of the other properties necessary
        in the input file for this component to run.

        To improve stability, this component can incorporate uplift into its
        internal mechanism. To use this, set *internal_uplift* to True, and . If you only have one module that requires this, do not add
        uplift manually in your loop; this method will include uplift
        automatically. If more than one of your components has this requirement,
        set *num_uplift_implicit_comps* to the total number of components that
        do.

        You can suppress this behaviour by setting *internal_uplift* to False.

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

        for i in range(repeats+1):
            # Calculate the gradients and sediment fluxes
            self.g[self.grid.active_links] = \
                self.grid.calculate_gradients_at_active_links(z)
            self.qs[self.grid.active_links] = (-self.kd *
                                               self.g[self.grid.active_links])

            # Calculate the net deposition/erosion rate at each node
            self.dqsds = self.grid.calculate_flux_divergence_at_nodes(
                self.qs[self.grid.active_links])

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
        """
        ::Updated docs go here::
        """
        self.diffuse(dt, **kwds)

    @property
    def time_step(self):
        """
        Returns time-step size (as a property).
        """
        return self.dt
