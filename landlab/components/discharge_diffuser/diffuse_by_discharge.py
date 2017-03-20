# -*- coding: utf-8 -*-
"""
This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router. Credit: Voller, Hobley, Paola.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

import numpy as np
from landlab import RasterModelGrid, Component, FieldError, INACTIVE_LINK, \
    CLOSED_BOUNDARY, CORE_NODE
import inspect
from landlab.utils.decorators import use_file_name_or_kwds


class DischargeDiffuser(Component):
    """Diffuse sediment proportional to an implicit water discharge value.

    This class implements Voller, Hobley, and Paola's scheme for sediment
    diffusion, where the diffusivity of the sediment is proportional to the
    local discharge of water. The method works by solving for a potential
    field describing the water discharge at all nodes on the grid, which
    enforces both mass conservation and flow downhill along topographic
    gradients. This routine is designed to construct sediment fans.

    Note that both the water and sediment discharges are calculated together
    within the component.

    The algorithm uses a rule that looks like:

        q_sed = q_water * (S - S_crit)

    where S_crit is a critical slope threshold. [MODIFY THIS]

    It is VITAL you initialize this component AFTER setting boundary
    conditions.

    The primary method of this class is :func:`run_one_step`.

    Construction::

        DischargeDiffuser(grid, ...)

    Notes
    -----
    This is a "research grade" component, and is subject to dramatic change
    with little warning. No guarantees are made regarding its accuracy or
    utility. It is not recommended for user use yet!

    Parameters
    ----------
    grid : ModelGrid
        A grid.


    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> import numpy as np
    >>> mg = HexModelGrid(4, 6, dx=2., shape='rect', orientation='vertical')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> Q_in = mg.add_ones('node', 'water__unit_flux_in')
    >>> z += mg.node_y.copy()
    >>> potfr = PotentialityFlowRouter(mg)
    >>> potfr.run_one_step()
    >>> Q_at_core_nodes = np.array(
    ...     [ 17.02012846,  16.88791903,  13.65746194,  14.85578934,
    ...       11.41908145,  11.43630865,   8.95902559,  10.04348075,
    ...        6.28696459,   6.44316089,   4.62478522,   5.29145188])
    >>> np.allclose(mg.at_node['surface_water__discharge'][mg.core_nodes],
    ...             Q_at_core_nodes)
    True
    """
    _name = 'DischargeDiffuser'

    _input_var_names = ('topographic__elevation',
                        'water__discharge_in',
                        'sediment__discharge_in',
                        )

    _output_var_names = ('topographic__elevation',
                         'surface_water__discharge',
                         'flow__potential',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'water__unit_flux_in': 'm/s',
                  'surface_water__discharge': 'm**3/s',
                  'flow__potential': 'm**3/s',
                  'surface_water__depth': 'm',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'water__unit_flux_in': 'node',
                    'surface_water__discharge': 'node',
                    'flow__potential': 'node',
                    'surface_water__depth': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'water__unit_flux_in': (
            'External volume water per area per time input to each node ' +
            '(e.g., rainfall rate)'),
        'surface_water__discharge': (
            'Magnitude of volumetric water flux out of each node'),
        'flow__potential': (
            'Value of the hypothetical field "K", used to force water flux ' +
            'to flow downhill'),
        'surface_water__depth': (
            'If Manning or Chezy specified, the depth of flow in the cell, ' +
            'calculated assuming flow occurs over the whole surface'),
                  }

    _min_slope_thresh = 1.e-24
    # if your flow isn't connecting up, this probably needs to be reduced

    @use_file_name_or_kwds
    def __init__(self, grid, slope, **kwds):
        """Initialize flow router.
        """
        if RasterModelGrid in inspect.getmro(grid.__class__):
            assert grid.number_of_node_rows >= 3
            assert grid.number_of_node_columns >= 3
            self._raster = True
        else:
            self._raster = False

        assert self._raster = True  # ...for now

        self._grid = grid
        self._slope = slope

        # hacky fix because water__discharge is defined on both links and nodes
        for out_field in self._output_var_names:
            if self._var_mapping[out_field] == 'node':
                try:
                    self.grid.add_zeros(self._var_mapping[out_field],
                                        out_field, dtype=float)
                except FieldError:
                    pass
            else:
                pass
            try:
                self.grid.add_zeros('node', 'surface_water__discharge',
                                    dtype=float)
            except FieldError:
                pass
        ni = grid.number_of_node_rows
        nj = grid.number_of_node_columns

        self._K = grid.zeros('node', dtype=float)
        self._prevK = grid.zeros('node', dtype=float)
        self._znew = grid.zeros('node', dtype=float)
        # discharge across north, south, west, and east face of control volume
        self._Qn = np.zeros((ni, nj), dtype='float')
        self._Qs = np.zeros((ni, nj), dtype='float')
        self._Qw = np.zeros((ni, nj), dtype='float')
        self._Qe = np.zeros((ni, nj), dtype='float')

        # coefficenst used in solition of flow conductivity K
        self._app = np.zeros((ni, nj), dtype='float')
        self._apz = np.zeros((ni, nj), dtype='float')
        self._aww = np.zeros((ni, nj), dtype='float')
        self._awp = np.zeros((ni, nj), dtype='float')
        self._awz = np.zeros((ni, nj), dtype='float')
        self._aee = np.zeros((ni, nj), dtype='float')
        self._aep = np.zeros((ni, nj), dtype='float')
        self._aez = np.zeros((ni, nj), dtype='float')
        self._ass = np.zeros((ni, nj), dtype='float')
        self._asp = np.zeros((ni, nj), dtype='float')
        self._asz = np.zeros((ni, nj), dtype='float')
        self._ann = np.zeros((ni, nj), dtype='float')
        self._anp = np.zeros((ni, nj), dtype='float')
        self._anz = np.zeros((ni, nj), dtype='float')

        self._slx = np.empty((ni, nj), dtype=float)
        self._sly = np.empty((ni, nj), dtype=float)
        self._Qsed_w = np.empty((ni, nj), dtype=float)
        self._Qsed_e = np.empty((ni, nj), dtype=float)
        self._Qsed_n = np.empty((ni, nj), dtype=float)
        self._Qsed_s = np.empty((ni, nj), dtype=float)

    def run_one_step(self, dt, Qsource, **kwds):
        """
        """
        grid = self.grid
        ni = grid.number_of_node_rows
        nj = grid.number_of_node_columns
        z = grid.at_node['topographic__elevation']
        Qsp = grid.at_node['water__discharge_in'].reshape(
            (grid.number_of_node_rows, grid.number_of_node_columns))
        Qsource = grid.at_node['sediment__discharge_in'].reshape(
            (grid.number_of_node_rows, grid.number_of_node_columns))
        mismatch = 10000.

        ######STABILITY ANALYSIS GOES HERE
        dt_stab = dt

        # elevation at current and new time
        # Note a horizonal surface is the initial condition
        eta = z.reshape((ni, nj))
        # etan = self._znew.reshape((grid.number_of_node_rows,
        #                            grid.number_of_node_columns))

        # do the sediment diffusion
        for dir in ('W', 'E', 'S', 'N'):
            self._grad_on_link(dir)
            Cslope = np.sqrt(slx**2 + sly**2)
            self._link_sed_flux_from_slope(Cslope, self._slope, dir)

        try:
            Qin = Qsource.reshape((ni, nj))
        except AttributeError:
            Qin = float(Qsource)  # if both fail, we're in trouble
        eta[:] += dt_stab*(
            self._Qsed_e + self._Qsed_n + self._Qsed_w + self._Qsed_s +
            Qin)/self.grid.dx/self.grid.dy

        # do the water routing on links
        # These calculations are based on assuming that the flow is a sheet
        # flow that can be characterized with a potential equation. If this
        # flow is isotropic (an assumption we should revisit) with this model
        # the flow discharge in the x-direction (for example) can be calculated
        # as a constant (K the 'flow conductivity') times the component of the
        # sediment slope in that direction. It helps to define a 'slope
        # velocity' u, with components ustar=-deta/dx and vstar=-deta/dx which
        # allows us to write down the following advection like gov. eq. for
        # the flow  ----div(Ku)+Q=0---where Q represents external flow inputs

        # Since we can readily determine u from the current sediment topography
        # We solve this equation for K using an upwind scheme







        # do the ortho nodes first, in isolation
        g = grid.calc_grad_at_link(z)
        if self.equation != 'default':
            g = np.sign(g) * np.sqrt(np.fabs(g))
            # ^...because both Manning and Chezy actually follow sqrt
            # slope, not slope
        # weight by face width - NO, because diags
        # g *= grid.width_of_face[grid.face_at_link]
        link_grad_at_node_w_dir = (
            g[grid.links_at_node] * grid.active_link_dirs_at_node)
        # active_link_dirs strips "wrong" face widths

        # now outgoing link grad sum
        outgoing_sum = (np.sum((link_grad_at_node_w_dir).clip(0.), axis=1) +
                        self._min_slope_thresh)
        pos_incoming_link_grads = (-link_grad_at_node_w_dir).clip(0.)

        if not self.route_on_diagonals or not self._raster:
            while mismatch > 1.e-6:
                K_link_ends = self._K[grid.neighbors_at_node]
                incoming_K_sum = (pos_incoming_link_grads*K_link_ends
                                  ).sum(axis=1) + self._min_slope_thresh
                self._K[:] = (incoming_K_sum + Qwater_in)/outgoing_sum
                mismatch = np.sum(np.square(self._K-prev_K))
                prev_K = self._K.copy()

            upwind_K = grid.map_value_at_max_node_to_link(z, self._K)
            self._discharges_at_link[:] = upwind_K * g
            self._discharges_at_link[grid.status_at_link == INACTIVE_LINK] = 0.
        else:
            # grad on diags:
            gwd = np.empty(grid._number_of_d8_links, dtype=float)
            gd = gwd[grid.number_of_links:]
            gd[:] = (z[grid._diag_link_tonode] - z[grid._diag_link_fromnode])
            gd /= (grid._length_of_link_with_diagonals[grid.number_of_links:])
            if self.equation != 'default':
                gd[:] = np.sign(gd)*np.sqrt(np.fabs(gd))
            diag_grad_at_node_w_dir = (gwd[grid._diagonal_links_at_node] *
                                       grid._diag_active_link_dirs_at_node)

            outgoing_sum += np.sum(diag_grad_at_node_w_dir.clip(0.), axis=1)
            pos_incoming_diag_grads = (-diag_grad_at_node_w_dir).clip(0.)
            while mismatch > 1.e-6:
                K_link_ends = self._K[grid.neighbors_at_node]
                K_diag_ends = self._K[grid._diagonal_neighbors_at_node]
                incoming_K_sum = ((pos_incoming_link_grads * K_link_ends
                                   ).sum(axis=1) +
                                  (pos_incoming_diag_grads * K_diag_ends
                                   ).sum(axis=1) + self._min_slope_thresh)
                self._K[:] = (incoming_K_sum + Qwater_in) / outgoing_sum
                mismatch = np.sum(np.square(self._K - prev_K))
                prev_K = self._K.copy()

            # ^this is necessary to suppress stupid apparent link Qs at flow
            # edges, if present.
            upwind_K = grid.map_value_at_max_node_to_link(z, self._K)
            upwind_diag_K = np.where(
                z[grid._diag_link_tonode] > z[grid._diag_link_fromnode],
                self._K[grid._diag_link_tonode],
                self._K[grid._diag_link_fromnode])
            self._discharges_at_link[:grid.number_of_links] = upwind_K * g
            self._discharges_at_link[grid.number_of_links:] = (
                upwind_diag_K * gd)
            self._discharges_at_link[grid._all_d8_inactive_links] = 0.

        np.multiply(self._K, outgoing_sum, out=self._Qw)
        # there is no sensible way to save discharges at links, if we route
        # on diagonals.
        # for now, let's make a property

        # now process uval and vval to give the depths, if Chezy or Manning:
        if self.equation == 'Chezy':
            # Chezy: Q = C*Area*sqrt(depth*slope)
            grid.at_node['surface_water__depth'][:] = (
                grid.at_node['flow__potential'] / self.chezy_C /
                self.equiv_circ_diam) ** (2. / 3.)
        elif self.equation == 'Manning':
            # Manning: Q = w/n*depth**(5/3)
            grid.at_node['surface_water__depth'][:] = (
                grid.at_node['flow__potential'] * self.manning_n /
                self.equiv_circ_diam) ** 0.6
        else:
            pass

    def run_one_step(self, **kwds):
        """Route surface-water flow over a landscape.

        Both convergent and divergent flow can occur.
        """
        self.route_flow(**kwds)

    @property
    def discharges_at_links(self):
        """Return the discharges at links.

        Note that if diagonal routing, this will return number_of_d8_links.
        Otherwise, it will be number_of_links.
        """
        return self._discharges_at_link

    def _grad_on_link(self, direction):
        """
        Updates slx and sly with link gradient values according to `direction`.

        direction = {'E', 'N', 'S', 'W'}
        """
        slice_e = (slice(0, ni, 1), slice(1, nj, 1))
        slice_n = (slice(1, ni, 1), slice(0, nj, 1))
        slice_w = (slice(0, ni, 1), slice(0, nj-1, 1))
        slice_s = (slice(0, ni-1, 1), slice(0, nj, 1))
        slice_ne = (slice(1, ni, 1), slice(1, nj, 1))
        slice_sw = (slice(0, ni-1, 1), slice(0, nj-1, 1))
        slice_nw = (slice(1, ni, 1), slice(0, nj-1, 1))
        slice_se = (slice(0, ni-1, 1), slice(1, nj, 1))

        if direction == 'W':
            self._slx[slice_e] = (eta[slice_w] - eta[slice_e])/self.grid.dx
            self._slx[:, 0] = 0.  # W col gets 0

            self._sly.fill(0.)
            self._sly[slice_ne] += eta[slice_sw]
            self._sly[slice_se] -= eta[slice_nw]
            self._sly[slice_n] += eta[slice_s]
            self._sly[slice_s] -= eta[slice_n]

            self._sly[0, 1:] += eta[0, :-1]  # S row add node to W not SW
            self._sly[0, :] += eta[0, :]  # S row add self not S
            self._sly[-1, 1:] -= eta[-1, :-1]  # N row less node to W not NW
            self._sly[-1, :] -= eta[-1, :]  # N row less self not N
            self._sly[1:, 0] += eta[:-1, 0]  # W col add node to S not SW
            self._sly[:-1, 0] -= eta[1:, 0]  # W col less node to N not NW
            self._sly[0, 0] += eta[0, 0]  # SW corner add self, not SW
            self._sly[-1, 0] -= eta[-1, 0]  # NW corner less self, not NW
            self._sly *= 0.25
            self._sly /= self.grid.dy

        elif direction == 'E':
            self._slx[slice_w] = (eta[slice_e] - eta[slice_w])/self.grid.dx
            self._slx[:, -1] = 0.  # E col gets 0

            self._sly.fill(0.)
            self._sly[slice_nw] += eta[slice_se]
            self._sly[slice_sw] -= eta[slice_ne]
            self._sly[slice_n] += eta[slice_s]
            self._sly[slice_s] -= eta[slice_n]

            self._sly[0, :-1] += eta[0, 1:]  # S row add node to E not SE
            self._sly[0, :] += eta[0, :]  # S row add self not S
            self._sly[-1, :-1] -= eta[-1, 1:]  # N row less node to E not NE
            self._sly[-1, :] -= eta[-1, :]  # N row less self not N
            self._sly[1:, -1] += eta[:-1, -1]  # E col add node to S not SE
            self._sly[:-1, -1] -= eta[1:, -1]  # E col less node to N not NE
            self._sly[0, -1] += eta[0, -1]  # SE corner add self, not SE
            self._sly[-1, -1] -= eta[-1, -1]  # NE corner less self, not NE
            self._sly *= 0.25
            self._sly /= self.grid.dy

        elif direction == 'S':
            self._sly[slice_n] = (eta[slice_s] - eta[slice_n])/self.grid.dy
            self._sly[0, :] = 0.  # S col gets 0

            self._slx.fill(0.)
            self._slx[slice_ne] += eta[slice_sw]
            self._slx[slice_nw] -= eta[slice_se]
            self._slx[slice_e] += eta[slice_w]
            self._slx[slice_w] -= eta[slice_e]

            self._slx[1:, 0] += eta[:-1, 0]  # W col add node to S not SW
            self._slx[:, 0] += eta[:, 0]  # W col add self not W
            self._slx[1:, -1] -= eta[:-1, -1]  # E col less node to S not SE
            self._slx[:, -1] -= eta[:, -1]  # E col less self not E
            self._slx[0, 1:] += eta[0, :-1]  # S row add node to W not SW
            self._slx[0, :-1] -= eta[0, 1:]  # S row less node to E not SE
            self._slx[0, 0] += eta[0, 0]  # SW corner add self, not SW
            self._slx[0, -1] -= eta[0, -1]  # SE corner less self, not SE
            self._slx *= 0.25
            self._slx /= self.grid.dx

        elif direction == 'N':
            self._sly[slice_s] = (eta[slice_n] - eta[slice_s])/self.grid.dy
            self._sly[-1, :] = 0.  # N col gets 0

            self._slx.fill(0.)
            self._slx[slice_se] += eta[slice_nw]
            self._slx[slice_sw] -= eta[slice_ne]
            self._slx[slice_e] += eta[slice_w]
            self._slx[slice_w] -= eta[slice_e]

            self._slx[:-1, 0] += eta[1:, 0]  # W col add node to N not NW
            self._slx[:, 0] += eta[:, 0]  # W col add self not W
            self._slx[:-1, -1] -= eta[1:, -1]  # E col less node to N not NE
            self._slx[:, -1] -= eta[:, -1]  # E col less self not E
            self._slx[-1, 1:] += eta[-1, :-1]  # N row add node to W not NW
            self._slx[-1, :-1] -= eta[-1, 1:]  # N row less node to E not NE
            self._slx[-1, 0] += eta[-1, 0]  # NW corner add self, not NW
            self._slx[-1, -1] -= eta[-1, -1]  # NE corner less self, not NE
            self._slx *= 0.25
            self._slx /= self.grid.dx

        else:
            raise NameError("direction must be {'E', 'N', 'S', 'W'}")

    def _link_sed_flux_from_slope(self, S_val, S_thresh, direction):
        """
        Update the sed flux array for a given link dir, assuming a critical S.
        """
        if direction == 'W':
            dir_sed_flux = self._Qsed_w
            dir_water_flux = self._Qw
            thisslice = (slice(0, -1, 1), slice(1, -1, 1))
            deadedge = (slice(0, -1, 1), slice(0, 1, 1))
        elif direction == 'E':
            dir_sed_flux = self._Qsed_e
            dir_water_flux = self._Qe
            thisslice = (slice(0, -1, 1), slice(1, -2, 1))
            deadedge = (slice(0, -1, 1), slice(-2, -1, 1))
        elif direction == 'N':
            dir_sed_flux = self._Qsed_n
            dir_water_flux = self._Qn
            thisslice = (slice(0, -2, 1), slice(0, -1, 1))
            deadedge = (slice(-2, -1, 1), slice(0, -1, 1))
        elif direction == 'S':
            dir_sed_flux = self._Qsed_s
            dir_water_flux = self._Qs
            thisslice = (slice(1, -1, 1), slice(0, -1, 1))
            deadedge = (slice(0, 1, 1), slice(0, -1, 1))
        else:
            raise NameError("direction must be {'E', 'N', 'S', 'W'}")
        slope_diff = (S_val - S_thresh).clip(0.)
        dir_sed_flux[thisslice] = (dir_water_flux[thisslice] *
                                   slope_diff[thisslice])
        dir_sed_flux[deadedge] = 0.
