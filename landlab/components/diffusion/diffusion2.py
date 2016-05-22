# -*- coding: utf-8 -*-
"""
This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router. Credit: Voller, Hobley, Paola.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

import numpy as np
from landlab import RasterModelGrid, Component, FieldError
import inspect
from landlab.utils.decorators import use_file_name_or_kwds

_ALPHA = 0.15   # time-step stability factor
# ^0.25 not restrictive enough at meter scales w S~1 (possible cases)


class LinearDiffuser2(Component):
    """
    This class implements a linear diffusion scheme on a raster grid. It is
    more precise than the LinearDiffuser, and less prone to grid artifacts,
    but requires a raster and is slower.

    The router takes account of the elevation values at diagonal nodes when
    performing the diffusion (i.e., it resolves slopes on patches, not links.)
    """
    _name = 'LinearDiffuser2'

    _input_var_names = ('topographic__elevation',
                        )

    _output_var_names = ('topographic__elevation',
                         )

    _var_units = {'topographic__elevation': 'm',
                  }

    _var_mapping = {
        'topographic__elevation': 'node',
        }

    _var_doc = {
        'topographic__elevation':
            'Land surface topographic elevation',
        }

    _min_slope_thresh = 1.e-24
    # ^if your flow isn't connecting up, this probably needs to be reduced

    @use_file_name_or_kwds
    def __init__(self, grid, linear_diffusivity=None,
                 critical_slope=0., **kwds):
        assert RasterModelGrid in inspect.getmro(grid.__class__)
        assert grid.number_of_node_rows >= 3
        assert grid.number_of_node_columns >= 3
        assert linear_diffusivity is not None, "Supply a diffusivity!"

        self._grid = grid
        if type(linear_diffusivity) in (float, int):
            self._kd = self.grid.empty('node')
            self._kd[:] = linear_diffusivity
        elif type(linear_diffusivity) is str:
            self._kd = mg.at_node[linear_diffusivity]
        else:
            assert linear_diffusivity.size == self.grid.number_of_nodes
            self._kd = linear_diffusivity
        if type(critical_slope) is np.ndarray:
            self._critS = self.grid.empty('node')
            self._critS[:] = critical_slope
        else:
            self._critS = critical_slope

        ncols = grid.number_of_node_columns
        nrows = grid.number_of_node_rows
        # create the blank node maps; assign the topo to an internally held
        # 2D map with dummy edges:
        self.elev_raster = np.empty((nrows+2, ncols+2), dtype=float)
        self.hR = np.zeros((nrows+2, ncols+2), dtype=float)
        hR = self.hR

        self.uE = np.zeros_like(hR)
        self.uW = np.zeros_like(hR)
        self.uN = np.zeros_like(hR)
        self.uS = np.zeros_like(hR)

        self.qsedE = np.zeros_like(hR)
        self.qsedW = np.zeros_like(hR)
        self.qsedN = np.zeros_like(hR)
        self.qsedS = np.zeros_like(hR)

        # setup slices for use in IDing the neighbors
        self._Es = (slice(1, -1), slice(2, ncols+2))
        self._NEs = (slice(2, nrows+2), slice(2, ncols+2))
        self._Ns = (slice(2, nrows+2), slice(1, -1))
        self._NWs = (slice(2, nrows+2), slice(0, -2))
        self._Ws = (slice(1, -1), slice(0, -2))
        self._SWs = (slice(0, -2), slice(0, -2))
        self._Ss = (slice(0, -2), slice(1, -1))
        self._SEs = (slice(0, -2), slice(2, ncols+2))
        self._core = (slice(1, -1), slice(1, -1))
        self._corecore = (slice(2, -2), slice(2, -2))
        # ^the actual, LL-sense core (interior) nodes of the grid

        # assert CFL condition:
        CFL_prefactor = _ALPHA * self.grid.link_length[
            :self.grid.number_of_links] ** 2.
        # ^ link_length can include diags, if not careful...
        self._CFL_actives_prefactor = CFL_prefactor[self.grid.active_links]
        # ^note we can do this as topology shouldn't be changing
        if type(linear_diffusivity) in (float, int):
            self.const_dt = self._CFL_actives_prefactor/linear_diffusivity
        else:
            self.const_dt = None

    def diffuse(self, dt, **kwds):
        """
        """
        upwind_kd_links = self.grid.map_max_of_link_nodes_to_link(self._kd)
        if self.const_dt is None:
            all_poss_dts = self._CFL_actives_prefactor/upwind_kd_links[
                self.grid.active_links]
            dt_internal = all_poss_dts.min()
        else:
            dt_internal = self.const_dt
        if dt > dt_internal:
            repeats = int(dt//dt_internal) + 1
            dt_internal = dt/repeats  # cleaner this way
        else:
            repeats = 1
            dt_internal = dt

        # create aliases for ease
        hR = self.hR
        hR[1:-1, 1:-1].flat = self._grid.at_node['topographic__elevation']

        Es = self._Es
        NEs = self._NEs
        Ns = self._Ns
        NWs = self._NWs
        Ws = self._Ws
        SWs = self._SWs
        Ss = self._Ss
        SEs = self._SEs
        core = self._core
        slope = self._critS

        qsedE = self.qsedE
        qsedW = self.qsedW
        qsedN = self.qsedN
        qsedS = self.qsedS
        uE = self.uE
        uN = self.uN
        uW = self.uW
        uS = self.uS

        for i in range(repeats):
            # update the dummy edges of our variables - these all act as closed
            # nodes (the inner, true boundaries are handled elsewhere... or
            # they should be closed anyway!):
            # note this isn't sufficient of we have diagonals turned on, as
            # flow can still occur on them
            hR[0, 1:-1] = hR[1, 1:-1]
            hR[-1, 1:-1] = hR[-2, 1:-1]
            hR[1:-1, 0] = hR[1:-1, 1]
            hR[1:-1, -1] = hR[1:-1, -2]
            hR[(0, -1, 0, -1), (0, -1, -1, 0)] = hR[(1, -2, 1, -2),
                                                    (1, -2, -2, 1)]
            qsedE.fill(0.)
            qsedW.fill(0.)
            qsedN.fill(0.)
            qsedS.fill(0.)
            # define the fluxes:
            grads_on_links = self.grid.calculate_gradients_at_links(hR[core])
            # apply an upwind scheme per VV:
            grads_on_links *= upwind_kd_links
            uE[core].flat = grads_on_links[self.grid.links_at_node[
                :, 0]]*-self.grid.link_dirs_at_node[:, 0]
            uN[core].flat = grads_on_links[self.grid.links_at_node[
                :, 1]]*-self.grid.link_dirs_at_node[:, 1]
            uW[core].flat = grads_on_links[self.grid.links_at_node[
                :, 2]]*-self.grid.link_dirs_at_node[:, 2]
            uS[core].flat = grads_on_links[self.grid.links_at_node[
                :, 3]]*-self.grid.link_dirs_at_node[:, 3]
            # note the sign switch as VV used out is +ve convention, LL uses
            # out is -ve. Process should automatically zero the grad of any
            # nonexistent link.
            # update the u BCs
            for BC in (uW, uE, uN, uS):
                BC[0, 1:-1] = BC[1, 1:-1]
                BC[-1, 1:-1] = BC[-2, 1:-1]
                BC[1:-1, 0] = BC[1:-1, 1]
                BC[1:-1, -1] = BC[1:-1, -2]
                BC[(0, -1, 0, -1), (0, -1, -1, 0)] = BC[(1, -2, 1, -2),
                                                        (1, -2, -2, 1)]

            hgradEx = (hR[core]-hR[Es])/self.grid.dx
            hgradEy = hR[SEs]-hR[NEs]+hR[Ss]-hR[Ns]
            hgradEy *= 0.25
            hgradEy /= self.grid.dy
            CslopeE = np.sqrt(np.square(hgradEx)+np.square(hgradEy))
            thetaE = np.arctan(np.fabs(hgradEy)/(np.fabs(hgradEx)+1.e-10))
            pgradEx = uE[core]  # pgrad is VV's vv, a velocity
            pgradEy = uN[core]+uS[core]+uN[Es]+uS[Es]
            pgradEy *= 0.25
            vmagE = np.sqrt(np.square(pgradEx)+np.square(pgradEy))
            # now resolve the effective flow magnitudes to downhill
            theta_vE = np.arctan(np.fabs(pgradEy)/(np.fabs(pgradEx)+1.e-10))
            vmagE *= np.cos(np.fabs(thetaE-theta_vE))
            qsedE = np.sign(hgradEx)*vmagE*(CslopeE-slope).clip(0.)*np.cos(
                thetaE)

            hgradWx = (hR[Ws]-hR[core])/self.grid.dx
            hgradWy = hR[SWs]-hR[NWs]+hR[Ss]-hR[Ns]
            hgradWy *= 0.25
            hgradWy /= self.grid.dy
            CslopeW = np.sqrt(np.square(hgradWx)+np.square(hgradWy))
            thetaW = np.arctan(np.fabs(hgradWy)/(np.fabs(hgradWx)+1.e-10))
            pgradWx = uW[core]
            pgradWy = uN[core]+uS[core]+uN[Ws]+uS[Ws]
            pgradWy *= 0.25
            vmagW = np.sqrt(np.square(pgradWx)+np.square(pgradWy))
            theta_vW = np.arctan(np.fabs(pgradWy)/(np.fabs(pgradWx)+1.e-10))
            vmagW *= np.cos(np.fabs(thetaW-theta_vW))
            qsedW = np.sign(hgradWx)*vmagW*(CslopeW-slope).clip(0.)*np.cos(
                thetaW)

            hgradNx = hR[NWs]-hR[NEs]+hR[Ws]-hR[Es]
            hgradNx *= 0.25
            hgradNx /= self.grid.dx
            hgradNy = (hR[core]-hR[Ns])/self.grid.dy
            CslopeN = np.sqrt(np.square(hgradNx)+np.square(hgradNy))
            thetaN = np.arctan(np.fabs(hgradNy)/(np.fabs(hgradNx)+1.e-10))
            pgradNx = uE[core]+uW[core]+uE[Ns]+uW[Ns]
            pgradNx *= 0.25
            pgradNy = uN[core]
            vmagN = np.sqrt(np.square(pgradNx)+np.square(pgradNy))
            theta_vN = np.arctan(np.fabs(pgradNy)/(np.fabs(pgradNx)+1.e-10))
            vmagN *= np.cos(np.fabs(thetaN-theta_vN))
            qsedN = np.sign(hgradNy)*vmagN*(CslopeN-slope).clip(0.)*np.sin(
                thetaN)

            hgradSx = hR[SWs]-hR[SEs]+hR[Ws]-hR[Es]
            hgradSx *= 0.25
            hgradSx /= self.grid.dx
            hgradSy = (hR[Ss]-hR[core])/self.grid.dy
            CslopeS = np.sqrt(np.square(hgradSx)+np.square(hgradSy))
            thetaS = np.arctan(np.fabs(hgradSy)/(np.fabs(hgradSx)+1.e-10))
            pgradSx = uE[core]+uW[core]+uE[Ss]+uW[Ss]
            pgradSx *= 0.25
            pgradSy = uS[core]
            vmagS = np.sqrt(np.square(pgradSx)+np.square(pgradSy))
            theta_vS = np.arctan(np.fabs(pgradSy)/(np.fabs(pgradSx)+1.e-10))
            vmagS *= np.cos(np.fabs(thetaS-theta_vS))
            qsedS = np.sign(hgradSy)*vmagS*(CslopeS-slope).clip(0.)*np.sin(
                thetaS)

            hR[core] += dt_internal*(qsedS+qsedW-qsedN-qsedE)

        self.grid.at_node['topographic__elevation'][:] = hR[core].flat

    def run_one_step(self, dt, **kwds):
        self.diffuse(dt, **kwds)
