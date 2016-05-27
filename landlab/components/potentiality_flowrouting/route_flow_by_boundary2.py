# -*- coding: utf-8 -*-
"""
This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router. Credit: Voller, Hobley, Paola.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

###in the diagonal case, even closed edges can produce "drag". Is this right?
#Could suppress by mirroring the diagonals

import numpy as np
from landlab import RasterModelGrid, Component, FieldError, INACTIVE_LINK
import inspect
from landlab.utils.decorators import use_file_name_or_kwds


class PotentialityFlowRouter(Component):
    """
    This class implements Voller, Hobley, and Paola's experimental matrix
    solutions for flow routing.
    Options are permitted to allow "abstract" routing (flow enforced downslope,
    but no particular assumptions are made about the governing equations), or
    routing according to the Chezy or Manning equations. This routine assumes
    that water is distributed evenly over the surface of the cell in deriving
    the depth, and does not assume channelization. You will need to back-
    calculate channel depths for yourself using known widths at each node
    if that is what you want.

    It is VITAL you initialize this component AFTER setting boundary
    conditions.

    Note
    ----
    This is a "research grade" component, and is subject to dramatic change
    with little warning. No guarantees are made regarding its accuracy or
    utility. It is not recommended for user use yet!
    """
    _name = 'PotentialityFlowRouter'

    _input_var_names = ('topographic__elevation',
                        'water__unit_flux_in',
                        )

    _output_var_names = ('water__discharge',
                         'water__discharge_x_component',
                         'water__discharge_y_component',
                         'flow__potential',
                         'water__depth',
                         )

    _var_units = {'topographic__elevation' : 'm',
                  'water__unit_flux_in' : 'm/s',
                  'water__discharge' : 'm**3/s',
                  'water__discharge_x_component' : 'm**3/s',
                  'water__discharge_y_component' : 'm**3/s',
                  'flow__potential' : 'm**3/s',
                  'water__depth': 'm',
                  }

    _var_mapping = {'topographic__elevation' : 'node',
                  'water__unit_flux_in' : 'node',
                  'water__discharge' : 'node',
                  'water__discharge_x_component' : 'node',
                  'water__discharge_y_component' : 'node',
                  'flow__potential' : 'node',
                  'water__depth': 'node',
                  }

    _var_doc = {'topographic__elevation' : 'Land surface topographic elevation',
                  'water__unit_flux_in' : 'External volume water per area per time input to each node (e.g., rainfall rate)',
                  'water__discharge' : 'Magnitude of volumetric water flux through each node',
                  'water__discharge_x_component' : 'x component of resolved water flux through node',
                  'water__discharge_y_component' : 'y component of resolved water flux through node',
                  'flow__potential' : 'Value of the hypothetical field "K", used to force water flux to flow downhill',
                  'water__depth': 'If Manning or Chezy specified, the depth of flow in the cell, calculated assuming flow occurs over the whole surface',
                  }

    _min_slope_thresh = 1.e-24 #if your flow isn't connecting up, this probably needs to be reduced


    @use_file_name_or_kwds
    def __init__(self, grid, method='D8', flow_equation='default',
                 Chezys_C=None, Mannings_n=0.03, return_components=False,
                 **kwds):
#### give val for Chezy's C!!!
        """Initialize flow router.

        Parameters
        ----------
        grid : ModelGrid
            A grid
        params : dict
            Input parameters. Optional parameters are:

            *  `flow_equation` : options are ``default`` (default),
               ``Manning``, or ``Chezy``.  If Equation is ``Manning`` or
               ``Chezy``, you must also specify:

               *  ``Mannings_n`` (if ``Manning``) : float
               *  ``Chezys_C`` (if ``Chezy``) : float
        """
        assert RasterModelGrid in inspect.getmro(grid.__class__)
        assert grid.number_of_node_rows >= 3
        assert grid.number_of_node_columns >= 3
        assert type(return_components) is bool

        self._grid = grid
        self.equation = flow_equation
        assert self.equation in ('default', 'Chezy', 'Manning')
        if self.equation == 'Chezy':
            self.chezy_C = Chezys_C
        elif self.equation == 'Manning':
            self.manning_n = Mannings_n
        assert method in ('D8', 'D4')
        if method == 'D8':
            self.route_on_diagonals = True
        else:
            self.route_on_diagonals = False
        self.return_components = return_components

        ncols = grid.number_of_node_columns
        nrows = grid.number_of_node_rows
        # create the blank node maps; assign the topo to an internally held
        # 2D map with dummy edges:
        self.elev_raster = np.empty((nrows+2, ncols+2), dtype=float)
        self._aPP = np.zeros_like(self.elev_raster)
        self._aWW = np.zeros_like(self.elev_raster)
        self._aWP = np.zeros_like(self.elev_raster)
        self._aEE = np.zeros_like(self.elev_raster)
        self._aEP = np.zeros_like(self.elev_raster)
        self._aNN = np.zeros_like(self.elev_raster)
        self._aNP = np.zeros_like(self.elev_raster)
        self._aSS = np.zeros_like(self.elev_raster)
        self._aSP = np.zeros_like(self.elev_raster)
        self._uE = np.zeros_like(self.elev_raster)
        self._uW = np.zeros_like(self.elev_raster)
        self._uN = np.zeros_like(self.elev_raster)
        self._uS = np.zeros_like(self.elev_raster)
        self._uNE = np.zeros_like(self.elev_raster)
        self._uNW = np.zeros_like(self.elev_raster)
        self._uSE = np.zeros_like(self.elev_raster)
        self._uSW = np.zeros_like(self.elev_raster)
        self._K = np.zeros_like(self.elev_raster)

        # extras for diagonal routing:
        self._aNWNW = np.zeros_like(self.elev_raster)
        self._aNWP = np.zeros_like(self.elev_raster)
        self._aNENE = np.zeros_like(self.elev_raster)
        self._aNEP = np.zeros_like(self.elev_raster)
        self._aSESE = np.zeros_like(self.elev_raster)
        self._aSEP = np.zeros_like(self.elev_raster)
        self._aSWSW = np.zeros_like(self.elev_raster)
        self._aSWP = np.zeros_like(self.elev_raster)
        self._totalfluxout = np.empty((nrows, ncols), dtype=float)
        self._meanflux = np.zeros_like(self._totalfluxout)
        self._xdirfluxout = np.zeros_like(self._totalfluxout)
        self._ydirfluxout = np.zeros_like(self._totalfluxout)
        self._xdirfluxin = np.zeros_like(self._totalfluxout)
        self._ydirfluxin = np.zeros_like(self._totalfluxout)

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
                self.grid.add_zeros('node', 'water__discharge', dtype=float)
            except FieldError:
                pass

        # make and store a 2d reference for node BCs
        self._BCs = 4*np.ones_like(self.elev_raster)
        self._BCs[self._core].flat = self._grid.status_at_node
        BCR = self._BCs  # for conciseness below
        # these are conditions for boundary->boundary contacts AND
        # anything->closed contacts w/i the grid, both of which forbid flow
        self.boundaryboundaryN = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._Ns] > 0), np.logical_or(
                BCR[self._Ns] == 4, BCR[self._core] == 4))
        self.boundaryboundaryS = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._Ss] > 0), np.logical_or(
                BCR[self._Ss] == 4, BCR[self._core] == 4))
        self.boundaryboundaryE = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._Es] > 0), np.logical_or(
                BCR[self._Es] == 4, BCR[self._core] == 4))
        self.boundaryboundaryW = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._Ws] > 0), np.logical_or(
                BCR[self._Ws] == 4, BCR[self._core] == 4))
        self.boundaryboundaryNE = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._NEs] > 0), np.logical_or(
                BCR[self._NEs] == 4, BCR[self._core] == 4))
        self.boundaryboundaryNW = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._NWs] > 0), np.logical_or(
                BCR[self._NWs] == 4, BCR[self._core] == 4))
        self.boundaryboundarySE = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._SEs] > 0), np.logical_or(
                BCR[self._SEs] == 4, BCR[self._core] == 4))
        self.boundaryboundarySW = np.logical_or(np.logical_and(
            BCR[self._core] > 0, BCR[self._SWs] > 0), np.logical_or(
                BCR[self._SWs] == 4, BCR[self._core] == 4))
        self.notboundaries = BCR[self._core] == 0

        self.equiv_circ_diam = 2.*np.sqrt(grid.dx*grid.dy/np.pi)
        # ^this is the equivalent seen CSWidth of a cell for a flow in a
        # generic 360 direction

    def route_flow(self, **kwds):
        """
        """
        grid = self.grid
        self._K = grid.at_node['flow__potential']
        self._Qw = grid.at_node['water__discharge']
        z = grid.at_node['topographic__elevation']
        qwater_in = grid.at_node['water__unit_flux_in'].copy()
        qwater_in[grid.node_at_cell] *= grid.area_of_cell
        prev_K = self._K.copy()
        mismatch = 10000.
        # do the ortho nodes first, in isolation
        g = grid.calc_grad_at_link(z)
        if self.equation != 'default':
            g = np.sign(g)*np.sqrt(np.fabs(g))
            # ^...because both Manning and Chezy actually follow sqrt
            # slope, not slope
        # weight by face width - NO, because diags
        # g *= grid.width_of_face[grid.face_at_link]
        link_grad_at_node_w_dir = (
            g[grid.links_at_node]*grid.active_link_dirs_at_node)
        # active_link_dirs strips "wrong" face widths

        # now outgoing link grad sum
        outgoing_sum = (np.sum(link_grad_at_node_w_dir.clip(0.), axis=1) +
                        self._min_slope_thresh)
        pos_incoming_link_grads = (-link_grad_at_node_w_dir).clip(0.)

        if not self.route_on_diagonals:
            while mismatch > 1.e-6:
                K_link_ends = self._K[grid.neighbors_at_node]
                incoming_K_sum = (pos_incoming_link_grads*K_link_ends
                                  ).sum(axis=1) + self._min_slope_thresh
                self._K[:] = (incoming_K_sum + qwater_in)/outgoing_sum
                mismatch = np.sum(np.square(self._K-prev_K))
                prev_K = self._K.copy()

            upwind_K = grid.map_max_of_link_nodes_to_link(self._K)
            self._discharges_at_link = upwind_K * g
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
                incoming_K_sum = ((pos_incoming_link_grads*K_link_ends
                                   ).sum(axis=1) +
                                  (pos_incoming_diag_grads*K_diag_ends
                                   ).sum(axis=1) + self._min_slope_thresh)
                self._K[:] = (incoming_K_sum + qwater_in)/outgoing_sum
                mismatch = np.sum(np.square(self._K-prev_K))
                prev_K = self._K.copy()

            upwind_K = grid.map_max_of_link_nodes_to_link(self._K)
            upwind_diag_K = np.amax(
                (self._K[grid._diag_link_tonode],
                 self._K[grid._diag_link_fromnode]), axis=0)
            self._discharges_at_link = np.empty(grid._number_of_d8_links)
            self._discharges_at_link[:grid.number_of_links] = upwind_K * g
            self._discharges_at_link[grid.number_of_links:] = (
                upwind_diag_K * gd)
            self._discharges_at_link[grid._all_d8_inactive_links] = 0.

        np.multiply(self._K, outgoing_sum, out=self._Qw)
        # there is no sensible way to save discharges at links, if we route
        # on diagonals.
        # for now, let's make a property

        if self.return_components:
            pass

        # now process uval and vval to give the depths, if Chezy or Manning:
        if self.equation == 'Chezy':
            # Chezy: Q = C*Area*sqrt(depth*slope)
            grid.at_node['water__depth'][:] = (
                grid.at_node['flow__potential']/self.chezy_C /
                self.equiv_circ_diam)**(2./3.)
        elif self.equation == 'Manning':
            # Manning: Q = w/n*depth**(5/3)
            grid.at_node['water__depth'][:] = (
                grid.at_node['flow__potential']*self.manning_n /
                self.equiv_circ_diam)**0.6
        else:
            pass

    def run_one_step(self, **kwds):
        self.route_flow(**kwds)

    @property
    def discharges_at_links(self):
        """Return the discharges at links.

        Note that if diagonal routing, this will return number_of_d8_links.
        Otherwise, it will be number_of_links.
        """
        return self._discharges_at_link
