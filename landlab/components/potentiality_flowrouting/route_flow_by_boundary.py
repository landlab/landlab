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
from landlab import RasterModelGrid, ModelParameterDictionary, Component, FieldError
import inspect

##########################
#THIS DOESN'T ACCOUNT FOR CELL AREA YET!!!!!!
##########################

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
                         'water__discharge',
                         )

    _var_units = {'topographic__elevation' : 'm',
                  'water__unit_flux_in' : 'm/s',
                  'water__discharge' : 'm**3/s',
                  'water__discharge_x_component' : 'm**3/s',
                  'water__discharge_y_component' : 'm**3/s',
                  'flow__potential' : 'm**3/s',
                  'water__depth': 'm',
                  'water__discharge' : 'm**3/s',
                  }

    _var_mapping = {'topographic__elevation' : 'node',
                  'water__unit_flux_in' : 'node',
                  'water__discharge' : 'node',
                  'water__discharge_x_component' : 'node',
                  'water__discharge_y_component' : 'node',
                  'flow__potential' : 'node',
                  'water__depth': 'node',
                  'water__discharge' : 'link',
                  }

    _var_doc = {'topographic__elevation' : 'Land surface topographic elevation',
                  'water__unit_flux_in' : 'External volume water per area per time input to each node (e.g., rainfall rate)',
                  'water__discharge' : 'Magnitude of volumetric water flux through each node',
                  'water__discharge_x_component' : 'x component of resolved water flux through node',
                  'water__discharge_y_component' : 'y component of resolved water flux through node',
                  'flow__potential' : 'Value of the hypothetical field "K", used to force water flux to flow downhill',
                  'water__depth': 'If Manning or Chezy specified, the depth of flow in the cell, calculated assuming flow occurs over the whole surface',
                  'water__discharge' : 'Water fluxes on links',
                  }

    _min_slope_thresh = 1.e-24 #if your flow isn't connecting up, this probably needs to be reduced


    def __init__(self, grid, params):
        self.initialize(grid, params)


    def initialize(self, grid, params):
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

        self._grid = grid

        if type(params) == str:
            input_dict = ModelParameterDictionary(params)
        else:
            assert type(params) == dict
            input_dict = params

        #ingest the inputs
        try:
            self.equation = input_dict['flow_equation']
        except KeyError:
            self.equation = 'default'
        assert self.equation in ('default', 'Chezy', 'Manning')
        if self.equation == 'Chezy':
            self.chezy_C = float(input_dict["Chezys_C"])
        if self.equation == 'Manning':
            self.manning_n = float(input_dict["Mannings_n"])

        ncols = grid.number_of_node_columns
        nrows = grid.number_of_node_rows
        #create the blank node maps; assign the topo to an internally held 2D map with dummy edges:
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

        #extras for diagonal routing:
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

        #setup slices for use in IDing the neighbors
        self._Es = (slice(1,-1),slice(2,ncols+2))
        self._NEs = (slice(2,nrows+2),slice(2,ncols+2))
        self._Ns = (slice(2,nrows+2),slice(1,-1))
        self._NWs = (slice(2,nrows+2),slice(0,-2))
        self._Ws = (slice(1,-1),slice(0,-2))
        self._SWs = (slice(0,-2),slice(0,-2))
        self._Ss = (slice(0,-2),slice(1,-1))
        self._SEs = (slice(0,-2),slice(2,ncols+2))
        self._core = (slice(1,-1),slice(1,-1))
        self._corecore = (slice(2,-2),slice(2,-2)) #the actual, LL-sense core (interior) nodes of the grid

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
            try:
                self.grid.add_zeros('link', 'water__discharge', dtype=float)
            except FieldError:
                pass

        #make and store a 2d reference for node BCs
        self._BCs = 4*np.ones_like(self.elev_raster)
        self._BCs[self._core].flat = self._grid.status_at_node
        BCR = self._BCs #for conciseness below
        #these are conditions for boundary->boundary contacts AND anything->closed contacts w/i the grid, both of which forbid flow
        self.boundaryboundaryN = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Ns]>0), np.logical_or(BCR[self._Ns]==4, BCR[self._core]==4))
        self.boundaryboundaryS = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Ss]>0), np.logical_or(BCR[self._Ss]==4, BCR[self._core]==4))
        self.boundaryboundaryE = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Es]>0), np.logical_or(BCR[self._Es]==4, BCR[self._core]==4))
        self.boundaryboundaryW = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Ws]>0), np.logical_or(BCR[self._Ws]==4, BCR[self._core]==4))
        self.boundaryboundaryNE = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._NEs]>0), np.logical_or(BCR[self._NEs]==4, BCR[self._core]==4))
        self.boundaryboundaryNW = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._NWs]>0), np.logical_or(BCR[self._NWs]==4, BCR[self._core]==4))
        self.boundaryboundarySE = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._SEs]>0), np.logical_or(BCR[self._SEs]==4, BCR[self._core]==4))
        self.boundaryboundarySW = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._SWs]>0), np.logical_or(BCR[self._SWs]==4, BCR[self._core]==4))
        self.notboundaries = BCR[self._core] == 0

        self.equiv_circ_diam = 2.*np.sqrt(grid.dx*grid.dy/np.pi)
        #this is the equivalent seen CSWidth of a cell for a flow in a generic 360 direction


    def route_flow(self, route_on_diagonals=True, return_components=False):
        """
        """
        #create aliases for ease
        hR = self.elev_raster
        aPP = self._aPP
        aWW = self._aWW
        aWP = self._aWP
        aEE = self._aEE
        aEP = self._aEP
        aNN = self._aNN
        aNP = self._aNP
        aSS = self._aSS
        aSP = self._aSP
        totalfluxout = self._totalfluxout
        meanflux = self._meanflux
        uE = self._uE
        uW = self._uW
        uN = self._uN
        uS = self._uS
        K = self._K
        Es = self._Es
        NEs = self._NEs
        Ns = self._Ns
        NWs = self._NWs
        Ws = self._Ws
        SWs = self._SWs
        Ss = self._Ss
        SEs = self._SEs
        core = self._core
        one_over_dx = 1./self._grid.dx
        one_over_dy = 1./self._grid.dy
        one_over_diagonal = 1./np.sqrt(self._grid.dx**2+self._grid.dy**2)
        qwater_in = self._grid.at_node['water__unit_flux_in'].reshape((self._grid.number_of_node_rows, self._grid.number_of_node_columns))
        prev_K = K.copy()
        bbN = self.boundaryboundaryN
        bbS = self.boundaryboundaryS
        bbE = self.boundaryboundaryE
        bbW = self.boundaryboundaryW
        if return_components:
            flux_error = np.zeros_like(self._totalfluxout)

        if route_on_diagonals:
            #extras for diagonal routing:
            aNWNW = self._aNWNW
            aNWP = self._aNWP
            aNENE = self._aNENE
            aNEP = self._aNEP
            aSESE = self._aSESE
            aSEP = self._aSEP
            aSWSW = self._aSWSW
            aSWP = self._aSWP

            bbNE = self.boundaryboundaryNE
            bbNW = self.boundaryboundaryNW
            bbSE = self.boundaryboundarySE
            bbSW = self.boundaryboundarySW

        #paste in the elevs
        hR[1:-1,1:-1].flat = self._grid.at_node['topographic__elevation']

        #update the dummy edges of our variables - these all act as closed nodes (the inner, true boundaries are handled elsewhere... or they should be closed anyway!):
        #note this isn't sufficient of we have diagonals turned on, as flow can still occur on them
        hR[0,1:-1] = hR[1,1:-1]
        hR[-1,1:-1] = hR[-2,1:-1]
        hR[1:-1,0] = hR[1:-1,1]
        hR[1:-1,-1] = hR[1:-1,-2]
        hR[(0,-1,0,-1),(0,-1,-1,0)] = hR[(1,-2,1,-2),(1,-2,-2,1)]

        aNN[core] = (-hR[core]+hR[Ns]).clip(0.)
        aNP[core] =  (hR[core]-hR[Ns]).clip(0.)
        aSS[core] = (-hR[core]+hR[Ss]).clip(0.)
        aSP[core] =  (hR[core]-hR[Ss]).clip(0.)
        aEE[core] = (-hR[core]+hR[Es]).clip(0.)
        aEP[core] =  (hR[core]-hR[Es]).clip(0.)
        aWW[core] = (-hR[core]+hR[Ws]).clip(0.)
        aWP[core] =  (hR[core]-hR[Ws]).clip(0.)
        aNN[core] *= one_over_dy
        aNP[core] *= one_over_dy
        aSS[core] *= one_over_dy
        aSP[core] *= one_over_dy
        aEE[core] *= one_over_dx
        aEP[core] *= one_over_dx
        aWW[core] *= one_over_dx
        aWP[core] *= one_over_dx

        #disable lateral flow betw boundary nodes:
        aEE[core][bbE] = 0.
        aEP[core][bbE] = 0.
        aWW[core][bbW] = 0.
        aWP[core][bbW] = 0.
        aNN[core][bbN] = 0.
        aNP[core][bbN] = 0.
        aSS[core][bbS] = 0.
        aSP[core][bbS] = 0.
        #diag correction happens below

        if self.equation != 'default':
            for grad in (aEE,aEP,aWW,aWP,aNN,aNP,aSS,aSP):
                np.sqrt(grad[core], out=grad[core]) #...because both Manning and Chezy actually follow sqrt slope, not slope

        if route_on_diagonals:
            #adding diagonals
            aNENE[core] = (-hR[core]+hR[NEs]).clip(0.)
            aNEP[core] =   (hR[core]-hR[NEs]).clip(0.)
            aSESE[core] = (-hR[core]+hR[SEs]).clip(0.)
            aSEP[core] =   (hR[core]-hR[SEs]).clip(0.)
            aSWSW[core] = (-hR[core]+hR[SWs]).clip(0.)
            aSWP[core] =   (hR[core]-hR[SWs]).clip(0.)
            aNWNW[core] = (-hR[core]+hR[NWs]).clip(0.)
            aNWP[core] =   (hR[core]-hR[NWs]).clip(0.)
            aNENE[core] *= one_over_diagonal
            aNEP[core] *= one_over_diagonal
            aSESE[core] *= one_over_diagonal
            aSEP[core] *= one_over_diagonal
            aSWSW[core] *= one_over_diagonal
            aSWP[core] *= one_over_diagonal
            aNWNW[core] *= one_over_diagonal
            aNWP[core] *= one_over_diagonal

            #disable lateral flow betw boundaries...
            aNENE[core][bbNE] = 0.
            aNEP[core][bbNE] = 0.
            aNWNW[core][bbNW] = 0.
            aNWP[core][bbNW] = 0.
            aSESE[core][bbSE] = 0.
            aSEP[core][bbSE] = 0.
            aSWSW[core][bbSW] = 0.
            aSWP[core][bbSW] = 0.

            if self.equation != 'default':
                for grad in (aNENE,aNEP,aNWNW,aNWP,aSESE,aSEP,aSWSW,aSWP):
                    np.sqrt(grad[core], out=grad[core]) #...because both Manning and Chezy actually follow sqrt slope, not slope

        if not route_on_diagonals:
            aPP[core] = aWP[core]+aEP[core]+aSP[core]+aNP[core]+self._min_slope_thresh
        else:
            aPP[core] = (aWP[core]+aEP[core]+aSP[core]+aNP[core]
                        +aNEP[core]+aSEP[core]+aSWP[core]+aNWP[core])+self._min_slope_thresh

        mismatch = 10000.
        self.loops_needed = 0

        #this explicit solution could easily be replaced with a better solver
        while mismatch>1.e-3:
            if not route_on_diagonals:
                K[core] = (aWW[core]*K[Ws]+aEE[core]*K[Es]+aSS[core]*K[Ss]+aNN[core]*K[Ns]
                                    +qwater_in)/aPP[core]
            else:
                K[core] = (aWW[core]*K[Ws]+aEE[core]*K[Es]+aSS[core]*K[Ss]+aNN[core]*K[Ns]
                          +aNENE[core]*K[NEs]+aSESE[core]*K[SEs]
                          +aSWSW[core]*K[SWs]+aNWNW[core]*K[NWs]
                          +qwater_in)/aPP[core]

            mismatch = np.sum(np.square(K[core]-prev_K[core]))
            self.loops_needed += 1
            prev_K = K.copy()
            #print mismatch

            for BC in (K,):
                BC[0,1:-1] = BC[1,1:-1]
                BC[-1,1:-1] = BC[-2,1:-1]
                BC[1:-1,0] = BC[1:-1,1]
                BC[1:-1,-1] = BC[1:-1,-2]
                BC[(0,-1,0,-1),(0,-1,-1,0)] = BC[(1,-2,1,-2),(1,-2,-2,1)]

        if route_on_diagonals:
            outdirs = (aNP,aSP,aEP,aWP,aNWP,aNEP,aSWP,aSEP)
            indirs = ((aNN,K[Ns]),(aSS,K[Ss]),(aEE,K[Es]),(aWW,K[Ws]),(aNWNW,K[NWs]),(aNENE,K[NEs]),(aSWSW,K[SWs]),(aSESE,K[SEs]))
        else:
            outdirs = (aNP,aSP,aEP,aWP)
            indirs = ((aNN,K[Ns]),(aSS,K[Ss]),(aEE,K[Es]),(aWW,K[Ws]))

        totalfluxout.fill(0.)
        for array in outdirs:
            totalfluxout += array[core]
        totalfluxout *= K[core]
        no_outs = np.equal(totalfluxout,0.)
        #this WON'T work if there is no out flux, i.e., a BC, so -
        if np.any(no_outs):
            for (inarray, inK) in indirs:
                totalfluxout[no_outs] += inarray[core][no_outs]*inK[no_outs]
        yes_outs = np.logical_not(no_outs)
        meanflux[:] = totalfluxout
        meanflux[yes_outs] -= 0.5*qwater_in[yes_outs]
        #so note the BC nodes get the value of the fluxes they RECEIVE, without any addition of flux
        #& also, nodes at the TOP of the network, with no explicit ins, get averaged values, not just their out values


        if return_components:
            #this takes a nontrivial number of additional calculations, so is optional
            #we aim to return the MEAN flow direction.
            if route_on_diagonals:
                prefactor_to_x = np.cos(np.arctan(self._grid.dy/self._grid.dx))
                prefactor_to_y = np.sin(np.arctan(self._grid.dy/self._grid.dx))
                xdir_mod_out = (0.,0.,1.,-1.,-prefactor_to_x,prefactor_to_x,-prefactor_to_x,prefactor_to_x)
                ydir_mod_out = (1.,-1.,0.,0.,prefactor_to_y,prefactor_to_y,-prefactor_to_y,-prefactor_to_y)
                xdir_mod_in = (0.,0.,-1.,1.,prefactor_to_x,-prefactor_to_x,prefactor_to_x,-prefactor_to_x)
                ydir_mod_in = (-1.,1.,0.,0.,-prefactor_to_y,-prefactor_to_y,prefactor_to_y,prefactor_to_y)
            else:
                xdir_mod_out = (0.,0.,1.,-1.)
                ydir_mod_out = (1.,-1.,0.,0.)
                xdir_mod_in = (0.,0.,-1.,1.)
                ydir_mod_in = (-1.,1.,0.,0.)
            #out is easier, so do it first:
            self._xdirfluxout.fill(0.)
            self._ydirfluxout.fill(0.)
            for (array,mod) in zip(outdirs,xdir_mod_out):
                self._xdirfluxout += mod*array[core]
            for (array,mod) in zip(outdirs,ydir_mod_out):
                self._ydirfluxout += mod*array[core]
            #now in
            self._xdirfluxin.fill(0.)
            self._ydirfluxin.fill(0.)
            for ((array,Kin),mod) in zip(indirs,xdir_mod_in):
                self._xdirfluxin += mod*array[core]*Kin
            for ((array,Kin),mod) in zip(indirs,ydir_mod_in):
                self._ydirfluxin += mod*array[core]*Kin
            #modify for doing the means:
            self._xdirfluxin[yes_outs] *= 0.5
            self._ydirfluxin[yes_outs] *= 0.5
            #flux outs all get adjusted:
            self._xdirfluxout *= 0.5
            self._ydirfluxout *= 0.5
            #NB: handling of degenerate cases where there is no angle defined (though rare) could be problematic.
            # => assign an arbitary zero angle in these cases
            #so...
            mean_x = self._xdirfluxin + self._xdirfluxout
            mean_y = self._ydirfluxin + self._ydirfluxout
            apparent_flux = np.sqrt(mean_x*mean_x + mean_y*mean_y)
            degenerate_fluxes = np.equal(apparent_flux,0.)
            if np.any(degenerate_fluxes):
                #THIS PART HAS NOT BEEN TESTED
                defined_fluxes = np.logical_not(degenerate_fluxes)
                flux_error[defined_fluxes] = meanflux[defined_fluxes]/apparent_flux[defined_fluxes]
                #zeros in the degenerate slots
                self._grid.at_node['water__discharge_x_component'][:] = (mean_x*flux_error).flat
                #now modify so we get the right answer for total flux, pointing arbitrarily N
                self._grid.at_node['water__discharge_y_component'][defined_fluxes] = (mean_y[defined_fluxes]*flux_error[defined_fluxes]).flat
                self._grid.at_node['water__discharge_y_component'][degenerate_fluxes] = meanflux[degenerate_fluxes].flat
            else:
                flux_error[:] = meanflux/apparent_flux
                self._grid.at_node['water__discharge_y_component'][:] = (mean_y*flux_error).flat
                self._grid.at_node['water__discharge_x_component'][:] = (mean_x*flux_error).flat

        #save the output
        self._grid.at_link['water__discharge'][self._grid.links_at_node[:, 3]] = uS[core].flat
        self._grid.at_link['water__discharge'][self._grid.links_at_node[:, 2]] = uW[core].flat
        self._grid.at_link['water__discharge'][self._grid.links_at_node[:, 1]] = uN[core].flat
        self._grid.at_link['water__discharge'][self._grid.links_at_node[:, 0]] = uE[core].flat
        self._grid.at_node['flow__potential'][:] = K[core].flat
        self._grid.at_node['water__discharge'][:] = meanflux.flat
        #the x,y components are created above, in the if statement

        #now process uval and vval to give the depths, if Chezy or Manning:
        if self.equation == 'Chezy':
            #Chezy: Q = C*Area*sqrt(depth*slope)
            self._grid.at_node['water__depth'][:] = (self._grid.at_node['flow__potential']/self.chezy_C/self.equiv_circ_diam)**(2./3.)
        elif self.equation == 'Manning':
            #Manning: Q = w/n*depth**(5/3)
            self._grid.at_node['water__depth'][:] = (self._grid.at_node['flow__potential']*self.manning_n/self.equiv_circ_diam)**0.6
        else:
            pass
