# -*- coding: utf-8 -*-
"""
This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router and diffuser. It needs to be broken into routing
and sed mobility methods post hoc.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

import numpy as np
from landlab import RasterModelGrid, ModelParameterDictionary, Component, FieldError
import inspect
from pylab import show, figure
from landlab.plot.imshow import imshow_node_grid

class PotentialityFlowRouter(Component):
    """
    """
    _name = 'PotentialityFlowRouter'
    
    _input_var_names = set(['topographic_elevation',
                            'water__volume_flux_in',
                            ])
    
    _output_var_names = set(['water__volume_flux_magnitude',
                             'water__volume_flux_xcomponent',
                             'water__volume_flux_ycomponent',
                             'potentiality_field',
                             'water__volume_flux',
                             ])
                             
    _var_units = {'topographic_elevation' : 'm',
                  'water__volume_flux_in' : 'm**3/s',
                  'water__volume_flux_magnitude' : 'm**3/s',
                  'water__volume_flux_xcomponent' : 'm**3/s',
                  'water__volume_flux_ycomponent' : 'm**3/s',
                  'potentiality_field' : 'm**3/s',
                  'water__volume_flux' : 'm**3/s',
                  }
    
    _var_mapping = {'topographic_elevation' : 'node',
                  'water__volume_flux_in' : 'node',
                  'water__volume_flux_magnitude' : 'node',
                  'water__volume_flux_xcomponent' : 'node',
                  'water__volume_flux_ycomponent' : 'node',
                  'potentiality_field' : 'node',
                  'water__volume_flux' : 'link',
                  }
    
    _var_defs = {'topographic_elevation' : 'Land surface topographic elevation',
                  'water__volume_flux_in' : 'External volume water input to each node (e.g., rainfall)',
                  'water__volume_flux_magnitude' : 'Magnitude of volumetric water flux through each node',
                  'water__volume_flux_xcomponent' : 'x component of resolved water flux through node',
                  'water__volume_flux_ycomponent' : 'y component of resolved water flux through node',
                  'potentiality_field' : 'Value of the hypothetical field "K", used to force water flux to flow downhill',
                  'water__volume_flux' : 'Water fluxes on links',
                  }

    
    def __init__(self, grid, params):
        self.initialize(grid, params)
    
    
    def initialize(self, grid, params):
        assert RasterModelGrid in inspect.getmro(grid.__class__)
        assert grid.number_of_node_rows >= 3
        assert grid.number_of_node_columns >= 3
        
        self._grid = grid
        
        if type(params) == str:
            input_dict = ModelParameterDictionary(params)
        else:
            assert type(params) == dict
            input_dict = params
        
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
        
        for out_field in self._output_var_names:
            if self._var_mapping[out_field]=='node':
                try:
                    self._grid.at_node[out_field]
                except FieldError:
                    self._grid.at_node[out_field] = np.empty(self._grid.number_of_nodes, dtype=float)
            elif self._var_mapping[out_field]=='link':
                try:
                    self._grid.at_link[out_field]
                except FieldError:
                    self._grid.at_link[out_field] = np.empty(self._grid.number_of_links, dtype=float)
        
        #make and store a 2d reference for node BCs
        self._BCs = 4*np.ones_like(self.elev_raster)
        self._BCs[self._core].flat = self._grid.get_node_status()
        BCR = self._BCs #for conciseness below
        #these are conditions for boundary-boundary contacts AND core-closed contacts w/i the grid, both of which forbid flow
        self.boundaryboundaryN = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Ns]>0), np.logical_and(BCR[self._core]==0, BCR[self._Ns]==4))
        self.boundaryboundaryS = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Ss]>0), np.logical_and(BCR[self._core]==0, BCR[self._Ss]==4))
        self.boundaryboundaryE = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Es]>0), np.logical_and(BCR[self._core]==0, BCR[self._Es]==4))
        self.boundaryboundaryW = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._Ws]>0), np.logical_and(BCR[self._core]==0, BCR[self._Ws]==4))
        self.boundaryboundaryNE = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._NEs]>0), np.logical_and(BCR[self._core]==0, BCR[self._NEs]==4))
        self.boundaryboundaryNW = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._NWs]>0), np.logical_and(BCR[self._core]==0, BCR[self._NWs]==4))
        self.boundaryboundarySE = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._SEs]>0), np.logical_and(BCR[self._core]==0, BCR[self._SEs]==4))
        self.boundaryboundarySW = np.logical_or(np.logical_and(BCR[self._core]>0, BCR[self._SWs]>0), np.logical_and(BCR[self._core]==0, BCR[self._SWs]==4))
        self.notboundaries = BCR[self._core] == 0
        
    
    def route_flow(self, route_on_diagonals=True):
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
        corecore = self._corecore
        one_over_dx = 1./self._grid.dx
        one_over_dy = 1./self._grid.dy
        one_over_diagonal = 1./np.sqrt(self._grid.dx**2+self._grid.dy**2)
        qwater_in = self._grid.at_node['water__volume_flux_in'].reshape((self._grid.number_of_node_rows, self._grid.number_of_node_columns))
        prev_K = K.copy()
        BCR = self._BCs
        bbN = self.boundaryboundaryN
        bbS = self.boundaryboundaryS
        bbE = self.boundaryboundaryE
        bbW = self.boundaryboundaryW
        
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
        hR[1:-1,1:-1].flat = self._grid.at_node['topographic_elevation']
        
        #update the dummy edges of our variables - these all act as closed nodes (the inner, true boundaries are handled elsewhere... or they should be closed anyway!):
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
        
        
        
#        #the routing routine should not "see" a flat - flow should get routed across it unimpeded
#        #find the flat nodes; have to iterate
#        flatE = np.isclose(aEE[corecore],aEP[corecore])
#        flatW = np.isclose(aWW[corecore],aWP[corecore])
#        flatN = np.isclose(aNN[corecore],aNP[corecore])
#        flatS = np.isclose(aSS[corecore],aSP[corecore])
#        sumflatE = np.sum(flatE)
#        sumflatW = np.sum(flatW)
#        sumflatN = np.sum(flatN)
#        sumflatS = np.sum(flatS)
#        oldsumflatE = sumflatE
#        oldsumflatW = sumflatW
#        oldsumflatN = sumflatN
#        oldsumflatS = sumflatS
#        while sum((sumflatE,sumflatN,sumflatS,sumflatW))>0:
#            print 'Cardinals: ', sum((sumflatE,sumflatN,sumflatS,sumflatW))
#            #figure(1)
#            #imshow_node_grid(self._grid, aNN[core])
#            #figure(2)
#            #imshow_node_grid(self._grid, aNP[core])
#            #figure(3)
#            #imshow_node_grid(self._grid, aSS[core])
#            #figure(4)
#            #imshow_node_grid(self._grid, aSP[core])
#            #show()
#            #oldEP=aEP[core][flatE].copy()
#            #oldWP=aWP[core][flatW].copy()
#            #oldNP=aNP[core][flatN].copy()
#            #oldSP=aSP[core][flatS].copy()
#            #E
#            flow_is_inE = np.logical_and(aWW[corecore]>0., flatE)
#            flow_is_outE = np.logical_and(aWP[corecore]>0., flatE) #note we won't be switching if it's flat
#            flow_is_inW = np.logical_and(aEE[corecore]>0., flatW)
#            flow_is_outW = np.logical_and(aEP[corecore]>0., flatW)
#            flow_is_inN = np.logical_and(aSS[corecore]>0., flatN)
#            flow_is_outN = np.logical_and(aSP[corecore]>0., flatN)
#            flow_is_inS = np.logical_and(aNN[corecore]>0., flatS)
#            flow_is_outS = np.logical_and(aNP[corecore]>0., flatS)
#            
#            aEP[corecore][flow_is_inE] = aWW[corecore][flow_is_inE] #if it's flat, aEE must already be 0...
#            aEE[corecore][flow_is_outE] = aWP[corecore][flow_is_outE]
#            #propagate to its counterpart to the E:
#            aWW[Es][core][flow_is_inE] = aWW[corecore][flow_is_inE]
#            aWP[Es][core][flow_is_outE] = aWP[corecore][flow_is_outE]
#            
#            aWP[corecore][flow_is_inW] = aEE[corecore][flow_is_inW]
#            aWW[corecore][flow_is_outW] = aEP[corecore][flow_is_outW]
#            aEE[Ws][core][flow_is_inW] = aEE[corecore][flow_is_inW]
#            aEP[Ws][core][flow_is_outW] = aEP[corecore][flow_is_outW]
#            
#            aNP[corecore][flow_is_inN] = aSS[corecore][flow_is_inN]
#            aNN[corecore][flow_is_outN] = aSP[corecore][flow_is_outN]
#            aSS[Ns][core][flow_is_inN] = aSS[corecore][flow_is_inN]
#            aSP[Ns][core][flow_is_outN] = aSP[corecore][flow_is_outN]
#            
#            aSP[corecore][flow_is_inS] = aNN[corecore][flow_is_inS]
#            aSS[corecore][flow_is_outS] = aNP[corecore][flow_is_outS]
#            aNN[Ss][core][flow_is_inS] = aNN[corecore][flow_is_inS]
#            aNP[Ss][core][flow_is_outS] = aNP[corecore][flow_is_outS]
#            #figure(1)
#            #imshow_node_grid(self._grid, aNN[core])
#            #figure(2)
#            #imshow_node_grid(self._grid, aNP[core])
#            #figure(3)
#            #imshow_node_grid(self._grid, aSS[core])
#            #figure(4)
#            #imshow_node_grid(self._grid, aSP[core])
#            #show()
#            
#            ###not convinced the following is the best way to go about this... Could leave occasional blanks?
#            flatE = np.isclose(aEE[corecore],aEP[corecore])
#            sumflatE = np.sum(flatE)
#            if sumflatE==oldsumflatE:
#                sumflatE=0 #no further change is happening; at least one fully blank row. Stop iterating
#            oldsumflatE = sumflatE
#            #oldsumflatE = min((sumflatE, oldsumflatE))
#            flatW = np.isclose(aWW[corecore],aWP[corecore])
#            sumflatW = np.sum(flatW)
#            if sumflatW==oldsumflatW:
#                sumflatW=0
#            oldsumflatW = sumflatW
#            #oldsumflatW = min((sumflatW, oldsumflatW))
#            flatN = np.isclose(aNN[corecore],aNP[corecore])
#            sumflatN = np.sum(flatN)
#            if sumflatN==oldsumflatN:
#                sumflatN=0
#            oldsumflatN = sumflatN
#            #oldsumflatN = min((sumflatN,oldsumflatN))
#            flatS = np.isclose(aSS[corecore],aSP[corecore])
#            sumflatS = np.sum(flatS)
#            if sumflatS==oldsumflatS:
#                sumflatS=0
#            oldsumflatS = sumflatS
#            #oldsumflatS = min((sumflatS, oldsumflatS))

        #figure(1)
        #imshow_node_grid(self._grid, aNN[core])
        #figure(2)
        #imshow_node_grid(self._grid, aNP[core])
        #figure(3)
        #imshow_node_grid(self._grid, aSS[core])
        #figure(4)
        #imshow_node_grid(self._grid, aSP[core])
        #show()
        #figure(1)
        #imshow_node_grid(self._grid, aEE[core])
        #figure(2)
        #imshow_node_grid(self._grid, aEP[core])
        #figure(3)
        #imshow_node_grid(self._grid, aWW[core])
        #figure(4)
        #imshow_node_grid(self._grid, aWP[core])
        #show()
        #this DOES NOT work
        
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
            
        
#            flatNE = np.isclose(aNENE[corecore],aNEP[corecore])
#            flatSE = np.isclose(aSESE[corecore],aSEP[corecore])
#            flatSW = np.isclose(aSWSW[corecore],aSWP[corecore])
#            flatNW = np.isclose(aNWNW[corecore],aNWP[corecore])
#            sumflatNE = np.sum(flatNE)
#            sumflatSE = np.sum(flatSE)
#            sumflatSW = np.sum(flatSW)
#            sumflatNW = np.sum(flatNW)
#            oldsumflatNE = sumflatNE
#            oldsumflatSE = sumflatSE
#            oldsumflatSW = sumflatSW
#            oldsumflatNW = sumflatNW
#            while sum((sumflatNE,sumflatNW,sumflatSE,sumflatSW))>0:
#                print 'diags: ', sum((sumflatNE,sumflatNW,sumflatSE,sumflatSW))
#                flow_is_inNE = np.logical_and(aNENE[corecore]>0., flatNE)
#                flow_is_outNE = np.logical_and(aNEP[corecore]>0., flatNE) #note we won't be switching if it's flat
#                flow_is_inSE = np.logical_and(aSESE[corecore]>0., flatSE)
#                flow_is_outSE = np.logical_and(aSEP[corecore]>0., flatSE)
#                flow_is_inSW = np.logical_and(aSWSW[corecore]>0., flatSW)
#                flow_is_outSW = np.logical_and(aSWP[corecore]>0., flatSW)
#                flow_is_inNW = np.logical_and(aNWNW[corecore]>0., flatNW)
#                flow_is_outNW = np.logical_and(aNWP[corecore]>0., flatNW)
#                
#                aNEP[corecore][flow_is_inNE] = aSWSW[corecore][flow_is_inNE] #if it's flat, aEE must already be 0...
#                aNENE[corecore][flow_is_outNE] = aSWP[corecore][flow_is_outNE]
#                #propagate to its counterpart to the E:
#                aSWSW[NEs][core][flow_is_inNE] = aSWSW[corecore][flow_is_inNE]
#                aSWP[NEs][core][flow_is_outNE] = aSWP[corecore][flow_is_outNE]
#                
#                aNWP[corecore][flow_is_inNW] = aSESE[corecore][flow_is_inNW]
#                aNWNW[corecore][flow_is_outNW] = aSEP[corecore][flow_is_outNW]
#                aSESE[NWs][core][flow_is_inNW] = aSESE[corecore][flow_is_inNW]
#                aSEP[NWs][core][flow_is_outNW] = aSEP[corecore][flow_is_outNW]
#                
#                aSEP[corecore][flow_is_inSE] = aSWSW[corecore][flow_is_inSE]
#                aSESE[corecore][flow_is_outSE] = aSWP[corecore][flow_is_outSE]
#                aSWSW[SEs][core][flow_is_inSE] = aSWSW[corecore][flow_is_inSE]
#                aSWP[SEs][core][flow_is_outSE] = aSWP[corecore][flow_is_outSE]
#                
#                aSWP[corecore][flow_is_inSW] = aSESE[corecore][flow_is_inSW]
#                aSWSW[corecore][flow_is_outSW] = aSEP[corecore][flow_is_outSW]
#                aSESE[SWs][core][flow_is_inSW] = aSESE[corecore][flow_is_inSW]
#                aSEP[SWs][core][flow_is_outSW] = aSEP[corecore][flow_is_outSW]
#                
#                flatNE = np.isclose(aNENE[corecore],aNEP[corecore])
#                sumflatNE = np.sum(flatNE)
#                if sumflatNE>=oldsumflatNE:
#                    sumflatNE=0 #no further change is happening; at least one fully blank row. Stop iterating
#                oldsumflatNE = min((sumflatNE, oldsumflatNE))
#                flatNW = np.isclose(aNWNW[corecore],aNWP[corecore])
#                sumflatNW = np.sum(flatNW)
#                if sumflatNW>=oldsumflatNW:
#                    sumflatNW=0
#                oldsumflatNW = min((sumflatNW, oldsumflatNW))
#                flatSE = np.isclose(aSESE[corecore],aSEP[corecore])
#                sumflatSE = np.sum(flatSE)
#                if sumflatSE>=oldsumflatSE:
#                    sumflatSE=0
#                oldsumflatSE = min((sumflatSE, oldsumflatSE))
#                flatSW = np.isclose(aSWSW[corecore],aSWP[corecore])
#                sumflatSW = np.sum(flatSW)
#                if sumflatSW>=oldsumflatSW:
#                    sumflatSW=0
#                oldsumflatSW = min((sumflatSW, oldsumflatSW))
        
        
        if not route_on_diagonals:
            aPP[core] = aWP[core]+aEP[core]+aSP[core]+aNP[core]+1.e-12
        else:
            aPP[core] = (aWP[core]+aEP[core]+aSP[core]+aNP[core]
                        +aNEP[core]+aSEP[core]+aSWP[core]+aNWP[core])+1.e-12
                        
        mismatch = 10000.
        self.loops_needed = 0

        #this explicit solution could easily be replaced with a better solver        
        while mismatch>1.e-3:
            if not route_on_diagonals:
                K[core] = (aWW[core]*K[Ws]+aEE[core]*K[Es]+aSS[core]*K[Ss]+aNN[core]*K[Ns]
                                    +qwater_in)/aPP[core]
#####This needs to be qw/2 because...?
            else:
                K[core] = (aWW[core]*K[Ws]+aEE[core]*K[Es]+aSS[core]*K[Ss]+aNN[core]*K[Ns]
                          +aNENE[core]*K[NEs]+aSESE[core]*K[SEs]
                          +aSWSW[core]*K[SWs]+aNWNW[core]*K[NWs]
                          +qwater_in)/aPP[core]

            mismatch = np.sum(np.square(K[core]-prev_K[core]))
            self.loops_needed += 1
            prev_K = K.copy()
            print mismatch
            
            for BC in (K,):
                BC[0,1:-1] = BC[1,1:-1]
                BC[-1,1:-1] = BC[-2,1:-1]
                BC[1:-1,0] = BC[1:-1,1]
                BC[1:-1,-1] = BC[1:-1,-2]
                BC[(0,-1,0,-1),(0,-1,-1,0)] = BC[(1,-2,1,-2),(1,-2,-2,1)]
        
        if not route_on_diagonals:
            uW[core] = aWW[core]*K[Ws]-aWP[core]*K[core]
            uE[core] = -aEE[core]*K[Es]+aEP[core]*K[core]
            uN[core] = -aNN[core]*K[Ns]+aNP[core]*K[core]
            uS[core] = aSS[core]*K[Ss]-aSP[core]*K[core]
        else:
            prefactor_to_x = np.arctan(self._grid.dy/self._grid.dx)
            prefactor_to_y = np.arctan(self._grid.dx/self._grid.dy)
            uW[core] = aWW[core]*K[Ws]-aWP[core]*K[core] + (aNWNW[core]*K[NWs]-aNWP[core]*K[core]+aSWSW[core]*K[SWs]-aSWP[core]*K[core])*prefactor_to_x
            uE[core] = -aEE[core]*K[Es]+aEP[core]*K[core] + (-aNENE[core]*K[NEs]+aNEP[core]*K[core]-aSESE[core]*K[SEs]+aSEP[core]*K[core])*prefactor_to_x
            uN[core] = -aNN[core]*K[Ns]+aNP[core]*K[core] + (-aNWNW[core]*K[NWs]+aNWP[core]*K[core]-aNENE[core]*K[NEs]+aNEP[core]*K[core])*prefactor_to_y
            uS[core] = aSS[core]*K[Ss]-aSP[core]*K[core] + (aSWSW[core]*K[SWs]-aSWP[core]*K[core]+aSESE[core]*K[SEs]-aSEP[core]*K[core])*prefactor_to_y
        
        uval = uW[core]+uE[core]
        vval = uN[core]+uS[core]
        uval[self.notboundaries] /= 2.
        vval[self.notboundaries] /= 2. #because we want mean values, but not where one of the values is 0 because it's a BC
        
        #save the output
        self._grid.at_link['water__volume_flux'][self._grid.node_links()[0]] = uS[core].flat #[S,W,N,E], (4,nnodes)
        self._grid.at_link['water__volume_flux'][self._grid.node_links()[1]] = uW[core].flat
        self._grid.at_link['water__volume_flux'][self._grid.node_links()[2]] = uN[core].flat
        self._grid.at_link['water__volume_flux'][self._grid.node_links()[3]] = uE[core].flat
        self._grid.at_node['potentiality_field'][:] = K[core].flat
        self._grid.at_node['water__volume_flux_xcomponent'][:] = uval.flat
        self._grid.at_node['water__volume_flux_ycomponent'][:] = vval.flat
        self._grid.at_node['water__volume_flux_magnitude'][:] = np.sqrt(uval*uval+vval*vval).flat
        
####Outstanding issues - 1. BC handling (flow comes back in from edges ATM); 2. flow routing on flats (?)